(ns encode-manager.core
  (:gen-class)
  (:use clojure.tools.cli)
 
(:use [clojure.java.shell :only [sh]])
(:import java.io.File)
(:use [clojure.core.match :rename {compile match-compile}]))

(def infile "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/TfPeakSeq-X-lncRnaRegions.tab")
(def scratch-out  "/home/wespisea/scratch/encodeTFcomb/")
(def outfile-smash "/home/wespisea/work/research/researchProjects/encode/data/combinedTF-LncRegionFiles.tab")
(def outfile-cmds "/home/wespisea/work/research/researchProjects/encode/scripts/joinTF-autogen.sh")



;;(def h (parseCsvIntoArrayHash infile  "\t"))  
;;(printOutFileList h outfile-smash)
;; (printOutCmdList h outfile-cmds)
;;  (map wrapRunJob (getJoinCmds h))

(defn runAll[] (do (load "joinTfSites")(def h (parseCsvIntoArrayHash infile  "\t")) (printOutFileList h outfile-smash) ))


;;;; start ....
(defn wrapRunJob[job-cmd]
 (let [c-opt 2 m-opt 4  ]
 (str "/home/wespisea/bin/runJob.pl -c " c-opt " -m " m-opt " -i \"" job-cmd "\" & " )))


(defn filterByLncRegion[regionTerm col]
  (map (fn[x](vector (first x) (filter  (fn[y](= regionTerm (y "lncRegion")))   (second x)) )   )  
   col))

(defn getCellLabVec[h]
(reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x  (distinct (sort(map #(% "cell") h)))  y  (rest (distinct (sort(map #(% "lab") h))))] ["cell" x "lab" y])))

(defn getTfLabVec[h]
(reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x  (distinct (sort(map #(% "antibody") h)))  y  (rest (distinct (sort(map #(% "lab") h))))] ["antibody" x "lab" y])))

;;(let [op (apply str (interpose "\n" (map #(getEntryInfo % ",") (filterByLncRegion "body" lab-tf))))](spit outfile-smash  op))


(defn printOutCmdList[h file]
   (let [op (apply str (interpose "\n" (map wrapRunJob (getJoinCmds h))))]
  (spit file op)))



(defn getJoinCmds[h]
(let 
  [labs   (rest (distinct (sort(map #(% "lab") h))))
   cells  (distinct (sort(map #(% "cell") h)))
   TFs   (distinct (sort(map #(% "antibody") h)))
   col {}
   lab-cell (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x cells y labs] ["cell" x "lab" y]))
   lab-tf (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x TFs y labs] ["antibody" x "lab" y]))
   lncRegions ["body" "promoter_proximal" "promoter_distal"]]
   (flatten (for[x [[lab-cell "antibody"][lab-tf "cell"]] y lncRegions] (map (fn[entry](joinTogetherExpr entry (second x) y)) (filterByLncRegion y (first x)))))))

(defn getFileReport[h]
(let 
  [labs   (rest (distinct (sort(map #(% "lab") h))))
   cells  (distinct (sort(map #(% "cell") h)))
   TFs   (distinct (sort(map #(% "antibody") h)))
   col {}
   lab-cell (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x cells y labs] ["cell" x "lab" y]))
   lab-tf (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x TFs y labs] ["antibody" x "lab" y]))
   lncRegions ["body" "promoter_proximal" "promoter_distal"]]
   (flatten (for[x [[lab-cell "TF"][lab-tf "cell"]] y lncRegions] (map (fn[entry](joinTogetherExpr entry (second x) y)) (filterByLncRegion y (first x)))))))





(defn printOutFileList[h file]
 (let [op (apply str (interpose "\n" (map #(getEntryInfo % "," scratch-out) (reduce #(cons %2 %1) (filterByLncRegion "body" (getCellLabVec h))(filterByLncRegion "body" (getTfLabVec h))))))]
  (spit file op)))

(defn getEntryInfo[entry delim dir]
"return information about each grouped entry,print to file..."
(let [[label1 val1 label2 val2] (first entry)
		acc (apply str(interpose "::" (map (fn[x](x "dccAccession")) (second entry))))
		filename (apply str(interpose "::" (map (fn[x](x "filename")) (second entry))))
		cell (apply str(interpose "::" (distinct (map (fn[x](x "cell")) (second entry)))))
		antibody (apply str(interpose "::" (distinct (map (fn[x](x "antibody")) (second entry)))))
		file (apply str (interpose "/" [dir val2 val1]))
]
 (apply str (interpose delim (list (str label1 "=" val1 "-"label2 "=" val2) file filename cell antibody acc )))))



;;(map (fn[x](joinTogetherExpr x "cell" "body"))(filterByLncRegion "body" lab-tf))






(defn cmdfnHelperList[ filelist tmp curr]
  (if (= 0 (count filelist))
	curr
	(cmdfnHelperList (rest filelist) tmp (list curr "join " tmp " " (first filelist) " > " (str tmp ".1") "; mv " (str tmp ".1 ") " " tmp "; "))))


(defn cmdfnList[filelist headers target] ;;headers need "lncRNA in their name, for col=1
  (apply str (flatten (list " echo " (interpose " " (cons "transcript_id" headers)) "  > " target "; "
  (cond 
     (= 1 (count filelist))  ( list "cat " (first filelist) " >> " target "; ")
     (= 2 (count filelist))  ( list "join " (first filelist) " " (second filelist) " >> " target "; ")
     (< 2 (count filelist))  (let [ tmp (str "/home/wespisea/tmp" (rand-int 1000000000))
				    tmp (str "/home/wespisea/tmp" (rand-int 1000000000))]
				(list 
						(cmdfnHelperList 
									(-> filelist rest rest) 
									tmp
									(list "join " (first filelist) " " (second filelist) " > " tmp"; "))
						"cat " tmp " >> " target "; rm " tmp "; ")))))))


(defn joinTogetherExpr[entry headCol lncRegionFileName]
"create the join expression to combined TF lnc region crosses"
   (let [head (first entry)
         body (second entry)
         outDir (str  scratch-out (head 3) "/" (head 1) "/") 
	files (map (fn[x](.replace (x "crossFileLoc") ".bed" ".summary" )) (second entry))
	headers (map (fn[x](x headCol)) (second entry))
	target (str outDir lncRegionFileName)]
	(str "mkdir -p " outDir  "; " (cmdfnList files headers target))))



  (defn getJoinCmds[h]
    (let
        [labs   (rest (distinct (sort(map #(% "lab") h))))
            cells  (distinct (sort(map #(% "cell") h)))
            TFs   (distinct (sort(map #(% "antibody") h)))
            col {}
            lab-cell (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x cells     y labs] ["cell" x "lab" y]))
            lab-tf (reduce #(let[res (searchFAnotNil %2 h)] (if (not (nil? res)) (cons [%2 res] %1) %1 ))  []  (for [x TFs y l    abs] ["antibody" x "lab" y]))
            lncRegions ["body" "promoter_proximal" "promoter_distal"]]
        (flatten (for[x [[lab-cell "antibody"][lab-tf "cell"]] y lncRegions] (map (fn[entry](joinTogetherExpr entry (second x    ) y)) (filterByLncRegion y (first x)))))))
 
    (defn printOutCmdList[h file]
         (let [op (apply str (interpose "\n" (map wrapRunJob (getJoinCmds h))))]
             (spit file op)))
 
 
 (defn extractRepIdr[file cellType rnaExtract delim fileEnding]
   (let [tmp (str file (rand-int 10000))
         src (.replace file ".gz" "")
         trg (.replace file ".gz" fileEnding)
         ctype (repeat 3 (str cellType "-" rnaExtract))
         head (apply str(interpose delim (cons "transcript" (map #(str (.toUpperCase cellType) "-" rnaExtract "-" % )(list "C    OMB" "RPKM1" "RPKM2" "IDR")))))]
     (str "echo \"" head "\" > " trg "; cat " src " |awk '{print $12," delim "$6,", delim "$14," delim "$16," delim "$18" del    im "}'|sed 's/[\";]//g' >> " trg ";")))
 
 
 (defn prepGtfFiles [fileHash dir delim fileEnding]
   (map
     (fn [entry]
       (extractRepIdr (entry "localfile") (entry "cell") (entry "rnaExtract") delim fileEnding))
     fileHash))
 




