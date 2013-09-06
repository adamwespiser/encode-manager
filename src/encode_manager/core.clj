 (ns encode-manager.core
  (:gen-class)
  (:use clojure.tools.cli)
 
  (:use [clojure.java.shell :only [sh]])
  (:import java.io.File)
  (:use [clojure.core.match :rename {compile match-compile}]));;end of namespace dec

;; EXPR Data encode rnaseq
;; (def ld "/home/wespisea/scratch/encodeRnaSeq/caltech")
;; (def h (parseFileHash (str ld "/files.txt")))
;; (outputCsvRecord h (str ld "/files.tab") "\t")
;; (def h1 (parseCsvIntoArrayHash (str ld "/files.tab") "\t"))
;; def rnas  (filter #(contains? (set our-cell-types) (% "cell")) (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "RnaSeq" "view" "TranscriptGencV3c" "localization" "cell"] h1))))
;;(def rna-hash (reduce #(assoc %1 (%2 "filename") (merge %2 {"localfile" (entryToLocalFileInDirExpr %2 ld )}  {"url" (entryToUrl %2)} )) {} rnas))
;; (outputCsvRecord rna-hash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr_caltech.tab" "\t")
;;(def h1 (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr_caltech.tab" "\t"))



;; Expr data rna-seq hg19 dashboard(cshl)
;;create encode peak seq list
;;  (def h (parseCsvIntoArrayHash encode-integration-tab "\t"))
;; EXPR Data
;; (def h (parseFileHash "/home/wespisea/work/research/researchProjects/encode/data/hg19_RNA_dashboard_files.txt"))
;; (outputCsvRecord h "/home/wespisea/scratch/encodeRnaSeq/cshl/files.txt" "\t")
;; (def h1 (parseCsvIntoArrayHash "/home/wespisea/scratch/encodeRnaSeq/cshl/files.txt" "\t"))
;; def rnas  (filter #(contains? (set our-cell-types) (% "cell")) (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "RnaSeq" "view" "TranscriptGencV7" "localization" "cell"] h1))))
;;(def rna-hash (reduce #(assoc %1 (%2 "filename") (merge %2 {"localfile" (entryToLocalFileInDirExpr %2 "/home/wespisea/scratch/encodeRnaSeq/cshl")}  {"url" (entryToUrl %2)} )) {} rnas))
;; (outputCsvRecord rna-hash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr.tab" "\t")
;;(def h1 (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr.tab" "\t"))
;; $cat /home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr.tab|awk '{print "wget",$1,"-O",$11}'>tmp.sh; chmod u+x tmp.sh; ./tmp.sh
;; 
;; (cmdfnList (map #(.replace (% "localfile") ".gz" "tmp" ) h1) (map (fn[x](str (.toUpperCase (x "cell")) "-" (x "rnaExtract")))h1) "/home/wespisea/scratch/encodeRnaSeq/cshl/allCellsCombined.space")

;;create encode peak seq list
;;  (def h (parseCsvIntoArrayHash encode-integration-tab "\t"))
;; (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "ChipSeq" "softwareVersion" "PeakSeq"] h));;(def TFs  (filter #(contains? (set our-cell-types) (% "cell")) (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "ChipSeq" "softwareVersion" "PeakSeq"] h))));; (def TF-hash (reduce #(assoc %1 (%2 "filename") (merge %2 {"localfile" (entryToLocalFileInDir %2 "/home/wespisea/scratch/encodePeakSeq")}  {"url" (entryToUrl %2)} )) {} TFs))
;; ( outputCsvRecord TF-hash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/TfPeakSeq.tab" "\t")

;;download peakseq and confirm;;  (def h1 (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/TfPeakSeq.tab" "\t"))
;; (downloadEntries h1)


;; RNA Expression WITH RPKM 1 and 2 & IDR
;; (def h0 (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/RnaSeqExpr.tab" "\t"))
;; (def h01 (filter #(not (= "NHEK" (% "cell"))) h0))
;; (def h1 (filter #(not (= "HMEC" (.toUpperCase (% "cell")))) h01))
;; (def seqDir "/home/wespisea/scratch/encodeRnaSeq/cshl/")
;; (def tmpEnd "tmp5")
;; (spit "/home/wespisea/sandbox/tmp.sh" (apply str (prepGtfFiles h1 seqDir " " tmpEnd)))
;;  $ chmod u+x /home/wespisea/sandbox/tmp.sh; /home/wespisea/sandbox/tmp.sh
;; for Rep5, change NHEK cell name to NHEK1;; (spit "/home/wespisea/sandbox/tmpCompile.sh" (cmdfnList (map #(.replace (% "localfile") ".gz" tmpEnd ) h1) (flatten (map (fn[exprName](list (str exprName "-RPKM1") (str  exprName "-RPMK2") (str exprName "-IDR")))    (map (fn[x](str (.toUpperCase (x "cell")) "-" (x "rnaExtract")))h1))) "/home/wespisea/scratch/encodeRnaSeq/cshl/allCellsCombined_2reps_idr.space"))



;;create encode peak seq list
;;  (def h (parseCsvIntoArrayHash encode-integration-tab "\t"))
;; (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "ChipSeq" "softwareVersion" "PeakSeq"] h))
;;(def TFs  (filter #(contains? (set our-cell-types) (% "cell")) (filter #(not (nil? %)) (searchFileArrayHash ["dataType" "ChipSeq" "softwareVersion" "PeakSeq"] h))))
;; (def TF-hash (reduce #(assoc %1 (%2 "filename") (merge %2 {"localfile" (entryToLocalFileInDir %2 "/home/wespisea/scratch/encodePeakSeq")}  {"url" (entryToUrl %2)} )) {} TFs))
;; ( outputCsvRecord TF-hash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/TfPeakSeq.tab" "\t")

;;download peakseq and confirm
;;  (def h1 (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/TfPeakSeq.tab" "\t"))
;; (downloadEntries h1)


;; Tf peak seq combos
;; (def h (parseCsvIntoArrayHash "/home/wespisea/work/research/researchProjects/encode/data/TfPeakSeq-X-lncRnaRegions.tab" "\t"))  


(def info " This file consists of functions desgined to parse encode \"file.txt\" descriptions of experiments for the encode II dataset and return them in hashmap.\n\n(def ds(parseCsvIntoArrayHash all-encode.tab \"\\t\"))")
(def cell-file "/home/adam/Downloads/cellTypeInformationEncode.tab")
(def url-string "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhRnaSeq/files.txt")
(def files-string "/home/adam/Downloads/encodeFiles.txt")
(def all-encode "/home/adam/programming/clojure/encode-manager/data/encodeAll_native.clj")
(def all-encode-tab "/home/adam/programming/clojure/encode-manager/data/encodeFilesAll.tab")
(def encode-integration "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/encode-integration-data-files.txt")
(def encode-integration-tab "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/encode-integration-data-files.tab")

(def our-cell-types '("GM12878" "HeLa-S3" "HepG2" "K562" "H1-hESC"))

(defmacro home-dir[] "/home/wespisea/scratch")


(def makeStripCmd[file]
  (let [in (.replace file ".gz" "")
        out (.replace file ".gz" "tmp")]
  (str "cat " in "| awk \'{print $12, \"\t\",\$6}\'|sed \'s/[\";]//g\' > " out ";")))



(defn entryToLocalFileInDirExpr [entry dir]
  "determine the local file name of an encode experiment"
  (let [exprdir (entry "exprdir")
        filename (entry "filename")]
   (str dir "/" (last (.split filename "/" )))))


(defn entryToLocalFileInDir [entry dir]
  "determine the local file name of an encode experiment"
  (let [exprdir (entry "exprdir")
        filename (entry "filename")]
   (str dir "/" filename)))


(defn entryToLocalFile [entry]
  "determine the local file name of an encode experiment"
  (let [exprdir (entry "exprdir")
        filename (entry "filename")]
   (str (home-dir) "/" exprdir "::" filename)))



(defn md5sumEntryLocal[entry]
  (first (.split ((clojure.java.shell/sh "md5sum" (entryToLocalFile entry)) :out) " ")))

(defn entryFileExists[entry]
  (let [exprdir (entry "exprdir")
        filename (entry "filename")]
   (.exists (File. (entryToLocalFile entry)))))


(defn entryToUrl [entry]
  "convert an entry to a url"
  (if (entry "url")
    (entry "url")
   
  (let [urldir (entry "urldir")
        exprdir  (entry "exprdir")
        filename (entry "filename")]
  (str urldir exprdir  "/" filename ))))

(defmacro wget[url target]
  `(clojure.java.shell/sh "wget" ~url "-O" ~target "--quiet"))

(defn confirmMd5Entry[entry]
  (if (entryFileExists entry)
      (= (entry "md5sum") (md5sumEntryLocal entry)) 
       false))

(defn searchFileArrayHash [inVec arrayRec]
"accespts the file-array-hash and a vector of keys and values"
(let [y (partition 2 inVec)
      f (map first y)
      s  (vec (map second y))]
      ;;s  (vec (map #(if (list? %) (cons :or %) % ) (map second y)))]
  (map (fn [hash-line]
         (let [x (vec (map hash-line f))]
             (match [x] 
                   [s] hash-line))) 
  arrayRec  )))


(defn searchFAnotNil [inVec arrayRec]
  (let [result (filter (fn[x](not (nil? x))) (searchFileArrayHash inVec arrayRec))]
    (if (= 0 (count result))
        nil 
        result)))


(defn downloadEntries[entries]
 (letfn [(download[entry]
         (if-not (entryFileExists entry)
            (wget (entryToUrl entry) (entryToLocalFile entry))))]
  (pmap #(future (download %)) entries)))

 
(defn getCellTypeArray[file]
  (map #(.split % "\t") (.split (slurp file) "\n")))

(defn selectTypes[types rec]
  "applied to a filtered, or unfiltered rec structure"
  (map #(cons (first %) (map (second %) types)) rec))

(defn getAllKeys[rec]
  "returns a list of all the keys in the file"
  (distinct (flatten (map #(keys (second %)) rec))))


(defn getWebsite[url]
" download a website, or file using the slurp function"
  (slurp url))

(defn getExpr[line]
  "for a file.txt as a string, get the experiment name"
  (first (.split line "\t")))

(defn getHashEntry[line]
  "for a line of file.txt, load the name=value pairs into a hash"
  (let [arr (.split (second (.split line "\t")) "; ")]
    (reduce #(assoc %1 (first (.split %2 "=")) (second (.split %2 "=")) ) {} arr )))

(defn processLine[line]
  "determine the experiment, and name=value hash for a line, returns a hash ds"
  (list (getExpr line) (getHashEntry line)))

(defn processLine-integration [line]
  (let [f (getExpr line)
        s-hash (getHashEntry line)
        tokens (.split f "/")
        file-name (last tokens)
        url-dir (nth (reverse (apply list tokens)) 1) ]
   (hash-map file-name (merge s-hash {"url" f "exprdir" url-dir}))))

(defn parseEncodeIntegrationFiles[sn]
  "returns a hash-map"
(reduce merge (map processLine-integration sn)))


(defn parseFilesList[url]
"returns a list of ( {exp1->{name1->val1, name2->val2...}, {exp2->...}}) "
  (let [hm {}]
    (map #(processLine %) (.split (getWebsite url) "\n"))))

(defn parseFileHash[url]
" returns a hash of {expr -> {name -> val, name2 -> val2 ... }, ... }"
  (let [hmm {}]
    (reduce #(assoc %1 (getExpr %2)(getHashEntry %2)) {} (.split (getWebsite url) "\n" ))  ))

(defn parseFileHashWithDir[url expr-dir url-entry]
"returns a from parseFileHash, execpt value for name=url directory,experiment dir are added"
  (reduce #(assoc %1 (getExpr %2)(merge (getHashEntry %2) (into {} (map vector ["exprdir" "urldir"] [expr-dir url-entry]))))
          {}
          (.split (getWebsite url) "\n")))


(defn getHashEntryCsv[line k]
  "reads the file.txt saved as .csv from file, needs a file line and list of keys, IN ORDER"
  (into {}
    (filter
        #(not (= "" (second %)))
	(map vector k line) )))

(defn inputFileNoComment[file comChar]
  " returns a .split'd array with commented out files removed"
   (filter #(not (= (.substring % 0 1) comChar))
	   (.split (slurp file) "\n")))




(defn parseCsvIntoHash[file delim]
  "load the lines of .csv into a hash, works on all-encode-tab"
  (let [arr (map #(.split % delim) (inputFileNoComment file "#"))
	k (rest (first arr))]
    (reduce #(assoc %1 (first %2) (getHashEntryCsv (rest %2) k)) {} (rest arr))))

(defn parseCsvIntoArrayHash[file delim]
  "load the lines of .csv into a hash works on all-encode-tab"
  (let [arr (map #(.split % delim) (inputFileNoComment file "#"))
	k (first arr)]
    (reduce #(cons (getHashEntryCsv  %2 k) %1) '() (rest arr))))


(defn loadDsFromFile[file]
  (read-string (slurp file)))

(defn getDs [] (loadDsFromFile all-encode))

(defn help-message
    "banner help for encode description file parser"
    [opts args banner]
    (let [{:keys [help input field delim]} opts]
        (println (info))
	(println (str "input: " input))
	(println (str "field: " field))
	(println (str "delim: " delim))
	(println banner)
	(System/exit 0)))

(defn getCellByAttr[rec name value]
  "returns the cell-lines that match for a given term/type"
  (let [n (.toUpperCase name)
  v (.toUpperCase value)]
    (map #( first %)
       (filter #(= (get (second %) name) value)
       rec))))



(defn getCellByAttrFromFile[name value]
  "returns the cell-lines that match for a given term/type"
  (let [n (.toUpperCase name)
        v (.toUpperCase value)]
    (map #( first %) 
       (filter #(= (get (second %) name) value) 
       (parseCsvIntoHash cell-file "\t")))))

(defn collectCellVals[name1]
"returns all possible valus for cell-line name from name=val hash"
  (let [search-fn (if (list? name) #(first (val %)) #(val %))
        name-query (if-not (list? name1) (list name1) name1)]
  (reverse (sort-by #(val %)
    (reduce #(if (contains? %1 %2) 
                 (assoc %1 %2 (+ 1 (%1 %2)))
                 (assoc %1 %2 1)) 
                 {} (map #(map % name-query) (vals (parseCsvIntoHash cell-file "\t")))
)))))
(defn collectCellVals1[rec name1]
"returns all possible valus for cell-line name from name=val hash"
  (let [search-fn (if (list? name) #(first (val %)) #(val %))
        name-query (if-not (list? name1) (list name1) name1)]
  (reverse (sort-by #(val %)
    (reduce #(if (contains? %1 %2)
                 (assoc %1 %2 (+ 1 (%1 %2)))
                 (assoc %1 %2 1))
                 {} (map #(map % name-query) (vals rec))
)))))
;;getAllKeys
(defn reportValuesForEncodeHash [rec]
 "reports all the values associated with an encode hash entry"
  (let [keys (getAllKeys rec)]
   (map #(take 10 (collectCellVals1 rec %)) keys)))



(defn cellTypeHeader[file]
  (set (first (getCellTypeArray cell-file))))


(defn generateUrlForExprDat[entryHash]
 "generates the file url for a given experiment"
  (let[{:strs [urldir exprdir]} (val entryHash)] (str urldir exprdir "/" (key entryHash))))



(defn convertValsToCsv[outer-k delim]
    "creates a partial function, currying the supplied keys, into a function that accepts a hash and prints out only the values corresponding to the curried keys in a concat'd .csv string"
    (fn [inner-x]
       (apply str
        (interpose delim (cons (first inner-x)  (map (second inner-x) outer-k))
	            ))))

(defn convertValsToCsvT[outer-k]
(fn [inner-x](count (cons (first inner-x) (map (second inner-x) outer-k)))))
    
(defn seqToILstring[seq delim]
  "return a string of interleaved sequence values with delim"
  (apply str
	 (interpose delim seq )))
    

(defn hashToCsvString[rec delim]
  "takes a loaded in hash and prints out a csv addr"
  (let [k (getAllKeys rec)]
   (seqToILstring
	    (cons
              (seqToILstring (cons "filename" k) delim)
	      (map (convertValsToCsv k delim) rec))
	   "\n")))
   
(defn outputCsvRecord[rec file delim]
"prints a file record as a .csv with missing vals if value not found"
  (spit file (hashToCsvString rec delim)))


(defn nameValueIntoString[h]
  "takes a key value hash and appends a string with \"key=\"value entries"
  (let [klist (keys h)
        vlist (vals h)]
    (apply str (interleave
           (map str (map #(apply str % "=") klist) vlist)
	  (repeat "; ")))))


(defn getCellInfo[celltype]
    "search the cell-file for info"
      (map  #(println (first %)(nameValueIntoString (second %)))
	         (filter #(= (first %) celltype )
			        (parseCsvIntoHash cell-file "\t"))))


(defn processHashEntryIntoString[h]
  "reconstruct the line of original file"
  (let [expr (first h)
	nameVal (second h)]
    (apply str expr "\t" (nameValueIntoString nameVal) "\n")))

(defn convertExprHashStr[rec]
"reconstruct the original file format, no data loss"
  (apply str (map processHashEntryIntoString rec)))

(defn outputAsOriginal[rec file]
  "print out the hash record in the original format"
  (spit file (convertExprHashStr rec)))



(defn getAllEncodeFiles[files-str]
"with the list of files.txt, download those and export a hash of all encode files"
(let [cinfo (parseCsvIntoHash files-string ",")]
(reduce #(merge %1 (parseFileHashWithDir ((val %2) "url") (key %2) ((val %2) "urldir"))) {} cinfo)))


(def banner "the purpose of this script is to handle the encode data descriptor (file.txt, that have there information in \"key=\"value format.\nTwo options are available, --convert and --select.")
(defn -main [& args]
  (let [[options args banner] (cli args
				 ["-h" "--help" "Show help" :default false :flag true]
				 ["-i" "--input" "The input file" :default "notSpecified"]
				 ["-d" "--delim" "The delimiter for own files" :default ","]
                                 ["-c" "--convert" "Convert key=val into csv form" :default false]
				 ["-s" "--select" "Which cols to display example \"cell,localization,lab\",this is case sensitive" :default "all"])
		  pi 2.35
		 {:keys [delim field input select convert]} options]
        (if (:help options)

        (help-message options args banner)
       
	(println (hashToCsvString (parseFileHash  input))))

	)

	)
