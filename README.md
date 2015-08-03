# encode-manager
A large suite of exploratory and predictive analysis to better understand the features of long non-coding RNA.

## Objective
Find the features of lncRNA that are associated with biological function and predict the set of lncRNA most likely to directly impact human physiology. Answering this question is important, as it will provide us with insights into what makes parts of the genome biologically important, as well as characterizing the background level of noise responsible spurious transcripts to form 


## PageRank Analysis
[pageRank white paper](https://www.rose-hulman.edu/~bryan/googleFinalVersionFixed.pdf)
After collection a list of over 100 lncRNA with demonstrated function, we wanted to compare their RNA expression, or abundance in different cell types, to lncRNA without known function. One hypothesis was that functional lncRNA are tightly clustered in expression space with the remainder generally dispersed. To measure this closeness between lncRNA, we decided to use the page rank algorithm: using euclidian distance between points in expression space instead of incoming page links.Ideally, if functional lncRNA are tightly grouped, they will have a much higher page rank value.     
[Here's an example comparing pagerank of closely grouped, versus dispersed points in a distribution](./plots/eigenVal_intuition.png)
From the analysis, we realized that non-functional lncRNA are generally clustered around the origin of expression space, consistent with their low expression, and functional lncRNA are much more dispersed.    
[pagerank importance values for functional lncRNA](./plots/rnaSeq-eigenRank/functionalTypes/bothPullDowns-all_biotypes/bothPullDowns-all_biotypesrank-bars.pdf)
To implement the pagerank algorithm I decided to extend the my R code base with C++ after realizing that calculating eigenvectors for dataset was too slow in base R.  Simply, the Rcpp package allowed me to add a C++ function as an inline string(see `./analysis/rnaSeq/eigenRank_08272013.R`, `srcCalcEigenCppD`), wrap the function, then use it as a regular R function.     

## Implementing logistic regression with custom decision boundaries
Using N-dimensional data, logistic regression can be used to apply a circular decision boundary. Increasing the degree polynomial of the objective function, the decision boundary can can therefore be more complicated. To implement logistic regression, the first and second derivatives to the objective function are used to optimize a weight matrix, which is then applied to the test set.  Although many packages exist for this purpose, implementing the logistic regression algorithm gives greater control, and is a great learning exercise.     
These two plots should be the same, which means that the formula used in the implemented logistic regression, is equivalent to applying a polynomial of degree 2 to the data, then passing the data into a regular logistic regression.     
[ellipse implemented into the algorithm](./plots/fullAnalysisExperiment/test/logReg/mlclass/ellipse_ex2data2_lambda=10_degree=2.pdf)     
[polynomial of degree=2 applied to data, then passed to log reg](./plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data2_degree=2.pdf)     
Code:     
[cost function](./analysis/rnaSeq/logisticRegPcaExprData.R#L139-149)    
[gradient function](./analysis/rnaSeq/logisticRegPcaExprData.R#L202-217)
[using optimization to find a log. reg. parameter](./analysis/rnaSeq/logisticRegPcaExprData.R#L482-509)    
(note: the optimization doesn't use the 2nd differential, as the optimization works better this way)     


## Extending the learning problem to account for unlabeled data
####Problem:
In most training sets, we some set of features, labeled as 'positive' or 'negative'. These two labels are then used to train an objective function. However, there may not be examples that are classified as either positive or negative.  This is exactly the case with predicting functional lncRNA. We have published reports of positive examples, but no one publishes papers with negative results. (This should change, ask me about the problems with academic publishing...).  
#####Soltion:
[Elkan and Noto](./docs/elkan_posonly.pdf)
In a brilliant paper, Elkan determined a way to adjust the probability outcome of any algorithm to account for unlabeled examples. Elkan's `c` is multiplied by the probability of a example being labeled, to get the probability of a example being positive.    
[code to estimate c.](/analysis/rnaSeq/fullAnalysis-Jan2014.R#L1532-1536)     
[plots when applied to ridge regression](./plots/fullAnalysisExperiment/test/logReg/ridgeRegression/PULearning/)     



## Using entropy based approaches to quantify tissue specificity
[Cabilli, 2011](http://genesdev.cshlp.org/content/25/18/1915.full.pdf)     
[code, to calculate tissue specificity](./analysis/rnaSeq/calcTissueSpecificity.R)    




