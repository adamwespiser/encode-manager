# encode-manager
A large suite of exploratory and predictive analysis to better understand the features of long non-coding RNA.

## Objective
Find the features of lncRNA that are associated with biological function and predict the set of lncRNA most likely to directly impact human physiology. Answering this question is important, as it will provide us with insights into what makes parts of the genome biologically important, as well as characterizing the background level of noise responsible spurious transcripts to form 


## PageRank Analysis
[pageRank white paper](https://www.rose-hulman.edu/~bryan/googleFinalVersionFixed.pdf)
After collection a list of over 100 lncRNA with demonstrated function, we wanted to compare their RNA expression, or abundance in different cell types, to lncRNA without known function. One hypothesis was that functional lncRNA are tightly clustered in expression space with the remainder generally dispersed. To measure this closeness between lncRNA, we decided to use the page rank algorithm: using euclidian distance between points in expression space instead of incoming page links.Ideally, if functional lncRNA are tightly grouped, they will have a much higher page rank value.     
[Here's an example comparing pagerank of closely grouped, versus dispesed points in a distribution](./plots/eigenVal_intuition.png)
From the analysis, we realized that non-functional lncRNA are generally clusted around the origin of expression space, cosistent with their low expression, and functional lncRNA are much more dispersed.    
[pagerank importance values for functional lncRNA](./plots/rnaSeq-eigenRank/functionalTypes/bothPullDowns-all_biotypes/bothPullDowns-all_biotypesrank-bars.pdf)
To implement the pagerank algorithm I decided to extend the my R code base with C++ after realizing that calculating eigenvectors for dataset was too slow in base R.  Simply, the Rcpp package allowed me to add a C++ function as an inline string(see `./analysis/rnaSeq/eigenRank_08272013.R`, `srcCalcEigenCppD`), wrap the function, then use it as a regular R function.     

## Implementing logistic regression with custom decision boundaries
Using N-dimensional data, logistic regression can be used to apply a circular decision boundary. Increasing the degree polynomial of the objective function, the decision boundary can can therefore be more complicated. To implement logistic regression, the first and second derivatives to the objective function are used to optimize a weight matrix, which is then applied to the test set.  Although many packages exist for this purpose, implementing the logistic regression algorithm gives greater control, and is a great learning exercise.






## Extending the learning problem to account for unlabled data

## Using entropy based approaches to quantify tissue specificity
[Cabilli, 2011](http://genesdev.cshlp.org/content/25/18/1915.full.pdf)



