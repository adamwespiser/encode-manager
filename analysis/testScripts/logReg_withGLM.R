# for reproducibility (data generation only) 
set.seed(123)
# sample size
Nsims <- 1e3
# generate data
df <- data.frame(x1=rnorm(Nsims, 0, 2), 
                 x2=runif(Nsims, -2, 2), 
                 x3=rnorm(Nsims))
# generate response
df[, "y"] <- rbinom(Nsims, 1, with(df, exp(x1+3*x2+2*x3+x2*x3+0.5*x3^2)/(1+exp(x1+3*x2+2*x3+x2*x3+0.5*x3^2))))

# glm (your version)
glm(y ~ x1*x1 + x1*x2 + x1*x3 + x2*x2 + x2*x3 + x3*x3, 
    data =df, family=binomial) 
# equivalent but simpliried version
glm(y ~  x1*x2 + x1*x3 + x2*x3, 
    data = df, family=binomial) 
# another equivalent version 
# I think this is what you mean by your input vector
vars <- names(df)[names(df) != "y"] 
# write down the formula
form <- paste("y ~", do.call(paste, c(as.list(do.call(paste, c(expand.grid(vars, vars), sep=":"))), sep=" + ")))
# use this formula
glm(form, data = df, family=binomial) 

# different version that includes x1^2, x2^2 and x3^2 
# you had these terms in the original version, but they weren't used 
# because you didn't put them in the I(...) notation
glm(y ~ x1*x2 + x1*x3 + x2*x3 + I(x1*x1) + I(x2*x2)  + I(x3*x3), 
    data = df, family=binomial) 
# and here's the same thing using the poly function (as @Ben Bolker suggested)
# this used the option raw=TRUE for the output to be comparable with the above, 
# however usually you should use ortogonal polynomials instead. 
glm(y ~ poly(x1, x2, x3, degree=2, raw=TRUE), data=df, family=binomial)

# yet another version without the non-squared terms (as @ Aaron suggested)
glm(y ~ x1:x2 + x1:x3 + x2:x3 + I(x1*x1) + I(x2*x2)  + I(x3*x3), 
    data = df, family=binomial) 
# you can also use define a formula similar to the one I suggested in the first version. 