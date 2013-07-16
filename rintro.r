# comments
x <- c(10.3, 4.2, 5.4, 6.1, 7.50)
length(x);
min(x);
max(x);
mean(x);
fivenum(x); # min, lower hinge, median, upper hinge, max
sd(x);
var(x);
y = 2 *x + 1;
seq(start, end, by=.2, length = 10);
rep(x, times=5);
rep(x, each=5);
ls() # to list all the objects
rm() # to remove the objects specified
objects() # similar to ls()
paste(c("X","Y"), 1:10, sep="") # to paste character vectors


# index vectors; selecting and modifying subsets of data
y <- (x+1)[!is.na(x) & x>0] # which contains non missing values of x and positive x

x[6] # 6th element in the vecor x
c("x","y")[rep(c(1,2,2,1), times=4)]
y<-x[-(1:5)] # 1 to 5 values in the vector x are excluded

fruit <- c(5,10,1,20); # to index using character strings
names(fruit)<-c('oranges', 'banana', 'apple', 'peach')
lunch<-fruit[c('oranges', 'apple'); # prints 5 1

x[is.na(x)] <- 0

# matrices are multi dimen generalizations of vector
# factors provide a commpact way to handle categotical data
# lists are genral form of vector in which variois elements need not be same type. in general they proide a convenient way to return results of statistical data
# data frames are matrix like, columns can be different types. one row per observational unit
# functions

as.something() # used as a coersion of modes from one to another
attributes(object);

# class, for simple objeccts this is just the mode.. but this helps for object oriented programing in R. 

unclass(object) # o remove temporary effects of a class

#factors, ordered and unordred
state<-c("tas", "sa", "qld", "nsw", "nsw", "nt", "wa", "wa","qld", "vic", "nsw", "vic", "qld", "qld", "sa", "tas","sa", "nt", "wa", "vic", "qld", "nsw", "nsw", "wa","sa", "act", "nsw", "vic", "vic", "act");
statef<-factor(state);
levels(statef);

tapply();
incomes <- c(60, 49, 40, 61, 64, 60, 59, 54, 62, 69, 70, 42, 56,61, 61, 61, 58, 51, 48, 65, 49, 49, 41, 48, 52, 46,59, 46, 58, 43);
tapply(incomes, statef, mean); # gives mean income of each state, the factor level

# uses ordered() to have the levels in a specific order

# arrays or matrices
z<-1:15000
dim(z)<-c(3,5,100);
z[1,3,45] # indexing of a matrix

x<-array(1:20, dim=c(4,5)); # generate a 4x5 array

xb<-matrix(0,n,b);
xv<-matrix(0,n,v);
ib<-cbind(1:n,blocks);
iv<-cbind(1:n,varieties);

xb[ib]<-1
xv[iv]<-1
x,-cbind(xb,xv)

N <- crossprod(xb,xv);
N <- table(blocks, varieties);


z<-array(data_vector, dim_vector);

#outer product.. a and b are two numeric vectors.. aand outer prod is wose dimensions is obtained by concatenating their two dimenstions and data vector is got by forming all possssssible products %o%

ab <- a %o% b
ab <- outer(a,b,"*");
f<-function(x,y) cos(y)/(1+x*);
z<-outer(x,y,f)

t() # transpose of an array
nrow(); ncol() # they give no of rows and no of cols
diag() # diagnol matrix 

solve(A,b) # solves A*x = b

ev <- eigen(Sm);
svd()<- # gives singular value dicomposition


#least squares fit

ans <- lsfit(X,y);

#frequency tables from factors

statefr <- table(statef)
statefr <- tapply(statef, statef, length)
cut()
factor(cuut(incomes, breaks = 35+10*(0:7))) -> incomef
table(incomef, statef)



# lists and data frames
Lst <- list(name="Fred", wife="Mary", no.children = 3, child.ages=c(4,7,9))
Lst$name = "anil"
Lst[[1]] = "anil"
Lst[5] = list(matrix=Mat);
list.ABC = c(list.A, list.B, list.C) # concatenating lists

# data frames
accountants <- data.frame(home=statef, loot=incomes, shot=incomef)
attach()
search() # gives the search path for objects
detach()
head(dataframe) #first entries of the dataframe
tail(dataframe) #last entries of the dataframe
names(dataframe) #column names of the dataframe

cut()

# when you attach a data.frame the columns will be available as global variables and you can operate on them by directly mentioning them without using $ symbol.. and you can see them ls(2) or objects(2) and if you change any of the attached variables, it will get copied into a new variable with the changes and gets into level 1, so now you can see the changed variable in ls(1) and the unchanged variable is still in ls(2) and when you detach the it the data.frame remains unchanged.

# working with data frames, useful convention that allows you to work with many differnet problems comfortable in the same working directory
# 1. gather together all variables for any well defined and separate problem in a data frame under a suitably informative name
# 2. when working with a problem attach the appropriate data frame at position 2, and use the workingn directory at level 1 for operational quantities and temporary variables;
# 3. before leaving a problem, add any variables you wish to keep for future reference to the data frame using the $ form of assignment, and then detach()
# 4.finally remove all unwanted variables from the working directory and keep it as clean of left-over temporary variables as possible.


# reading data from files
# i normally use rstudio so importing data from files is done through gui
read.table("houses.data", header=TRUE);
scan()

data(); data("iris"); # use data to load available data sets into r env
View(object); 
edit(object); # for viewing and editing in inbuild data viewer and editor of r

xnew <- edit(data.frame()); # to enter new data via spreadsheet interface

# probability distributions

# functions for cumulative dist, probability density and the quantile function and to simulate the form of the distribution

beta    #beta
binom   #binomial
cauchy  #Cauchy
chisq   #chi-squared
exp     #exponential
f       #F
gamma   #gamma
geom    #geometric
hyper   #hypergeometric
lnorm   #log-normal
logis   #logistic
nbinom  #negative binomial
norm    #normal
pois    #Poisson
signrank#signed rank
unif    #uniform
weibull #weibull
wilcox  #Wilcoxon
t       #Students t

# prefix the names with 'd' for density and 'p' for the CDF, and 'q' for quantilee function and 'r' for simulation (random deviates)

# use package SuppDists for more distributions


attach(faithful);
summary(eruptions);
density(eruptions);
fivenum(eruptions);
stem(eruptions);
hist(eruptions)
hist(eruptions, seq(1.6, 5.2, 0.2), prob=TRUE);
lines(density(eruptions, bw=0.1))
rug(eruptions);
ecdf(eruptions);

long <- eruptions[eruptions>3]
plot(ecdf(long), do.points = FALSE, verticals=TRUE)
x <- seq(3, 5.4, 0.01)
lines(x, pnorm(x, mean = mean(long), sd = sd(var)), lty = 3);
par(pty = "s") # arrange for a square figure region
qqnorm(long); qqline(long);
shapiro.test(long);
boxplot(A, B);
t.test(A, B);
wilcox.test(A, B);




# Grouping, loops and conditional execution
ifelse(conditional, a, b) # returns vector of length longest among a, b.. a[i] if conditional[i] is true else b[i]

if (expr_1) expr_2 else expr_3

for (name in expr_1) expr_2

repeat expr
while (conditional) expr
break
next


# writing functions
name <- function(arg_1, arg_2, ...) expression

twosam <- function(y1, y2) {
    n1 <- length(y1); n2 <- length(y2)
    yb1 <- mean(y1); yb2 <- mean(y2)
    s1 <- var(y1); s2 <- var(y2)
    s <- ((n1-1)*s1 + (n-1)*s2)/(n1+n2-2)
    tst <- (yb1 -yb2)/sqrt(s*(1/n1 + 1/n2))
    tst

}

tstat <- twosam(data$male, data$female0; tstat


#customizing env, .First .Last


# statistical models in R
# suppose y, x, x0, x1, x2, .. are numerical variables, X is a matrix and A, B, C are factors.

y ~ x
y ~ 1+x # simple linear regression model of y on x.  both represent the same

y ~ 0 + x
y ~ -1 + x
y ~ x - 1 # linrear regresion without intercept

log(y) ~ x1 + x2 # multiple regression of transformation variable, log(y), on x1 and x2 wth an implicit intercept term

y ~ poly(x, 2)
y ~ 1 + x + I(x^2)
y ~ X + poly(x, 2)
y ~ A
y ~ A + x
y ~ A*B
 # and etc

fitted.model <- lm(formula, data = data.frame);
# generic functions for extracting model information

add1
deviance
formula
predict
step
alias
drop1
kappa
print
summary
anova
effects
labels
proj
vcov
coef
family
plot
residuals


new.model <- update(old.model, new.formula)

fm03 <- lm( y ~ x1 + x2 + x3, data = production)
fm4 <- update(fm03, . ~ . + x4)
smf4 <- update(fm4, sqrt(.) ~ .)

# learn generalized linear models
link functions etc





# plotting
par() # a function used to modify the list of graphics parameters for the current graphics device
plot(x,y)
plot(x)

pairs(X) # X is a data frame
coplot(a~b|c)
coplot(a~b|c+d)
qqnorm(x)
qqline(x)
qqplot(x,y)

hist(x)
hist(x, nclass=n)
hist(x, breaks=b,...)
dotchart(x,...)
image(x,y,z,...)
contour(x,y,z,...)
persp(x,y,z,...)

# plot arguments
add = TRUE
axes = FALSE
log = "x"
log = "xy"
type = "p" or "l" or "b" # points or lines or points connected by lines
type = "o" or "h" "s" "S" "n"

xlab = string
ylab = string
main = string
sub = string
abline(a,b)
abline(h=y)
abline(v=x)
abline(lm.obj)
polygon(x,y,...)
legend(x,y,legend,...)


help(Hershey)
locator()
text(locator(1), "Outlier", adj=0)
identify(x,y,labels)

# learn to use ggplot2 package for plotting



# packages
library()
library(randomForest)
search() # to see what packages are installed
help.start()


# getting a mode of a categorical feature
mymode <- function(x) as.numeric(table(x))[which.max(table(x))]

# numeric to factor and then back to numeric
# we cant just do as.numeric to the factor to get the actual value stored
# in the vector, instead we get the index of the level that is in the vector.

x = round(rnorm(100,sd=5))
f = as.factor(x)
y = as.numeric(f) # is not equal to x
z = as.numeric(levels(f)[f]) # which is equal to x


# must learn packages are ggplot2 and plyr





# common functions
with
within
which.max
which
round
subset
boxplot
transform
summarize
length
ddply
par # set or query graphical parameters
rep
seq
sample
unique
gl # to generate factors
interaction # to get a factor which represent the interaction of given factors
str # compactly display the structure of an arbitrary R object
head
tail
ctabs
table
cbind
rbind
c
image # display a color image
print
