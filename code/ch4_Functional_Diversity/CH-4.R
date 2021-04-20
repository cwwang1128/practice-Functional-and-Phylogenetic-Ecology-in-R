setwd("C://Users/user/Desktop")

########4.3 Quantifying the Functional Composition of Communities 
########Using the Moments of Trait Distributions
#(1)mean, (2)standard deviation, (3)skew, and (4)kurtosis are begin knowing the composition.
my.sample <- read.table("FD.example.sample.txt",header = T, row.names = 1)
my.sample
#sp. names as columan names, site names as row names.
traits <- read.table("FD.traits.txt",sep = "\t", header = T, row.names = 1)
traits
#sp. names as row names, trait names as columan names.

#contains the names of species present in first community (my.sample's first row)
spp <- names(my.sample[1, my.sample[1,] > 0])
#extracted from traits matrix by asking spp.
traits[spp, ]

#Also, can do it.
traits[names(my.sample[1, my.sample[1, ] > 0]), ]

#first, calculate the community mean for each trait 
mean(traits[names(my.sample[1, my.sample[1, ] > 0]), ], na.rm = T)
#but, it can't use "mean", can do the two way: "colMean" or "apply"
#way 1
colMeans(traits[names(my.sample[1, my.sample[1, ] > 0]), ], na.rm = T)
#"na.rm = T" use trait data matrix had missing trait values for some species. 
#way 2
apply(traits[names(my.sample[1, my.sample[1, ] > 0]), ] , MARGIN = 2 , mean)

#second, calculate the community standard deviation for each trait 
##but, it can't use "sd", can do "apply"
apply(traits[names(my.sample[1, my.sample[1, ] > 0]), ] , MARGIN = 2 , sd)

##High SD are indicative of more functional diversity,
##but they may be biased due to differences in the mean from community to community.


#To reduce this bias, can use "coefficient of variation": SD/mean
apply(traits[names(my.sample[1, my.sample[1, ] > 0]), ] , 2 , sd)/
  apply(traits[names(my.sample[1, my.sample[1, ] > 0]), ] , 2 , mean)


#The skew of the trait distribution for community 1
library("fBasics")
skewness(traits[names(my.sample[1, my.sample[1, ] > 0]), ], method = "moment",na.rm = T)
#High values of skewness may do not imply lower community functional diversity,
#but do indicate that most co-occurring species tend to have very similar trait values.

#The kurtosis of the traits in community 1
kurtosis(traits[names(my.sample[1, my.sample[1, ] > 0]), ], method = "moment",na.rm = T)
#Small kurtosis values indicate community trait distributions with ???fatter??? tails, 
#may indicate an increase in the average trait disparity between co-occurring species

#seting function for calculate community mean
mean.funk <- function(x){
  ##calculate mean for the column 1 in the trait matrix for sp. that are present in communtiy
  mean(traits[names(x[x> 0]), 1], na.rm = T)
}
#seting function for calculate community standard diviation
sd.funk <- function(x){
  sd(traits[names(x[x> 0]), 1], na.rm = T)
}
#seting function for calculate community skewness
skewness.funk <- function(x){
  skewness(traits[names(x[x> 0]), 1], method = "moment", na.rm = T)
}
#seting function for calculate community kurtosis
kurtosis.funk <- function(x){
  kurtosis(traits[names(x[x> 0]), 1], method = "moment", na.rm = T)
}

#apply those function
apply(my.sample, MARGIN = 1, mean.funk)
apply(my.sample, MARGIN = 1, sd.funk)
apply(my.sample, MARGIN = 1, skewness.funk)
apply(my.sample, MARGIN = 1, kurtosis.funk)



#"community-weighted mean (CWM)"
#mean trait value weighted by the relative abundance of each species.
#for more fitting the system

#relative abundances = abundance / the total abundance
my.ra.sample <- my.sample/ rowSums(my.sample)
my.ra.sample

#also, pakeage "vegan" is a way to calculate relative abundances
library(vegan)
#"decostand" is a variety code by its method.
my.ra.sample2 <- decostand(my.sample, method = "total", MARGIN = 1)
#method "total" divides each value in a row by the row sum or row total
##transforming the matrix to presence/absence (method = ???pa???),
##standardizing the values to a mean of zero and unit variance (method = ???standardize???),
##standardizing the values to range between zero and one (method = ???range???),
##dividing the values by the maximum value in the row or column (method = ???max???). 
my.ra.sample2
#transforming to relative abundance
rowSums(my.ra.sample2)
#checking data transforming to relative abundance. the code input should be 1


#For CWM, just weighted functional trait by the relative abundance of each species.
#they need have same length.
weighted.mean(traits[colnames(my.ra.sample), 1], my.sample[1, ])
weighted.mean(traits[colnames(my.ra.sample), 2], my.sample[1, ])
weighted.mean(traits[colnames(my.ra.sample), 2], my.sample[2, ])

weight.mean.funk <- function(x){
  weighted.mean(traits[names(x[x>0]), 1], x[x>0])
}
apply(my.sample, MARGIN = 1, weight.mean.funk)

#weighted standard deviation
library(SDMTools)
weight.sd.funk <- function(x){
  wt.sd(traits[names(x[x>0]), 1], x[x>0])
}
apply(my.sample, MARGIN = 1, weight.sd.funk)



########4.4 Dendrogram-based VS Euclidean distance-based meansures of FD
########4.4.1 Generating Trait Distance Matrices

#produce a distance matrix for the second trait, by using "dist" and "euclidean" method.
trait.2 <- as.matrix(traits[,2])
rownames(trait.2) = rownames(traits)
trait.2
my.dist.mat.2 <- dist(trait.2, method = "euclidean")
dist(traits, method = "euclidean")

#first, transform data to approximate a normal distribution for each trait. e.g., using a log()
traits.scaled <- apply(log(traits), MARGIN = 2, scale)
#all the trait values are scaled to approximately a mean of zero and unit variance.
#center: each values - mean
#scale: (each values - mean)/sd
#to reduce the differen size values impact the tests


#principal components (PC) analysis produce one PC axis for each input column (i.e., trait).
pc <- princomp(traits.scaled)

#To select the few and majority axes that explain the variance in the data,
#examine the proportion of the total variance explained by each axis.
summary(pc)
#A good rule: the PC axes that explain over 90 % or even 95 % of the variation.
#in the example, it's first three

#examine the trait loadings, to get most heavily weighted traits
print(pc$loadings, cutoff = 0.001)
#trait 1 and trait 2 most heavily influence the first PC axes
#traits 4 and 5 most heavily influence the second and third PC axes, respectively.

#each species lands on the first three PC axes and generate a distance matrix
print(pc$scores, cutoff = 0.001)
#genetate a distance matrix
pc.scores <- pc$scores[,1:3]
rownames(pc.scores) <- rownames(traits)
pc.scores

#calculate the multivariate Euclidean distance between all species 
pc.dist.mat <- dist(pc.scores, method = "euclidean")
pc.dist.mat
#the way less to be biased by the co-variation and differences measure scale. 
#it's most recommended to calculate distance matrix, generally.
#don't calculating a distance matrix using all trait simultaneously
#given that traits strongly covary

#since trait sometime is miss or mixed variable (i.e., continuous and categorical)
#"dist()" can't be used, alternative, "gowdis()" can be used.
#"gowdis()" is based on Gower distance.
library(FD)
gowdis(traits)

########4.4.2  Generating Trait Dendrograms
#1. generating a distance trait matrix
#2. generating a trait dendrograms, data by step 1

#1. generating a distance trait matrix
#if traits have covary, see last section, 
#trait values is miss or mixed variables, use "gowdis()"
my.dist <- dist(traits, method = "euclidean")
#2. generating a trait dendrograms, by step 1 matrix distance
#"hclust" based on "Unweighted Pair Group Method with Arithmetic Mean = UPGMA" 
my.dendro <- hclust(my.dist, method = "average")
plot(my.dendro)

#make a UPGMA dendrogram for the trait found in the second column of our trait matrix
trait.2 <- as.matrix(traits[,2])
rownames(trait.2) <- rownames(traits)
my.dendro.2 <- hclust(dist(trait.2, method = "euclidean"), method = "average")
plot(my.dendro.2)
#it is often important and interesting to perform all analyses on all traits at once 
#and each trait individually

########4.4.3  Pairwise and Nearest Neighbor Measures
#convert the dendrogram to a distance matrix 
dendro.dist.mat <- cophenetic(my.dendro)
#it will make less information.
#if not necessary, don't do it.

#so,here are using raw data output a square matrix
square.dist.mat <- as.matrix(pc.dist.mat)
#unweighted mean pairwise distance for row2 (community2).
mean(as.dist(square.dist.mat[names(my.sample[2, my.sample[2,] > 0]) ,names(my.sample[2, my.sample[2,] > 0])]))

#set function for unweighted mean pairwise distance.
trait.pw <- function(x){
  mean(as.dist(square.dist.mat[names(x[x > 0 ]),names(x[x > 0])]))
}
apply(my.sample, MARGIN = 1, trait.pw)

#also, "mpd()" can be used it. but, slower than apply (the function is "loop")
library("picante")
mpd(my.sample, square.dist.mat, abundance.weighted = F)


#The abundance-weighted version of the pairwise trait distance metric 
trait.pw.prime <- function(x){
  ##get the names of the sps in the coummunity
  com.names <- names(x[x > 0])
  ##make a matrix with one row containing abundances and names of all present sps
  com <- t(as.matrix(x[x > 0]))
  #Make functional distance matrix for taxa in community
  com.dist <- square.dist.mat[com.names,com.names]
  ##calculate the product of the abundances of all sps in the community
  abundance.products <-t(as.matrix(com[1,com[1, ] > 0, drop = F]))%*%as.matrix(com[1,com[1, ] > 0, drop = F])
  ##calculate a mean of the community functional distance matrix weighted by the products of all
  ##pairwise abundances
  weighted.mean(com.dist, abundance.products)
}
apply(my.sample, MARGIN = 1, trait.pw.prime)

#also, "mpd()" can be used it. but, slower than apply (the function is "loop")
mpd(my.sample, square.dist.mat, abundance.weighted = T)


#Extract or replace the diagonal of a matrix, or construct a diagonal matrix.
#to ensure we are not counting conspecifics in our nearest neighbor calculations.
new.square.dist <- square.dist.mat
diag(new.square.dist) <- NA

#calculate the mean nearest functional neighbor distance for the species in our third community
com.dist.mat <- new.square.dist[names(my.sample[3,my.sample[3,] > 0]), names(my.sample[3, my.sample[3,] >0])]
#This is the mean nearest neighbor distance for our community.
#i.e.,calculate the minimum value in each row
mean(apply(com.dist.mat, MARGIN = 1, min, na.rm =T), na.rm = T)
sd(apply(com.dist.mat, MARGIN = 1, min, na.rm =T), na.rm = T)
#A low variance indicates that species are relatively evenly placed in functional space.


#First the mean nearest functional neighbor distance
trait.mnn <- function(x){
  ##get the names of the sps present in the coummunity
  com.names <- names(x[x > 0])
  
  ##make the coummunity functional distance matrix by extracting those row and cloumns that have
  ##sps present in coummunity
  my.com.dist <- square.dist.mat[com.names, com.names]
  
  ##Set all diagonal values to NA so that the zeros for conspecific comparisons do not interfere
  ##with our calculation of nearest neighbors
  diag(my.com.dist) <- NA
  ##use apply() to calculate the minimun value in each functional distance matrix and take a mean
  ##of those values
  mean(apply(my.com.dist, MARGIN = 1, min), na.rm = T)
}

#the function for the standard deviation of nearest neighbor distances.
trait.sdnn <- function(x){
  ##get the names of the sps present in the coummunity
  com.names <- names(x[x > 0])
  
  ##make the coummunity functional distance matrix by extracting those row and cloumns that have
  ##sps present in coummunity
  my.com.dist <- square.dist.mat[com.names, com.names]
  
  ##Set all diagonal values to NA so that the zeros for conspecific comparisons do not interfere
  ##with our calculation of nearest neighbors
  diag(my.com.dist) <- NA
  ##use apply() to calculate the minimun value in each functional distance matrix and take a mean
  ##of those values
  sd(apply(my.com.dist, MARGIN = 1, min), na.rm = T)
}

apply(my.sample, MARGIN = 1, trait.mnn)
apply(my.sample, MARGIN = 1, trait.sdnn)
##above codes not work.


trait.mnn.prime <- function(x){
  ##get the names of the sps present in the coummunity
  com.names <- names(x[x > 0])
  
  ##make the coummunity functional distance matrix by extracting those row and cloumns that have
  ##sps present in coummunity
  my.com.dist <- square.dist.mat [com.names, com.names]
  
  ##Set all diagonal values to NA so that the zeros for conspecific comparisons do not interfere
  ##with our calculation of nearest neighbors
  diag(my.com.dist) <- NA
  ##use apply() to calculate the minimun value in each functional distance matrix and take a mean
  ##of those values
  weighted.mean(apply(my.com.dist, MARGIN = 1, min),x[x > 0], na.rm = T)
}

#the function for the standard deviation of nearest neighbor distances.
trait.sdnn.prime <- function(x){
  ##get the names of the sps present in the coummunity
  com.names <- names(x[x > 0])
  
  ##make the coummunity functional distance matrix by extracting those row and cloumns that have
  ##sps present in coummunity
  my.com.dist <- square.dist.mat[com.names, com.names]
  
  ##Set all diagonal values to NA so that the zeros for conspecific comparisons do not interfere
  ##with our calculation of nearest neighbors
  diag(my.com.dist) <- NA
  ##use apply() to calculate the minimun value in each functional distance matrix and take a mean
  ##of those values
  wt.sd(apply(my.com.dist, MARGIN = 1, min), x[x > 0], na.rm = T)
}

apply(my.sample, MARGIN = 1, trait.mnn.prime)
apply(my.sample, MARGIN = 1, trait.sdnn.prime)
##above codes not work.

#standard deviation of the nearest neighbor distance, 
#the unweighted and abundance-weighted mean
mntd(my.sample, square.dist.mat, abundance.weighted = F)
mntd(my.sample, square.dist.mat, abundance.weighted = F)
#low values of the mean nearest functional neighbor distance indicate functionally
#similar species co-occurring or low FD, whereas high values indicate functionally
#dissimilar species are co-occurring and therefore high FD


########4.4.4 Ranges and Convex Hulls
#The convex hull volume for a community is now commonly referred to as FRic or Functional Richness
#Strat from univariate equivalent for community 1
com.1.names <- names(my.sample[1, my.sample[1,] > 0])
#the range for each trait in the first community.
apply(traits[com.1.names, ], MARGIN = 2, max)- apply(traits[com.1.names, ], MARGIN = 2, min)

#write a function fot calculate reange of each trait(i.e., functional richness) simultaneously
range.function <- function(x){
  ##get the names of the sps present in our community
  com.names <- names(x[x > 0])
  ##Calculate the range for each trait
  apply(traits[com.names, ], MARGIN = 2, max)
  -apply(traits[com.names, ], MARGIN = 2, min)
}
apply(my.sample, MARGIN = 1, range.function)

library(geometry)
convhulln(traits[names(my.sample[1, my.sample[1,] > 0]), ])
#can't work deu to sample not enough
hull.output <- convhulln(traits[names(my.sample[2, my.sample[2, ] > 0]), ], options = "FA")
#but colum2 can do it

names(hull.output)
hull.output$vol

hull.function <- function(x){
  ##get the names of the sps present in community.
  com.names <- names(x[x>0])
  ##calculate the hull using all traits.
  convhulln(traits[com.names, ], options ="FA")$vol
}

apply(my.sample, MARGIN = 1, hull.function)

library(FD)
#use dbFD to calculate FRic
dbFD(traits, my.sample)$FRic

########4.4.5 other measures
#Functional Evenness metric utilizes a minimum spanning tree (MST) to connect all species 
#in functional trait space and measures the regularity of the species points along the branches of this
#tree and the regularity of their abundances
##Thus, we may expect it to be very similar to a nearest neighbor metric.
dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FEve

#Functional Divergence metric first measures the average distance of all species from the centroid
#of the trait space in a community and then sums the magnitude of the divergences from that mean.
##Thus higher values are supposed to indicate more dispersion
##towards the maximum and minimum of the range of traits.
dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FDiv

#The FDis metric calculates the distance of each species from the centroid of the community traits.
##it is not quite the same calculation as the pairwise distance
##between species we can expect that it might be highly correlated.
dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FDis


########4.5 comparing metrics of functional diversity

richness <- rowSums(decostand(my.sample, method = "pa", MARGIN = 1))
pw <- mpd(my.sample, as.matrix(dist(traits)), abundance.weighted = F)
pw.prime <- mpd(my.sample, as.matrix(dist(traits)), abundance.weighted = T)
nn <- mntd(my.sample, as.matrix(dist(traits)), abundance.weighted = F)
nn.prime <- mntd(my.sample, as.matrix(dist(traits)), abundance.weighted = T)
fric <- dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FRic
feve <- dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FEve
fdiv <- dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FDiv
fdis <- dbFD(traits[colnames(my.sample), ], my.sample, w.abun = T, stand.x = T)$FDis
outputs <- as.data.frame(cbind(richness, pw, pw.prime, nn, nn.prime, fric, feve, fdiv, fdis))
names(outputs) <- c("richness", "pw", "pw.prime", "nn", "nn.prime", "FRic", "FEve", "FDiv", "Fdis")
plot(outputs, pch = 16)
cor(outputs)

#We see that the PW and FDis metrics are almost identical as has been previously
#pointed out and that FDis and FDiv, while being
#correlated, are not identical indices. Lastly, we see that the hull/FRic metric and the
#mean nearest neighbor metric are related to species richness, indicating that they
#should perhaps be considered in the context of a null model for comparative analyses