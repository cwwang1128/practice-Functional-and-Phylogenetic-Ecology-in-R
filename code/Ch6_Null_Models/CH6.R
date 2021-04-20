##Chapter 6 Null Models
setwd("C://Users/bill/Desktop/ch6_null_models")
#6.2 Background============================================================
#6.2.1 Why Use Null Models for Phylogenetic and Functional Analyses?
library(picante)
pd.sample <- readsample("null.pd.example.sample.txt")
pd.phylo <- read.tree("null.pd.example.phylo.txt")
faith.output <- pd(pd.sample, pd.phylo)
plot(faith.output[,2], faith.output[,1], xlab = "Species Richness",
     ylab = "Faith's PD", pch = 16)

mpd.output <- mpd(pd.sample, cophenetic(pd.phylo), abundance.weighted = F)
plot(faith.output[,2], mpd.output, xlab = "Species Richness", ylab = "MPD",
     pch = 16)

#6.4 Randomizing Community Data Matrices in R===============================
#6.4.1 Unconstrained Randomizations=========================================
##This typically results in an inflation of type I error and 
##therefore they are not recommended for most analyses.

##Say we have a meta-community containing species six species (A, B, C, D,
##E, and F) and in the three communities in our study system we know their 
##presence or absences.
com.data <- matrix(c(1,0,1,0,1,0,1,1,1,1,1,0,1,1,1,0,0,1), nrow = 3)
rownames(com.data) = c("com1", "com2", "com3")
colnames(com.data) = c("A", "B", "C", "D", "E", "F")
##sp richness for each community
rowSums(com.data)

##If we wanted to randomly assemble community one (i.e., com1), we would then
##randomly draw four species from the list of A through F without replacement
sample(colnames(com.data), size = 4, replace = F)

##This could be repeated many times to produce many null or random communities
##for com1. Here we will only replicate the process ten times.
t(replicate(10, sample(colnames(com.data), size = 4, replace = F)))

##To perform this randomization for all communities to output a community data
##matrix we can use the randomizeMatrix() function in the picante package.
library(picante)
reps <- replicate(5, randomizeMatrix(com.data, null.model = "richness"))
apply(reps, MARGIN = 3, rowSums)
#the row sums do not change (i.e., the community species richness is fixed)
apply(reps, MARGIN = 3, colSums)
#the column sums are not consistent for species indicating that the occupancy 
#rate, the number of communities a species is in, is not maintained.


#6.4.2 Constrained Randomizations================================================
##"Independent Swap" null model
reps.is <- replicate(5, randomizeMatrix(com.data, null.model = "independentswap"))

##We can now see that the row and column sums (i.e., the community species
##richness and species occupancy rate) are now consistent across randomizations.
apply(reps.is, MARGIN = 3, rowSums)
apply(reps.is, MARGIN = 3, colSums)


##apply in abundance
com.data.a <-  matrix(c(3,0,2,0,10,0,4,2,3,4,5,0,6,1,2,0,0,7), nrow = 3)
rownames(com.data.a) = c("com1", "com2", "com3")
colnames(com.data.a) = c("A", "B", "C", "D", "E", "F")

##Now run the independent swab null model and sum the rows and columns.
reps.a.is <- replicate(5, randomizeMatrix(com.data.a, null.model = "independentswap"))
apply(reps.a.is, MARGIN = 3, rowSums)
##We see that the row sums are in some instances different indicating that the total
##abundance of a community varies from randomization to randomization.
apply(reps.a.is, MARGIN = 3, colSums)
##The column sums on the other hand are consistent indicating the total abundance of
##a species in the system is fi xed in the randomizations.

##Though, unfortunately, when we consider all of this evidence we see that adding
##abundance to this problem that the total abundance in a community is no longer fixed.



#6.5 Randomizing Phylogenetic Data=================================================
#6.5.1 Unconstrained Randomizations================================================
my.phylo <- read.tree("null.example.phylo.txt")
plot(my.phylo)
##The goal of the unconstrained phylogenetic randomization is to simply randomize
##the taxa names on the tips or terminal branches on the phylogeny.
my.phylo.rand <- my.phylo
my.phylo.rand$tip.label <- sample(my.phylo$tip.label, length(my.phylo$tip.label)
                                  ,replace = F)
plot(my.phylo.rand)


##Now use the function tipShuffle() to randomize the names of taxa on the phylogeny.
my.phylo.rand.2 <- tipShuffle(my.phylo)
plot(my.phylo.rand.2)


#The MPD metric
##When coding any null model in R, there are a couple of options: for-loop
##and replicate. And use replicate is faster than for-loop
my.sample <- read.table("null.example.sample.txt", sep = "\t", row.names = 1, header = T)
rand.mpd.fun <- function(x){
  tmp.phylo <- tipShuffle(x)
  
  ##we use mpd() here, but this function itself
  ##is slow due to its use of for() loops through
  ##rows in the community data matrix. See Chapter
  ##3 for an alternative and faster approach.
  mpd(my.sample, cophenetic(tmp.phylo))
}

##Next we can replicate this function ten times to produce ten random MPD values
##(columns of the output) for each community (rows of the output).
null.output <- replicate(10, rand.mpd.fun(my.phylo))
null.output
##When coding null models, it is generally an idea to fi rst utilize a small number of
##iterations and check the output of the null model. If there is an error in the code or
##it is slow, you do not want to spend time waiting for thousands of iterations to run
##to find out. When looking at the null model output matrix of any null that you have
##coded, one should immediately check if the results are indeed random.

hist(null.output[1,])
abline(v = mpd(my.sample, cophenetic(my.phylo))[1], col = "red", lwd = 2)
##It is possible that your observed value lands far outside the null distribution and you
##should replot the histogram with an adjusted x -axis that could incorporate your
##observed value.

##standardize effect size (S.E.S .)
ses.1 <- (mpd(my.sample, cophenetic(my.phylo))[1] - mean(null.output[1,])
          ) / sd(null.output[1,])
ses.1
ses.all <- (mpd(my.sample, cophenetic(my.phylo)) - apply(null.output, MARGIN = 1, mean)
            ) / apply(null.output, MARGIN = 1, sd)
##Positive S.E.S. values indicate an observed diversity value, MPD in this case, than
##average null value and negative values indicate an observed diversity value that is
##lower than average.
##That said, comparing S.E.S. values across communities or studies may be hindered
##by differences in species richness causing null distributions to take different shapes.
##For this reason and others, it is preferable to also compute where the observed value 
##ranks in the null distribution.

p.val.all <- apply(cbind(mpd(my.sample, cophenetic(my.phylo)), null.output), MARGIN = 1,
                   rank)[1,]/11
##This value tells us where the observed value ranks in the overall distribution (i.e.,
##the quantile score). This value can be utilized to calculate a P -value by dividing it
##by the number of null model iterations plus one.
##The value reported is the P -value for where the observed lands in the null probability
##distribution. In this example we have only 10 iterations of the null model and one
##observed value making 11 total values and therefore a denominator of 11.


##The PD metric
observed.pd <- pd(my.sample, my.phylo, include.root = F)

rand.pd.fun <- function(x){
  tmp.phylo <- tipShuffle(x)
  pd(my.sample, tmp.phylo)[,1]
}
null.output <- replicate(10, rand.pd.fun(my.phylo))

ses.all <- (observed.pd - apply(null.output, MARGIN = 1, mean)
            ) / apply(null.output, MARGIN = 1, sd)
p.val.all <- apply(cbind(observed.pd, null.output), MARGIN = 1, rank)[1,]/11


#6.5.2 Constrained Randomizations===============================================
library("Rsundials")
hardy.K <- 3
##Rsundials can't get...

#6.6.1 Unconstrained Randomizations=============================================
traits <- read.table("null.example.traits.txt",sep = "\t", header = T, row.names = 1)
rand.traits <- traits
rand.traits
replicate(5, sample(rownames(rand.traits), length(rownames(rand.traits)),replace = F))
trait.shuffle.funk <- function(x){
  x <- x[colnames(my.sample), ]
  rownames(x) <- sample(rownames(x), length(rownames(x)),replace = F)
  mpd(my.sample, as.matrix(dist(x[colnames(my.sample), ])))
}
replicate(10, trait.shuffle.funk(traits))

#6.6.2 Constrained Randomizations=============================================
max.sam.1 <- max(traits[names(my.sample[my.sample > 0]), 1])

##Calculate the maximum for trait 1 in community 1

max.nam.1 <-colnames(my.sample[1,which(my.sample[1,]>0)])
max(traits[max.nam.1, 1])

function(x){
  max.nam.1 <-colnames(x[1,which(my.sample[1,]>0)])
  max(traits[max.nam.1, 1])
  
}


str(my.sample)
shuff.constrain.nn <- function(x){
  ##Calculate the maximum and minimum value for
  ##trait 1 in community
  max.sam.1 <- NULL
  for(i in 1:3){
    max.sam.1 <-cbind(max.sam.1,
                      max(traits[colnames(my.sample[1,which(my.sample[i,]>0)]), 1]))
  }
  min.sam.1 <- NULL
  for(i in 1:3){
    min.sam.1 <-cbind(min.sam.1,
                      min(traits[colnames(my.sample[1,which(my.sample[i,]>0)]), 1]))
  }
    
  ##Calculate the maximum and minimum value for
  ##trait 2 in community
  max.sam.2 <- NULL
  for(i in 1:3){
    max.sam.2 <-cbind(max.sam.2,
                      max(traits[colnames(my.sample[1,which(my.sample[i,]>0)]), 2]))
  }
  min.sam.2 <- NULL
  for(i in 1:3){
    min.sam.2 <-cbind(min.sam.2,
                      min(traits[colnames(my.sample[2,which(my.sample[1,]>0)]), 2]))
  }
  ##Calculate the maximum and minimum value for
  ##trait 3 in community
  max.sam.3 <- max(traits[names(x[x > 0]), 3])
  min.sam.3 <- min(traits[names(x[x > 0]), 3])
  ##Get the names of all the species in your trait
  ##data that ate less or equal to the maximum
  ##trait 1 value observed in community and names
  ##of species greater than or equal to the
  ##observed minimum trait 1 value in community.
  tr.1.a <- names(traits[traits[ ,1] <= max.sam.1, 1])
  tr.1.b <- names(traits[traits[ ,1] >= min.sam.1, 1])
  ##Repeat the above for trait 2
  tr.2.a <- names(traits[traits[ ,2] <= max.sam.2, 2])
  tr.2.b <- names(traits[traits[ ,2] >= min.sam.2, 2])
  ##Repeat the above for trait 3
  tr.3.a <- names(traits[traits[ ,3] <= max.sam.3, 3])
  tr.3.b <- names(traits[traits[ ,3] >= min.sam.3, 3])
  ##Intersect all of the names from above such
  ##that the only remaining names are species that
  ##fall within the three dimensional range for
  ##the community. this is the new species pool.
  pruned.names <- Reduce(intersect, list(tr.1.a, tr.1.b, tr.2.a, tr.2.b,
                                         tr.3.a, tr.3.b))
  ##Prune smaple to only include species in new 
  ##species pool.
  pruned.sam <- my.sample[, pruned.names]
  ##prune trait matrix to only include species in new
  ##species pool
  pruned.matrix <- traits[pruned.names, ]
  ##Shuffle the species names on the pruned trait
  ##matrix
  rownames(pruned.matrix)<- sample(rownames(pruned.matrix),
                                   length(rownames(pruned.matrix)), replace = F)
  ##Calculate mean nearest neighbor distance.
  mntd(pruned.sam, as.matrix(dist(pruned.matrix)), abundance.weighted = F)
}
apply(replicate(999, apply(my.sample, 1, shuff.constrain.nn)),3,function(x){diag(x)})
##Error! should check it!