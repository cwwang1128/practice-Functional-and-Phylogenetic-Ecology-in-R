setwd("C://Users/user/Desktop")
#5.3 Tree-Based Measures of Phylogenetic Beta Diversity

#5.3.1 UniFrac
##UniFrac: (pd.total - pd.shared) / pd.total
my.sample <- read.table("beta.example.sample.txt",sep="\t",row.names = 1,header = T)
library("picante")
my.phylo <- read.tree("beta.example.phylo.txt")
#unweighted UniFrac
unifrac(my.sample, my.phylo)
##it is work!!! 

##including the relative abundances of species into the calculation
##1 step is to calculate a raw weighted UniFrac value
##Given the calculation utilizes the relative abundances of species in each community we first
##transform our community data matrix to transform raw counts of species to relative abundances.
my.ra.sample <- my.sample/rowSums(my.sample)

#2 step is to loop through each branch to calculate the product of its length and the absolute
#difference of the relative abundances of the species in the two communities subtended by that branch.
node.bls <- matrix(NA, nrow = length(my.phylo$edge.length), ncol = 6)
node.bls[, 1:2] <- my.phylo$edge

plot.phylo(my.phylo, show.tip.label = F)
nodelabels(bg = "white")
tiplabels(bg = "white")

node.bls[,3] <- my.phylo$edge.length
head(node.bls)

## In book, he use for-loop to introduce the analysis
##there is unsing the package to do.
library(GUniFrac)
GUniFrac(my.sample, my.phylo, alpha = c(0, 0.5, 1))
##The GUniFrac() function in the GUniFrac package can run the analyses with multiple
##different levels of alpha at the same time to produce an array of dissimilarity matrices.


#5.3.2 Phylogenetic Sorenson's Index
##The PhyloSor metric is a metric of similarity such that as the phylogenetic similarity of
##two communities increases the PhyloSor value increases.
phylosor(my.sample, my.phylo)


##5.4 Distance-Based Measures of Phylogenetic and Functional Beta Diversity

#5.4.1 Pairwise Measures
##In a distance-based framework this overall similarity can be calculated as the distance between
##all species or individuals in one community to all species or individuals in a second community
#first step: make a distance matrix for our phylogeny or trait dataset.
p.dist.mat <- cophenetic(my.phylo)
#we will calculate a Euclidean distance matrix from the raw trait values.
traits <- read.table("beta.example.traits.txt", sep= "\t", row.names = 1, header = T)
f.dist.mat = as.matrix(dist(traits, method = "euclidean"))
#second step:  to extract the names of the species present in our community
com.spp.1 <-names(my.sample[1, my.sample[1, ] > 0])
com.spp.2 <-names(my.sample[2, my.sample[2, ] > 0])
##Here we are just using the phylogenetic distance matrix, but the exact same procedure could
##be done using the trait distance matrix to measure functional beta diversity
p.dist.mat[com.spp.1, com.spp.2]
mean(p.dist.mat[com.spp.1, com.spp.2])

##If the use of conspecifics was undesirable for some reason they could easily be removed 
##from the analysis by replacing all zero values in the community phylogenetic
##or functional distance matrix with NA values.
new.p.dist.mat <- p.dist.mat
new.p.dist.mat[new.p.dist.mat == 0 ] = NA
mean(p.dist.mat[com.spp.1, com.spp.2], na.rm = T)

##using function from picante package
comdist(my.sample, cophenetic(my.phylo), abundance.weighted = F)

#loop is not efficiency, so we creat a function
get.presents <- function(x){
  names(x[x > 0])
}
list.of.names <- apply(my.sample, 1, get.presents)
Dpw.apply.function <- function(x){
  ##write a function within our main function.
  ##This function calculates the mean phylogenetic
  ##distance of all sps found in community z
  ##and community x. the communoty x comes from
  ##our main dunction and the community z comes
  ##from this function. This is, in essence, how
  ##we replace two nested loops.
  tmp.function<- function(z){
    mean(p.dist.mat[x, z])
  }
  ##Apply our tmp.function from above to all
  ##levels in our list.of.names list. In other
  ##words, for community x, we will calculate the
  ##Dpw between that community and all other
  ##communities including itself.
  lapply(list.of.names, FUN = tmp.function)
}
dpw.output <- lapply(list.of.names, Dpw.apply.function)
do.call(cbind, dpw.output)


#5.4.2 Nearest Neighbor Measures
##This is the question that nearest neighbor distance-based beta diversity analyses ask
##and they are conceptually similar to beta diversity indices that simply ask how many
##species or genera are shared between two communities.
com.spp.1 <-names(my.sample[1, my.sample[1, ] > 0])
com.spp.2 <-names(my.sample[2, my.sample[2, ] > 0])
p.dist.mat[com.spp.1, com.spp.2]
#Because the ultimate goal is to calculate the nearest neighbor values for each species
#in both communities we can use two apply() functions, one for the rows and
#one for the columns, and calculate a mean value to produce the Dnn result
mean(c(apply(p.dist.mat[com.spp.1, com.spp.2], MARGIN = 1, min, na.rm=T),
       apply(p.dist.mat[com.spp.1, com.spp.2], MARGIN = 2, min, na.rm=T)
       ))
##Again note that this will assign NA
##values to non-conspecifics that have no phylogenetic branch length separating them
##or that share identical trait values and as such you will want to take care before
##applying this approach.
p.dist.mat.noc <- p.dist.mat
p.dist.mat.noc[p.dist.mat.noc==0]=NA
mean(c(apply(p.dist.mat.noc[com.spp.1, com.spp.2], MARGIN = 1, min, na.rm=T),
       apply(p.dist.mat.noc[com.spp.1, com.spp.2], MARGIN = 2, min, na.rm=T)
))
##using function from picante package
comdistnt(my.sample, cophenetic(my.phylo), abundance.weighted = F, exclude.conspecifics = F)

## weighted by abundance
comdistnt(my.sample, cophenetic(my.phylo), abundance.weighted = T, exclude.conspecifics = F)

#loop is not efficiency, so we creat a function
get.presents <- function(x){
  names(x[x > 0])
}
list.of.names <- apply(my.sample, 1, get.presents)
Dnw.apply.function <- function(x){
  ##here is the only change from the code in the
  ##previous subsection for Dpw
  tmp.function<- function(z){
    mean(c(apply(p.dist.mat.noc[x, z], MARGIN = 1, min, na.rm=T),
           apply(p.dist.mat.noc[x, z], MARGIN = 2, min, na.rm=T)))
  }
  lapply(list.of.names, FUN = tmp.function)
}
do.call(cbind, lapply(list.of.names,Dnw.apply.function))


#5.5 Other Metrics
##Rao's Dbeta
##The equation for this metric is defined simply as the pairwise phylogenetic
##(or functional) distance between all individuals in two communities
rao.output <- raoD(my.sample, my.phylo)
rao.output$Dkl
rao.output$H

#5.6 Comparing Metrics
unifrac.output <- c(unifrac(my.sample, my.phylo))
p.sor.output <- c(phylosor(my.sample, my.phylo))
Dpw <- c(comdist(my.sample, cophenetic(my.phylo), abundance.weighted = F))
Dpw.prime <- c(comdist(my.sample, cophenetic(my.phylo), abundance.weighted = T))
Dnn <- c(comdistnt(my.sample, cophenetic(my.phylo), abundance.weighted = F, exclude.conspecifics = T))
Dnn.prime <- c(comdistnt(my.sample, cophenetic(my.phylo), abundance.weighted = T, exclude.conspecifics = T))
Dnn.prime.noc <- c(comdistnt(my.sample, cophenetic(my.phylo), abundance.weighted = T, exclude.conspecifics = F))
rao.D <- c(as.dist(raoD(my.sample, my.phylo)$Dkl))
rao.H <- c(as.dist(raoD(my.sample, my.phylo)$H))
outputs <- as.data.frame(cbind( c(unifrac.output), c(p.sor.output), c(Dpw), c(Dpw.prime), 
                                c(Dnn), c(Dnn.prime), c(Dnn.prime.noc), c(rao.D), c(rao.H)))
names(outputs) <- c("UniFrac", "Psor", "Dpw", "Dpw'","Dnn", "Dnn'", "Dnn'.noc", "RaoD", "RaoH")
plot(outputs, pch = 16)
cor(outputs)
