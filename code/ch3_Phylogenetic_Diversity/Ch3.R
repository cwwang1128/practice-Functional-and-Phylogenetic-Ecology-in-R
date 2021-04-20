setwd("C://Users/bill/Desktop/Functional and Phylogenetic Ecology in R/supplementary material/ch3_Phylogenetic Diversity")

##########3.3 ?€œCommunity?€? Datasets
#Binary indicating
##The example table is space with column headers and the row names in column one.
pa.matrix<-read.table("pa.matrix.txt",sep=" ",header = T, row.names = 1)
##Summing the rows = richness for each site.
rowSums(pa.matrix)
##Summing the columns = the number of sites a taxon occupies.
##This value could be divided by the number of sites (i.e., the number of rows)
##to calculate the occupancy rate for a taxon in the study system
colSums(pa.matrix)

#indicating the number of individuals
abund.matrix<-read.table("abund.matrix.txt",sep=" ",header = T, row.names = 1)
##row sums are the total abundance of all taxa per site
rowSums(abund.matrix)
##column sums are the total abundances of taxa in the study system.
colSums(abund.matrix)

#indicating the relative abundance of taxa
ra.matrix<-read.table("ra.matrix.txt",sep=" ",header = T, row.names = 1)
##row sums should equal one
rowSums(ra.matrix)
##column sums equal the sum of the relative abundances for a taxon across sites
colSums(ra.matrix)

library(picante)
#read in your three-column text file with no column headers and examine the output.
my.3.sample<-readsample("matrix.3col.txt")
#put ones in the second column to simply indicate presences of an individual. 
#The values are the numbers of individuals for a taxon per site.
decostand(my.3.sample, method = "pa", MARGIN = 1)
#Similarly, if we wanted to transform our original matrix into relative abundance we
#can use the same function and the ?€œtotal?€? method, which divides each value in a row
#by the total of all values in a row if the MARGIN is set to one.
my.ra.3.sample<-decostand(my.3.sample, method = "total", MARGIN = 1)
my.ra.3.sample

#save a copy of your community data matrix
write.table(my.ra.3.sample, "my.ra.3.matrix.txt", sep="/t", 
            row.names = T, col.names = T, quote = F)
#save your community data matrix in the three column format
writesample(my.ra.3.sample, "my.new.3col.data.txt")



#######3.4 Tree-Based Measures of Phylogenetic Diversity
#Faith?€™s Index is described as the sum of the branch lengths 
#connecting all species in an assemblage.
##reading our community data matrix
#first row contained the taxa names
#the first column contained the name of the assemblage.
my.sample<-read.table("PD.example.sample.txt",sep = "\t", row.names = 1, header = T)

#The example phylogenetic tree containing all of the species in our community data
my.phylo<- read.tree("PD.example.phylo.txt")
my.phylo

#Potential problems that could easily be detected are nonsensical taxa names 
#on the tips of the tree that do not match the names of taxa in your community data matrix
plot(my.phylo)
#also can use is.ultrametric()


#calculating Faith? Index for a single community, and not all communities at once,
#we must extract the data from the community for all species that are present
#(i.e., have an abundance or relative abundance greater than zero).
#the first row in our example community data matrix for only those columns (i.e., species)
#where the values (i.e., the abundance or relative abundance) are greater than zero.
my.sample[1, my.sample[1, ] > 0]

##Because only the species present in the community are in this output, we can simply
#calculate the species richness of the community by asking R for the length of this vector.
species.richness.com.1 <- length(my.sample[1, my.sample[1, ]>0])
species.richness.com.1


library(geiger)
#use the treedata() function to prune our original phylogeny.
#for comparative analyses, it can match and organize phylogenetic and trait data
#the function requires (1) a phylogenetic tree
#and (2)trait data's a matrix with the row names or a vector being the taxa names
#Although not using actual trait data, utilize the vector of abundances and names 
#of taxa extracted from community data matrix as the ?€œtrait?€? data 
#pruning original phylogeny to only contain the species present in first community.
treedata(my.phylo, data=t(my.sample[1, my.sample[1, ]>0]))
#modify the code to avoid the warning message that several taxa (i.e., tips)
#in the phylogeny were not found in the trait dataset
#and were therefore pruned from the phylogeny.
pruned.tree = treedata(my.phylo, data = t(my.sample[1, my.sample[1, ] > 0]),
                       warnings = F)$phy
#Faith?€™s Index can now be easily computed using the pruned phylogenetic tree containing
#only species found in the community by summing all branch (i.e., edge) lengths.
sum(pruned.tree$edge.length) 
#This code only provides Index for a single community.

##to calculate the Faith?€™s Index for all communities at once.
#(1) write a simple function named prune.sum.function() for pruning
#(2) called my.phylo using row x from the community data matrix.
#We then sum the branch lengths in the pruned phylogeny and output the result.
prune.sum.function<-function(x){
  tmp.tree<-treedata(my.phylo, x[x > 0], warnings = F)$phy
  sum(tmp.tree$edge.length)
}
apply(my.sample, MARGIN = 1,prune.sum.function)
#same above code
apply(my.sample, 1, function(x) {
  sum(treedata(my.phylo, x[x > 0], warnings = F)$phy$edge.length)
})

#package is often rely on loop(). But apply() is faster than loop(), especially large data.
#so, use apply() is good at fast deal with large data.
#use package's function is good at easy deal with small data.

###use picante package calculate Faith's index
#The original Faith?€™s Index did not include the root of the pruned phylogenetic tree
pd(my.sample, my.phylo, include.root = F)

#EH (?€œEvolutionary History?€? or ?€œEvolutionary Heritage?€?, modified Faith's index) include root.
#The rationale for including the root is that it provides more information 
#regarding the long evolutionary history leading up to the species found in the community.
pd(my.sample, my.phylo, include.root = T)

#but,PD does not includ relative abundance of individual sp. in the assemblages.
#it is fine for conservation assessments, but terribly for community structure and diversity.
# so, modified Faith?€™s index weighted by abundance that is called "Weighted Faith?€™s Index".

#To calculate the Weighted Faith?€™s Index for a single community,
#(1) by extracting from a larger phylogeny the phylogeny containing only the species
com.1.phylo<- treedata(my.phylo, data=t(my.sample[1, my.sample[1, ]>0]))$phy

#(2)generate an empty matrix filled with NA values
#that hold the output for the variables quantify for each individual branch.
#The number of rows corresponds to the parameter "n" in the above equation.
branches <- matrix(NA, nrow(com.1.phylo$edge), ncol = 4 )

#In the first two columns placed the beginning and ending (i.e., basal and terminal) node
#for each branch by asking for the edges stored in our community phylogeny object.
#We will use this information to ask R what species are subtended by each branch on the phylogeny.
branches[,1:2] <- com.1.phylo$edge
#In the third column of the matrix we will place the length of all edges (i.e., branches)
#This column therefore contains the parameter "l" in the above equation.
branches[,3] <- com.1.phylo$edge.length
head(branches)
#we have the nodes defining each edge in the phylogeny and the length of that edge.

#To visualize what the node numbers in the first two columns correspond to in our phylogeny
plot.phylo(com.1.phylo, show.tip.label=F)
nodelabels(bg = "white")
tiplabels(bg = "white")
#visualize what the information in the matrix we have created indicates.

#write a simple for() loop to go through each row in our matrix,
#ask for the species descendent from the node in column two, and average their abundance.
#Summing this information would provide the denominator in the above equation.

for(i in 1:nrow(branches)){
  ##using the terminal node number for each branch
  ##in colum two of the node.bls object we pull
  ##out the taxa names subtended by the branch
  ##described in each row of the node.bls object.
  leaves.node <- tips(com.1.phylo, branches[i,2])
  
  ##Average the abundance of the sp. found in
  ##community 1, which is the first row in our
  ##my.sample object. this is equal to Ai in the
  ##raw weighted.Faith equation.
  branches[i, 4] = mean(t(my.sample[1, leaves.node]))
  
  ##close loop
}

#the number of branches = the "n" parameter in the equation = the number of rows
number.of.branches <- nrow(com.1.phylo$edge)
#The denominator(??†æ??) is calculated by summing column four.
#each row in column four is an "A" value (average abundance).
denominator <- sum(branches[,4])
#The numerator = the sum of (individual branch lengths (li) branches[,3]
#* average abundance of the species (A) branches[,4])
numerator <- sum(branches[,3] * branches[,4])
#We now have the values to calculate the Weighted Faith?€™s Index for a single community.
weight.faith.output <- number.of.branches * (numerator / denominator)
weight.faith.output
#This is obviously suboptimal and was only used to demonstrate the general approach.

#A better method!!!!!
#(1)write a function to calculate the Weighted Faith?€™s Index 
#(2)apply this approach to the rows in the community data matrix.

#the function will take an input community data matrix and trim a larger phylogeny called 
#my.phylo to an individual community phylogeny for each row in the community data matrix. 
weighted.faith.function <- function(x){
  ##extract the names of sp. in a community with an abundance greater than zero and use
  ##that information to make a pruned phylogeny for that community
  tmp.tree = treedata(my.phylo, data = x[x > 0], warnings = F)$phy
  ##Create empty branches matrix
  branches = matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
  ##Fill first two columns of the matrix with node numbers defining each edge
  branches[,1:2] = tmp.tree$edge
  ##Fill the third column with the length of each branch
  branches[,3] = tmp.tree$edge.length
  
  ##loop through each set of leaves to calculate the mean abundance (Ai) for those sp.
  for(i in 1:nrow(branches)){
    leaves.node=tips(tmp.tree, branches[i,2])
    branches[i,4]=mean(t(x[leaves.node]))
  }
  ##Lastly calculated te Weighted Faith's Index
  n.of.branches=nrow(tmp.tree$edge)
  denominator=sum(branches[,4])
  numerator=sum(branches[,3]*branches[,4])
  n.of.branches*(numerator/denominator)
}

all.weighted.output <- apply(my.sample, MARGIN = 1, weighted.faith.function)
all.weighted.output
#We have now completed the calculation of the Weighted Faith?€™s Index.


##########3.5 Distance-Based Measures of Phylogenetic Diversity
#phylogenetic diversity metrics are most generally fall into one of two
#categories?€? pairwise or nearest neighbor.
##########3.5.1 Pairwise Measures
#Distance-based measures of phylogenetic diversity utilize a phylogenetic distance
#matrix or phylogenetic variance?€“covariance matrix to quantify a metric of relatedness
#between species or taxa in a community or sample.
#A phylogenetic distance matrix

plot.phylo(my.phylo)
#get phylogenetic distance matrix
dist.mat<-cophenetic(my.phylo)
dist.mat[1:4,1:4] 

#Get phylogenetic variance?€“covariance (VCV) matrix

vcv(my.phylo)
com.1 <- my.sample[1, my.sample[1, ]>0]
com.1
names(com.1)
dist.mat.com.1 <- dist.mat[names(com.1), names(com.1)]
mean(as.dist(dist.mat.com.1))
new.mpd.function <- function(x){
  ##Get the names of sp. present in the community
  com.names <- names(x [ x > 0 ])
  ##Calculate mpd by extracting the lower triangle
  ##of a phylogenetic distance matrix comprised of only
  ##the sp.in our communtiy
  mean(as.dist(dist.mat[com.names, com.names]))
}
apply(my.sample, MARGIN = 1, new.mpd.function)
mpd(my.sample, cophenetic(my.phylo), abundance.weighted = F)
new.mpd.a.function <- function(x){
  ##Get the names of sp. present in the community
  com.names <- names(x[x > 0])
  ##Make a matrix with one row containing aboundances 
  ##and names of all present sp.
  com <- t(as.matrix(x[x > 0]))
  
  ##Make phylogenetic distance matrix for taxa in community
  com.dist <- dist.mat[com.names, com.names]
  
  ##Calculate the product of the abundances of all sp. in the community
  abundance.products <- t(as.matrix(com[1, com[1, ] > 0, drop = F]))
  %*% as.matrix(com[1, com[1, ] > 0, drop = F])
  
  ##calculate a mean of the community phylogenetic distance matrix weighted
  ##by the products of all pairwise abundances.
  weighted.mean(com.dist, abundance.products)
}
apply(my.sample, MARGIN = 1, new.mpd.a.function)
mntd(my.sample, cophenetic(my.phylo), abundance.weighted = T)
psv(my.sample, my.phylo)
pse(my.sample, my.phylo)
psr(my.sample, my.phylo)


############3.5.2
new.mntd.function <- function(x){
  ##Get the names of the sp. present in a community
  com.names <- names(x[x > 0])
  ##Make the communtiy phylogenetic distacne matrix by extracting those rows
  ##and columns that have sp. present in our community.
  my.com.dist <- dist.mat[com.names, com.names]
  ##set all diagonal values to NA so that the zero for conspecific comparisons
  #do not interfere with our calculation of nearest neighbors.
  diag(my.com.dist) <- NA
  
  #use apply() to calculate the minimun value in each row of the community
  #phylogenetic distance matrix and take a mean of those values
  mean(apply(my.com.dist, MARGIN = 1, min), na.rm = T)
}

apply(my.sample, MARGIN = 1, new.mntd.function)

mntd(my.sample, cophenetic(my.phylo), abundance.weighted = F)

new.mntd.a.function<-function(x){
#Get the names of sp. present in community.
com.names <- names(x[x > 0])  
##Make the community phylogenetic distance matrix using the
##names of the sp. present in the community.
my.com.dist <- dist.mat[com.names, com.names]

##Place NA values in the diagonals
diag(my.com.dist) <- NA
##calculate a mean of the minimum values in each row of the community
##phylogenetic distance matrix weighted by the abundances of the sp.
##present in the community
weighted.mean(apply(my.com.dist, 1, min, na.rm = T), x[x > 0])
}
apply(my.sample, MARGIN = 1, new.mntd.a.function)
mntd(my.sample, cophenetic(my.phylo), abundance.weighted = T)

sntd.function <- function(x){
  
  #Get the names of the sp. present in a community.
  com.names <- names(x[x > 0])  
  ##Make the communtiy phylogenetic distacne matrix by extracting those rows
  ##and columns that have sp. present in our community.
  my.com.dist <- dist.mat[com.names, com.names]
  
  ##set all diagonal values to NA so that the zero for conspecific comparisons
  #do not interfere with our calculation of nearest neighbors.
  diag(my.com.dist) <- NA

  ##Use apply() to calculate the minimun value in each row of the community phylogentic
  ##distance matrix and take a mean of those values.
  sd(apply(my.com.dist, MARGIN = 1, min), na.rm = T)
}
apply(my.sample, MARGIN = 1, sntd.function)
library(SDMTools)
sntd.a.function <- function(x){
  #Get the names of sp. present in community.
  com.names <- names(x[x > 0])  
  ##Make the community phylogenetic distance matrix using the
  ##names of the sp. present in the community.
  my.com.dist <- dist.mat[com.names, com.names]
  
  ##Place NA values in the diagonals
  diag(my.com.dist) <- NA
  ##calculate a mean of the minimum values in each row of the community
  ##phylogenetic distance matrix weighted by the abundances of the sp.
  ##present in the community
  wt.sd(apply(my.com.dist, 1, min, na.rm =T), x[x > 0])
}
apply(my.sample, MARGIN = 1, sntd.a.function)

#########3.6
richness <- rowSums(decostand(my.sample, method = "pa", MARGIN = 1 ))
pd.no.root <- pd(my.sample, my.phylo, include.root =  F)[,1]
pd.root <- pd(my.sample, my.phylo, include.root =  F)[,1]
pw <- mpd(my.sample, cophenetic(my.phylo), abundance.weighted = F)
pw.prime <- mpd(my.sample, cophenetic(my.phylo), abundance.weighted = T)
Rao.pw <- raoD(my.sample, my.phylo)$Dkk
helmus.psv <- psv(my.sample, my.phylo)[,1]
helmus.pse <- pse(my.sample, my.phylo)[,1]
helmus.psr <- psr(my.sample, my.phylo)[,1]
nn <- mntd(my.sample, cophenetic(my.phylo), abundance.weighted = F)
nn.prime <- mntd(my.sample, cophenetic(my.phylo), abundance.weighted = T)
output <- as.data.frame(cbind(richness, pd.no.root, pd.root, pw, pw.prime, Rao.pw,
                              helmus.psv, helmus.pse, helmus.psr, nn, nn.prime))
names(output) <-c("Richness", "Pd.No.Root", "Pd.Root", "MPD", "MPD.abund", "Rao", 
                  "Helmus.PSV", "Helmus.PSE", "Helmus.PSR", "MNTD", "MNTD.abund")
plot(output)
