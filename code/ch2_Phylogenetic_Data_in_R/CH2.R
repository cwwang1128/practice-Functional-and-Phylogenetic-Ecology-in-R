#2.2
library(ape)
setwd("C://Users/bill/Desktop/Functional and Phylogenetic Ecology in R/supplementary material/ch2_Phylogenetic Data in R")
#Newick format
my.phylo<-read.tree("my.phylo.newick.file.txt")
write.tree(my.phylo)
###Newick format
my.nexus.phylo<-read.nexus("my.phylo.nexus.file.txt")
write.nexus(my.phylo)

##see the object
my.phylo
#tips=terminal modes

##is.rooted
is.rooted(my.phylo)
##ll tips end at the same 
##point as is true for a phylogeny of extant sp. scaled to time)
is.ultrametric(my.phylo)

##data's names
names(my.phylo)

my.phylo$edge
my.phylo$tip.label
my.phylo$Nnode
my.phylo$node.label
my.phylo$edge.length

#clumn bind
cbind(my.phylo$edge, my.phylo$edge.length)

#2.3
plot(my.phylo)

plot(my.phylo, cex=0.8)

add.scale.bar( ,length = 0.1)


plot.phylo(my.phylo, type = "fan", show.tip.label = T, show.node.label = F, 
           edge.color = "blue", edge.width = 1, tip.color = "blue")

plot(my.phylo)
#place internal node number
nodelabels( ,col= "black", bg="gray")
#place tip node number
tiplabels( ,col= "black", bg="gray")

#2.4

my.subtrees<- subtrees(my.phylo)


my.subtrees[[15]]
plot(my.subtrees[[15]])

drop.tree<-drop.tip(my.phylo, c("e","j","s"))
plot(drop.tree)


##?‹¥ä¸¦é?žç‚ºbifurcatingï¼Œå¯ä»¥åˆ©?”¨"multi2di()"å°‡è»¸??†é??
my.poly.phylo<-read.tree("example.poly.txt")
my.poly.phylo
plot(my.poly.phylo)
my.di.phylo<-multi2di(my.poly.phylo)
my.di.phylo
plot(my.di.phylo)

#branching time for eaxh internal node for an ultrametric phylogeny
branching.times(my.phylo)
##a vector of values with a length equal to the number of internal nodes in the phylogeny.
##The order of the nodes is from the root of the phylogeny towards the tips


##phylogenetic distance matrix
p.dist.mat<- cophenetic(my.phylo)
p.dist.mat[1:4, 1:4]
##matrix has the species names as the row and column names
##the values are the sum of the branch lengths separating each pair of species
##two species are far apart on the phylogeny their distance will be larger than that for two closely related species


##phylogenetic variance-covariance(VCV) matrix
vcv(my.phylo)
##The diagonal values of the matrix are the distances from the root to the tip that contains that species. 
##The off diagonal values indicate the amount of shared branch length between two species.
##Assuming a Brownian Motion model of trait evolution, 
#the diagonal values are used to estimate the expected variance in a trait, and the off diagonal values are used to estimate the expected covariance in the trait values between two species.
##Essentially high off diagonal values mean species that are more closely related and are expected to have more similar trait values.




#2.5
##generate a random phylogeny with random splitting containing 40 species.
new.tree<- rtree(40, rooted = T, tip.label = NULL)
plot(new.tree)
##This generated a random coalescent tree.
new.ultra.tree<-rcoal(40)
plot(new.ultra.tree)


for (i in 1:9){
  ##Make a random tree with 10 terminal taxa
  temp.random.tree<-rtree(10)
  ##Write the random tree to your working
  ##directory with each new file having a new number
  write.tree(temp.random.tree, paste(i,".txt",seq=""))
##close the loop
}

#Exercises
##Exercises2.1
##Make a Newick file for five species (speciesA, speciesB, speciesC, speciesD, and speciesE) 
##where speciesA and speciesE are most closely related to one another,
##speciesB and speciesD are most closely related to one another,
##and speciesC is most closely related to speciesA and speciesE.
##The file should have no branch lengths.
##Read this file into R and plot it with blue branches.
Exercises2.1.rawdata<-"((C,(A,E)),(D,B));"
Exercises2.1.data<-read.tree(text=Exercises2.1.rawdata)
write.tree(Exercises2.1.data)
plot(Exercises2.1.data,edge.color = "blue")
title("Exercises2.1")

##Exercises2.2 #unsure answer
Exercises2.2.rawdata<-"((C:3,(A:3,E:3):3):3,(D:3,B:3):3);"
Exercises2.2.data<-read.tree(text=Exercises2.2.rawdata)
write.tree(Exercises2.2.data)
summary.phylo(Exercises2.2.data)
plot(Exercises2.2.data)
title("Exercises2.2")

##Exercises2.3

##Exercises2.4

##Exercises2.5

##Exercises2.6

##Exercises2.7

##Exercises2.8
