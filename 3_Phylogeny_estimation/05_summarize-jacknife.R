#############
## IQ TREE JACKNIFING
############

#The R packages needed
library(Biostrings)
library(data.table)
library(ape)
library(phangorn)
library(phytools)
options(stringsAsFactors = FALSE)

###### PARAMETER SETUP ###########

###### Directory Setup ###########
jack.reps = "/Volumes/Armored/Mantellidae_Subfamily/Sequence_capture/Jacknife_Trees/IQTrees_all"

jack.reps = "/Volumes/Armored/Mantellidae_Subfamily/Transcriptome-full/Jacknife_Trees"

#Outgroup selection 
outgroup.taxa<-c("Rana_kukunoris", "Pelophylax_nigromaculatus")


#########################################
# 1. Create majority rule consensus thing 
######################################### 

#Creates majority rule consensus tree from jacknifed trees
setwd(jack.reps)
tree.files = list.files(pattern = ".", full.names = F, recursive = F)

#Loop through trees and collapse poorly supported nodes into polytomies
system("rm ../alltrees.tre")
for (i in 1:length(tree.files)){
  #read in tree file
  temp.tree<-read.tree(tree.files[i])
  root.tree<-root(temp.tree, outgroup = outgroup.taxa, resolve.root = T)
  write.tree(root.tree, file = "../alltrees.tre", append = T)
}#end i loop

#Load back in the trees
all.trees<-read.tree(file = "../alltrees.tre")

#Computes the tree that most likely represents the topology of everything. 
rf.tree = averageTree(all.trees,method="symmetric.difference")
plotTree(root(rf.tree,outgroup=outgroup.taxa,resolve.root=TRUE))

#Consensus edges
t1<-consensus.edges(all.trees)
plotTree(t1,fsize=0.4)

t2<-consensus.edges(all.trees,if.absent="ignore")
plotTree(t2,fsize=0.4)

t3<-consensus.edges(all.trees,method="least.squares")
plotTree(t3,fsize=0.4)


### Older


#Much bette
plotTree(compute.brlen(root(all.trees[[ii]],outgroup=outgroup.taxa,
                            resolve.root=TRUE)))


densiTree(all.trees[1:500], type = "cladogram")

#Makes a majority rule consensus tree from the trees
cons.tree<-consensus(all.trees, p = 0.5, check.labels = TRUE)

write.tree(cons.tree, file = paste0(out.dir, "/majority_consensus.tre"))
unlink(paste0(out.dir, "/alltrees.tre"))
#
# ##### END SCRIPT
#
