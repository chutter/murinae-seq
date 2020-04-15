library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

options(stringsAsFactors = FALSE)

###### PARAMETER SETUP ALL ####
#work.dir = "/Volumes/Armored/Mantellidae_Subfamily/Transcriptomes_full"
#tree.dir = "/Volumes/Armored/Mantellidae_Subfamily/Transcriptomes_full/IQTrees_trans"
#astral.path = "/usr/local/bin/Astral/astral.5.14.2.jar" #if R can't find path

work.dir = "/home/c111h652/scratch/Mantellidae_All/Trees/Gene_Trees"
tree.dir = "/home/c111h652/scratch/Mantellidae_All/Trees/Gene_Trees/IQTrees_all-markers"
astral.path = "/home/c111h652/programs/Astral/Astral" #if R can't find path
astral.version = "astral.5.14.2.jar"

#Settings
output.file = "all-markers" #File name

mem = 40 #in GB
poly.lim = 10 #Bootstrap when to collapse polytomies
max.tax = 118 #The total number of samples
miss.data = 1 #Missing data allowed to keep tree, 1 to not considering missing data

#Drop taxa from the tree, do not drop outgroups
drop.taxa = c()


###############################################################################
###############################################################################
###################### Functions                      #########################
###############################################################################
###############################################################################

# Adapted di2multi function from the ape package to plot polytomies
# based on numeric node support values
# (di2multi does this based on edge lengths)
# Needs adjustment for unrooted trees as currently skips the first edge
di2multi4node <- function (phy, tol = 0.5) {
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  if (is.na(as.numeric(phy$node.label[2])))
    stop("node labels can't be converted to numeric values")
  if (is.null(phy$node.label))
    stop("the tree has no node labels")
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) 
    return(phy)
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) 
      next
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}

###############################################################################
###############################################################################
###################### Do not edit below              #########################
###############################################################################
###############################################################################

setwd(tree.dir)
tree.files = list.files(".")
tree.files = tree.files[grep(".tre", tree.files)]
#Loop through trees and collapse poorly supported nodes into polytomies
for (i in 1:length(tree.files)){
  #read in tree file
  temp.tree<-read.tree(tree.files[i])
  
  #Drops the outgroups if you want
  if (length(drop.taxa) != 0){ temp.tree = drop.tip(temp.tree, drop.taxa) }
  
  if (length(temp.tree$tip.label) < 4 ){ next }
  #If missing data isn't 1
  if (miss.data != 1){
    if (length(temp.tree$tip.label) <= floor(max.tax*miss.data) ){ next }
  } #end outer if
    
  #Removes node labels and polytomies
  if(length(temp.tree$node.label) == 0){
    #stop()
    write.tree(temp.tree, file = paste0(work.dir, "/", output.file, "_genetrees.tre"), append = T)
  } else {
    #replace empty with 0s
    temp.tree$node.label[temp.tree$node.label == ""]<-"0"
    #Find and collapse nodes with bl close to 0 from above 
    new.tree<-di2multi4node(temp.tree, tol = poly.lim) 
    write.tree(new.tree, file = paste0(work.dir, "/", output.file, "_genetrees.tre"), append = T)
  }#end else
}#end if

#Save final tree file into astral folder
setwd(work.dir)
#Run Astral 
system(paste0("java -Xmx", mem, "g -D", paste0('"java.library.path=', astral.path, '/lib/"'), 
              " -jar ", astral.path, "/", astral.version, 
              " -i ",  output.file, "_genetrees.tre -t 2 -o ", output.file, "_astral.tre"))



###### END SCRIPT


