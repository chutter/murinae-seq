library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(data.table)

#Parallelization
library(parallel)
library(foreach)
library(doParallel)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warnings in loops

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#General directory meanings
work.dir = "/Volumes/Rodents/Rodents/Selection-test/Site-models"
align.dir = "/Volumes/Rodents/Rodents/Alignments/Genes"

#General directory meanings
work.dir = "/home/c111h652/scratch/Rodents/Selection-test/Site-models"
align.dir = "/home/c111h652/scratch/Rodents/Genes"

threads = 16
mem.val = 160

################################################################################################
################################################################################################
############## Step 1. Set up the tree                              ############################
################################################################################################
################################################################################################
################################################################################################

#Starting stuff
setwd(work.dir)

#Load in tree and convert to zero length branches
tree = read.tree("starting_tree.tre")
tree$edge.length = NULL
tree$node.label = NULL
write.tree(tree, file = "rodents_nobl.tre")

################################################################################################
################################################################################################
############## Step 2. SITE MODELS CODEML                           ############################
################################################################################################
################################################################################################
################################################################################################

#Go to working directory
setwd(paste0(align.dir))

#Finds loci from the file names
file.names = list.files(".", recursive = T)

#Goes to the working directory
setwd(paste0(work.dir))

if (file.exists("output") != T) { dir.create("output")} else {
  done.gene = list.dirs("output/.", full.names = F, recursive = F)
  file.names = file.names[!file.names %in% paste0(done.gene, ".phy")]
}#end resume check

#Sets up multiprocessing
cl = makeCluster(threads)
registerDoParallel(cl)
mem.val = floor(mem.val/threads)

#Loops through each locus and does operations on them
foreach(i=1:length(file.names), .packages = c("ape", "stringr", "Biostrings","data.table", "Rsamtools")) %dopar% {
#for (i in 1:length(file.names)){
  
  ##############
  #STEP 1: Get the taxa in the alignments for these loci
  ##############
  
  #Copies alignment
  gene.name = gsub(".phy$", "", file.names[i])
  setwd(work.dir)
  out.dir = paste0(work.dir, "/output/", gene.name)
  dir.create(out.dir)
   
  #Fixes tree
  tree.file = read.tree("rodents_nobl.tre")
  align = DNAStringSet(readAAMultipleAlignment(file = paste0(align.dir, "/", file.names[i]), format = "phylip"))
  
  align.names = names(align)
  not.tree = tree.file$tip.label[!tree.file$tip.label %in% align.names]
  new.tree = drop.tip(tree.file, not.tree)
  write.tree(new.tree, file = paste0(out.dir, "/rodents.tre") )
  
  #The m1 models
  setwd(out.dir)
  dir.create("model1")
  setwd(paste0(out.dir, "/model1"))
  system(paste0("cp ", work.dir, "/codeml-control-m1.ctl ", out.dir, "/model1"))
  system(paste0("cp ", align.dir, "/", file.names[i], " ", out.dir, "/model1"))
  system(paste0("mv ", out.dir, "/model1/", file.names[i], " ", out.dir, "/model1/marker.phy"))
  system(paste0("cp ", out.dir, "/rodents.tre ", out.dir, "/model1"))
  system(paste0("codeml codeml-control-m1.ctl"))  
  
  #The m2 models
  setwd(out.dir)
  dir.create("model2")
  setwd(paste0(out.dir, "/model2"))
  system(paste0("cp ", work.dir, "/codeml-control-m2.ctl ", out.dir, "/model2"))
  system(paste0("cp ", align.dir, "/", file.names[i], " ", out.dir, "/model2"))
  system(paste0("mv ", out.dir, "/model2/", file.names[i], " ", out.dir, "/model2/marker.phy"))
  system(paste0("cp ", out.dir, "/rodents.tre ", out.dir, "/model2"))
  system(paste0("cp ", work.dir, "/codeml-control-m2.ctl ", out.dir, "/model2"))
  system(paste0("codeml codeml-control-m2.ctl"))  
  
}#end i loop
  
stopCluster(cl)


#### End Script

