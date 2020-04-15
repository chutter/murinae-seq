#############
## IQ TREE JACKNIFING
############

#The R packages needed
library(Biostrings)
library(data.table)
library(ape)
library(phangorn)
options(stringsAsFactors = FALSE)

###### PARAMETER SETUP ###########

#Set up how you want to decide the size of each alignment. Can only choose one! 
selection.method<-"basepairs" #Number of BASE PAIRS to include in each jacknife run
#selection.method<-"genes" #OR
jacknife.size<-50000 #Number of GENES or BASEPAIRS
jacknife.reps<-1000 #Number of total jacknife runs

#Sets up the running parameters and other details
locus.completeness<-50 #50% complete matrix, filters out loci with less than 50% sampled taxa
min.locus.size<-100 #min individual locus size
codon.partition<-F #partition by codon position or not 
merge.parts<-T #acts like parition finder and merges partitions if they have similar rates
threads<-8

###### Directory Setup ###########
work.dir<-"/home/c111h652/scratch/Mantellidae/Jacknife_Trees"
out.dir<-"/home/c111h652/scratch/Mantellidae/Jacknife_Trees/IQTrees_exon"
align.dir<-"/home/c111h652/scratch/Mantellidae/Alignments/exon-only_trimmed"

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

setwd(align.dir)
locus.names<-list.files(pattern = ".", full.names = F, recursive = F)

#Go to wd and grab files and get sample names
dir.create(out.dir)
setwd(out.dir)
dir.create("jacknife_trees")

#########################################
# 1. Create percent completeness list of loci
######################################### 

#Get number of runs based off number of files
#Loops through each locus and writes each species to end of file
bad.loci<-c()
taxa.count<-data.table(Locus = locus.names, Count = as.numeric(0), Length = as.numeric(0))
for (i in 1:length(locus.names)){
  #Reads in files
  align<-readAAMultipleAlignment(file = paste(align.dir, "/", locus.names[i], sep =""), format = "phylip")
  
  #Use taxa remove
  tax.names<-rownames(align)

  #########################################
  # A. Excludes specific loci based on selected parameters
  #########################################
  #removes too short loci
  if (ncol(align) <= min.locus.size){ bad.loci<-append(bad.loci, locus.names[i]) }
  
  #########################################
  # C. Generates loci stats for percent completeness matrices 
  #########################################
  
  taxa.count[i,2]<-length(tax.names)
  taxa.count[i,3]<-ncol(align)
  
} #end i loop

fin.tax<-taxa.count[!taxa.count$Locus %in% bad.loci,]  
max.tax<-max(fin.tax$Count)
dataset<-fin.tax[(fin.tax$Count/max.tax)*100 > locus.completeness]

#########################################
# 2. Create sampling datasets and concatenate and run
######################################### 

#Gets number of runs
locus.list<-dataset$Locus
for (i in 1:jacknife.reps){

  #samples without replacement for loci 
  setwd(out.dir)
  
  #chooses loci based on number of genes
  if (selection.method == "genes"){
    sample.loci<-locus.list[sample.int(n = length(locus.list), size = jacknife.size, replace = F)]
    if (jacknife.size >= length(locus.list)){ stop("ERROR. Jacknife size greater than number of available genes.")}
  }#end if
  
  if (selection.method == "basepairs"){
    #Gets the sampling order
    sample.order<-locus.list[sample.int(n = length(locus.list), size = length(locus.list), replace = F)]
    #Uses a while loop to randomly add loci
    loci.count<-1
    total.bp<-0
    sample.loci<-c()
    while (total.bp <= jacknife.size){
      #create list 
      #add loci to list
      sample.loci<-append(sample.loci, sample.order[loci.count])
      #Gets size
      temp.dataset<-dataset[dataset$Locus %in% sample.loci]
      total.bp<-sum(temp.dataset$Length)
      loci.count<-loci.count+1
    }#end WHILE
  }#end IF
  
  #Copies files to folder
  dir.create(paste0(out.dir, "/temp_trees"))
  system(paste0("cp ", paste0(align.dir, "/", sample.loci, collapse = " "),  " ", out.dir, "/temp_trees"))
  #concatenates files
  setwd(paste0(out.dir, "/temp_trees"))
  system(paste0("AMAS.py concat -f phylip -d dna -i *phy -u phylip --part-format raxml"))
  
  #Copies and pastes the concatenated files to output directory, deletes and renames
  system(paste0("cp concatenated.out partitions.txt ", out.dir))
  setwd(out.dir)
  system(paste0("mv concatenated.out iqtree_jacknife.phy"))
  system(paste0("mv partitions.txt iqtree_jacknife_parts.txt"))
  system(paste0("rm -r ", out.dir, "/temp_trees"))
  
  #Sets up parameter type selections from above
  part.scheme<-"MFP"
  part.file<-""
  if (merge.parts == T){ part.scheme<-paste0(part.scheme, "+MERGE") }
  if (merge.parts == F){ part.file<-" -spp iqtree_jacknife_parts.txt" }
  if (codon.partition == T){ codon.st<-" -st CODON" } else { codon.st<-"" }
  
  #Runs IQTree
  system(paste0("iqtree -s ", out.dir, "/iqtree_jacknife.phy", part.file,
                " -bb 1000 -nt ", threads, " -m ", part.scheme, codon.st, " -rcluster 10 -msub nuclear"))
  
  #Delete extra files
  del.files<-dir(path = out.dir)
  to.delete<-del.files[grep(pattern = "treefile$", x = del.files, invert=T)]
  unlink(paste0(out.dir, "/", to.delete))
  
  #moves tree file to save place
  if (merge.parts == F){ system(paste0("mv iqtree_jacknife_parts.txt.treefile jacknife_", i, ".tre")) } else {
    system(paste0("mv iqtree_jacknife.phy.treefile jacknife_", i, ".tre"))}
  system(paste0("cp jacknife_", i, ".tre ", out.dir, "/jacknife_trees"))
  unlink(paste0("jacknife_", i, ".tre"))

  #Removes loci already done to fullfill the jacknifing
  #locus.list<-locus.list[!locus.list %in% sample.loci]
  
}#end i loop


##### END SCRIPT: MOVE ONTO 05 TO SUMMARIZE RESULTS