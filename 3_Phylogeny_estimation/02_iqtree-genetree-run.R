#############
## IQ TREE
############

###### PARAMETER SETUP ####
#work.dir = "/Users/chutter/Dropbox/Research/WIP/Mantellidae_SeqCap/Gene_Trees"
#out.dir = "/Users/chutter/Dropbox/Research/WIP/Mantellidae_SeqCap/Gene_Trees/IQTrees_exon"
#align.dir = "/Volumes/Armored/Mantellidae_All/Alignments/exon-only_trimmed"

work.dir<-"/home/c111h652/scratch/Mantellidae_All/Trees/Gene_Trees"
out.dir<-"/home/c111h652/scratch/Mantellidae_All/Trees/Gene_Trees/IQTrees_all-markers"
align.dir<-"/home/c111h652/scratch/Mantellidae_All/Alignments/all-markers_trimmed"
threads = 8

#########################

#Go to wd and grab files and get sample names
setwd(align.dir)
file.names = list.files(pattern = ".", full.names = F, recursive = F)

#Checks for output directory; makes it if it doesn't exist
if (file.exists(out.dir) == F){ dir.create(out.dir) }
setwd(out.dir)

########################################
# First pass using hard-coded threads
########################################

done.loci = list.files(pattern = ".", full.names = F, recursive = F)
done.loci = gsub(".treefile", "", done.loci)
todo.loci = file.names[!file.names %in% done.loci]
if (length(todo.loci) == 0){ quit() }

#Loops through each locus and writes each species to end of file
for (i in 1:length(todo.loci)){

  #create output place
  system(paste0("cp ", align.dir, "/", todo.loci[i], " ", out.dir, "/", todo.loci[i]))
  
  #Run IQ tree on tree
  system(paste0("iqtree -s ", out.dir,"/", todo.loci[i], " -bb 1000 -nt ", threads, " -m MFP+MERGE -rcluster 20 -msub nuclear"))

  #Delete extra files
  del.files = dir(path = out.dir)
  temp.delete = del.files[grep(pattern = todo.loci[i], x = del.files, invert=F)]
  to.delete = temp.delete[grep(pattern = "treefile$", x = temp.delete, invert=T)]
  if (length(to.delete) != 0){ unlink(paste0(out.dir, "/", to.delete)) }
  
} #i loopo end
  

########################################
# Second pass using AUTO threads to pick up failed loci
########################################
  
done.loci = list.files(pattern = ".", full.names = F, recursive = F)
done.loci = gsub(".treefile", "", done.loci)
todo.loci = file.names[!file.names %in% done.loci]
if (length(todo.loci) == 0){ quit() }

#Loops through each locus and writes each species to end of file
for (i in 1:length(todo.loci)){
  
  #create output place
  system(paste0("cp ", align.dir, "/", todo.loci[i], " ", out.dir, "/", todo.loci[i]))
  
  #Run IQ tree on tree
  system(paste0("iqtree -s ", out.dir,"/", todo.loci[i], " -bb 1000 -nt AUTO -m MFP+MERGE -rcluster 20 -msub nuclear"))
  
  #Delete extra files
  del.files = dir(path = out.dir)
  temp.delete = del.files[grep(pattern = todo.loci[i], x = del.files, invert=F)]
  to.delete = temp.delete[grep(pattern = "treefile$", x = temp.delete, invert=T)]
  if (length(to.delete) != 0){ unlink(paste0(out.dir, "/", to.delete)) }
  
} #i loopo end

  
### END SCRIPT
  