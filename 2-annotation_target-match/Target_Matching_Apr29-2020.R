#REQUIRED PACKAGES
library(ape)
library(stringr)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

#install.packages("seqinr", dependencies = T, repos = "https://cloud.r-project.org")
library(seqinr)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warnings in loops

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#This script does the following:
#1. Matches the loci to the contigs, saves them to a new file
#2. Also finds the potential paralogs, removes them, and saves them to a separate file
  
#Set up directories
threads = "8" #threads, keep quotes
contig.save = "Rodent_Match"  #This is your save name for the big contig match file

# #CLUSTER directories
# work.dir = "/Volumes/Rodents/Rodents" #Your main project directory
# proc.dir = "/Volumes/Rodents/Rodents/Processed_Samples_MERGE"
# contig.dir = "/Volumes/Rodents/Rodents/Assembled_Contigs"
# genome.file = "/Volumes/Rodents/Rodents/GCF_000001635.26_GRCm38.p6_genomic.fna"
# #prot.file = "/Volumes/Armored/Rodents/GCF_000001635.26_GRCm38.p6_protein.faa"
# prot.file = "/Volumes/Rodents/Rodents/Genomes/Mus/protein-mouse_exons.fa"

#CLUSTER directories
work.dir = "/Volumes/Rodents/Rodents" #Your main project directory 
proc.dir = "/Volumes/Rodents/Rodents/Processed_Samples_MERGE"
contig.dir = "/Volumes/Rodents/Rodents/Assembled_Contigs"
prot.file = "/Volumes/Rodents/Selected_Transcripts/Mus_best_prot.fa"


##############################################################################################
#####################  1.Match loci to contigs                   #############################
#####################                                            #############################
##############################################################################################

########################################################################  
# Database setup
########################################################################

#Creates the directory if it doesn't exist
if (file.exists(proc.dir) == F){ dir.create(proc.dir) }

#headers for the blast db
headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
           "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

#gets lists of directories and files with sample names
setwd(contig.dir)
file.names = list.files(pattern = "", full.names = F, recursive = T)
out.dir = "assembled-contigs" 

#Matching and processing for each sample
for (i in 2:length(file.names)){
  
  #Sets up working directories for each species
  sample = gsub(pattern = "_contigs.fa$", replacement = "", x = file.names[i])
  species.dir = paste(proc.dir, "/", sample, sep = "")
  #Creates species directory if none exists
  if (file.exists(species.dir) == F){ dir.create(species.dir) }
  if (file.exists(paste(species.dir, "/", out.dir, sep = "")) == F) {dir.create(paste(species.dir, "/", out.dir, sep = "")) }
  setwd(paste(species.dir, "/", out.dir, sep = ""))
  
  #if (file.exists(paste0(sample, "_orthologs.fa")) == T){  next }
  
  # # DEDUPE almost exact duplicate removal
  system(paste("dedupe.sh in=",contig.dir, "/", file.names[i], " ordered=t overwrite=true ",
                " out=", sample, "_dd.fa", " minidentity=97", sep = ""), ignore.stderr = T)
  
  #Make blast database for the probe loci
  system(paste0("makeblastdb -in ", sample, "_dd.fa",
                " -parse_seqids -dbtype nucl -out ", sample, "_nucl-blast_db"))
  
  #Matches samples to proteins
  system(paste0("tblastn -task tblastn-fast -db ", sample, "_nucl-blast_db -evalue 0.001 -seg no",
                " -query ", prot.file, " -out ", sample, "_prot-match.txt", 
                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                " -num_threads ", threads))
   
  match.data = fread(paste(sample, "_prot-match.txt", sep =""), sep = "\t", header = F, stringsAsFactors = FALSE)
  setnames(match.data, headers)
  
  #Matches need to be greater than 12
  filt.data = match.data[match.data$matches > 12,]
  #Percent identitiy must match 50% or greater
  filt.data = filt.data[filt.data$pident >= 50,]  

  #Sorting: exon name, contig name, bitscore higher first, evalue
  setorder(filt.data, qName, tName, -pident, -bitscore, evalue)
  
  #Make sure the hit is greater than 50% of the reference length
  filt.data = filt.data[filt.data$matches >= (0.5 * filt.data$qLen),]
  #Mus: 100,741
  #Rat: 105,751
  
  #########################################################################
  #Part A: Fixes duplicate same contigs matching to the same locus (2 entries)
  #########################################################################
  
  contig.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)
  
  #Saves non duplicated data
  good.data = filt.data[!filt.data$qName %in% contig.names,]
  
  new.data = c()  
  save.paralog = c()
  for (j in 1:length(contig.names)) {
    sub.match = filt.data[filt.data$qName %in% contig.names[j],]
    #Skips if 1 row
    if (nrow(sub.match) == 1){ next }
    
    #Saves highest bitscore
    save.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),]
    
    #Saves longest if equal bitscores
    save.match = save.match[abs(save.match$qStart-save.match$qEnd) == max(abs(save.match$qStart-save.match$qEnd)),]
    
    #saves top match here  
    if (nrow(save.match) >= 2){ 
      save.match = save.match[1,]
      }
  
    new.data = rbind(new.data, save.match)
  
  } #end j
  
  #Saves final dataset  
  save.data = rbind(good.data, new.data)
  
  #Reads in contigs
  contigs = scanFa(FaFile(paste(sample, "_dd.fa", sep = "")))
  red.contigs = contigs[names(contigs) %in% save.data$tName]
  dup.contigs = save.data$tName[duplicated(save.data$tName)]
  dup.match = save.data[save.data$tName %in% dup.contigs,]
  dup.data = dup.match[order(dup.match$tName)]
  
  #Loops through each potential duplicate
  dup.loci = unique(dup.data$tName)
  fix.seq = DNAStringSet()
  for (j in 1:length(dup.loci)){
    #pulls out data that matches to multiple contigs
    sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
    sub.data = sub.data[order(sub.data$tStart)]
    
    #Fixes direction and adds into data
    #Finds out if they are overlapping
    for (k in 1:nrow(sub.data)){
        new.start = min(sub.data$tStart[k], sub.data$tEnd[k])
        new.end = max(sub.data$tStart[k], sub.data$tEnd[k])
        sub.data$tStart[k] = new.start
        sub.data$tEnd[k] = new.end
    }#end k loop
    
    #Saves them if it is split up across the same locus
    if (length(unique(sub.data$tName)) == 1 && length(unique(sub.data$qName)) == 1){
      spp.seq = contigs[names(contigs) %in% sub.data$qName]
      names(spp.seq) = paste(sub.data$tName[1], "_|_", sample, sep ="")
      fix.seq = append(fix.seq, spp.seq)
      next
    }
    
    #Cuts the node apart and saves separately
    sub.data$tStart = sub.data$tStart-(sub.data$qStart-1)
    #If it ends up with a negative start
    sub.data$tStart[sub.data$tStart <= 0] = 1
    #Fixes ends
    sub.data$tEnd = sub.data$tEnd+(sub.data$qLen-sub.data$qEnd)
    
    #Fixes if the contig is smaller than the full target locus
    sub.data$tEnd[sub.data$tEnd >= sub.data$tLen] = sub.data$tLen[1]
    
    starts = c()
    ends = c()
    starts[1] = 1
    for (k in 1:(nrow(sub.data)-1)){
      ends[k] = sub.data$tEnd[k]+floor((sub.data$tStart[k+1]-sub.data$tEnd[k])/2)
      starts[k+1] = ends[k]+1
    } #end k loop
    ends = append(ends, sub.data$tLen[1])
    
    #Looks for overlapping contigs
    tmp = ends-starts
    if(length(tmp[tmp < 0 ]) != 0){ 
      save.paralog = append(save.paralog, sub.data$qName)
      next
    }
    #Collects new sequence fragments
    spp.seq = contigs[names(contigs) %in% sub.data$tName]
    new.seq = DNAStringSet()
    for (k in 1:length(starts)){ new.seq = append(new.seq, subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }  
    
    #renames and saves
    names(new.seq) = paste(sub.data$qName, "_|_", sample, sep ="")
    fix.seq = append(fix.seq, new.seq)
  } #end j loop
  
  #Writes the base loci
  temp.data = save.data[!save.data$qName %in% gsub("_\\|_.*", "", names(fix.seq)),]
  base.data = temp.data[!temp.data$tName %in% save.paralog,]
  base.loci = contigs[names(contigs) %in% base.data$tName]
  sort.data = base.data[match(names(base.loci), base.data$tName),]
  #Name and finalize
  names(base.loci) = paste(sort.data$qName, "_|_", sample, sep ="")
  fin.loci = append(base.loci, fix.seq)
  fin.loci = fin.loci[width(fin.loci) >= 60]
  
  #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
  temp = fin.loci[duplicated(names(fin.loci)) == T]
  if(length(temp) != 0){ stop("DUPLICATE FOUND") }
  
  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(fin.loci))
  write.fasta(sequences = final.loci, names = names(final.loci), 
              paste(sample, "_orthologs.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  paralog.names = new.data[new.data$qName %in% save.paralog,]
  paralog.contigs = contigs[names(contigs) %in% paralog.names$tName]
  write.loci = as.list(as.character(paralog.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(sample, "_paralogs.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  system(paste("rm ", sample, "_dd.fa ", sep = ""))  
  system(paste0("rm *nucl-blast_db*"))
  print(paste(sample, " probe matching complete. ", length(final.loci), " found!"))

}#end i loop
  
########################################################################  
# STEP 2
# Output a file of summary stats
########################################################################

#Sets up data summary
setwd(work.dir)

header.data<-c("Sample", "noContigs", "orthoLoci", "minLen", "maxLen", "meanLen", "medianLen")   
samples<-gsub("_contigs.fa$", "", file.names)
prelim.data<-data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
setnames(prelim.data, header.data)
prelim.data[, Sample:=as.character(samples)]
merge.contigs<-DNAStringSet()

#Cycles through each assembly run and assesses each
for (i in 1:length(samples)){
  
  #Gets raw contigs to count them
  setwd(paste(proc.dir, "/", samples[i], "/", out.dir, sep = ""))
 
  #Gets length of raw contigs
  contigs<-scanFa(FaFile(paste(contig.dir, "/", samples[i], "_contigs.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("noContigs", header.data), value = length(contigs) )
  
  #Gets length of raw contigs
  og<-scanFa(FaFile(paste(samples[i], "_orthologs.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("orthoLoci", header.data), value = length(og) )
  
  set(prelim.data, i =  match(samples[i], samples), j = match("minLen", header.data), value = min(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("maxLen", header.data), value = max(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("meanLen", header.data), value = mean(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("medianLen", header.data), value = median(width(og)) )
  
  merge.contigs<-append(merge.contigs, og)
  
}#End loop for things

#Saves combined, final dataset
setwd(work.dir)
write.csv(prelim.data, file = paste(contig.save, "-sample-assessment.csv", sep = ""))

final.loci<-as.list(as.character(merge.contigs))
write.fasta(sequences = final.loci, names = names(final.loci), paste(contig.save, "_contigs.fa", sep = ""), nbchar = 1000000, as.string = T)

#### END SCRIPT 
# New files will be located in each species folder within 'Processed_Samples', in the folder assembled_loci
