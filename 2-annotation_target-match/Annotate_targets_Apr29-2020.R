############################################################
# For rodent exomes, 04.2020
# Retrieves the exons from each selected transcript and
# concatenates them together as CDS's.
############################################################

#REQUIRED PACKAGES
library(ape)
library(stringr)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(seqinr)
library(BSgenome)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warning


find.orf = function(input.seq, codons = F, min.size = 80){
  #Sets up data
  # input.seq<-trimmed[j]
  codon.table<-data.frame(Start = rep(0,6), End = rep(0,6), Frame = c("F1", "F2", "F3", "R1", "R2", "R3"))
  for.seq<-as.character(input.seq)
  
  #Gets codon stuff
  TAA<-matchPattern("TAA", for.seq)
  TGA<-matchPattern("TGA", for.seq)
  TAG<-matchPattern("TAG", for.seq)
  
  #Forward Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F1")
    codon.table<-rbind(codon.table, temp.table)
  } 
  
  #Forward Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Forward Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F3")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Sets up data
  rev.seq<-as.character(reverseComplement(input.seq))
  
  #Gets codon stuff
  TAA<-matchPattern("TAA", rev.seq)
  TGA<-matchPattern("TGA", rev.seq)
  TAG<-matchPattern("TAG", rev.seq)
  
  #Rev Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R1")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R3")
    codon.table<-rbind(codon.table, temp.table)
  } #end if 
  
  if (codons == T) { return(codon.table) }
  
  if (codons == F) {  
    frames<-unique(codon.table$Frame)
    orf.frame<-data.frame()
    for (x in 1:length(frames)){
      temp.codon<-codon.table[codon.table$Frame == frames[x],]
      temp.codon<-temp.codon[order(temp.codon$Start),]
      
      if (temp.codon$Start[1] == 0){
        temp.start<-as.numeric(gsub("F|R", "", temp.codon$Frame))
        add.frame<-data.frame(FrameStart = temp.start, FrameEnd = width(input.seq), 
                              Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
        orf.frame<-rbind(orf.frame, add.frame)
        next
      }
      #Goes through each of the given directions codons and converts to frame ranges
      temp.frame<-data.frame()
      for (y in 1:(nrow(temp.codon)+1)){
        #First y the start is 1, otherwise take from previous end
        if (y == 1){ frame.start<-as.numeric(gsub("F|R", "", temp.codon$Frame[y])) } else { frame.start<-temp.frame$FrameEnd[y-1]+4 }
        
        #Gets end by subtracting from the codon start
        frame.end<-temp.codon$Start[y]-1
        temp.frame<-rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
      } # end y loop
      
      temp.frame$FrameEnd[nrow(temp.frame)]<-width(input.seq)
      
      #Adds all the data together
      add.frame<-cbind(temp.frame, Size = (temp.frame$FrameEnd-temp.frame$FrameStart)+1, Frame = frames[x])
      orf.frame<-rbind(orf.frame, add.frame)
      
    } #end x loop
    
    orf.frame<-orf.frame[orf.frame$Size >= min.size,]
    return(orf.frame)
  } # end else
  
}# END FUNCTION

##################################################################################################
##################################################################################################
# Step 1: Coding exons 
##################################################################################################
##################################################################################################

# Hardcoded file name locations
ref = "/Volumes/Rodents/Rodents/Genomes/Mus/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
gtf.file = "/Volumes/Rodents/Rodents/Genomes/Mus/Mus_musculus.GRCm38.99.gtf"
transcript.file = "/Volumes/Rodents/Selected_Transcripts/selected-transcripts.txt"
out.dir = "/Volumes/Rodents/Selected_Transcripts"
logfilename = "get_selected_seqs.log"

#Read in genome
gen.seq = scanFa(FaFile(ref))
names(gen.seq) = gsub(" .*", "", names(gen.seq))

#Read in selected transcripts
transcript.data = data.table(read.table(transcript.file, sep = "\t", header = T))
header.names = colnames(transcript.data)
header.names = gsub("\\.", "_", header.names)
header.names = gsub("__", "_", header.names)
setnames(transcript.data, header.names) 

#Read in GTF
gtf.data = data.table(read.table(gtf.file, sep = "\t"))
header.names = c("chrom", "method", "feature_type", "start", "end", "un1", "strand", "un2", "info")
setnames(gtf.data, header.names) 
gtf.data[, un1 := NULL]
gtf.data[, un2 := NULL]
gtf.data[, method := NULL]
gtf.data[, gene_id := str_match(gtf.data$info, "gene_id (.*?);")[,2] ]
gtf.data[, transcript_id := str_match(gtf.data$info, "transcript_id (.*?);")[,2] ]
gtf.data[, prot_id := str_match(gtf.data$info, "protein_id (.*?);")[,2] ]
gtf.data[, exon_id := str_match(gtf.data$info, "exon_id (.*?);")[,2] ]
gtf.data[, exon_no := str_match(gtf.data$info, "exon_number (.*?);")[,2] ]

gtf.data[, info := NULL]
gtf.data = gtf.data[(gtf.data$end - gtf.data$start) >= 50,]
gtf.data = gtf.data[gtf.data$feature_type != "",]
gtf.data$feature_type = gsub("_", "-", gtf.data$feature_type)

#Next loop through each gene from the selected transcript
#And use the gtf file to obtain the exon delimitation
#And separate and save the exons 

header.data = c("target_name", "feature_type", "coding", "gene_id", "transcript_id", 
                "protein_id", "exon_id", "exon_no", "chrom", "start", "end", "length", "strand")   
collect.data = data.table(matrix(as.character(0), nrow = nrow(gtf.data), ncol = length(header.data)))
setnames(collect.data, header.data)
collect.data[, exon_no:=as.integer(exon_no)]
collect.data[, start:=as.integer(start)]
collect.data[, end:=as.integer(end)]

final.seqs = DNAStringSet()
index.val = as.integer(1)
for (i in 1:nrow(transcript.data)){
  
  temp.gtf = gtf.data[gtf.data$transcript_id %in% transcript.data$Transcript_stable_ID[i], ]
  temp.gtf = temp.gtf[order(start)]
  
  if (nrow(temp.gtf) == 0){ next }
  
  #Combines the exon and cds data together to more easily give ID to cds data
  gtf.cds = temp.gtf[temp.gtf$feature_type == "CDS",]
  gtf.cds[, exon_id:=NULL]
  gtf.exon = temp.gtf[temp.gtf$feature_type == "exon",]
  gtf.exon = data.frame(exon_id = gtf.exon$exon_id, exon_no = gtf.exon$exon_no)
  gtf.done = merge(gtf.cds, gtf.exon, by = "exon_no", all = F)
  
  if (nrow(gtf.done) == 0){ next }

  gtf.ranges = makeGRangesFromDataFrame(gtf.done,
                           keep.extra.columns=T,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field="chrom",
                           start.field="start",
                           end.field="end",
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
  
  exon.seqs = BSgenome::getSeq(gen.seq, gtf.ranges)
  names(exon.seqs) = paste0("chr", gtf.done$chrom, "_trid", 
                            gtf.ranges$transcript_id, "_exid", 
                            gtf.ranges$exon_id)
  
  for (j in 1:length(exon.seqs)){
    orf.data = find.orf(exon.seqs[j], codons = F, min.size = width(exon.seqs[j])*0.90)
    orf.data = orf.data[orf.data$Frame != "R1",]
    orf.data = orf.data[orf.data$Frame != "R2",]
    orf.data = orf.data[orf.data$Frame != "R3",]
    
    if (nrow(orf.data) == 0){ next }
        
    best.orf = orf.data[orf.data$Size == max(orf.data$Size),][1,]

    if (best.orf$Frame == "F1"){ exon.seqs[j] = subseq(exon.seqs[j], start = 1, end = width(exon.seqs[j])) }
    if (best.orf$Frame == "F2"){ exon.seqs[j] = subseq(exon.seqs[j], start = 2, end = width(exon.seqs[j])) }
    if (best.orf$Frame == "F3"){ exon.seqs[j] = subseq(exon.seqs[j], start = 3, end = width(exon.seqs[j])) }
    
    if (round(width(exon.seqs[j])/3) != width(exon.seqs[j])/3){
      exon.seqs[j] = subseq(exon.seqs[j], start = 1, end = width(exon.seqs[j])-1)
      if (round(width(exon.seqs[j])/3) != width(exon.seqs[j])/3){
        exon.seqs[j] = subseq(exon.seqs[j], start = 1, end = width(exon.seqs[j])-1)
      }#end 2nd if
    }#end first if
    
    #Obtains relevant metadata
    exon.ranges = gtf.done[gtf.done$exon_id == gsub(".*_exid", "", names(exon.seqs)[j])]
    
    #Collect the data
    set(collect.data, i =  index.val, j = match("target_name", header.data), value = names(exon.seqs)[j] )
    set(collect.data, i =  index.val, j = match("feature_type", header.data), value = exon.ranges$feature_type)
    set(collect.data, i =  index.val, j = match("coding", header.data), value = "Yes")
    set(collect.data, i =  index.val, j = match("gene_id", header.data), value = exon.ranges$gene_id)
    set(collect.data, i =  index.val, j = match("transcript_id", header.data), value = exon.ranges$transcript_id)
    set(collect.data, i =  index.val, j = match("protein_id", header.data), value = exon.ranges$prot_id)
    set(collect.data, i =  index.val, j = match("exon_id", header.data), value = exon.ranges$exon_id)
    set(collect.data, i =  index.val, j = match("exon_no", header.data), value = exon.ranges$exon_no)
    set(collect.data, i =  index.val, j = match("chrom", header.data), value = exon.ranges$chrom)
    set(collect.data, i =  index.val, j = match("start", header.data), value = exon.ranges$start)
    set(collect.data, i =  index.val, j = match("end", header.data), value = exon.ranges$end)
    set(collect.data, i =  index.val, j = match("length", header.data), value = exon.ranges$end-exon.ranges$start)
    set(collect.data, i =  index.val, j = match("strand", header.data), value = exon.ranges$strand)
    index.val = as.integer(index.val + 1)
  }#end j loop
  
  #Save final sequences 
  final.seqs = append(final.seqs, exon.seqs)
  
}# end i loop

#Saves the sequences 
setwd(out.dir)

#Mus filtered exons
save.seqs = final.seqs[width(final.seqs) >= 60]
save.seqs = save.seqs[names(save.seqs) %in% collect.data$target_name]
final.loci = as.list(as.character(save.seqs))
write.fasta(sequences = final.loci, names = names(final.loci), 
            "Mus_best_cds.fa", nbchar = 1000000, as.string = T)

#Translate to proteins
trans.seqs = Biostrings::translate(save.seqs)
final.loci = as.list(as.character(trans.seqs))
write.fasta(sequences = final.loci, names = names(final.loci), 
            "Mus_best_prot.fa", nbchar = 1000000, as.string = T)

#Writes metadata for sequences
write.data = collect.data[collect.data$target_name != "0",]
write.data = write.data[write.data$target_name %in% names(save.seqs),]
write.csv(write.data, file = "Mus_best_cds_metadata.csv")

##################################################################################################
##################################################################################################
# Step 2: Non-coding exons and UTRs
##################################################################################################
##################################################################################################

### Gotta deal with UTR and non-coding exons differently 

header.data = c("target_name", "feature_type", "coding", "gene_id", "transcript_id", 
                "protein_id", "exon_id", "exon_no", "chrom", "start", "end", "length", "strand")   
collect.data = data.table(matrix(as.character(0), nrow = nrow(gtf.data), ncol = length(header.data)))
setnames(collect.data, header.data)
collect.data[, start:=as.integer(start)]
collect.data[, end:=as.integer(end)]

final.seqs = DNAStringSet()
index.val = as.integer(1)
for (i in 1:nrow(transcript.data)){
  
  temp.gtf = gtf.data[gtf.data$transcript_id %in% transcript.data$Transcript_stable_ID[i], ]
  temp.gtf = temp.gtf[order(start)]
  
  if (nrow(temp.gtf) == 0){ next }

  #Combines the exon and cds data together to more easily give ID to cds data
  gtf.cds = temp.gtf[temp.gtf$feature_type == "CDS",]
  gtf.cds[, exon_id:=NULL]
  gtf.exon = temp.gtf[temp.gtf$feature_type == "exon",]
  gtf.exon = data.frame(exon_id = gtf.exon$exon_id, exon_no = gtf.exon$exon_no)
  gtf.comb = merge(gtf.cds, gtf.exon, by = "exon_no", all = T)
  
  non.code.exon = gtf.comb[is.na(gtf.comb$chrom) == T,]$exon_no
  non.code = temp.gtf[temp.gtf$exon_no %in% non.code.exon,]
  gtf.utr = temp.gtf[grep("_utr", temp.gtf$feature_type)]
  gtf.done = rbind(non.code, gtf.utr)
  
  if (nrow(gtf.done) == 0){ next }
  
  gtf.ranges = makeGRangesFromDataFrame(gtf.done,
                                        keep.extra.columns=T,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field="chrom",
                                        start.field="start",
                                        end.field="end",
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE)
  
  # exon.range = gtf.ranges[gtf.ranges$feature_type == "exon",]
  # utr.range = gtf.ranges[gtf.ranges$feature_type != "exon",]
  # 
  # exon.seqs = BSgenome::getSeq(gen.seq, exon.range)
  # utr.seqs = BSgenome::getSeq(gen.seq, utr.range)
  # 
  # names(exon.seqs) = paste0("Noncoding-exon_chr", exon.range$chrom, "_trid", 
  #                           exon.range$transcript_id, "_exid", 
  #                           exon.range$exon_id)
  # 
  # names(utr.seqs) = paste0("Noncoding-utr_chr", utr.range$chrom, "_trid", 
  #                           utr.range$transcript_id, "_", 
  #                           utr.range$feature_type, "-", rep(1:length(utr.seqs)))
  # 
  # nc.seqs = append(exon.seqs, utr.seqs)
  # 
  save.seqs = DNAStringSet()
  utr.count = 0
  for (j in 1:length(gtf.ranges)){
 
    #Gets the sequence
    nc.seq = BSgenome::getSeq(gen.seq, gtf.ranges[j])
    
    #Obtains relevant metadata
    exon.range = gtf.done[j,]
 
    if (exon.range$feature_type == "exon"){
      names(nc.seq) = paste0("Noncoding-exon_chr", exon.range$chrom, "_trid", 
                                exon.range$transcript_id, "_exid", 
                                exon.range$exon_id)
    }#end if
      
    if (exon.range$feature_type != "exon"){
      utr.count = utr.count + 1
      names(nc.seq) = paste0("Noncoding-utr_chr", exon.range$chrom, "_trid", 
                               exon.range$transcript_id, "_", 
                               exon.range$feature_type, "-", utr.count)
    }#end if 
    
    #Collect the data
    set(collect.data, i =  index.val, j = match("target_name", header.data), value = names(nc.seq) )
    set(collect.data, i =  index.val, j = match("feature_type", header.data), value = exon.range$feature_type)
    set(collect.data, i =  index.val, j = match("coding", header.data), value = "Yes")
    set(collect.data, i =  index.val, j = match("gene_id", header.data), value = exon.range$gene_id)
    set(collect.data, i =  index.val, j = match("transcript_id", header.data), value = exon.range$transcript_id)
    set(collect.data, i =  index.val, j = match("protein_id", header.data), value = exon.range$prot_id)
    set(collect.data, i =  index.val, j = match("exon_id", header.data), value = exon.range$exon_id)
    set(collect.data, i =  index.val, j = match("exon_no", header.data), value = exon.range$exon_no)
    set(collect.data, i =  index.val, j = match("chrom", header.data), value = exon.range$chrom)
    set(collect.data, i =  index.val, j = match("start", header.data), value = exon.range$start)
    set(collect.data, i =  index.val, j = match("end", header.data), value = exon.range$end)
    set(collect.data, i =  index.val, j = match("length", header.data), value = exon.range$end-exon.range$start)
    set(collect.data, i =  index.val, j = match("strand", header.data), value = exon.range$strand)
    index.val = as.integer(index.val + 1)
    
    save.seqs = append(save.seqs, nc.seq)
  }#end j loop
  
  #Save final sequences 
  final.seqs = append(final.seqs, save.seqs)
  
}# end i loop

#Saves the sequences 
setwd(out.dir)

#Mus filtered exons
save.seqs = final.seqs[width(final.seqs) >= 60]
final.loci = as.list(as.character(save.seqs))
write.fasta(sequences = final.loci, names = names(final.loci), 
            "Mus_best_noncode.fa", nbchar = 1000000, as.string = T)

#Write the table data
write.data = collect.data[collect.data$target_name != "0",]
write.data = write.data[write.data$target_name %in% names(save.seqs),]
write.csv(write.data, file = "Mus_best_noncode_metadata.csv")


#### ENd script