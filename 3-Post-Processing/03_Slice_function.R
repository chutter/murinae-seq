library(ape)
library(GenomicRanges)
library(Biostrings)

#Consensus function
make.consensus<- function (input.alignment, method = c("majority", "threshold", "IUPAC", 
                                                     "profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = c("DNA", "RNA")) {
  
  #input.alignment<-trimmed
  #Converts alignment to matrix of characters to be used
  new.align<-strsplit(as.character(input.alignment), "")
  align.in<-matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)
  
  #Does based on method
  method <- match.arg(method)
  
  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC, 
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold, 
                  res, NA)
    names(res) <- NULL
  }
  
  out.consensus<-DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}
#Slice function
slice.trim<-function(input.align, slice.size.bp = 100, threshold = 0.45){
  
  #makes consensus sequence for comparison
  #input.align<-trimal.align
  input.con<-make.consensus(input.align, method = "majority")
  names(input.con)<-"Reference_Locus"
  
  comb.align<-append(input.align, input.con)
  
  #Gets slice information ready
  slice.no<-ceiling(max(width(input.align))/slice.size.bp)
  slice.start<-1
  slice.end<-slice.size.bp
  
  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){ 
    slice.end<-max(width(input.align))
  }#end if check
  output.align<-DNAStringSet()
  for (x in 1:slice.no){
    
    #Slice alignment into number of slices 
    sliced.align<-subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences 
    bad.align<-pairwise.inf.sites(sliced.align, "Reference_Locus")
    #Remove bad sequence chunks
    rem.seqs<-bad.align[bad.align >= threshold]
    good.align<-sliced.align[!names(sliced.align) %in% names(rem.seqs)]
    #Makes replacement gap seqs for the bad ones
    blank.align<-DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if
    
    #Saves the slices and cats
    save.slice<-append(good.align, blank.align)
    save.slice<-save.slice[order(names(save.slice))]
    save.names<-names(save.slice)
    output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align)<-save.names
    
    #Gets new start and stop
    slice.start<-slice.start+100
    slice.end<-slice.end+100
    #checks to see if its out of bounds
    if (slice.end > max(width(input.align))){ 
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){ break } else {
        save.slice<-subseq(comb.align, start = slice.start, end = slice.end)
        save.slice<-save.slice[order(names(save.slice))]
        save.names<-names(save.slice)
        output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
        names(output.align)<-save.names
        break
        }
    }#end if
  }#end x loop
  
  #Removes reference
  output.align<-output.align[names(output.align) != "Reference_Locus"]  
  #removes gap only taxa
  str.splitted<-strsplit(as.character(output.align), "")
  x.align<-as.matrix(as.DNAbin(str.splitted) )
  len.temp<-as.character(as.list(x.align))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= 20]  
  return.align<-output.align[!names(output.align) %in% names(spp.rem)]
  return(return.align)      
}#end FUNCTION

pairwise.inf.sites<-function(x, y) {
  #Alignment should be DNAStringSet
  # x<-m.align
  # y<-"Reference_Locus"
  temp.align<-strsplit(as.character(x), "")
  mat.align<-lapply(temp.align, tolower)
  m.align<-as.matrix(as.DNAbin(mat.align))
  
  #Filters out weirdly divergent sequences
  new.align<-as.character(m.align)
  new.align[new.align == "n"]<-"-"
  new.align[is.na(new.align) == T]<-"-"
  ref<-new.align[rownames(new.align) == y,]
  summary.data<-c()
  all.pars<-c()
  all.over<-c()
  for (z in 1:nrow(new.align)) {
    #Site counter
    pars<-0
    overlap<-0
    tar<-new.align[z,]
    combined<-matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,]<-ref
    combined[2,]<-tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col<-vector("character", length = nrow(combined))
      seq.col<-combined[,k]
      #not equal to -
      f.char<-seq.col[seq.col != '-'] 
      #don't count missing seq
      if (length(f.char) <= 1) { next }
      
      if (length(f.char) >= 2){
        overlap<-overlap+1
        if (f.char[1] != f.char [2]) { pars<-pars+1 }
      }#end if
    }#ends informative sites loop
    all.pars<-append(all.pars, pars)
    all.over<-append(all.over, overlap)
  }# ends seq loop
  #Summarizes and returns data
  summary.data<-all.pars/all.over
  summary.data[is.nan(summary.data)]<-0
  names(summary.data)<-rownames(new.align)
  return(summary.data)
}

write.phy<-function (x, file = "", interleave = FALSE, strict = FALSE){
  str2cha <- function(x) {
    unlist(strsplit(x, ""))
  }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""), 
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i], 
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames))) 
      cat("WARNING: Truncation of taxon names created", 
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ", 
                                                                  diff[i]), collapse = ""), sep = "")
  }
  if (!interleave) 
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1) 
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }
  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm))) 
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1) 
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1]))) 
      sm <- c(sm, "")
    phy <- c(phy, sm)
  }
  if (file == "") {
    cat(phy, sep = "\n")
  }
  else {
    write(phy, file = file)
  }
}


################################################################################################
################################################################################################
############## Step 1. Trim the exon-intron alignments combined     ############################
################################################################################################
################################################################################################
################################################################################################

#Out directory
input.dir = "input_directory"
out.dir = "output_directory"

#Filtering
min.len = 80 # Minimum alignment length to keep (after slicing)
min.taxa = 4 #Deletes if drops below this taxa number
min.cov = 0.5 #Minimum sample coverage (breadth) across alignment

#You can play with these settings if its too conservative
min.to.slice = 100 #Minimum size of alignment to slice 
slice.size = 80 #The slice size that it cuts up alignments
threshold = 0.40 #The threshold distance from consensus for throwing out a sample slice

#Creates direcotries
dir.create(out.dir)
setwd(input.dir)
locus.names = list.files(".")

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
 
  #Reads in alignment
  align = DNAStringSet(readAAMultipleAlignment(file = paste0(input.dir, "/", locus.names[i]), format = "phylip"))

  #Slice trims it, only above 100bp other it 
  if (width(align)[1] < min.to.slice){
  #Removes samples that too short individually
  write.temp = strsplit(as.character(align), "")
  aligned.set = as.matrix(as.DNAbin(write.temp) )
  #Writes alignment in output folder if not large enough
  write.phy(aligned.set, file= paste0(out.dir, "/", locus.names[i]), interleave = F)
  next
  }
  
  #Applies slice function
  red.align = slice.trim(align, slice.size.bp = slice.size, threshold = threshold)
  
  #Removes samples that too short individually
  write.temp = strsplit(as.character(red.align), "")
  aligned.set = as.matrix(as.DNAbin(write.temp) )
  len.temp = as.character(as.list(aligned.set))
  len.loci = lapply(len.temp, function (x) x[x != "-"])
  spp.len = unlist(lapply(len.loci, function (x) length(x)))
  spp.rem = spp.len[spp.len < (max(spp.len) * as.numeric(min.cov))]  
  spp.rem = append(spp.rem, spp.len[spp.len < min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(out.dir, "/", locus.names[i]), interleave = F)

}# end i looop
  