library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(data.table)
library(treeio)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warnings in loops

read.mlc = function (mlcfile){
  mlc = readLines("output.txt")
  mlc = readLines(mlcfile)
  x = grep("dN & dS for each branch", mlc)
  y = grep("Time used", mlc)

  #Checks if all have data
  if (length(x) == 0 || length(y) == 0) { return(NULL) }
  
  #Formats it
  mlc = mlc[x:y-2]
  hi = grep("dN/dS", mlc)
  cn = unlist(strsplit(mlc[hi], " "))
  cn = cn[cn != ""]
  
  ii = grep("\\d+\\.\\.\\d+", mlc)
  info = mlc[ii]
  info = sub("^\\s+", "", info)
  info = sub("\\s+$", "", info)
  res <- lapply(info, function(x) {
    y <- unlist(strsplit(x, "\\s+"))
    edge <- unlist(strsplit(y[1], "\\.\\."))
    yy <- c(edge, y[-1])
    as.numeric(yy)
  }) %>% do.call("rbind", .)
  
  row.names(res) <- NULL
  colnames(res) <- c("parent", "node", cn[-1])
  colnames(res) <- gsub("\\*", "_x_", colnames(res))
  colnames(res) <- gsub("\\/", "_vs_", colnames(res))
  return(res)
} 

find.BEB = function (mlcfile){
  mlc = readLines("output.txt")
  mlc = readLines(mlcfile)
  x = grep("dN & dS for each branch", mlc)
  y = grep("Time used", mlc)

  #Checks if all have data
  if (length(x) == 0 || length(y) == 0) { return(NULL) }
  
  #Formats it
  mlc = mlc[x:y-2]
  hi = grep("Bayes Empirical Bayes", mlc)+4
  cn = unlist(strsplit(mlc[hi], " "))
  cn = cn[cn != ""]
  
  start = grep("Bayes Empirical Bayes", mlc)+6
  
  data.lines = c()
  ind.val = 0
  stop = T
  while (stop == T){
    curr.line = start + ind.val
    curr.data = unlist(strsplit(mlc[curr.line], " "))
    curr.data = curr.data[curr.data != ""]
    curr.data = curr.data[curr.data != "+-"]
    
    if (length(curr.data) == 0){ 
      stop = F
      break
      }
    
    curr.final = data.frame(Codon = curr.data[1], AA = curr.data[2], Post = curr.data[3],
                            postMean = curr.data[4], SEw = curr.data[5])

    data.lines = rbind(data.lines, curr.final)
    ind.val = ind.val + 1
  }# end while loop
    
  return(data.lines)
  
}#end function
  

find.lnl = function (mlcfile){
  mlc = readLines(mlcfile)
  lnl.line = mlc[grep("lnL", mlc)]
  lnl.line = gsub("):", "", lnl.line)
  
  cn = unlist(strsplit(lnl.line, " "))
  cn = cn[cn != ""]
  
  res = data.frame(ntime = as.numeric(cn[2]), np = as.numeric(cn[4]), lnl = as.numeric(cn[5]))
  return(res)
}#end function 

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#General directory meanings
work.dir = "/Volumes/Rodents/Selection-test/Site-models/output" #zipped files
align.dir = "/Volumes/Rodents/Rodents/Alignments/Genes"
threads = 2
mem.val = 8

################################################################################################
################################################################################################
############## Step 1. Parse the shit out of it                     ############################
################################################################################################
################################################################################################
################################################################################################

#Starting stuff
setwd(work.dir)

#Finds loci from the file names
file.names = list.files(".", recursive = F)

#Loops through each locus and does operations on them
save.data = c()
save.beb = c()

#Loops through each locus and does operations on them
#foreach(i=1:length(file.names), .packages = c("ape", "stringr", "Biostrings","data.table", "Rsamtools")) %dopar% {
  
for (i in 1:length(file.names)){
  
  ##############
  #STEP 1: Get the taxa in the alignments for these loci
  ##############
  
  #Copies alignment
  setwd(work.dir)
  gene.name = gsub(".zip$","", file.names[i])
  if (file.exists(gene.name) == T ) { system(paste0("rm -r ", work.dir, "/", gene.name)) }
  system(paste0("unzip ", file.names[i]), ignore.stdout = T)
  
  #Fixes tree
  #tree.file = read.tree("rodents_nobl.tre")
  #align = DNAStringSet(readAAMultipleAlignment(file = paste0(align.dir, "/", file.names[i]), format = "phylip"))
  
  #align.names = names(align)
  #not.tree = tree.file$tip.label[!tree.file$tip.label %in% align.names]
  #new.tree = drop.tip(tree.file, not.tree)
  #write.tree(new.tree, file = paste0(out.dir, "/rodents.tre") )
  
  #Copies working files
  
  #The m1 models
  setwd(paste0(work.dir, "/", gene.name))
  setwd(paste0("model1"))
  
  #Get stats
  lik1 = find.lnl("output.txt")
  if (nrow(lik1) == 0){ next }
  lik1 = cbind(Model = "m0", lik1)
  mlc1 = read.mlc("output.txt")
  if (length(mlc1) == 0){ next }
  if (nrow(mlc1) == 0){ next }
  
  #The m2 models
  setwd(paste0(work.dir, "/", gene.name))
  setwd(paste0("model2"))
  
  #Get stats
  lik2 = find.lnl("output.txt")
  if (nrow(lik2) == 0){ next }
  lik2 = cbind(Model = "m1", lik2)
  mlc2 = read.mlc("output.txt")
  if (length(mlc2) == 0){ next }
  if (nrow(mlc2) == 0){ next }
  
  #Likelihood ratio test
  ratio = 2*(lik2$lnl - lik1$lnl)
  df = lik2$np-lik1$np
  pval = pchisq(ratio, df=df, lower.tail=FALSE)

  #Save the data
  temp.data = data.frame(Gene = gene.name, m0_lnl = lik1$lnl, m0_np = lik1$np,
                         m1_lnl = lik2$lnl, m1_np = lik2$np, m0_dnds = as.numeric(mlc1[1,6]),
                         m1_dnds = as.numeric(mlc2[1,6]), df = df, pval = pval)

  if (is.nan(temp.data$m0_dnds) == T){ next }
  if (is.nan(temp.data$m1_dnds) == T){ next }
  if (temp.data$m1_dnds == 0){ selection = "neutral" }
  if (temp.data$m1_dnds > 1){ selection = "positive" }
  if (temp.data$m1_dnds < 1){ selection = "negative" }

  beb2 = find.BEB("output.txt")
  
  if (is.null(beb2) == T){ 
    n.beb = 0
    s.beb = 0
  } else{ 
    n.beb = nrow(beb2)
    temp.beb = beb2[beb2$Post >= 0.95,]
    s.beb = nrow(temp.beb)
    #Aggregate beb significane data
    beb.data = cbind(Gene = gene.name, beb2) 
    save.beb = rbind(save.beb, beb.data)
  }#end if
  
  #test for beb significance
  temp.data = cbind(temp.data, Selection = selection, nBEB = n.beb, sigBEB = s.beb)
  save.data = rbind(save.data, temp.data)
  system(paste0("rm -r ", work.dir, "/", gene.name))
  
  #Read in results data
  #rst.data = read.paml_rst("rst")
}#end i loop
  

setwd(work.dir)
setwd("..")
write.data = save.data
write.csv(write.data, file = "Site-selection_results.csv", row.names = F)

write.data = save.beb
write.csv(write.data, file = "Site-selection_results_BEB.csv", row.names = F)

#Write go-termed list of positive selected genes only 
#Make blast database for the probe loci


# system(paste0("makeblastdb -in GCF_000001635.26_GRCm38.p6_protein.faa",
#               " -parse_seqids -dbtype prot -out gene_nucl-blast_db"))
# 
# #Matches samples to proteins
# system(paste0("blastp -db gene_nucl-blast_db -evalue 0.001 -seg no",
#               " -query protein_mouse_genes.fa -out prot-prot-match.txt", 
#               " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
#               " -num_threads ", threads))
# 
# #headers for the blast db
# headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
#            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")


genome.dir = "/Volumes/Rodents/Genomes/Mus_Genome"
setwd(genome.dir)

ensembl.data = fread("ensembl_data.txt")
genome.temp = read.gff("Mus_GRCm38.p6_genomic.gff")
genome.data = data.table(genome.temp)
genome.data[, MatchID:=as.character("0")]

#genome.data = genome.data[genome.data$type == "gene",]
genome.data$MatchID = gsub(".*Dbxref=", "", genome.data$attributes)
genome.data$MatchID = gsub("Name=.*", "", genome.data$MatchID)
genome.data$MatchID = gsub("RFAM=.*", "", genome.data$MatchID)
genome.data$MatchID = gsub("Note=.*", "", genome.data$MatchID)
genome.data$MatchID = gsub("\\..*", "", genome.data$MatchID)
genome.data = genome.data[genome.data$type == "CDS",]

#Read in data to match 
sel.results = fread("/Volumes/Rodents/Selection-test/Site-models/Site-selection_results.csv")
sel.results = sel.results[sel.results$pval <= 0.001,]
sel.results = sel.results[sel.results$sigBEB >= 1,]
sel.results = sel.results[sel.results$Selection == "positive",]
sel.results[, Ensembl:=as.character("0")]

gene.names = gsub("Mm.*", "", sel.results$Gene)
gene.names = gsub(".{1}$", "", gene.names)

for (i in 1:length(gene.names)){
  
  gene.line = genome.data[grep(gene.names[i], genome.data$MatchID),]
  
  if (nrow(gene.line) == 0){ next }
  
  temp.id = gsub(".*GeneID:", "", gene.line$attributes[1])
  gene.id = gsub(",.*", "", temp.id)
  
  #Gets ensembl data
  temp.ensembl = ensembl.data[ensembl.data$`NCBI gene ID` == gene.id,]
  if (nrow(temp.ensembl) == 0){ next }
  
  sel.results$Ensembl[i] = unique(temp.ensembl$`Gene stable ID`)[1]
}#end i loop


setwd("/Volumes/Rodents/Selection-test/Site-models/")
write.data = sel.results
write.csv(write.data, file = "Site-selection_results_Ensembl.csv", row.names = F)

sel.write = sel.results[sel.results$Ensembl != "0",]
sel.write = sel.write[sel.write$pval <= 0.001,]
sel.write = sel.write[sel.write$sigBEB >= 1,]
sel.write = sel.write[sel.write$selection == "positive",]

fileConn = file("Ensemble_codes_filtered.txt")
writeLines(sel.write$Ensembl, fileConn)
close(fileConn)




####### END SCRIPT
###########################
