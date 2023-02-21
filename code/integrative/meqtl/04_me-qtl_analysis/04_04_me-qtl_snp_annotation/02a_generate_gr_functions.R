library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(IRanges)
library(arules)
library(reshape2)

#' Function to generate an eQTL GRanges object (with MAF bins)
#'
#' @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
#' @param mdata Merging data (e.g. Background SNPs)
#' @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObjecteQTL <- function(input, mdata, ofile){
  
  # merge with MAF info
  input <- merge(input, mdata, by = "SNP")
  
  # generate GRanges object
  input_gr <- GenomicRanges::GRanges(seqnames = input$CHR, 
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)), 
                                                               end = as.numeric(as.character(input$POS))),
                                     snp_id = input$SNP, 
                                     bin = input$bin)
  # remove duplicates
  input_gr <- input_gr[!duplicated(input_gr)] 
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "UCSC"
  
  # save GRanges object to an rds file
  saveRDS(input_gr, file =  ofile)
  
  # return finalized GRanges object
  return(input_gr)
}

#' Function to generate a GRanges object (for all datasets with Chromosome, base pair position and rs SNP ID information)
#'
#' @param input Input data of which an GRanges object should be made of (needs CHR, BP, P and SNP as column names)
#' @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObject <- function(input, ofile){
  
  # filter for p-value (P) < 0.05
  input <- input[input$P < 0.05, ]
  
  # remove all rows including NA in any column
  input <- input[complete.cases(input), ]
  
  # generate GRanges object
  input_gr <- GenomicRanges::GRanges(seqnames = input$CHR,
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)),
                                                               end = as.numeric(as.character(input$POS))),
                                     snp_id = input$SNP,
                                     p_value = input$P)  # generate GRanges object only with chr, position, rs SNP ID and p-value
  # remove duplicates
  input_gr <- input_gr[!duplicated(input_gr)] 
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "UCSC"
  
  # save GRanges object to an rds file
  saveRDS(input_gr, file = ofile)
  
  return(input_gr)
}