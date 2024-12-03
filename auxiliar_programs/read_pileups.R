# Read pileup step of the cfChIP pipeline
# Description: Make figures of specified regions of the genome for one sample
# Output: new folder in working directory with a plot for each important region

library(Gviz)
library(rtracklayer)
library(org.Hs.eg.db)

########### USER INPUT (only change variables in this section) ##############
buffer <- 1000 # default = 500, can change with user input
args <- commandArgs( trailingOnly = TRUE )
sample_path = args[1] # path to sample bigwig
bed_path = args[2]  # bed file of important regions to plot
chrom_sizes_path = args[3] # chromosome sizes for this genome model
gen = args[4] # genome model
#############################################################################

# read in files
bed_con <- file(bed_path, open="r") # bed file of specific regions to snapshot
chrom_sizes <- read.table(chrom_sizes_path, header = FALSE, sep="\t")

gtrack <- GenomeAxisTrack() 

# go through each region in bed file
while (TRUE) {
  # get region info
  line = readLines(bed_con, n=1)
  if (identical(line, character(0))) { # no more lines
    break
  }
  tokens <- strsplit(line, "\t")
  chr <- tokens[[1]][1]
  start <- strtoi(tokens[[1]][2])
  end <- strtoi(tokens[[1]][3])
  
  # ensure buffer within chromosome bounds
  start_buffer <- start - buffer
  end_buffer <- end + buffer
  end_chrom <- chrom_sizes[which(chrom_sizes$V1 == chr), "V2"]  
  if (start_buffer < 0) {
    start_buffer <- 0
  } else if (end_buffer > end_chrom) {
    end_buffer <- end_chrom
  }
  
  # set up tracks
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  sample_dt <- DataTrack(range = sample_path, genome = gen, chromosome = chr)
  
  # NCBI RefSeq track
  ucscGenes <- UcscTrack(genome=gen, chromosome = chr, table="ncbiRefSeq", 
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack", 
                         rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
                         symbol = 'name', transcript = "name", strand = "strand", 
                         name = "RefSeq", stacking = 'pack', showID = T, 
                         geneSymbol = T, transcriptAnnotation="symbol")
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes_named <- ucscGenes
  ranges(ucscGenes_named) <- z

  # highlight track for region without buffer
  htrack <- HighlightTrack(trackList = list(gtrack, ucscGenes_named, sample_dt), 
                           start = start, width = end-start, chromosome = chr, 
                           col = "azure3", fill = "azure2")
  
  # plot region and save
  plot_name <- paste(chr,"_",start,"_",end,".pdf", sep="")
  pdf(plot_name)
  plotTracks(list(itrack,htrack), from=start_buffer,to=end_buffer, chromosome=chr,type="hist") 
  dev.off()
  
}
close(bed_con)

