require(optparse) 
library(glue)
library(tidyverse)

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

setwd(opt$outdir)

# -----------------------------------------------------------------------------
# Merge all individual .RDS files into one matrix
# -----------------------------------------------------------------------------

# Combine .RDS files in the output directory specified by user
countdf <- list.files(getwd(), pattern=glob2rx("counts_df_*.rds"), full.name=TRUE) %>%
  map_dfc(readRDS)

# Write as tsv
write.table(countdf, file='all.counts.tsv', sep='\t', quote=FALSE, row.names=FALSE)

# -----------------------------------------------------------------------------
# END report
# -----------------------------------------------------------------------------

cat(paste("Combined counts matrix written to: all.counts.tsv", "\n"), file=stdout())
