# Creating mutpos files for unique mutations
# Input: mutpos file with total mutations
# Output: mutpos file with unique mutations

args <- commandArgs(TRUE)
in_file <- as.character(args[1])

mutpos <- read.table(in_file, header=TRUE) # Depending
colnames(mutpos) <- c('chr','ref','pos','depth','muts','T','C','G','A','ins','del','N')
nts <- c('T','C','G','A')
for (i in 1:nrow(mutpos)) {
  if (mutpos[i,'muts'] > 0) {
    for (n in nts) {
      if (mutpos[i,n] > 1) {
        mutpos[i,n] <- 1
      }
    }
    mutpos[i,'muts'] <- sum(mutpos[i,nts])
  }
}

out_file <- paste(substr(in_file,1,nchar(in_file)-7),'-unique.mutpos',sep='')
write.table(mutpos, sep="\t", file=out_file, row.names=FALSE, col.names=FALSE, quote=FALSE)