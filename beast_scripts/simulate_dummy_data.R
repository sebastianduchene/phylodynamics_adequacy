library(phangorn)

tr <- rtree(100)

tr$tip.label <- paste0('s_', 1:100, '_', round(runif(100, 0, 5), 2))

plot(tr)

s <- as.DNAbin(simSeq(tr))

write.dna(s, file = 'dummy_sequences.fasta', format = 'fasta', nbcol = -1, colsep = '')

