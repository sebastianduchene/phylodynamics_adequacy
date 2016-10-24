library(phangorn)

tr <- read.tree('basel_code/bd_sim_100.trees')

s1 <- as.DNAbin(simSeq(tr[[1]], l = 100))

write.dna(s1, '~/Desktop/simulated_seq.fasta', format = 'fasta', nbcol = -1, colsep ='')


