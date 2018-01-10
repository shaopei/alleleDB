# R --vanilla --slave --args $(pwd) Fisher_excat_test.txt < FisherExactTest_Motif_SNPs_vocano.R 

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2]

data = read.table(input_f, header = T,sep = "\t")
num_of_TF = dim(data)[1]

colnames(data)[2] = 'oddsratio'
colnames(data)[3] = 'pvalue'

# Make a basic volcano plot
pdf(paste(input_f,'volcanoPlot.pdf', sep='_'))
with(data, plot(log2(oddsratio), -log10(pvalue), 
                   pch=20, frame.plot=FALSE, cex=-log10(pvalue),ylim=c(0,max(ceiling(-log10(0.05/num_of_TF)),max(-log10(pvalue)))),
                   #xlim=c(-3,3)
                   #xlim=c(min(-3,min(log2(oddsratio)), na.rm = T),max(3,max(log2(oddsratio)),na.rm = T))
                   ))
abline(v=0)
abline(h=-log10(0.05), col='blue',lty=2, lwd=2)
abline(h=-log10(0.05/num_of_TF), col='red',lty=2, lwd=2)

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(data, pvalue<.1 & oddsratio > 1 ), points(log2(oddsratio), -log10(pvalue), pch=20, col='dark orange', cex=-log10(pvalue)))
#with(subset(data, pvalue<.1 & oddsratio < 1 ), points(log2(oddsratio), -log10(pvalue), pch=20, col='dark green', cex=-log10(pvalue)))

with(subset(data,  oddsratio > 1 ), points(log2(oddsratio), -log10(pvalue), pch=20, col='dark orange', cex=-log10(pvalue)))
with(subset(data,  oddsratio < 1 ), points(log2(oddsratio), -log10(pvalue), pch=20, col='dark green', cex=-log10(pvalue)))
legend("topright", legend = c("Enriched in Concordant", "Enriched in Discordant"), col = c('dark green', 'dark orange'),
       pch=20)

# Label points with the textxy function from the calibrate plot
with(subset(data, pvalue<.05), 
     text(log2(oddsratio), -log10(pvalue), labels= TF_name, cex=1, pos= 3))
dev.off()