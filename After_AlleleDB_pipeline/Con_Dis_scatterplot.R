# R --vanilla --slave --args $(pwd) GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1.txt < Con_Dis_scatterplot.R 

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2] #GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1.txt
dummy = 1

# code start
#setwd('/Users/shaopei/Box Sync/temp')
Strand_count = read.table(input_f,header=T,sep = "\t")
original_col_number = dim(Strand_count)[2]
summary(Strand_count)
Con_col_1 = 'dark green'#'springgreen4'
#Con_col_2 = 'blue'
Dis_col_1 ='dark orange'
#Dis_col_2 = 'dark red'
filter_out_col = 'gray77'

Strand_count$log_mat_over_pat_plus = with(Strand_count,log2((mat_count_plus+dummy)/(pat_count_plus+dummy)))
Strand_count$log_mat_over_pat_minus = with(Strand_count, log2((mat_count_minus+dummy)/(pat_count_minus+dummy)))
# the ratio between reads count on minus and plus strand
Strand_count$log2_minus_over_plus = log2(Strand_count$log_mat_over_pat_minus / Strand_count$log_mat_over_pat_plus)

# Stats about the number of regions
# Concordant regions are filtered by log2_minus_over_plus
d = sum(Strand_count$Groseq.Discordance=='Discordant')
c = sum(Strand_count$Groseq.Discordance=='Concordant')
f = sum((Strand_count$log2_minus_over_plus[Strand_count$Groseq.Discordance == 'Concordant'] < -1.5) | (Strand_count$log2_minus_over_plus[Strand_count$Groseq.Discordance == 'Concordant'] > 1.5))
cat('# of Concordant region before filter = ', c,'\n')
cat('# of Concordant region after filter = ', c-f,'\n')
cat('# of Discordant region = ', d,'\n')

# hist of log2_minus_over_plus
pdf(paste('Concordant','dummy',dummy,'minus_over_plus_hist.pdf', sep='_'))
#hist(Strand_count$log2_minus_over_plus[sub])
hist(Strand_count$log2_minus_over_plus[Strand_count$Groseq.Discordance=='Concordant'])
dev.off()

# scatter plot 
pdf(paste(input_f,'dummy',dummy,'scatterPlot.pdf', sep='_'))
with(Strand_count, plot(log_mat_over_pat_plus, log_mat_over_pat_minus, 
                             col=ifelse(Groseq.Discordance == 'Concordant', Con_col_1 , Dis_col_1 ),
                             xlim=c(-6,6),ylim=c(-6,6),
                             #xlim=c(-600,600), ylim=c(-600,600), 
                             pch=20,frame.plot=FALSE, main = paste('dummy=',dummy)))
abline(v=0)
abline(h=0)
abline(a=0, b=2^(1.5), col=Con_col_1,lty=2)
abline(a=0, b=1/2^(1.5), col=Con_col_1, lty=2)
# color the filter out points
sub = (Strand_count$log2_minus_over_plus < -1.5) | (Strand_count$log2_minus_over_plus > 1.5) & (Strand_count$Groseq.Discordance == 'Concordant')
with(Strand_count, points(log_mat_over_pat_plus[sub], log_mat_over_pat_minus[sub],  
                               col= filter_out_col, pch=20))
legend("topleft", legend = c(paste("Concordant (n=",c-f,")",sep=""), paste("Discordant  (n=",d,")",sep="")), col = c(Con_col_1, Dis_col_1),
       pch=20)
dev.off()

# export the input file with region pass filter
new = Strand_count[ (Strand_count$Groseq.Discordance == 'Discordant') | (Strand_count$log2_minus_over_plus >= -1.5) & (Strand_count$log2_minus_over_plus <= 1.5),1:original_col_number]
colnames(new)[1] = paste('#',colnames(new)[1],sep='')
#View(new)
write.table(new,paste(substr(input_f, 1, nchar(input_f)-4),'_ConFiltered.bed',sep=''), row.names=FALSE, sep="\t", quote = FALSE)
