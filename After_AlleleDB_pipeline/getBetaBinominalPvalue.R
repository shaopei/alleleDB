# R --vanilla --slave --args $(pwd) counts.txt < test.R 
# R --vanilla --slave --args $(pwd) counts.txt_plus counts_minus.txt < test.R 

args=(commandArgs(TRUE))
setwd(args[1])

if(length(args)==3){
  count_plus_bin = read.table(args[2], header=T, sep = "\t")
  count_minus_bin = read.table(args[3], header=T, sep = "\t")
  count_plus_bin$strand='plus'
  count_minus_bin$strand='minus'
  count_SSCombined_bin = rbind.data.frame(count_plus_bin, count_minus_bin)
}else{
  count_SSCombined_bin = read.table(args[2],header=T,sep = "\t")
}

####functions
get.mat_allele_count <- function(count){
  mat_allele_count = apply(count, 1,function(x){
    if(x['mat_all']=='T') return (as.numeric(x['cT']));
    if(x['mat_all']=='C') return (as.numeric(x['cC']));
    if(x['mat_all']=='G') return (as.numeric(x['cG']));
    if(x['mat_all']=='A') return (as.numeric(x['cA']))
    else return (NA)})
  return (mat_allele_count)
}
get.pat_allele_count <- function(count){
  pat_allele_count = apply(count, 1,function(x){
    if(x['pat_all']=='T') return (as.numeric(x['cT']));
    if(x['pat_all']=='C') return (as.numeric(x['cC']));
    if(x['pat_all']=='G') return (as.numeric(x['cG']));
    if(x['pat_all']=='A') return (as.numeric(x['cA']))
    else return (NA)})
  return (pat_allele_count)
}
rowMins <- function(input_matrix){
  rowMin=c()
  for (i in 1 :dim(input_matrix)[1]){
    rowMin=c(rowMin, min(input_matrix[i,]))
  }
  return(rowMin)
}

betabinom.test.rho <- function(hap_count, total_count, rho){
  temp = pbetabinom (hap_count, total_count, 0.5, rho)
  temp = rowMins(cbind(1,temp))
  temp_min = rowMins(cbind(temp, 1-temp))
  return(2*temp_min)
}
###main
count <- count_SSCombined_bin[count_SSCombined_bin$SymCls !='Weird', ]
count$mat_allele_count <- get.mat_allele_count(count)
count$pat_allele_count <- get.pat_allele_count(count)
count$total.reads.count <-count$mat_allele_count + count$pat_allele_count



library('VGAM')
fit_all.rho <- vglm(formula = cbind(count$pat_allele_count,count$mat_allele_count) ~ 1, family = "betabinomial")
cat("use all snps\n")
Coef(fit_all.rho)
fit_Asym.rho <- vglm(formula = cbind(count$pat_allele_count[count$SymCls =='Asym'],count$mat_allele_count[count$SymCls =='Asym']) ~ 1, family = "betabinomial")
cat("use Asym snps\n")
Coef(fit_Asym.rho)
fit_Sym.rho <- vglm(formula = cbind(count$pat_allele_count[count$SymCls =='Sym'],count$mat_allele_count[count$SymCls =='Sym']) ~ 1, family = "betabinomial")
cat("use Sym snps\n")
Coef(fit_Sym.rho)

## make plot
# parameters for plot 
Rho_all = Coef(fit_all.rho)[2]
Rho_Asym = Coef(fit_Asym.rho)[2]
Rho_Sym = Coef(fit_Sym.rho)[2]
#cat("Rho_all =", Rho_all, "\n")
#cat("Rho_Asym =", Rho_Asym, "\n")
#cat("Rho_Sym =", Rho_Sym, "\n")


dir.create("betabinomial")
setwd(paste(getwd(),"betabinomial",sep='/'))

##plot1
pdf('counts_hist_and_simmulated_counts.pdf',width=17, height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
# observed data
empirical.reads.ratio.hist = hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), plot = F)
hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), ylim=c(0,2*max(empirical.reads.ratio.hist$counts)), col='grey')

# simulated data using Rho_all
count$expected.reads.count.all <- rbetabinom(rep(1,length(count$total.reads.count)), count$total.reads.count, prob = 0.5,rho = Rho_all)
expected.reads.ratio.hist <- hist(count$expected.reads.count.all/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark green')
cat("SSE to simulated data using Rho_all\n")
sum((empirical.reads.ratio.hist$counts - expected.reads.ratio.hist$counts)^2)

# simulated data using Rho_Asym and Rho_Sym
count$expected.reads.count.SymAsym <- count$expected.reads.count.all
count$expected.reads.count.SymAsym[count$SymCls =='Asym'] <-  rbetabinom(rep(1,length(count$total.reads.count[count$SymCls =='Asym'])), count$total.reads.count[count$SymCls =='Asym'],prob = 0.5,rho = Rho_Asym)
count$expected.reads.count.SymAsym[count$SymCls =='Sym'] <-  rbetabinom(rep(1,length(count$total.reads.count[count$SymCls =='Sym'])), count$total.reads.count[count$SymCls =='Sym'],prob = 0.5,rho = Rho_Sym)
expected.reads.ratio.hist <- hist(count$expected.reads.count.SymAsym/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark red')
cat("SSE to simulated data using Rho_Asym and Rho_Sym\n")
sum((empirical.reads.ratio.hist$counts - expected.reads.ratio.hist$counts)^2)
#251,776,786

# simulated data using Binomial 
count$expected.reads.count.bin <-  rbetabinom(rep(1,length(count$total.reads.count)), count$total.reads.count, 0.5, rho = 0)
expected.reads.ratio.hist <- hist(count$expected.reads.count.bin/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='black')

#Legend
legend("topleft", legend = c('Binomial',paste('Betabinomial (rho=',round(Rho_all,digits = 2),')',sep=''),
                             paste('Betabinomial (Asym.rho=',round(Rho_Asym,digits = 2),', Sym.rho=',round(Rho_Sym,digits = 2),')',sep='')), 
       col=c('black','dark green','dark red'), pch=1, cex=1.5, pt.cex=2) # optional legend
dev.off()

##plot2
pdf('counts_hist_and_bin_significant.pdf',width=17, height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), col='red', main="Red is Binomial Pvalue < 0.05")
hist(count$mat_allele_count[count$SymPval>0.05]/count$total.reads.count[count$SymPval>0.05], breaks=seq(0,1,0.05), col='dark grey', add=T)
dev.off()

##plot3 
#using rho estimated from the snps classified as Sym in Binomial model 
#use the lower count of mat or pat
count$SymPval.beta.SymModel.lower = betabinom.test.rho(rowMins(cbind(count$pat_allele_count,count$mat_allele_count)), count$total.reads.count, rho = Rho_Sym )
pdf('counts_hist_and_betabin_significant_lowercount_RhoSym.pdf',width=17, height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), col='red', main="Red is Beta-Binomial Pvalue < 0.05")
hist(count$mat_allele_count[count$SymPval.beta.SymModel.lower>0.05]/count$total.reads.count[count$SymPval.beta.SymModel.lower>0.05], breaks=seq(0,1,0.05), col='dark grey', add=T)
dev.off()

##plot4
#using all snps
#use the lower count of mat or pat
count$SymPval.beta.OneModel.lower = betabinom.test.rho(rowMins(cbind(count$pat_allele_count,count$mat_allele_count)), count$total.reads.count, rho = Rho_all )
pdf('counts_hist_and_betabin_significant_lowercount_RhoAll.pdf',width=17, height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), col='red', main="Red is Beta-Binomial Pvalue < 0.05")
hist(count$mat_allele_count[count$SymPval.beta.OneModel.lower>0.05]/count$total.reads.count[count$SymPval.beta.OneModel.lower>0.05], breaks=seq(0,1,0.05), col='dark grey', add=T)
dev.off()
