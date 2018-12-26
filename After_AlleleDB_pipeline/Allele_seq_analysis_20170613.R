args=(commandArgs(TRUE))

setwd(args[1])


setwd("/Users/shaopei/Box Sync/Danko_lab_work/Transcrition_directionality_project/fdr_Test/")

count_plus_beta = read.table("counts_plus_beta_mat.txt",header=T,sep = "\t")
count_minus_beta = read.table("counts_minus_beta_mat.txt",header=T,sep = "\t")
count_SSCombined_beta = rbind.data.frame(count_plus_beta, count_minus_beta)

count_plus_bin = read.table("counts_plus_bin.txt",header=T,sep = "\t")
count_minus_bin = read.table("counts_minus_bin.txt",header=T,sep = "\t")
count_plus_bin$strand='plus'
count_minus_bin$strand='minus'
count_SSCombined_bin = rbind.data.frame(count_plus_bin, count_minus_bin)

  
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
count_SSCombined_beta$mat_allele_count <- get.mat_allele_count(count_SSCombined_beta)
count_SSCombined_beta$pat_allele_count <- get.pat_allele_count(count_SSCombined_beta)
count_SSCombined_beta$total.reads.count <-count_SSCombined_beta$mat_allele_count + count_SSCombined_beta$pat_allele_count

count_SSCombined_bin$mat_allele_count <- get.mat_allele_count(count_SSCombined_bin)
count_SSCombined_bin$pat_allele_count <- get.pat_allele_count(count_SSCombined_bin)
count_SSCombined_bin$total.reads.count <-count_SSCombined_bin$mat_allele_count + count_SSCombined_bin$pat_allele_count


count <- count_SSCombined_bin[count_SSCombined_bin$SymCls !='Weird', ]
#count$strand <- count_SSCombined_bin$strand[count_SSCombined_bin$SymCls !='Weird']
library('VGAM')
fit_NoW.ab <- vglm(formula = cbind(count$pat_allele_count,count$mat_allele_count) ~ 1, family = "betabinomialff")
Coef(fit_NoW.ab)

fit_NoW.rho <- vglm(formula = cbind(count$pat_allele_count,count$mat_allele_count) ~ 1, family = "betabinomial")
Coef(fit_NoW.rho)
#mu       rho 
#0.5011688 0.2409120 



fit_Asym.ab <- vglm(formula = cbind(count$pat_allele_count[count$SymCls =='Asym'],count$mat_allele_count[count$SymCls =='Asym']) ~ 1, family = "betabinomialff")
Coef(fit_Asym.ab)
fit_Sym.ab <- vglm(formula = cbind(count$pat_allele_count[count$SymCls =='Sym'],count$mat_allele_count[count$SymCls =='Sym']) ~ 1, family = "betabinomialff")
Coef(fit_Sym.ab)
fit_Sym.rho <- vglm(formula = cbind(count$pat_allele_count[count$SymCls =='Sym'],count$mat_allele_count[count$SymCls =='Sym']) ~ 1, family = "betabinomial")
Coef(fit_Sym.rho)
#mu        rho 
#0.50032003 0.01722301 


#library('ibb')
#bb?binom.test Pvalue <- bb.test(count$mat_allele_count, count$total.reads.count,count$SymCls)
#bbtestPvalue <- bb.test(count$mat_allele_count, count$total.reads.count,count$SymCls,n.threads = 8)

hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), ylim=c(0,50000), col='red')
hist(count$mat_allele_count[count$SymPval.beta.SymModel>0.02]/count$total.reads.count[count$SymPval.beta.SymModel>0.02], breaks=seq(0,1,0.05), ylim=c(0,30000), col='dark grey', add=T)
#hist(count$mat_allele_count[count$SymPval.beta.SymModel.rho>0.05]/count$total.reads.count[count$SymPval.beta.SymModel.rho>0.05], breaks=seq(0,1,0.05), ylim=c(0,30000), col='dark grey', add=T)
empirical.reads.ratio.hist = hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), plot = F)
s1_all = 1.571751 
s2_all = 1.579117
s1_Asym = 0.3774592 
s2_Asym = 0.3816091
s1_Sym = 28.51298 
s2_Sym = 28.54951
Rho_Sym = 0.01722301 



count$expected.reads.count.all <- rbetabinom.ab(rep(1,length(count$total.reads.count)), count$total.reads.count, shape1 = s1_all, shape2 = s2_all)
expected.reads.ratio.hist <- hist(count$expected.reads.count.all/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark green')
sum((empirical.reads.ratio.hist$counts - expected.reads.ratio.hist$counts)^2)
#331,938,240
count$expected.reads.count.Asym <- count$expected.reads.count.all
count$expected.reads.count.Asym[count$SymCls =='Asym'] <-  rbetabinom.ab(rep(1,length(count$total.reads.count[count$SymCls =='Asym'])), count$total.reads.count[count$SymCls =='Asym'], shape1 = s1_Asym, shape2 = s2_Asym)
#expected.reads.ratio.hist <- hist(count$expected.reads.count.Asym/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
#lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark blue')
count$expected.reads.count.Asym[count$SymCls =='Sym'] <-  rbetabinom.ab(rep(1,length(count$total.reads.count[count$SymCls =='Sym'])), count$total.reads.count[count$SymCls =='Sym'], shape1 = s1_Sym, shape2 = s2_Sym)
expected.reads.ratio.hist <- hist(count$expected.reads.count.Asym/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark red')
sum((empirical.reads.ratio.hist$counts - expected.reads.ratio.hist$counts)^2)
#251,776,786
count$expected.reads.count.bin <-  rbetabinom(rep(1,length(count$total.reads.count)), count$total.reads.count, 0.5, rho = 0)
expected.reads.ratio.hist <- hist(count$expected.reads.count.bin/count$total.reads.count, breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='black')


legend("topleft", legend = c('Binomial','Betabinomial(one model)','Betabinomial(two models)'), col=c('black','dark green','dark red'), pch=1) # optional legend

count$SymPval.beta.SymModel = betabinom.test.ab(count$mat_allele_count, count$total.reads.count, s1_Sym, s2_Sym)
#count$SymPval.beta.SymModel.pat = betabinom.test.ab(count$pat_allele_count, count$total.reads.count, s1_Sym, s2_Sym)
#count$SymPval.beta.SymModel.rho = betabinom.test.rho(count$mat_allele_count, count$total.reads.count,Rho_Sym)

count$SymPval.beta.AllModel = betabinom.test.ab(count$mat_allele_count, count$total.reads.count, s1_all, s2_all)
#count$SymPval.beta.AllModel.pat = betabinom.test.ab(count$pat_allele_count, count$total.reads.count, s1_all, s2_all)

N=500
plot(0:N, dbetabinom.ab(0:N,N,s1_all,s2_all))
plot(0:N, dbetabinom.ab(0:N,N,s1_Asym, s2_Asym))
plot(0:N, dbetabinom.ab(0:N,N,s1_Sym, s2_Sym))
hist(rbetabinom.ab(0:N,N,s1_all,s2_all)/N)
hist(rbetabinom.ab(0:N,N,s1_Asym,s2_Asym)/N)
hist(rbetabinom.ab(0:N,N,s1_Sym,s2_Sym)/N)

hist(count$mat_allele_count[count$SymCls=='Asym']/count$total.reads.count[count$SymCls=='Asym'], breaks=seq(0,1,0.05), ylim=c(0,30000), col='red')
expected.reads.ratio.hist <- hist(count$expected.reads.count.Asym[count$SymCls=='Asym']/count$total.reads.count[count$SymCls=='Asym'], breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark red')

hist(count$mat_allele_count[count$SymCls=='Sym']/count$total.reads.count[count$SymCls=='Sym'], breaks=seq(0,1,0.05), ylim=c(0,30000), col='blue', add=T)
expected.reads.ratio.hist <- hist(count$expected.reads.count.Asym[count$SymCls=='Sym']/count$total.reads.count[count$SymCls=='Sym'], breaks=seq(0,1,0.05), plot=F)
lines(expected.reads.ratio.hist$mids, expected.reads.ratio.hist$counts, type='o', col='dark blue')

par(mfrow=c(1,1))
hist(count$SymPval)
hist(count$SymPval.beta.SymModel)
hist(count$SymPval.beta.AllModel)
#### FDR
#sim_n = 5
#n = dim(count)[1]
#act_pvals=count$SymPval.beta.SymModel.rho # pval as reported in counts file
#cnt_sums=count$total.reads.count
getFDR.binomial<-function(cnt_sums, act_pvals, sim_n){
    n = length(act_pvals)
    act_pvals = sort(act_pvals)
    sim_pvals_m=c()
    for (j in 1 :sim_n){
      sim_pvals=c()
      sim_counts = rbetabinom(rep(1,length(cnt_sums)), cnt_sums, 0.5, rho = 0)
      for (i in 1 :length(sim_counts)){
        sim_pvals=c(sim_pvals, binom.test(sim_counts[i], cnt_sums[i], p=0.5)$p.value)
      }
      sim_pvals_m=cbind(sim_pvals_m, sort(sim_pvals))
    }
    pvs = c(seq(0,0.009,0.001), seq(0.01, 0.09, 0.01), seq(0.1, 0.9, 0.1))
    for (pv in pvs){
      Nact = sum(act_pvals<=pv)
      mean_Nsims =  sum(sim_pvals_m <= pv) / sim_n
      FDR=mean_Nsims/(Nact+1)
      cat(paste(pv, Nact, mean_Nsims, FDR, sep = '\t'))
      cat('\n')
    }
}

getFDR.betabinomial<-function(cnt_sums, act_pvals, sim_n, shape1, shape2){
  n = length(act_pvals)
  act_pvals = sort(act_pvals)
  sim_pvals_m=c()
  for (j in 1 :sim_n){
    sim_pvals=c()
    sim_counts = rbetabinom.ab(rep(1,length(cnt_sums)), cnt_sums, shape1, shape2)
    for (i in 1 :length(sim_counts)){
      sim_pvals=c(sim_pvals, betabinom.test.ab(sim_counts[i], cnt_sums[i], shape1, shape2))
    }
    sim_pvals_m=cbind(sim_pvals_m, sort(sim_pvals))
  }
  pvs = c(seq(0,0.009,0.001), seq(0.01, 0.09, 0.01), seq(0.1, 0.9, 0.1))
  for (pv in pvs){
    Nact = sum(act_pvals<=pv)
    mean_Nsims =  sum(sim_pvals_m <= pv) / sim_n
    FDR=mean_Nsims/(Nact+1)
    cat(paste(pv, Nact, mean_Nsims, FDR, sep = '\t'))
    cat('\n')
  }
}

#act_pvals=count$SymPval.beta.SymModel.rho # pval as reported in counts file
#cnt_sums=count$total.reads.count
getFDR.binomial(count$total.reads.count, count$SymPval, 5)
getFDR.betabinomial(count$total.reads.count, count$SymPval.beta.SymModel.rho, 5, s1_Sym, s2_Sym)

####Need to use function ###
rowMins <- function(input_matrix){
  rowMin=c()
  for (i in 1 :dim(input_matrix)[1]){
    rowMin=c(rowMin, min(input_matrix[i,]))
  }
  return(rowMin)
}
betabinom.test.ab <- function(hap_count, total_count, shape1, shape2){
  temp = pbetabinom.ab (hap_count, total_count, shape1,shape2)
  temp = rowMins(cbind(1,temp))
  temp_min = rowMins(cbind(temp, 1-temp))
  return(2*temp_min)
}
betabinom.test.rho <- function(hap_count, total_count, rho){
  temp = pbetabinom (hap_count, total_count, 0.5, rho)
  temp = rowMins(cbind(1,temp))
  temp_min = rowMins(cbind(temp, 1-temp))
  return(2*temp_min)
}

##### temp

count_a1 <- count_plus_beta_a1
count_a1$mat_allele_count <- get.mat_allele_count(count_a1)
count_a1$pat_allele_count <- get.pat_allele_count(count_a1)
count_a1$total.reads.count <- count_a1$mat_allele_count + count_a1$pat_allele_count

library('VGAM')
fit.ab <- vglm(formula = cbind(count$mat_allele_count,count$total.reads.count-count$mat_allele_count) ~ 1, family = "betabinomialff")
Coef(fit.ab)
#shape1   shape2 
#1.572162 1.583614
fit.plus.ab <- vglm(formula = cbind(count_plus_beta$mat_allele_count,count_plus_beta$total.reads.count-count_plus_beta$mat_allele_count) ~ 1, family = "betabinomialff")
Coef(fit.plus.ab)




shape1=2.019315
shape2=2.031387
temp = pbetabinom.ab (count_a1$pat_allele_count, count_a1$total.reads.count, shape1,shape2)
temp_min = rowMins(cbind(temp, 1-temp))
temp_min= 2*temp_min
count_a1$SymPval.R.pat <- temp_min
temp = pbetabinom.ab (count_a1$mat_allele_count, count_a1$total.reads.count, shape1,shape2)
temp_min = rowMins(cbind(temp, 1-temp))
temp_min= 2*temp_min
count_a1$SymPval.R.mat <- temp_min


##
hist(count$mat_allele_count/count$total.reads.count, breaks=seq(0,1,0.05), ylim=c(0,20000), col='dark grey')
## tempt codes ####
#count[count$SymPval.R <=0.05,][count$SymPval >=0.05,]
count_a1$SymPval.smallerCount <- count_a1$SymPval
View(count_a1)
disgreeSet <- count_a1[count_a1$SymPval.R.mat <=0.05,]
disgreeSet <-disgreeSet[disgreeSet$SymCls != 'Weird',]
disgreeSet <-disgreeSet[disgreeSet$SymPval.R.pat >=0.05,]

disgreeSet$difference <- abs(disgreeSet$SymPval.R.pat - disgreeSet$SymPval.smallerCount)
disgreeSet <-disgreeSet[disgreeSet$SymPval.R.pat >=0.05,]

disgreeSet <-disgreeSet[disgreeSet$SymPval >=0.05,]
disgreeSet$difference <- abs(disgreeSet$SymPval.R.pat - disgreeSet$SymPval.smallerCount)

View(count_noW)
count_noW <- count_a1[count_a1$SymCls !=  'Weird',]
count_noW <- count_noW[count_noW$total.reads.count >= 12,]
plot(count_noW$SymPval.R.pat, count_noW$SymPval.R.mat)
plot(count_noW$SymPval, count_noW$SymPval.R.mat)
disgreeSet <- count_noW[count_noW$SymPval <= 0.05,]
disgreeSet <-disgreeSet[disgreeSet$SymPval.R.mat >=0.05,]