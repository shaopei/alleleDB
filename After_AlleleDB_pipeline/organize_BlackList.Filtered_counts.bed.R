print('R --vanilla --slave --args $(pwd) < script.R ')
# combine multiple *BlackList.Filtered_counts.bed into one SymCals files by "chrm" and "snppos"
# output is "d.SymCals.all.txt"

# setwd to where the TF_counts.txt locate
#setwd("~/Desktop/part_D/Danko_lab_work/Transcrition_directionality_project/ENCODE/batch1_interestingHets")
args=(commandArgs(TRUE))
setwd(args[1])

merge_file<-function(file_list, out.name){
 for (file in file_list){
   print (file)
  # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=TRUE, sep="\t")
    }
  
    # if the merged dataset does exist, append to it
    else {
      x = strsplit(file, "_")[[1]][1]
      temp_dataset <-read.table(file, header=TRUE, sep="\t")
      head(temp_dataset)
      dataset <- merge(dataset, temp_dataset, all=T, by = c("chrm", "snppos"),suffixes = c("",paste(".",x, sep="")))
      rm(temp_dataset)
    }  
  head(dataset)
  core.col = names(dataset)[names(dataset) %in% c("chrm", "snppos")]
  SymCls.col = names(dataset)[grep (glob2rx("*.SymCls"),names(dataset))]
  #SymPval.col = names(dataset)[grep (glob2rx("*.SymPval"),names(dataset))]
  #winning.col = names(dataset)[grep (glob2rx("*.winning"),names(dataset))]
  dataset = dataset [ , c(core.col,SymCls.col)]   
  #var.out = grep ( "^.*\\..*\\.d\\.SymCals\\..*\\.txt$",names(dataset))
  #print(var.out) 
  #if (length(var.out) > 0) {dataset_sub = subset(dataset, select = - var.out)
  #}else {dataset_sub = dataset}
 }
  write.table(dataset, out.name, sep="\t", row.names =F, quote=F)
  return (dataset)
}

file_list_c <- list.files(pattern=glob2rx("*BlackList.Filtered_counts.bed"))
print (file_list_c)
print (length(file_list_c))
for(i in 0:(floor(length(file_list_c)/10)-1)) {
  #for (j in seq(1+10*i, 10+10*i, by=1)){print (j)}
  merge_file(file_list_c[seq(1+10*i, 10+10*i, by=1)], paste("d.SymCals.",i,".txt", sep=""))
}
merge_file(file_list_c[seq(floor(length(file_list_c)/10)*10+1, length(file_list_c), by =1)], paste("d.SymCals.",floor(length(file_list_c)/10),".txt", sep=""))


d.SymCals.f_list = list.files(pattern="^d.SymCals.")    
d.SymCals.all = merge_file(d.SymCals.f_list, "d.SymCals.all.txt")