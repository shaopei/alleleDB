# R --vanilla --slave --args $(pwd) ../LEP_ZYG_ATGCA_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed < rtfbsdb_scan.R 

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2]




## Attach the package firstly
library(rtfbsdb);

##Create db from pre-installed dataset
db<- CisBP.extdata("Mus_musculus");
#  Zip file for Mus_musculus is downloaded on 2015_12_08

## Select all motifs from CisBP dataset
tfs <- tfbs.createFromCisBP(db);
#   1823 PWM(s) are defined in the CisBP dataset.
# ! 899 PWM file(s) are failed to be loaded ( Missing PWMs : 717 , Empty PWMs : 182 ).
# * 924 PWM(s) in the tfbs object.
# Warning message:
# In unzip(cisbp.db@zip.file, c(pwm.file), exdir = tmp.dir) :
#   requested file not found in the zip file

show(tfs)
# Species:  Mus_musculus
# TF number:  924
# Distance Matrix:  NULL
# Cluster Matrix:  NULL
# Expression:  NULL

file.twoBit <- "/local/storage/data/mm10/mm10.2bit"

# Example 1: Scan the whole genome
# Scan 2bit file within whole genome to find motif binding sites
#r1.scan <- tfbs.scanTFsite( tfs, file.twoBit, ncores = 7);

# Example 2: Scan a specified range
# Get a data frame from a plain-text bed file for your range of interest
Con_Dis_bed <- read.table(input_f, header=F);

# Scan 2bit file within all bed regions to find motif binding sites
r2.scan <- tfbs.scanTFsite( tfs, file.twoBit, Con_Dis_bed[,1:3], return.type="writedb", threshold=7, ncores = 20);
r3.scan <- tfbs.scanTFsite( tfs, file.twoBit, Con_Dis_bed[,1:3], threshold=7, ncores = 20);

#output
write.table(r3.scan$summary, "tfbs.scanTFsite.summary", sep="\t",quote = F)


