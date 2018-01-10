echo 'bash AfterAlleleDB_strandSpecific.sh Min_count sample_name'
# bash AfterAlleleDB_strandSpecific_F1_PACHY_TCGTA.sh 5 PACHY_TCGTA
# partameters
Min_count=$1
echo "Min_count         =$1" 
sample_name=$2
echo "sample_name         =$2" 
Max_Pvalus=1
d=150
dRegion=250

#locations of pipeline
PL=/workdir/sc2457/tools/After_AlleleDB_pipeline

# input files
plus_input_file=interestingHets_plus.txt
minus_input_file=interestingHets_minus.txt

# locations of dREG results
dREG_dir=/workdir/sc2457/mouse_AlleleSpecific/${sample_name}.dREG
zcat  ${dREG_dir}/out.dREG.peak.full.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $7, $3, $4,$5, "+"}' > tss_paired_${sample_name}_dREG_plus.bed
zcat  ${dREG_dir}/out.dREG.peak.full.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $2, $7, $4,$5, "-"}' > tss_paired_${sample_name}_dREG_minus.bed

tss_plus_f=tss_paired_${sample_name}_dREG_plus.bed
tss_minus_f=tss_paired_${sample_name}_dREG_minus.bed


# only anaysis autosome now
grep -v X ${plus_input_file} > ${plus_input_file:0:-4}_noX
grep -v X ${minus_input_file} > ${minus_input_file:0:-4}_noX



# filter input files based on Min reads count and Max P-value
R --vanilla --slave --args $(pwd) ${plus_input_file:0:-4}_noX ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 
R --vanilla --slave --args $(pwd) ${minus_input_file:0:-4}_noX ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 

###### identify Concordant and Discordant regions with defined d and dRegion ######
# d is the length to identify ASE SNPs from each TSS. For plus TSS, the beginning of TSS(+20nt?) to that +d
# dRegion (the length to keep in the Con_Dis regions), the legnth of the regions is (dRegion + distance_between_TSS + dRegion)
new_name=${sample_name}_d${d}_dRegion${dRegion}_withStrandSpecific_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt
python ${PL}/TSS_d_regions_CandD_directionality_index_20170612.py ${plus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt \
${minus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt ${tss_plus_f} ${tss_minus_f} \
${new_name} ${d} ${dRegion}
#python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/TSS_d_regions_CandD_directionality_index_20170612.py interestingHets_plusNoW.txt interestingHets_minusNoW.txt /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_plus.bed /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_minus.bed test.txt

# Con plus and minus change in the same direction, Did change in opposite direction.
# filter the Concordant region to remove the Regions with large difference between plus and minus strand (although same direction, but different in amplitude).
R --vanilla --slave --args $(pwd) ${new_name} < ${PL}/Con_Dis_scatterplot.R
# oupput is GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount${Min_count}_MaxPvalue${Max_Pvalus}_ConFiltered.bed
#GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed


##################
# Motif analysis #
##################
mkdir Motif_SNPs_Anaylsis
# identify SNPs overlap with Concordant or Discordant regions
intersectBed -wa -a /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.bed -b ${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed > Motif_SNPs_Anaylsis/${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.snp.bed  

cd Motif_SNPs_Anaylsis
grep Concordant ../${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed > Concordant.bed
grep Discordant ../${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed > Discordant.bed

# use rtfbsdb to scan the concrdant and discordant region for motif
R --vanilla --slave --args $(pwd) ../${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed < ${PL}/rtfbsdb_scan.R


for tf in `cat tfbs.scanTFsite.summary |awk '(NR>1){print $3}'|sort |uniq`
do
  # a specific TF
  for motif in  `grep $tf tfbs.scanTFsite.summary | awk '{print $2}' `
    do unstarch scan.db.db.starch |grep ${motif} >> tfbs.scanTFsite_${tf}
  done
  # identify motif sites overlap with SNPs
  # might have duplicates of binding sites due to more than one SNP is a region, but this will be taken care later in the pipeline
  intersectBed -wa -a tfbs.scanTFsite_${tf} -b ${sample_name}_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.snp.bed > tfbs.scanTFsite_${tf}_wSNP.bed
  
  # calculate Dis_SNP, Dis_Motif,Con_SNP, Con_Motif
  echo ${tf} >> TF_SNP_Motif_count
  # Dis_SNP, the number of Discordant regions with motif(s) that overlaps with SNP
  intersectBed -wa -wb -a tfbs.scanTFsite_${tf}_wSNP.bed -b Discordant.bed | awk '{print $9, $10, $11, $12}' |sort |uniq |wc -l >> TF_SNP_Motif_count
  # Dis_Motif, the number of Discordant regions with motif(s) (with or without SNP)
  intersectBed -wa -wb -a tfbs.scanTFsite_${tf} -b Discordant.bed | awk '{print $9, $10, $11, $12}' |sort |uniq |wc -l >> TF_SNP_Motif_count
  # Con_SNP, the number of Concordant regions with motif(s) that overlaps with SNP
  intersectBed -wa -wb -a tfbs.scanTFsite_${tf}_wSNP.bed -b Concordant.bed | awk '{print $9, $10, $11, $12}' |sort |uniq |wc -l >> TF_SNP_Motif_count
  # Con_Motif, the number of Concordant regions with motif(s) (with or without SNP)
  intersectBed -wa -wb -a tfbs.scanTFsite_${tf} -b Concordant.bed | awk '{print $9, $10, $11, $12}' |sort |uniq |wc -l >> TF_SNP_Motif_count
done

mkdir toremove
mv Concordant.bed Discordant.bed toremove/.
mv tfbs.scanTFsite_* toremove/.

python ${PL}/FisherExactTest_Motif_SNPs.py TF_SNP_Motif_count
R --vanilla --slave --args $(pwd) Fisher_excat_test.txt < ${PL}/FisherExactTest_Motif_SNPs_volcano.R 