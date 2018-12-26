echo 'bash AfterAlleleDB_strandSpecific.sh Min_count'
# partameters
Min_count=$1
echo "Min_count         =$1" 
Max_Pvalus=1
d=150
dRegion=250

#locations of pipeline
PL=/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/After_AlleleDB_pipeline
tss_plus_f=${PL}/tss_paired_gm12878_plus.bed
tss_minus_f=${PL}/tss_paired_gm12878_minus.bed


# input files
plus_input_file=interestingHets_plus.txt
minus_input_file=interestingHets_minus.txt

# filter input files based on Min reads count and Max P-value
R --vanilla --slave --args $(pwd) ${plus_input_file} ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 
R --vanilla --slave --args $(pwd) ${minus_input_file} ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 

# identify Concordant and Discordant regions with defined d and dRegion
# d is the length to identify ASE SNPs from each TSS. For plus TSS, the beginning of TSS(+20nt?) to that +d
# dRegion (the length to keep in the Con_Dis regions), the legnth of the regions is (dRegion + distance_between_TSS + dRegion)
new_name=GM12878_GroSeq_d${d}_dRegion${dRegion}_withStrandSpecific_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt
python ${PL}/TSS_d_regions_CandD_directionality_index_20170612.py ${plus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt \
${minus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt ${tss_plus_f} ${tss_minus_f} \
${new_name} ${d} ${dRegion}
#python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/TSS_d_regions_CandD_directionality_index_20170612.py interestingHets_plusNoW.txt interestingHets_minusNoW.txt /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_plus.bed /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_minus.bed test.txt

# Con plus and minus change in the same direction, Did change in opposite direction.
# filter the Concordant region to remove the Regions with large difference between plus and minus strand (although same direction, but different in amplitude).
R --vanilla --slave --args $(pwd) ${new_name} < ${PL}/Con_Dis_scatterplot.R
# oupput is GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount${Min_count}_MaxPvalue${Max_Pvalus}_ConFiltered.bed
#GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed