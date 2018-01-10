# need to be under the folder like 6_CNV.Peak.BlackList.Filtered_counts_bed_in_GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1
fastq_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_fastq.tsv' # for downloaded fastq
#fastq_metadata_fp ='/workdir/sc2457/ENCODE/reDoTF/metadata.tsv' # for old 
#Groseq_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1.txt'

print 'fastq_metadata_fp =', fastq_metadata_fp
#print 'Groseq_fp =', Groseq_fp



import copy
Con_regions={}
Dis_regions={}
with open('Concordant.bed') as Con_f:
    for l in Con_f:
        ll = l.strip().split('\t')
        region = (ll[0],ll[1],ll[2])
        Con_regions[region] = set(['NA'])

with open('Discordant.bed') as Dis_f:
    for l in Dis_f:
        ll = l.strip().split('\t')
        region = (ll[0],ll[1],ll[2])
        Dis_regions[region] = set(['NA'])




import glob
import scipy.stats as stats

file_names = glob.glob('*cordantRegions.Filtered_counts.bed')
#examine each Concordant or Docordant region, if there is one Asym --> Asym, if only Sym --> Sym
FastqFileAcc_Symcals_combined={}
FastqFileAcc_Symcals_combined_Asym={}
for f in file_names:
    fastq_acc = f.split('.')[0]
    Con_Dis = f.split('.')[-3]
    if fastq_acc not in FastqFileAcc_Symcals_combined:
        if Con_Dis == 'ConcordantRegions':
            FastqFileAcc_Symcals_combined[fastq_acc]={Con_Dis:copy.deepcopy(Con_regions)}
        elif Con_Dis == 'DiscordantRegions':
            FastqFileAcc_Symcals_combined[fastq_acc]={Con_Dis:copy.deepcopy(Dis_regions)}
        FastqFileAcc_Symcals_combined_Asym[fastq_acc]={Con_Dis:[]}
    else:
        if Con_Dis == 'ConcordantRegions':
            FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis]=copy.deepcopy(Con_regions)
        elif Con_Dis == 'DiscordantRegions':
            FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis]=copy.deepcopy(Dis_regions)
        FastqFileAcc_Symcals_combined_Asym[fastq_acc][Con_Dis]=[]
    with open(f, 'U') as f_in:
        for l in f_in.readlines()[1:]:
            ll = l.strip().split('\t')
            Region = (ll[6], ll[7], ll[8])
            #if Region not in FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis]:
            #    FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis][Region] = set()
            FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis][Region].add(ll[4])
    for regions in FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis]:
        if 'Asym' in FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis][regions]:
            FastqFileAcc_Symcals_combined_Asym[fastq_acc][Con_Dis].append('Asym')
        elif 'Sym' in FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis][regions]:
            FastqFileAcc_Symcals_combined_Asym[fastq_acc][Con_Dis].append('Sym')
        else:
            FastqFileAcc_Symcals_combined_Asym[fastq_acc][Con_Dis].append('NA')


#compile Asym and Sym counts
FastqFileAcc_Asym_count={}
for f in FastqFileAcc_Symcals_combined_Asym:
    Dis_Asym = FastqFileAcc_Symcals_combined_Asym[f]['DiscordantRegions'].count('Asym')
    Dis_Sym = FastqFileAcc_Symcals_combined_Asym[f]['DiscordantRegions'].count('Sym')
    Dis_NA = FastqFileAcc_Symcals_combined_Asym[f]['DiscordantRegions'].count('NA')
    Con_Asym = FastqFileAcc_Symcals_combined_Asym[f]['ConcordantRegions'].count('Asym')
    Con_Sym  = FastqFileAcc_Symcals_combined_Asym[f]['ConcordantRegions'].count('Sym')
    Con_NA = FastqFileAcc_Symcals_combined_Asym[f]['ConcordantRegions'].count('NA')
    FastqFileAcc_Asym_count[f] =  [Dis_Asym, Dis_Sym,Con_Asym, Con_Sym, Dis_NA, Con_NA]


#map fastq to target
fastq_metadata = [l.strip().split('\t') for l in open(fastq_metadata_fp, 'U').readlines()]
FastqFileAcc_ExperimentAcc_dic = {}
FastqFileAcc_Target_dic = {}
ExperimentAcc_Target_dic = {}
for l in fastq_metadata[1:]:
    experiment_accession = l[3]
    file_accession = l[0]
    target = l[16].split('-')[0] # for new meta
    #target = l[15].split('-')[0]  #for old
    #print target
    FastqFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
    FastqFileAcc_Target_dic[file_accession] = target
    ExperimentAcc_Target_dic[experiment_accession] = target


# perform Fisher's excat test
with open('Fisher_excat_test.txt' ,'w') as out:
    out.write('\t'.join(str(s) for s in ['FastqFileAcc','Target', 'ExperimentAcc',
                                                 's_oddsratio', 's_pvalue', 'Dis_Asym', 'Dis_Sym','Con_Asym', 'Con_Sym','Dis_NA', 'Con_NA'])+'\n')
    for f in FastqFileAcc_Asym_count:
        Dis_Asym, Dis_Sym,Con_Asym, Con_Sym,Dis_NA, Con_NA  = FastqFileAcc_Asym_count[f]
        s_oddsratio, s_pvalue = stats.fisher_exact([[Dis_Asym, Dis_Sym],[Con_Asym, Con_Sym]])
        try:
            out.write('\t'.join(str(s) for s in [f,FastqFileAcc_Target_dic[f], FastqFileAcc_ExperimentAcc_dic[f],
                                                 s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym, Dis_NA, Con_NA])+'\n')
        except KeyError:  #use ExpAcc
            try:
                out.write('\t'.join(str(s) for s in [f,ExperimentAcc_Target_dic[f], f, s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym, Dis_NA, Con_NA])+'\n')
            except KeyError:
                out.write('\t'.join(str(s) for s in [f,f, f, s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym, Dis_NA, Con_NA])+'\n')

