#python script.py path_To_metafata 
# use Assembly=='hg19', Output_type == 'peaks'
# pair up bed files with fastq files
# use the bed files from the same experiment session
# use NarrowPeak bed if available, else use BroadPeak bed


from sys import argv
fastq_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_fastq.tsv' # for downloaded fastq
NarrowPeakbed_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/ENCODE_bed/metadata_NarrowPeakbed.tsv' #argv[1]  # for donwloaded bed
additional_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_FastqAndAllBed_20170628.tsv' #to identify what else I need to download

bed_files_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/ENCODE_bed/' 
CNVFiltered_counts_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_counts/1_CNVFiltered_counts_bed/'
fastq_bed_pair_log = 'ENCODE_fastq_bed_pair_sameExpAcc.txt'

#make dictionary from NarrowPeakbed_metadata
NP_bed_metadata = [l.strip().split('\t') for l in open(NarrowPeakbed_metadata_fp, 'U').readlines()]
ExperimentAcc_BedFileAcc_dic = {}
BedFileAcc_ExperimentAcc_dic = {}
BedFileAcc_PeakType_dic = {}
Target_BedFileAcc_dic = {}
Target_ExperimentAcc_dic={}
ExperimentAcc_Target_dic = {}
for l in NP_bed_metadata[1:]:
    experiment_accession = l[3]
    Output_type = l[2]
    file_accession = l[0]+'.bed'
    Assembly=l[42]
    target = l[16].split('-')[0]
    peaktype = l[1].split(' ')[1]
    if Assembly=='hg19':
        if Output_type.count('peaks') >0:
            if experiment_accession not in ExperimentAcc_BedFileAcc_dic:
                ExperimentAcc_BedFileAcc_dic[experiment_accession] = []
            ExperimentAcc_BedFileAcc_dic[experiment_accession].append(file_accession)
            BedFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
            if target not in Target_BedFileAcc_dic:
                Target_BedFileAcc_dic[target] = []
            Target_BedFileAcc_dic[target].append(file_accession)
            if target not in Target_ExperimentAcc_dic:
                Target_ExperimentAcc_dic[target] =set()
            Target_ExperimentAcc_dic[target].add(experiment_accession)
            ExperimentAcc_Target_dic[experiment_accession] = target
            BedFileAcc_PeakType_dic[experiment_accession] = peaktype

#make dictionary from fastq_metadata
fastq_metadata = [l.strip().split('\t') for l in open(fastq_metadata_fp, 'U').readlines()]
ExperimentAcc_FastqFileAcc_dic = {}
FastqFileAcc_ExperimentAcc_dic = {}
Target_FastqFileAcc_dic = {}
FastqFileAcc_Target_dic = {}
for l in fastq_metadata[1:]:
    experiment_accession = l[3]
    file_accession = l[0]
    target = l[16].split('-')[0]
    if experiment_accession not in ExperimentAcc_FastqFileAcc_dic:
        ExperimentAcc_FastqFileAcc_dic[experiment_accession] = []
    ExperimentAcc_FastqFileAcc_dic[experiment_accession].append(file_accession)
    FastqFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
    if target not in Target_FastqFileAcc_dic:
        Target_FastqFileAcc_dic[target] = []
    Target_FastqFileAcc_dic[target].append(file_accession)
    FastqFileAcc_Target_dic[file_accession] = target
    ExperimentAcc_Target_dic[experiment_accession] = target
    if target not in Target_ExperimentAcc_dic:
        Target_ExperimentAcc_dic[target] =set()
        Target_ExperimentAcc_dic[target].add(experiment_accession)


# use Narrow Peak bed from the same Exp Acc as the fastq file
Exp_withNoBedFile_yet=set()
Target_name_of_Exp_withNoBedFile_yet=set()
with open(fastq_bed_pair_log, 'w') as out:
    out.write('\t'.join(['FastqFileAcc','Target','ExperimentAcc','ExperimentAcc share the same Target','PeakType','BedFileAcc from the same ExperimentAcc'])+'\n') #
    for f in FastqFileAcc_ExperimentAcc_dic:
        try:
            out.write('\t'.join([f,FastqFileAcc_Target_dic[f] ,FastqFileAcc_ExperimentAcc_dic[f],
                                 ','.join(Target_ExperimentAcc_dic[FastqFileAcc_Target_dic[f]]),
                                 BedFileAcc_PeakType_dic[FastqFileAcc_ExperimentAcc_dic[f]],
                                 ','.join(ExperimentAcc_BedFileAcc_dic[FastqFileAcc_ExperimentAcc_dic[f]])#,
                                 #','.join(Target_BedFileAcc_dic[FastqFileAcc_Target_dic[f]])
                                 ]) +'\n')
            ExperimentAcc_BedFileAcc_dic[FastqFileAcc_ExperimentAcc_dic[f]]
        except KeyError:
            Exp_withNoBedFile_yet.add(FastqFileAcc_ExperimentAcc_dic[f])
            Target_name_of_Exp_withNoBedFile_yet.add(FastqFileAcc_Target_dic[f])
        

len(Exp_withNoBedFile_yet)
#11
#Exp_withNoBedFile_yet
#set([('ENCSR000BUF', 'CREB1'), ('ENCSR904YPP', 'NR3C1'), ('ENCSR000AKE', 'H3K36me3'), ('ENCSR000AOX', 'H3K9me3'), ('ENCSR769ZTN', 'GTF2F1'), ('ENCSR679FAB', 'RNF2'), ('ENCSR974OFJ', 'KLF5'), ('ENCSR000BMI', 'SRF'), ('ENCSR000BMQ', 'EGR1'), ('ENCSR900XDB', 'ZFP36'), ('ENCSR009MBP', 'HSF1')])


# for those ExpAcc without NarrowPeak bed files, download Broadbeak bed files and use those. 
additional_metadata = [l.strip().split('\t') for l in open(additional_metadata_fp, 'U').readlines()]
NeedToDownload_BedFile_metadata=[]
for l in additional_metadata[1:]:
    experiment_accession = l[3]
    Output_type = l[2]
    file_format = l[1].split(' ')[0]
    peaktype = l[1].split(' ')[-1]
    file_accession = l[0]+'.'+file_format
    Assembly=l[42]
    target = l[16].split('-')[0]
    #if target in Target_withNoBedFile_yet:
    if experiment_accession in Exp_withNoBedFile_yet:
        if Assembly=='hg19':
            if Output_type.count('peaks') >0 and file_format=='bed':
                NeedToDownload_BedFile_metadata.append(l)
                if experiment_accession not in ExperimentAcc_BedFileAcc_dic:
                    ExperimentAcc_BedFileAcc_dic[experiment_accession] = []
                ExperimentAcc_BedFileAcc_dic[experiment_accession].append(file_accession)
                BedFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
                if target not in Target_BedFileAcc_dic:
                    Target_BedFileAcc_dic[target] = []
                Target_BedFileAcc_dic[target].append(file_accession)
                if target not in Target_ExperimentAcc_dic:
                    Target_ExperimentAcc_dic[target] =set()
                Target_ExperimentAcc_dic[target].add(experiment_accession)
                ExperimentAcc_Target_dic[experiment_accession] = target
                BedFileAcc_PeakType_dic[experiment_accession] = peaktype

with open('additonal_bed_to_download.txt', 'w') as out:
    out.write('#need to manually download this\n')
    out.write('#cd '+bed_files_fp+'\n')
    out.write('#xargs -n 1 curl -O -L < additonal_bed_to_download.txt \n')
    for l in NeedToDownload_BedFile_metadata:
        out.write(l[41]+'\n')


Exp_withNoBedFile_yet=set()
Target_name_of_Exp_withNoBedFile_yet=set()
with open(fastq_bed_pair_log, 'w') as out:
    out.write('\t'.join(['FastqFileAcc','Target','ExperimentAcc','ExperimentAcc share the same Target','PeakType','BedFileAcc from the same ExperimentAcc'])+'\n') #
    for f in FastqFileAcc_ExperimentAcc_dic:
        try:
            out.write('\t'.join([f,FastqFileAcc_Target_dic[f] ,FastqFileAcc_ExperimentAcc_dic[f],
                                 ','.join(Target_ExperimentAcc_dic[FastqFileAcc_Target_dic[f]]),
                                 BedFileAcc_PeakType_dic[FastqFileAcc_ExperimentAcc_dic[f]],
                                 ','.join(ExperimentAcc_BedFileAcc_dic[FastqFileAcc_ExperimentAcc_dic[f]])#,
                                 #','.join(Target_BedFileAcc_dic[FastqFileAcc_Target_dic[f]])
                                 ]) +'\n')
            #ExperimentAcc_BedFileAcc_dic[FastqFileAcc_ExperimentAcc_dic[f]]
        except KeyError:
            Exp_withNoBedFile_yet.add(FastqFileAcc_ExperimentAcc_dic[f])
            Target_name_of_Exp_withNoBedFile_yet.add(FastqFileAcc_Target_dic[f])


with open(fastq_bed_pair_log, 'U') as f:
    with open('make_2_CNV.PeakFiltered_counts.bed.sh','w') as out:
        out.write('mkdir 2_CNV.PeakFiltered_counts_bed\n')
        for l in f.readlines()[1:]:
            ll = l.strip().split('\t')
            bed_files = ' '.join([bed_files_fp+b for b in ll[5].split(',')])
            out.write('intersectBed -wa -a '+CNVFiltered_counts_fp+ll[0]+'.CNVFiltered_counts.bed -b '+bed_files+' > 2_CNV.PeakFiltered_counts_bed/'+ll[0]+'.CNV.PeakFiltered_counts.bed\n')
        for f in FastqFileAcc_ExperimentAcc_dic:
            if FastqFileAcc_ExperimentAcc_dic[f] in Exp_withNoBedFile_yet:
                out.write('cp '+CNVFiltered_counts_fp+f+'.CNVFiltered_counts.bed 2_CNV.PeakFiltered_counts_bed/.\n' )
