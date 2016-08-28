# replace -PATH-TO- with location of alleleDB_v2.0

bowtie-build --offrate 2 ../test_inp_HG00100_pgenome_hg19_chr20/20_HG00100_maternal.fa ../test_inp_HG00100_pgenome_hg19_chr20/HG00100_maternal > ../test_inp_HG00100_pgenome_hg19_chr20/bowtie_build.maternal.log
bowtie-build --offrate 2 ../test_inp_HG00100_pgenome_hg19_chr20/20_HG00100_paternal.fa ../test_inp_HG00100_pgenome_hg19_chr20/HG00100_paternal > ../test_inp_HG00100_pgenome_hg19_chr20/bowtie_build.paternal.log


/-PATH-TO-/alleleDB_v2.0/alleledb_pipeline/alleledb.sh \
	HG00100 test_inp_HG00100_pgenome_hg19_chr20 \
	/-PATH-TO-/alleleDB_v2.0/test/test_inp_HG00100_pgenome_hg19_chr20 \
	/-PATH-TO-/alleleDB_v2.0/test/test_inp_HG00100_pgenome_hg19_chr20/snp.calls.bed \
	/-PATH-TO-/alleleDB_v2.0/test/test_inp_HG00100_rnaseq_chr20/ \
	HG00100_ase_test \
	/-PATH-TO-/alleleDB_v2.0/alleledb_pipeline \
	/-PATH-TO-/alleleDB_v2.0/alleledb_pipeline/PIPELINE.mk \
	0.05 \
	ase \
	0	


