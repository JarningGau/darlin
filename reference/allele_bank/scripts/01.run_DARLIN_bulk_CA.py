import os

sample_dict = {
	'LL583_Gr_bulk_RNA_CA': {
		"fq1": "raw_fastq/LL583-Gr-3A_S10_L001_R1_001.fastq.gz",
		"fq2": "raw_fastq/LL583-Gr-3A_S10_L001_R2_001.fastq.gz",
	},
	'LL584_Gr_bulk_RNA_CA': {
		"fq1": "raw_fastq/LL584-Gr-3A_S11_L001_R1_001.fastq.gz",
		"fq2": "raw_fastq/LL584-Gr-3A_S11_L001_R2_001.fastq.gz",
	},
	'LL638_Gr_bulk_RNA_CA': {
		"fq1": "raw_fastq/LL638-Gr-3A_S12_L001_R1_001.fastq.gz",
		"fq2": "raw_fastq/LL638-Gr-3A_S12_L001_R2_001.fastq.gz",
	}
}

for sample_id, fq in sample_dict.items():
	cmd = f'darlin bulk run \
		--sample-id {sample_id} \
		--fq1 {fq["fq1"]} \
		--fq2 {fq["fq2"]} \
		--output-dir ./output \
		--threads 8 \
		--protocol pe250 \
		--locus Col1a1 \
		--reads-cutoff 1 2 \
		--umi-ld 1 \
		--lb-hd-relative 0.01 \
		--umi-len 12 \
		--keep-pear'
	os.system(cmd)
