# build the indexes (remember toggle for moni vcf or only fasta)

./columba_moni_build_pfp.sh \
-r ../data/simon_chrom19_index/columba_index \
-f ../data/simon_chrom19_references/chr19_only.fasta.gz \
-m ../data/simon_chrom19_references/chr19_only.fasta.gz \
-v ../data/simon_chrom19_references/ALL.chr19.subset.vcf.gz \
-s ../data/simon_chrom19_references/sample_list.txt \
-o ../data/simon_chrom19_index/moni_index simon_chrom19_output/

# run de allignent for the indexes (single ended)
./columba_moni_allign.sh \
  -L ../data/simon_chrom19_index/runs \
  -c ../data/simon_chrom19_index/columba_index \
  -m ../data/simon_chrom19_index/moni_index \
  -f ../data/simon_reads/reads_subset.fastq


# Sampling fasta and vcf files do it with samtools and use your brain to extract actual fcking proper data

# fastq file (reads)
# Extract 40,000 lines (10k reads) and compress
head -n 40000 ../data/simon_reads/SRR17981962.fastq > ../data/simon_reads/reads_subset.fastq
