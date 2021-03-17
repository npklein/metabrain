


module load GATK

for CHR in {1..22}
do

  java -Xmx200g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
    -V /groups/umcg-bios/tmp04/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.gg.vcf.gz \
    -o /groups/umcg-bios/tmp04/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.gg.vcf.gz \
    --restrictAllelesTo BIALLELIC \
    -selectType SNP \
    -L /groups/umcg-bios/tmp04/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.gg.vcf.gz \
    -nt 20

done



