OUTPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/"
JOBDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/jobs/filterRnaEditSites"

for CHR in {1..22}
do


echo "#!/bin/bash
#SBATCH --job-name=FilterRnaEditSites_chr$CHR
#SBATCH --output=FilterRnaEditSites_chr$CHR.out
#SBATCH --error=FilterRnaEditSites_chr$CHR.err
#SBATCH --time=5:59:59
#SBATCH --cpus-per-task 4
#SBATCH --mem 40gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --qos=dev

ENVIRONMENT_DIR=\".\"
set -e
set -u

echo \"## \"\$(date)\" Start \$0\"

#Load gatk module
module load GATK/3.7-Java-1.8.0_74
module list

mkdir -p $OUTPUTDIR


java -Xmx38g -Djava.io.tmpdir=\${TMPDIR} -jar \${EBROOTGATK}/GenomeAnalysisTK.jar \\
    -T SelectVariants \\
    -R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\
    -V /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.PASSonly.BiallelicSNVsOnly.noDiagnostics.gg.vcf.gz \\
    -o /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz \\
    --excludeIntervals /apps/data/rnaedit.com/Human_AG_all_hg19_v2.sorted.chr$CHR.intervals

cd $OUTPUTDIR
md5sum genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz > genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz.md5
md5sum genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz.tbi > genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz.tbi.md5
cd -

echo \"## \"\$(date)\" ##  \$0 Done \"

" > $JOBDIR/FilterRnaEditSites_chr$CHR.sh

done
