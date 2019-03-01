BATCH=exoma_Dravet_Mendelics_ago_18
THREADS=15
DATAPATH=/home/bioinf/exoma/exoma_Dravet_Mendelics_ago_18/raw

export TARGETS="/home/murilo/Dropbox/Nextera/Agillent/SureSelectXTAllexon_v6_HG38.list"
export BAITS="/home/murilo/Dropbox/Nextera/Agillent/SureSelectXTAllexon_v6_HG38.list"

export REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa"

export PICARD="java -jar -Xmx20G /opt/picard-tools-2.5.0/picard.jar"
export FASTQC=fastqc
export BWA="bwa"
export GATK="java -jar /opt/GenomeAnalysisTK.jar"

export KNOWN1=/home/bioinf/ref/GRCh38.p1/HG38_ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf
export KNOWN2=/home/bioinf/ref/GRCh38.p1/HG38_ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf

export VEP="/home/murilo/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl"

function picard_hsmetrics
{
    # arg1: input BAM
    # arg2: output dir
    OUTFILE=${2}/$(basename ${1}).metrics
    ${PICARD} CalculateHsMetrics VERBOSITY=WARNING BAIT_INTERVALS=${BAITS} TARGET_INTERVALS=${TARGETS} INPUT=${1} OUTPUT=${OUTFILE} PER_TARGET_COVERAGE=${OUTFILE}.unit REFERENCE_SEQUENCE=${REF}
}

function run_fastqc
{
    ## arg 1: it's the sample name
    OUTDIR=$3
    ${FASTQC} ${1} ${2} --outdir ${OUTDIR}
}

function bwa_align_md
{
    R2=`echo $1 | sed 's/_R1_/_R2_/'`

    sn=$(basename ${1})
    sn=${sn%_*}

    nomeAmostra=`echo $(basename ${1}) | cut -f 1 -d _`

    fcln=$(zcat ${1} | head -n 1 | cut -f3-4 -d:)
    fcln=${fcln/\:/\.}

    lane=$(zcat ${1} | head -n 1 | cut -f4 -d:)

    bamout="${3}/${sn}.$lane.bam"
    samout="${3}/${sn}.$lane.sam"

    header="@RG\\tID:${fcln}\\tSM:${nomeAmostra}\\tPL:ILLUMINA\\tLB:${sn}\\tCN:MENDELICS"

 
    ${BWA} mem -M -t ${2} -R ${header} ${REF} ${1} ${R2} > ${samout}
    ${PICARD} SortSam INPUT=${samout} OUTPUT=${bamout}.withdups.bam SO=coordinate
    rm ${samout}
    ${PICARD} MarkDuplicates INPUT=${bamout}.withdups.bam OUTPUT=${bamout} METRICS_FILE=${bamout}.metrics
    rm ${bamout}.withdups.bam
    ${PICARD} BuildBamIndex INPUT=${bamout}
    ${PICARD} ValidateSamFile INPUT=${bamout} OUTPUT=${bamout}.validation VALIDATE_INDEX=true

}

function run_gatk {
    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"
    RECALCSV="${INPUTBAM}.recal.csv"
    RECALBAM="${INPUTBAM}.recal.bam"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS}
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} --disable_bam_indexing -known $KNOWN1 -known $KNOWN2 && rm ${INPUTBAM}
    ${PICARD} BuildBamIndex INPUT=${REALNBAM}
    ${GATK} -T BaseRecalibrator -R ${REF} -I ${REALNBAM} -o ${RECALCSV} -knownSites $KNOWN1 -knownSites $KNOWN2 
    ${GATK} -T PrintReads -R ${REF} -I ${REALNBAM} -BQSR ${RECALCSV} -o ${RECALBAM} && rm ${REALNBAM}
    ${PICARD} BuildBamIndex INPUT=${RECALBAM}
}

function run_gatk_MERGED {
    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS} &&
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} --disable_bam_indexing &&
    rm ${INPUTBAM}
    #${PICARD} BuildBamIndex INPUT=${REALNBAM}
    #${GATK} -T HaplotypeCaller -R ${REF} -I ${REALNBAM} -o ${VCF} -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 
}

function HaplotypeCaller {
    REALNBAM="${1}"
    VCF="../results/$1.vcf"
    ${PICARD} BuildBamIndex INPUT=${REALNBAM} && 
    ${GATK} -T HaplotypeCaller -R ${REF} -I ${REALNBAM} -o ${VCF} -maxAltAlleles 10
}

function haplotypecaller {
    folder=$1                       
    myfiles=`ls $folder/*bam`
    INPUT=`echo $myfiles | sed 's/ / -I /g' | sed 's/^/-I /g'`
    VCF="results/$folder.vcf"
    ${GATK} -T HaplotypeCaller -R ${REF} $INPUT -o ${VCF} -maxAltAlleles 30
}


function mergegz
{
f=$1
cd $f
filename=`ls *gz | awk -F "_R" '{print $1}' | uniq`
for file in $filename
do
zcat $file*R1* | gzip > ../$file.R1.fastq.gz
zcat $file*R2* | gzip > ../$file.R2.fastq.gz
done
cd ..
}

function findadapters
{
     r1=$1
     r2=`echo $r1 | sed 's/R1/R2/g' | sed 's/trim\_1/trim\_2/'`
     filename=`echo $(basename ${r1}) | awk -F '_R' '{print $1}'`
     
     AdapterRemoval --identify-adapters --file1 $r1 --file2 $r2 > $2/$filename.AdapterRemoval
}

function trimadapters
{
     r1=$1
DATAPATH=$2
     r2=`echo $r1 | sed 's/R1/R2/g'`
     filename=`echo $(basename ${r1}) | awk -F '.R' '{print $1}'`

cutadapt -o $r1.cutadapt $r1 -u -140 -q 25,25
cutadapt -o $r2.cutadapt $r2 -u -140 -q 25,25
     
~murilo/trim_galore --gzip --paired $r1.cutadapt $r2.cutadapt --path_to_cutadapt /home/murilo/cutadapt-1.10/bin/cutadapt --clip_R1 10 --clip_R2 3 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --length 40 --max_n 3 -q 25 -o $DATAPATH

}

function trimadapters_ex4
{
     r1=$1
     r2=`echo $r1 | sed 's/R1/R2/'`
     filename=`echo $(basename ${r1}) | awk -F '.R' '{print $1}'`

#cutadapt -o $r1.cutadapt $r1 -u -150 -q 25,25
#cutadapt -o $r2.cutadapt $r2 -u -150 -q 25,25
     
~murilo/trim_galore --gzip --paired $r1 $r2 --path_to_cutadapt /home/murilo/cutadapt-1.10/bin/cutadapt --clip_R1 18 --clip_R2 18 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --length 40 --max_n 5 -q 25 -o $2

}

function mergeBAM {
#cd bams
sample=$1
bam=$sample.merged.bam
files=`ls $sample*.bam | awk '{print "INPUT=" $1}' | uniq`

$PICARD MergeSamFiles $files OUTPUT=$bam 2> $bam.MergeSamFiles.log &&
${PICARD} BuildBamIndex INPUT=${bam} && 
filesToRemove=`echo $files | sed 's/INPUT=//g'` &&
rm $filesToRemove
# cd ..
}

function removeSingles {
R2=`echo $1 | sed 's/.R1./.R2./'`
~murilo/Dropbox/bbmap/repair.sh in=$1 in2=$R2 out1=$1.fix.fastq.gz out2=$R2.fix.fastq.gz overwrite=true 2> $1.fix.report.txt 
}

export -f picard_hsmetrics
export -f run_fastqc
export -f bwa_align_md
export -f run_gatk
export -f findadapters
export -f trimadapters
export -f mergeBAM
export -f run_gatk_MERGED
export -f haplotypecaller
export -f trimadapters_ex4
export -f removeSingles

cd ${BATCH}

for fq in `ls ${DATAPATH}/*R1*.gz`
do
bwa_align_md $fq $THREADS bams
done

for bam in `ls bams/*bam`
do
run_gatk $bam '2>' $bam.log
done

parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln

for bam in `ls *bam | cut -f 1 -d . | sort | uniq`
do
mergeBAM $bam &
done

for bam in `ls bams/*bam`
do
run_gatk_MERGED $bam '2>' $bam.log
done


cd bams/ 
for bam in `ls *bam`
do
HaplotypeCaller $bam &
done

function GenotypeGVCFs {
f1=$1
f2=$2
NOME1=`echo $(basename ${f1}) | cut -f 1 -d _`
NOME2=`echo $(basename ${f2}) | cut -f 1 -d _`
$GATK -T GenotypeGVCFs -R $REF -o $NOME1.$NOME2.vcf -V $f1 -V $f2 --max_alternate_alleles 30
}

cd /home/bioinf/paineis/MTOR_jun_18/results


out='/home/bioinf/exoma/exoma_4_23_25/qc/aln/'
out='/home/bioinf/paineis/MTOR_jun_18/qc/aln'
for bam in `ls *bam`
do picard_hsmetrics $bam $out &
done

haplotypecaller bams

