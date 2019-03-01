BATCH=exoma_34_35_36_37
THREADS=9
DATAPATH=/home/bioinf/exoma/exoma_34_35_36_37/raw

export TARGETS=/home/murilo/Dropbox/Nextera/MedExome/MedExome_hg38_capture_targets.list 
export BAITS=/home/murilo/Dropbox/Nextera/MedExome/MedExome_hg38_empirical_targets.list 

export REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa"

export PICARD="java -jar -Xmx20G /opt/picard-tools-2.5.0/picard.jar"
export FASTQC=fastqc
export BWA="bwa"
export GATK="java -jar /opt/GenomeAnalysisTK.jar"
export KNOWN1="/home/bioinf/ref/var/clinvar-hg38_GATK.vcf"

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
    ## Arg 1: fq1
    ## Arg 2: fq2
    ## Arg 3: number of threads
    ## Arg 4: output dir

    sn=$(basename ${1})
    sn=${sn%_A*}
    sn=${sn%_T*}
    sn=${sn%_C*}
    sn=${sn%_G*}
    sn=${sn%_R*}

    fcln=$(zcat ${1} | head -n 1 | cut -f3-4 -d:)
    fcln=${fcln/\:/\.}

    lane=$(zcat ${1} | head -n 1 | cut -f4 -d:)

    bamout="${4}/${sn}.$lane.bam"
    samout="${4}/${sn}.$lane.sam"

    header="@RG\\tID:${fcln}\\tSM:${sn}\\tPL:ILLUMINA\\tLB:${sn}\\tCN:LaCTAD"

    date >> ${4}/${sn}.date
    
    ${BWA} mem -M -t ${3} -R ${header} ${REF} ${1} ${2} > ${samout}
    ${PICARD} SortSam INPUT=${samout} OUTPUT=${bamout}.withdups.bam SO=coordinate
    rm ${samout}
    ${PICARD} MarkDuplicates INPUT=${bamout}.withdups.bam OUTPUT=${bamout} METRICS_FILE=${bamout}.metrics
    rm ${bamout}.withdups.bam
    ${PICARD} BuildBamIndex INPUT=${bamout}
    ${PICARD} ValidateSamFile INPUT=${bamout} OUTPUT=${bamout}.validation VALIDATE_INDEX=true

    date >> ${4}/${sn}.date

}


function run_gatk {
    ## Arg 1: sample (this will be looked for at {sample}.bam                                                              

    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"
    RECALCSV="${INPUTBAM}.recal.csv"
    RECALBAM="${INPUTBAM}.recal.bam"
    VCF="results/$(basename ${INPUTBAM}).g.vcf"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS}
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} -compress 0 --disable_bam_indexing
    rm ${INPUTBAM}
    ${PICARD} BuildBamIndex INPUT=${REALNBAM}
    ${GATK} -T BaseRecalibrator -R ${REF} -I ${REALNBAM} -o ${RECALCSV} -knownSites ${KNOWN1}
    ${GATK} -T PrintReads -R ${REF} -I ${REALNBAM} -BQSR ${RECALCSV} -o ${RECALBAM}
    rm ${REALNBAM}
    ${PICARD} BuildBamIndex INPUT=${RECALBAM}
}

function run_gatk_MERGED {
    ## Arg 1: sample (this will be looked for at {sample}.bam                                                              

    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"
    VCF="results/$(basename ${INPUTBAM}).g.vcf"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS}
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} -compress 0 --disable_bam_indexing &&
    rm ${INPUTBAM}
    #${PICARD} BuildBamIndex INPUT=${REALNBAM}
    #${GATK} -T HaplotypeCaller -R ${REF} -I ${REALNBAM} -o ${VCF} -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 
}

function haplotypecaller {
    ## Arg 1: sample (this will be looked for at {sample}.bam                                                              

    INPUTBAM="${1}"
    VCF="results/$(basename ${INPUTBAM}).vcf"
    ${GATK} -T HaplotypeCaller -R ${REF} -I ${INPUTBAM} -o ${VCF}
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
     r2=`echo $r1 | sed 's/R1/R2/g'`
     filename=`echo $(basename ${r1}) | awk -F '.' '{print $1}'`
     exome=`echo $r1 | cut -f 7 -d /`

     cd $2
     
     ~murilo/TrimGalore-0.4.5/trim_galore --gzip --quality 20 --illumina --trim-n --length 40 --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --paired $r1 $r2 --path_to_cutadapt /home/murilo/cutadapt-1.10/bin/cutadapt

     mv $filename".R1_val_1.fq.gz" $filename".R1."$exome"_trim.fastq.gz"
     mv $filename".R2_val_2.fq.gz" $filename".R2."$exome"_trim.fastq.gz"

     mv $filename*_trimming_report.txt ../qc
     cd ..
}

function mergeBAM {
cd bams
sample=$1
bam=$sample.merged.bam
files=`ls $sample*.bam | awk '{print "INPUT=" $1}' | uniq`

$PICARD MergeSamFiles $files OUTPUT=$bam 2> $bam.MergeSamFiles.log &&
${PICARD} BuildBamIndex INPUT=${bam}
filesToRemove=`echo $files | sed 's/INPUT=//g'`
#rm $filesToRemove
cd ..
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

mkdir -p ${BATCH}/qc/aln ${BATCH}/bams ${BATCH}/results 

cd ${BATCH}

FQ1=$(ls ${DATAPATH}/*/*R1*gz)
FQ2=$(ls ${DATAPATH}/*/*R2*gz)

parallel --xapply -j ${THREADS} run_fastqc ::: ${FQ1} ::: ${FQ2} ::: qc
parallel -j ${THREADS} trimadapters ::: ${FQ1} ::: raw

FQ1=$(ls ${DATAPATH}/*R1*gz)
FQ2=$(ls ${DATAPATH}/*R2*gz)

parallel --xapply -j 1 bwa_align_md ::: ${FQ1} ::: ${FQ2} ::: ${THREADS} ::: bams

parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln

parallel -j ${THREADS} run_gatk ::: $(ls bams/*bam)

parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln

parallel -j ${THREADS} haplotypecaller ::: $(ls bams/*bam)

