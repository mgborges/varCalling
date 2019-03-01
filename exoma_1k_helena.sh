BATCH=helena_fastq
THREADS=15
DATAPATH=/iscsi/helena/helena_fastq/raw

export BIPMED=/home/murilo/bipmed_ann.HG38.vcf
export ABRAOM=/home/murilo/Dropbox/ABRAON/new_ABRaOM_60+_SABE_609_exomes_annotated.vSR.PASS_trimm_hg38.VEP.vcf
export BIPMED_ARRAY=/home/murilo/BIPMED_SNP_ARRAY_296_hg38.vcf.DECOMPOSE.DP.vcf

export KNOWN1=/home/bioinf/ref/GRCh38.p1/HG38_ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf
export KNOWN2=/home/bioinf/ref/GRCh38.p1/HG38_ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf

export TARGETS="/home/murilo/Dropbox/Nextera/Agillent/SureSelectXTAllexon_v6_HG38.list"
export BAITS="/home/murilo/Dropbox/Nextera/Agillent/SureSelectXTAllexon_v6_HG38.list"
export REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa"

export PICARD="java -jar -Xmx20G /opt/picard-tools-2.5.0/picard.jar"
export FASTQC=fastqc
export BWA="bwa"
export GATK="java -jar /opt/GenomeAnalysisTK.jar"
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
    sn=${sn%.R*}

    fcln=$(zcat ${1} | head -n 1 | cut -f3-4 -d: | cut -f 1 -d . | sed 's/@//')
    fcln=${fcln/\:/\.}

    lane=$(zcat ${1} | head -n 1 | cut -f4 -d: | cut -f 2 -d /)

    bamout="${4}/${sn}.$lane.bam"
    samout="${4}/${sn}.$lane.sam"

    header="@RG\\tID:${fcln}\\tSM:${sn}\\tPL:ILLUMINA\\tLB:${sn}\\tCN:BGI"
    
    ${BWA} mem -M -t ${3} -R ${header} ${REF} ${1} ${2} > ${samout}
    ${PICARD} SortSam INPUT=${samout} OUTPUT=${bamout}.withdups.bam SO=coordinate
    rm ${samout}
    ${PICARD} MarkDuplicates INPUT=${bamout}.withdups.bam OUTPUT=${bamout} METRICS_FILE=${bamout}.metrics
    rm ${bamout}.withdups.bam
    ${PICARD} BuildBamIndex INPUT=${bamout}
    ${PICARD} ValidateSamFile INPUT=${bamout} OUTPUT=${bamout}.validation VALIDATE_INDEX=true


}


function run_gatk_BIPMED {
    ## Arg 1: sample (this will be looked for at {sample}.bam                                                              

    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"
    RECALCSV="${INPUTBAM}.recal.csv"
    RECALBAM="${INPUTBAM}.recal.bam"
    VCF="results/$(basename ${INPUTBAM}).g.vcf"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS}
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} -compress 0 --disable_bam_indexing -known $BIPMED -known $ABRAOM -known $BIPMED_ARRAY
    rm ${INPUTBAM}
    ${PICARD} BuildBamIndex INPUT=${REALNBAM}
    ${GATK} -T BaseRecalibrator -R ${REF} -I ${REALNBAM} -o ${RECALCSV} -knownSites $BIPMED -knownSites $ABRAOM -knownSites $BIPMED_ARRAY
    ${GATK} -T PrintReads -R ${REF} -I ${REALNBAM} -BQSR ${RECALCSV} -o ${RECALBAM}
    rm ${REALNBAM}
    ${PICARD} BuildBamIndex INPUT=${RECALBAM}
}

function run_gatk_1k {
    ## Arg 1: sample (this will be looked for at {sample}.bam                                                              

    INPUTBAM="${1}"
    INTERVALS="${INPUTBAM}.intervals"
    REALNBAM="${INPUTBAM}.realn.bam"
    RECALCSV="${INPUTBAM}.recal.csv"
    RECALBAM="${INPUTBAM}.recal.bam"
    VCF="results/$(basename ${INPUTBAM}).g.vcf"

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS}
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} -compress 0 --disable_bam_indexing -known $KNOWN1 -known $KNOWN2
    rm ${INPUTBAM}
    ${PICARD} BuildBamIndex INPUT=${REALNBAM}
    ${GATK} -T BaseRecalibrator -R ${REF} -I ${REALNBAM} -o ${RECALCSV} -knownSites $KNOWN1 -knownSites $KNOWN2 
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
    folder=$1                       
    myfiles=`ls $folder/*bam`
    INPUT=`echo $myfiles | sed 's/ / -I /g' | sed 's/^/-I /g'`
    VCF="/iscsi/murilo/realignment_article/br_exom/$folder.vcf"
    java -jar -Xmx170G /opt/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REF} $INPUT -o ${VCF} -nct $THREADS -rf BadCigar
}

function haplotypecaller_GVCF {
    INPUT=$1
    SAMPLE=`echo $(basename ${INPUT}) | cut -f 1 -d .`
    FOLDER=`echo $INPUT | cut -f 1 -d '/'`
    VCF="/iscsi/murilo/realignment_article/br_exom/$SAMPLE.$FOLDER.g.vcf"
    java -jar /opt/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REF} -I $INPUT -o ${VCF} -rf BadCigar -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000  2> $VCF.$FOLDER.log
    bgzip $VCF
}


function mergegz
{
f=$1
cd $f
filename=`ls *gz | awk -F "_R" '{print $1}' | uniq`

for file in $filename
do
newfilename=`echo $file | sed 's/_//g'`
zcat $file*R1* | gzip > ../$newfilename.R1.fastq.gz &
zcat $file*R2* | gzip > ../$newfilename.R2.fastq.gz
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
     filename=`echo $(basename ${r1}) | awk -F '_R' '{print $1}'`

     cd $2
     
     ~murilo/trim_galore_zip/trim_galore --gzip --length 100 --paired $r1 $r2 --path_to_cutadapt /home/murilo/cutadapt-1.10/bin/cutadapt

     mv $filename"_R1_val_1.fq.gz" $filename"_R1_trim_1.fastq.gz"
     mv $filename"_R2_val_2.fq.gz" $filename"_R2_trim_2.fastq.gz"

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
export -f run_gatk_BIPMED
export -f run_gatk_1k
export -f findadapters
export -f trimadapters
export -f mergeBAM
export -f run_gatk_MERGED
export -f haplotypecaller

mkdir -p ${BATCH}/qc/aln ${BATCH}/bams ${BATCH}/results 

cd /iscsi/helena/helena_fastq

for bam in `ls *GBR*`; 
do 
samtools sort -n $bam -o $bam.sort &&
bamToFastq -i $bam.sort.bam -fq ../raw/originals/$bam.R1.fastq -fq2 ../raw/originals/$bam.R2.fastq & 
done

#for f in `ls *R1*`
#do 
#R1=$f
#R2=`echo $f | sed 's/R1/R2/g'`; 
#~murilo/Dropbox/bbmap/repair.sh in1=$R1 in2=$R2 out1=$R1.REPAIR.fastq.gz out2=$R2.REPAIR.fastq.gz outsingle=$R1.single.fastq.gz
#done

FQ1=$(ls ${DATAPATH}/*R1*gz)
FQ2=$(ls ${DATAPATH}/*R2*gz)

#parallel --xapply -j ${THREADS} run_fastqc ::: ${FQ1} ::: ${FQ2} ::: qc
#parallel -j ${THREADS} findadapters ::: ${FQ1} ::: qc

#parallel --xapply -j 1 bwa_align_md ::: ${FQ1} ::: ${FQ2} ::: ${THREADS} ::: bams

for R1 in $FQ1
do 
R2=`echo $R1 | sed 's/R1/R2/g'`
echo bwa_align_md $R1 $R2 ${THREADS} bams
done


parallel -j ${THREADS} run_gatk_1k ::: $(ls bams/*bam)


for file in `ls bams/*bam`
do
echo run_gatk_1k $file '&&'
done > tmp

run_gatk_1k bams/101.7.bam &&
run_gatk_1k bams/10.1.bam &&
run_gatk_1k bams/102.7.bam &&
run_gatk_1k bams/103.7.bam &&
run_gatk_1k bams/104.7.bam &&
run_gatk_1k bams/105.7.bam &&
run_gatk_1k bams/106.7.bam &&
run_gatk_1k bams/107.7.bam &&
run_gatk_1k bams/108.8.bam &&
run_gatk_1k bams/109.8.bam &&
run_gatk_1k bams/110.8.bam &
run_gatk_1k bams/111.8.bam &&
run_gatk_1k bams/11.1.bam &&
run_gatk_1k bams/112.8.bam &&
run_gatk_1k bams/113.8.bam &&
run_gatk_1k bams/114.8.bam &&
run_gatk_1k bams/115.8.bam &&
run_gatk_1k bams/116.8.bam &&
run_gatk_1k bams/117.8.bam &&
run_gatk_1k bams/118.8.bam &&
run_gatk_1k bams/119.8.bam &
run_gatk_1k bams/1.1.bam &&
run_gatk_1k bams/120.8.bam &&
run_gatk_1k bams/121.8.bam &&
run_gatk_1k bams/12.1.bam &&
run_gatk_1k bams/122.8.bam &&
run_gatk_1k bams/123.8.bam &&
run_gatk_1k bams/124.8.bam &&
run_gatk_1k bams/125.8.bam &&
run_gatk_1k bams/126.8.bam &&
run_gatk_1k bams/127.8.bam &
run_gatk_1k bams/128.8.bam &&
run_gatk_1k bams/129.8.bam &&
run_gatk_1k bams/130.8.bam &&
run_gatk_1k bams/131.2.bam &&
run_gatk_1k bams/13.1.bam &&
run_gatk_1k bams/132.2.bam &&
run_gatk_1k bams/14.1.bam &&
run_gatk_1k bams/15.1.bam &&
run_gatk_1k bams/16.1.bam &&
run_gatk_1k bams/17.1.bam &
run_gatk_1k bams/18.1.bam &&
run_gatk_1k bams/19.2.bam &&
run_gatk_1k bams/20.2.bam &&
run_gatk_1k bams/21.1.bam &&
run_gatk_1k bams/2.1.bam &&
run_gatk_1k bams/22.1.bam &&
run_gatk_1k bams/23.1.bam &&
run_gatk_1k bams/24.1.bam &&
run_gatk_1k bams/25.1.bam &&
run_gatk_1k bams/26.4.bam &
run_gatk_1k bams/27.4.bam &&
run_gatk_1k bams/28.4.bam &&
run_gatk_1k bams/29.4.bam &&
run_gatk_1k bams/296.3.bam &&
run_gatk_1k bams/30.4.bam &&
run_gatk_1k bams/31.4.bam &&
run_gatk_1k bams/3.1.bam &&
run_gatk_1k bams/32.5.bam &&
run_gatk_1k bams/33.5.bam &&
run_gatk_1k bams/34.5.bam &
run_gatk_1k bams/35.5.bam &&
run_gatk_1k bams/36.5.bam &&
run_gatk_1k bams/37.5.bam &&
run_gatk_1k bams/38.5.bam &&
run_gatk_1k bams/39.5.bam &&
run_gatk_1k bams/40.5.bam &&
run_gatk_1k bams/42.5.bam &&
run_gatk_1k bams/43.5.bam &&
run_gatk_1k bams/44.5.bam &&
run_gatk_1k bams/45.5.bam &
run_gatk_1k bams/46.5.bam &&
run_gatk_1k bams/47.5.bam &&
run_gatk_1k bams/48.5.bam &&
run_gatk_1k bams/49.5.bam &&
run_gatk_1k bams/50.5.bam &&
run_gatk_1k bams/51.5.bam &&
run_gatk_1k bams/5.1.bam &&
run_gatk_1k bams/52.5.bam &&
run_gatk_1k bams/53.5.bam &&
run_gatk_1k bams/54.5.bam &
run_gatk_1k bams/55.6.bam &&
run_gatk_1k bams/56.6.bam &&
run_gatk_1k bams/57.6.bam &&
run_gatk_1k bams/58.6.bam &&
run_gatk_1k bams/59.6.bam &&
run_gatk_1k bams/60.6.bam &&
run_gatk_1k bams/61.6.bam &&
run_gatk_1k bams/6.1.bam &&
run_gatk_1k bams/62.6.bam &&
run_gatk_1k bams/63.6.bam &
run_gatk_1k bams/64.6.bam &&
run_gatk_1k bams/65.6.bam &&
run_gatk_1k bams/66.6.bam &&
run_gatk_1k bams/67.6.bam &&
run_gatk_1k bams/68.6.bam &&
run_gatk_1k bams/69.6.bam &&
run_gatk_1k bams/70.6.bam &&
run_gatk_1k bams/71.6.bam &&
run_gatk_1k bams/7.1.bam &&
run_gatk_1k bams/72.6.bam &
run_gatk_1k bams/73.6.bam &&
run_gatk_1k bams/74.6.bam &&
run_gatk_1k bams/75.6.bam &&
run_gatk_1k bams/76.6.bam &&
run_gatk_1k bams/77.7.bam &&
run_gatk_1k bams/78.7.bam &&
run_gatk_1k bams/79.7.bam &&
run_gatk_1k bams/80.7.bam &&
run_gatk_1k bams/81.7.bam &&
run_gatk_1k bams/8.1.bam &
run_gatk_1k bams/82.7.bam &&
run_gatk_1k bams/83.7.bam &&
run_gatk_1k bams/84.7.bam &&
run_gatk_1k bams/85.7.bam &&
run_gatk_1k bams/9.1.bam &&
run_gatk_1k bams/94.7.bam &&
run_gatk_1k bams/95.7.bam &&
run_gatk_1k bams/96.7.bam &&
run_gatk_1k bams/97.7.bam &&
run_gatk_1k bams/98.7.bam &&
run_gatk_1k bams/99.7.bam &

cp -r bams BAMS_1K

run_gatk_BIPMED bams/101.7.bam &&
run_gatk_BIPMED bams/10.1.bam &&
run_gatk_BIPMED bams/102.7.bam &&
run_gatk_BIPMED bams/103.7.bam &&
run_gatk_BIPMED bams/104.7.bam &&
run_gatk_BIPMED bams/105.7.bam &&
run_gatk_BIPMED bams/106.7.bam &&
run_gatk_BIPMED bams/107.7.bam &&
run_gatk_BIPMED bams/108.8.bam &&
run_gatk_BIPMED bams/109.8.bam &&
run_gatk_BIPMED bams/110.8.bam &
run_gatk_BIPMED bams/111.8.bam &&
run_gatk_BIPMED bams/11.1.bam &&
run_gatk_BIPMED bams/112.8.bam &&
run_gatk_BIPMED bams/113.8.bam &&
run_gatk_BIPMED bams/114.8.bam &&
run_gatk_BIPMED bams/115.8.bam &&
run_gatk_BIPMED bams/116.8.bam &&
run_gatk_BIPMED bams/117.8.bam &&
run_gatk_BIPMED bams/118.8.bam &&
run_gatk_BIPMED bams/119.8.bam &
run_gatk_BIPMED bams/1.1.bam &&
run_gatk_BIPMED bams/120.8.bam &&
run_gatk_BIPMED bams/121.8.bam &&
run_gatk_BIPMED bams/12.1.bam &&
run_gatk_BIPMED bams/122.8.bam &&
run_gatk_BIPMED bams/123.8.bam &&
run_gatk_BIPMED bams/124.8.bam &&
run_gatk_BIPMED bams/125.8.bam &&
run_gatk_BIPMED bams/126.8.bam &&
run_gatk_BIPMED bams/127.8.bam &
run_gatk_BIPMED bams/128.8.bam &&
run_gatk_BIPMED bams/129.8.bam &&
run_gatk_BIPMED bams/130.8.bam &&
run_gatk_BIPMED bams/131.2.bam &&
run_gatk_BIPMED bams/13.1.bam &&
run_gatk_BIPMED bams/132.2.bam &&
run_gatk_BIPMED bams/14.1.bam &&
run_gatk_BIPMED bams/15.1.bam &&
run_gatk_BIPMED bams/16.1.bam &&
run_gatk_BIPMED bams/17.1.bam &
run_gatk_BIPMED bams/18.1.bam &&
run_gatk_BIPMED bams/19.2.bam &&
run_gatk_BIPMED bams/20.2.bam &&
run_gatk_BIPMED bams/21.1.bam &&
run_gatk_BIPMED bams/2.1.bam &&
run_gatk_BIPMED bams/22.1.bam &&
run_gatk_BIPMED bams/23.1.bam &&
run_gatk_BIPMED bams/24.1.bam &&
run_gatk_BIPMED bams/25.1.bam &&
run_gatk_BIPMED bams/26.4.bam &
run_gatk_BIPMED bams/27.4.bam &&
run_gatk_BIPMED bams/28.4.bam &&
run_gatk_BIPMED bams/29.4.bam &&
run_gatk_BIPMED bams/296.3.bam &&
run_gatk_BIPMED bams/30.4.bam &&
run_gatk_BIPMED bams/31.4.bam &&
run_gatk_BIPMED bams/3.1.bam &&
run_gatk_BIPMED bams/32.5.bam &&
run_gatk_BIPMED bams/33.5.bam &&
run_gatk_BIPMED bams/34.5.bam &
run_gatk_BIPMED bams/35.5.bam &&
run_gatk_BIPMED bams/36.5.bam &&
run_gatk_BIPMED bams/37.5.bam &&
run_gatk_BIPMED bams/38.5.bam &&
run_gatk_BIPMED bams/39.5.bam &&
run_gatk_BIPMED bams/40.5.bam &&
run_gatk_BIPMED bams/42.5.bam &&
run_gatk_BIPMED bams/43.5.bam &&
run_gatk_BIPMED bams/44.5.bam &&
run_gatk_BIPMED bams/45.5.bam &
run_gatk_BIPMED bams/46.5.bam &&
run_gatk_BIPMED bams/47.5.bam &&
run_gatk_BIPMED bams/48.5.bam &&
run_gatk_BIPMED bams/49.5.bam &&
run_gatk_BIPMED bams/50.5.bam &&
run_gatk_BIPMED bams/51.5.bam &&
run_gatk_BIPMED bams/5.1.bam &&
run_gatk_BIPMED bams/52.5.bam &&
run_gatk_BIPMED bams/53.5.bam &&
run_gatk_BIPMED bams/54.5.bam &
run_gatk_BIPMED bams/55.6.bam &&
run_gatk_BIPMED bams/56.6.bam &&
run_gatk_BIPMED bams/57.6.bam &&
run_gatk_BIPMED bams/58.6.bam &&
run_gatk_BIPMED bams/59.6.bam &&
run_gatk_BIPMED bams/60.6.bam &&
run_gatk_BIPMED bams/61.6.bam &&
run_gatk_BIPMED bams/6.1.bam &&
run_gatk_BIPMED bams/62.6.bam &&
run_gatk_BIPMED bams/63.6.bam &
run_gatk_BIPMED bams/64.6.bam &&
run_gatk_BIPMED bams/65.6.bam &&
run_gatk_BIPMED bams/66.6.bam &&
run_gatk_BIPMED bams/67.6.bam &&
run_gatk_BIPMED bams/68.6.bam &&
run_gatk_BIPMED bams/69.6.bam &&
run_gatk_BIPMED bams/70.6.bam &&
run_gatk_BIPMED bams/71.6.bam &&
run_gatk_BIPMED bams/7.1.bam &&
run_gatk_BIPMED bams/72.6.bam &
run_gatk_BIPMED bams/73.6.bam &&
run_gatk_BIPMED bams/74.6.bam &&
run_gatk_BIPMED bams/75.6.bam &&
run_gatk_BIPMED bams/76.6.bam &&
run_gatk_BIPMED bams/77.7.bam &&
run_gatk_BIPMED bams/78.7.bam &&
run_gatk_BIPMED bams/79.7.bam &&
run_gatk_BIPMED bams/80.7.bam &&
run_gatk_BIPMED bams/81.7.bam &&
run_gatk_BIPMED bams/8.1.bam &
run_gatk_BIPMED bams/82.7.bam &&
run_gatk_BIPMED bams/83.7.bam &&
run_gatk_BIPMED bams/84.7.bam &&
run_gatk_BIPMED bams/85.7.bam &&
run_gatk_BIPMED bams/9.1.bam &&
run_gatk_BIPMED bams/94.7.bam &&
run_gatk_BIPMED bams/95.7.bam &&
run_gatk_BIPMED bams/96.7.bam &&
run_gatk_BIPMED bams/97.7.bam &&
run_gatk_BIPMED bams/98.7.bam &&
run_gatk_BIPMED bams/99.7.bam &


#parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln
cp -r bams BAMS_NO_REALN

parallel -j ${THREADS} run_gatk_BIPMED ::: $(ls bams/*bam)

cp -r bams BAMS_BIPMED
cp -r BAMS_NO_REALN bams



parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln

cd /home/murilo/helena_WES_real

haplotypecaller_GVCF BAMS_1K/101.7.bam.recal.bam  &&
haplotypecaller_GVCF BAMS_1K/10.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/102.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/103.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/104.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/105.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/106.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/107.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/108.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/109.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/110.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/111.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/11.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/112.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/113.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/114.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/115.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/116.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/117.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/118.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/119.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/1.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/120.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/121.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/12.1.bam.recal.bam &
haplotypecaller_GVCF BAMS_1K/122.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/123.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/124.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/125.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/126.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/127.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/128.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/129.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/130.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/131.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/13.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/132.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/14.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/15.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/16.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/17.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/18.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/19.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/20.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/21.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/2.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/22.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/23.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/24.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/25.1.bam.recal.bam &
haplotypecaller_GVCF BAMS_1K/26.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/27.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/28.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/29.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/296.3.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/30.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/31.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/3.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/32.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/33.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/34.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/35.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/36.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/37.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/38.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/39.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/40.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/42.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/43.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/44.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/45.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/46.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/47.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/48.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/49.5.bam.recal.bam &
haplotypecaller_GVCF BAMS_1K/50.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/51.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/5.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/52.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/53.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/54.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/55.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/56.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/57.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/58.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/59.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/60.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/61.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/6.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/62.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/63.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/64.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/65.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/66.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/67.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/68.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/69.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/70.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/71.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/7.1.bam.recal.bam &
haplotypecaller_GVCF BAMS_1K/72.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/73.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/74.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/75.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/76.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/77.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/78.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/79.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/80.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/81.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/8.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/82.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/83.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/84.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/85.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/9.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/94.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/95.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/96.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/97.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/98.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_1K/99.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/101.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/10.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/102.7.bam.recal.bam &
haplotypecaller_GVCF BAMS_BIPMED/103.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/104.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/105.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/106.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/107.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/108.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/109.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/110.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/111.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/11.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/112.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/113.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/114.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/115.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/116.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/117.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/118.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/119.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/1.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/120.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/121.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/12.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/122.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/123.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/124.8.bam.recal.bam &
haplotypecaller_GVCF BAMS_BIPMED/125.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/126.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/127.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/128.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/129.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/130.8.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/131.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/13.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/132.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/14.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/15.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/16.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/17.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/18.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/19.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/20.2.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/21.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/2.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/22.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/23.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/24.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/25.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/26.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/27.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/28.4.bam.recal.bam &
haplotypecaller_GVCF BAMS_BIPMED/29.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/296.3.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/30.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/31.4.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/3.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/32.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/33.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/34.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/35.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/36.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/37.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/38.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/39.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/40.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/42.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/43.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/44.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/45.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/46.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/47.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/48.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/49.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/50.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/51.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/5.1.bam.recal.bam &
haplotypecaller_GVCF BAMS_BIPMED/52.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/53.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/54.5.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/55.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/56.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/57.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/58.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/59.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/60.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/61.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/6.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/62.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/63.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/64.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/65.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/66.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/67.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/68.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/69.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/70.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/71.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/7.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/72.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/73.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/74.6.bam.recal.bam &
haplotypecaller_GVCF BAMS_BIPMED/75.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/76.6.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/77.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/78.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/79.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/80.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/81.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/8.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/82.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/83.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/84.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/85.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/9.1.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/94.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/95.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/96.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/97.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/98.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_BIPMED/99.7.bam.recal.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/101.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/10.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/102.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/103.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/104.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/105.7.bam &
haplotypecaller_GVCF BAMS_NO_REALN/106.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/107.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/108.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/109.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/110.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/111.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/11.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/112.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/113.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/114.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/115.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/116.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/117.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/118.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/119.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/1.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/120.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/121.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/12.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/122.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/123.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/124.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/125.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/126.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/127.8.bam &
haplotypecaller_GVCF BAMS_NO_REALN/128.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/129.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/130.8.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/131.2.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/13.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/132.2.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/14.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/15.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/16.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/17.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/18.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/19.2.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/20.2.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/21.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/2.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/22.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/23.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/24.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/25.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/26.4.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/27.4.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/28.4.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/29.4.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/296.3.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/30.4.bam &
haplotypecaller_GVCF BAMS_NO_REALN/31.4.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/3.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/32.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/33.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/34.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/35.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/36.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/37.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/38.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/39.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/40.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/42.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/43.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/44.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/45.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/46.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/47.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/48.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/49.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/50.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/51.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/5.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/52.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/53.5.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/54.5.bam &
haplotypecaller_GVCF BAMS_NO_REALN/55.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/56.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/57.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/58.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/59.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/60.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/61.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/6.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/62.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/63.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/64.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/65.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/66.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/67.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/68.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/69.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/70.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/71.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/7.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/72.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/73.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/74.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/75.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/76.6.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/77.7.bam &
haplotypecaller_GVCF BAMS_NO_REALN/78.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/79.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/80.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/81.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/8.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/82.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/83.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/84.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/85.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/9.1.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/94.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/95.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/96.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/97.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/98.7.bam &&
haplotypecaller_GVCF BAMS_NO_REALN/99.7.bam &

amostras1k=`ls *1K*g.vcf.gz | sed 's/^/-V /g' | sed 's/\n//g'`
java -jar -Xmx50G /opt/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF --max_alternate_alleles 100 $amostras1k -o /iscsi/murilo/realignment_article/HELENA.1K.vcf &
amostrasNO_REALN=`ls *NO_REALN*g.vcf.gz | sed 's/^/-V /g' | sed 's/\n//g'`
java -jar -Xmx50G /opt/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF --max_alternate_alleles 100 $amostrasNO_REALN -o /iscsi/murilo/realignment_article/HELENA.NO_REALN.vcf &
amostrasBIPMED=`ls *BIPMED*g.vcf.gz | sed 's/^/-V /g' | sed 's/\n//g'`
java -jar -Xmx50G /opt/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF --max_alternate_alleles 100 $amostrasBIPMED -o /iscsi/murilo/realignment_article/HELENA.BIPMED.vcf &

cd /iscsi/murilo/realignment_article/HELENA/NOFILTER

for VCF in `ls *vcf`
do
java -jar /opt/snpEff/SnpSift.jar filter "(DP>10)" -a "(QUAL>30)" $VCF > ../$VCF.DP.vcf &&
vcftools --vcf ../$VCF.DP.vcf --recode --out ../$VCF.DP.GENE &&
mv ../$VCF.DP.GENE.recode.vcf ../$VCF.DP.QUAL.vcf && rm ../$VCF.DP.vcf &
done




cd /iscsi/murilo/realignment_article/HELENA/FILTERED

for VCF in `ls *vcf`
do
grep -v ^# $VCF | cut -f 1,2,4,5 | sed 's/\t/_/g' > $VCF.pos &
done

