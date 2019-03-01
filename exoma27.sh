## Processing for exome data
## This uses only single-subject calls
## By Benilton Carvalho - Nov/14
## Altered by Murilo - Jun/16
BATCH=exoma_27
THREADS=10


DATAPATH=/home/bioinf/exoma/exoma_27/raw
export TARGETS=/home/murilo/Dropbox/Nextera/Restrict/HG19/nexterarapidcapture_exome_targetedregions.list
export BAITS=/home/murilo/Dropbox/Nextera/Restrict/HG19/NexteraRapidCapture_Exome_Probes_v1.1.list


export REF="/home/bioinf/ref/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
export PICARD="java -jar -Xmx20G /opt/picard-tools-2.5.0/picard.jar"
export FASTQC=fastqc
export BWA="bwa"
export GATK="java -jar /opt/GenomeAnalysisTK.jar"
export KNOWN1="/home/benilton/gatkbundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
export KNOWN2="/home/benilton/gatkbundle/1000G_phase1.indels.hg19.sites.vcf"
export VEP="/home/murilo/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl"

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
    ${PICARD} SortSam INPUT=${samout} OUTPUT=${bamout}.withdups.bam SO=coordinate &&
    rm ${samout}
    ${PICARD} MarkDuplicates INPUT=${bamout}.withdups.bam OUTPUT=${bamout} METRICS_FILE=${bamout}.metrics &&
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

    ${GATK} -T RealignerTargetCreator -R ${REF} -I ${INPUTBAM} -o ${INTERVALS} -L ${TARGETS} -ip 200
    ${GATK} -T IndelRealigner -R ${REF} -I ${INPUTBAM} -targetIntervals ${INTERVALS} -o ${REALNBAM} -compress 0 --disable_bam_indexing -L ${TARGETS} -ip 200 &&
    rm ${INPUTBAM}
    ${PICARD} BuildBamIndex INPUT=${REALNBAM}
    ${GATK} -T BaseRecalibrator -R ${REF} -I ${REALNBAM} -o ${RECALCSV} -knownSites ${KNOWN1} -knownSites ${KNOWN2} -L ${TARGETS} -ip 200
    ${GATK} -T PrintReads -R ${REF} -I ${REALNBAM} -BQSR ${RECALCSV} -o ${RECALBAM} &&
    rm ${REALNBAM}
    ${PICARD} BuildBamIndex INPUT=${RECALBAM}
    ${GATK} -T HaplotypeCaller -R ${REF} -I ${RECALBAM} -o ${VCF} -L ${TARGETS} -ip 200 -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 
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
     filename=`echo $(basename ${r1}) | awk -F '.R' '{print $1}'`

     cd $2
     # Considerando que geralmente os transposomas tem 12pb, temos que se temos reads de 100pb temos que 100-12= 88pb
     ~murilo/trim_galore_zip/trim_galore --gzip --length 88 --paired $r1 $r2 --path_to_cutadapt /home/murilo/cutadapt-1.10/bin/cutadapt

     mv $filename".R1_val_1.fq.gz" $filename"_R1_trim_1.fastq.gz"
     mv $filename".R2_val_2.fq.gz" $filename"_R2_trim_2.fastq.gz"

     mv $filename*_trimming_report.txt ../qc
     cd ..
}

function VEP
{
	file=$1
	name=`echo $file | awk -F / '{print $2}' | awk -F . '{print $1}'`
	perl $VEP -i $file -o results/$name.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --gmaf --variant_class --vcf --port 3337 --force_overwrite --merged --offline --assembly GRCh37 --fork 7
}

export -f picard_hsmetrics
export -f run_fastqc
export -f bwa_align_md
export -f run_gatk
export -f findadapters
export -f trimadapters
export -f VEP

mkdir -p ${BATCH}/qc/aln ${BATCH}/bams ${BATCH}/results 

cd ${BATCH}

FQ1=$(ls ${DATAPATH}/*R1*gz)
FQ2=$(ls ${DATAPATH}/*R2*gz)

parallel --xapply -j ${THREADS} run_fastqc ::: ${FQ1} ::: ${FQ2} ::: qc
parallel -j ${THREADS} findadapters ::: ${FQ1} ::: qc

THREADS=6

parallel -j ${THREADS} trimadapters ::: ${FQ1} ::: raw

FQ1=$(ls ${DATAPATH}/*R1*trim*gz)
FQ2=$(ls ${DATAPATH}/*R2*trim*gz)

THREADS=10

parallel --xapply -j ${THREADS} run_fastqc ::: ${FQ1} ::: ${FQ2} ::: qc
parallel -j ${THREADS} findadapters ::: ${FQ1} ::: qc

THREADS=2
parallel --xapply -j 5 bwa_align_md ::: ${FQ1} ::: ${FQ2} ::: ${THREADS} ::: bams

THREADS=10
parallel -j ${THREADS} run_gatk ::: $(ls bams/*bam)

parallel -j ${THREADS} picard_hsmetrics ::: $(ls bams/*bam) ::: qc/aln



myfiles=`ls results/*-*.g.vcf`
gVCFlist=`echo $myfiles | sed 's/ / -V /g' | sed 's/^/-V /g'`

$GATK -T GenotypeGVCFs -R $REF -o results/SMS-RGA.vcf $gVCFlist --max_alternate_alleles 12 &


for file in `ls results/*_*.g.vcf`
do
	name=`echo $file | awk -F / '{print $2}' | awk -F . '{print $1}'`
	$GATK -T GenotypeGVCFs -R $REF -o results/$name.vcf -V $file --max_alternate_alleles 12 &
done

gzip results/*g.vcf

parallel -j 1 VEP ::: $(ls results/*.vcf)

