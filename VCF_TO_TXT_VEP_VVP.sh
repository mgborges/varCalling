function VEP_VVP
{
	VCF=$1
	POS=$2	
	/home/murilo/VVP-pub/vt/vt decompose -s $VCF -o $VCF.decomp.vcf

	vcftools --vcf $VCF.decomp.vcf --bed $POS --recode --min-meanDP 20 --minQ 30 --out $VCF
	mv $VCF.recode.vcf $VCF.decomp.vcf

	/home/murilo/ensembl-vep/vep -i $VCF.decomp.vcf -o $VCF.decomp.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite --vcf --merged

	/home/murilo/ensembl-vep/vep -i $VCF.decomp.vcf -o $VCF.decomp.VEP --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite  --merged

	cat $VCF.decomp.VEP.vcf | sed 's/^chr//' > $VCF.tmp
	rm $VCF.decomp.vcf

	/home/murilo/VVP-pub/VVP -i $VCF.tmp -d /home/murilo/VVP-pub/gnomad.062717.build -v CSQ,4,6,1,15 1> $VCF.vvp.out

	cat $VCF.vvp.out | sed 's/\./,/g' > $VCF.vvp.out.xls
	rm $VCF.tmp $VCF.vvp.out 
}

function VEPtoTable
{
	VCF=$1
	VEP=`echo $VCF | sed 's/.vcf$//g'`
	grep -v ^# $VEP > $VCF.CABECALHO
	VVP=`echo $VCF | sed 's/.decomp.VEP.vcf$/.vvp.out.xls/g'`

	java -jar /home/murilo/snpEff/SnpSift.jar extractFields $VCF CHROM POS ALT REF "GEN[*].GT" "GEN[*].AD" -e "." -s ";" > $VCF.tmp

	while read -r line
	do
	IMPACT=`echo $line | awk -F "IMPACT=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$IMPACT" == "" ]
	then
	IMPACT="."
	fi
	SYMBOL=`echo $line | awk -F "SYMBOL=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$SYMBOL" == "" ]
	then
	SYMBOL="."
	fi
	BIOTYPE=`echo $line | awk -F "BIOTYPE=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$BIOTYPE" == "" ]
	then
	BIOTYPE="."
	fi
	GMAF=`echo $line | awk -F "GMAF=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$GMAF" == "" ]
	then
	GMAF="."
	fi
	SIFT=`echo $line | awk -F "SIFT=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$SIFT" == "" ]
	then
	SIFT="."
	fi
	PolyPhen=`echo $line | awk -F "PolyPhen=" '{print $2}' | awk -F ";" '{print $1}'`
	if [ "$PolyPhen" == "" ]
	then
	PolyPhen="."
	fi
	otherinfo=`echo $line | cut -f 1-13 -d " "`
	TRANSCRITO=`echo $otherinfo | cut -f 5 -d " "`

	CHR=`echo $line | cut -f 2 -d " " | cut -f 1 -d ":"`
	POS=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | cut -f 1 -d "-"`
	ALL=`echo $line | cut -f 3 -d " "`


	verificator=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | grep -c "-"`

	if [ "$verificator" != "0" ] && [ "$ALL" != "-" ]
	then
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	elif [ "$verificator" != "0" ] && [ "$ALL" == "-" ]
	then
	POS=`bc <<< "$POS - 1"`
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	elif [ "$verificator" == "0" ] && [ "$ALL" == "-" ]
	then
	POS=`bc <<< "$POS - 1"`
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	else
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	fi

	REF_TECIDO=`echo $DP | cut -f 1 -d ';' | cut -f 1 -d ','`
	ALT_TECIDO=`echo $DP | cut -f 1 -d ';' | cut -f 2 -d ','`
	REF_SANGUE=`echo $DP | cut -f 2 -d ';' | cut -f 1 -d ','`
	ALT_SANGUE=`echo $DP | cut -f 2 -d ';' | cut -f 2 -d ','`

	MOSAICO_TECIDO=`bc <<< "$ALT_TECIDO * 100 / ($REF_TECIDO + $ALT_TECIDO)"`
	MOSAICO_SANGUE=`bc <<< "$ALT_SANGUE * 100 / ($REF_SANGUE + $ALT_SANGUE)"`

	cobertura_TECIDO=`bc <<< "$REF_TECIDO + $ALT_TECIDO"`
	cobertura_SANGUE=`bc <<< "$REF_SANGUE + $ALT_SANGUE"`

	VVP_SCORE=`grep $TRANSCRITO $VVP | grep $POS -m 1 | cut -f 8,13,18 | sed 's/\t/ /g'`

	echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $VVP_SCORE $GT $DP $cobertura_TECIDO $cobertura_SANGUE

	done < $VCF.CABECALHO
	rm $VCF.tmp $VCF.CABECALHO
}

function EXECUTE_VEP_VVP
{
	vcf=$1
	POS=$2	
	PATHWAY=$3
	VEP_VVP $vcf $POS && 
	VEPtoTable $vcf.decomp.VEP.vcf | grep -e shift -e missense_variant -e nonsense -e splice -e stop > $vcf.$PATHWAY.tmp && 
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_VVP.header $vcf.$PATHWAY.tmp > $vcf.$PATHWAY.filtered.xls && 
	rm $vcf.$PATHWAY.tmp $vcf*decomp* $vcf*vvp.out.xls &&
	echo DONE FOR $vcf
}

# Quanto aos filtros, pode aplicar para as consequencias: frame-shift, missense, nonsense, splicing site, stop-codon
POS=/home/murilo/Dropbox/Nextera/GATOR/gator.pos
PATHWAY="GATOR"

for vcf in `ls *vcf`
do
EXECUTE_VEP_VVP $vcf $POS $PATHWAY 
done

POS=/home/murilo/Dropbox/sh_scripts/MTOR_hg38.bed
PATHWAY="MTOR"

for vcf in `ls *vcf`
do
EXECUTE_VEP_VVP $vcf $POS $PATHWAY 
done

POS=/home/murilo/Dropbox/sh_scripts/NEUROG2.bed
PATHWAY="NEUROG2"

for vcf in `ls *vcf`
do
EXECUTE_VEP_VVP $vcf $POS $PATHWAY 
done

####### WHOLE EXOME

POS=/home/murilo/Dropbox/sh_scripts/GENOME.bed
PATHWAY="ALL"

for vcf in `ls *vcf`
do
EXECUTE_VEP_VVP $vcf $POS $PATHWAY 
done
