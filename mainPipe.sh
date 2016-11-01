#!/bin/bash

# Sequencing alignment pipeline

module load seq/cutadapt/1.11


# put all the scripts that pipeline uses into a folder and cd into it
# cd Codes/NETseqAlignment/
param=$1
adapter=$2

Samples=`sed -n "/<Sample Names:/,/>/p" $param | sed '$ d' | sed '/<.*/d'`
Notification=`sed -n "/<Notification Email/,/>/p" $param | sed '$ d' | sed '/<.*/d'`
User=`sed -n "/<Orchestra User ID/,/>/p" $param | sed '$ d' | sed '/<.*/d'`
indexDir=`sed -n "/<Index Directory/,/>/p" $param | sed '$ d' | sed '/<.*/d'`
initialFilesDir=`sed -n "/<Initial Files Directory/,/>/p" $param | sed '$ d' | sed '/<.*/d'`
projectName=`sed -n "/<Project Name:/,/>/p" $param | sed '$ d' | sed '/<.*/d'`

indexMatch=`sed -n "/<Index Match:/,/>/p" $param | sed '$ d' | sed '/<.*/d'`

baseDir="/groups/winston/ue4/${projectName}"



echo "baseDirectory: " $baseDir
echo "initial Dir: " $initialFilesDir

indexList=`echo -e $indexMatch`

#loopStart,g
for g in $indexList
do
    f=${g%-*}
    echo "Doing file "$f

    ### Cut the adapter sequence ###

    mkdir -p ${baseDir}/postCleaning/${f}/LogErr
    mkdir -p ${baseDir}/postCleaning/${f}/Removed

    outDir=${baseDir}/postCleaning/${f}
    preout1=${outDir}/${f}_noAdaptR_1.fastq
    bad1=${outDir}/Removed/${f}_removed_1

    Adapter=`less ${adapter}`
    echo "Adapter: " $Adapter
    prinseq=/groups/churchman/jd187/Program/PrInSeq/prinseq-lite-0.20.2/prinseq-lite.pl

    #@1,0,cutadapt: cut the adapter sequence
    cutadapt -f fastq -a $Adapter -O 3 -m 1 --error-rate=0.21 \
    --length-tag 'length=' -o ${preout1} ${initialFilesDir}/${f}.fastq

    # trim_right 1 to remove the A added by RT superscipt polymerase (when the fragment is smaller than the read length)

    #@2,1,clean: clean the fastq
    perl ${prinseq} -fastq ${preout1} \
    -out_good ${outDir}/${f}_cleaned -out_bad ${bad1} \
    -no_qual_header -min_len 7 -min_qual_mean 20 -trim_right 1 \
    -trim_ns_right 1 -trim_qual_right 20 -trim_qual_type mean -trim_qual_window 3 -trim_qual_step 1

    #@3,2,barcodeXtract: extract molecular barcode
    python /groups/churchman/jd187/NETseq/script/extractMolecularBarcode.py \
    ${outDir}/${f}_cleaned.fastq \
    ${outDir}/${f}_cleaned_noBarcode.fastq \
    ${outDir}/${f}_barcodeDistribution.txt \
    ${outDir}/${f}_ligationDistribution.txt


    ### Align reads without barcode ###

    mkdir -p ${baseDir}/TopHat2/${f}/LogErr
    outDir=${baseDir}/TopHat2/${f}

    index=${indexDir}/${g#*-}
    reads1=${baseDir}/postCleaning/${f}/${f}_cleaned_noBarcode.fastq

    seg=20

    #@4,3,tophat_no_barcode: align without barcode
    tophat2 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 -o $outDir --min-anchor-length 8 \
    --max-insertion-length 3 \
    --max-deletion-length 3 --num-threads 4 --max-multihits 100 \
    --library-type fr-firststrand --segment-mismatches 2 --no-coverage-search\
    --segment-length ${seg} \
    --b2-sensitive \
     $index ${reads1}


    ### Remove PCR duplicates ###
    BAMdir=${baseDir}/TopHat2/${f}

    #@5,4,removePCRdups: remove PCR duplicates
    python removePCRdupsFromBAM.py ${BAMdir}/accepted_hits.bam ${BAMdir}/accepted_hits_noPCRdup.bam

    #@6,5,sort_bam: sort bam file
    samtools sort ${BAMdir}/accepted_hits_noPCRdup.bam ${BAMdir}/accepted_hits_noPCRdup_sorted


    ### Calculate coverage ###
    mkdir -p ${baseDir}/Coverage/${f}

    script="/groups/churchman/ue4/Scripts/customCoverage.py"

    #@7,6,coverage: calculate coverage
    samtools view -q 50 ${baseDir}/TopHat2/${f}/accepted_hits_noPCRdup_sorted.bam | \
    python $script ${baseDir}/Coverage/${f}/${f}_TH

#loopEnd
done
