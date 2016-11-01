#!/bin/sh
xsub="bsub -q short -W 5:30 -u umuteser@gmail.com -N"
stamp=$(date -d "today" +"%Y%m%d%H%M")
mkdir -p flag
if [ -f flag/alljobs.jid ]; then
    checkJobs flag/alljobs.jid 
    [ $? == 1 ] && exit 0;
fi
cwd=`realpath ./flag`
rm flag/*.redo flag/*.jid flag/*.killed 2>/dev/null
declare -A nmap
declare -A mmap
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

    echo step: 1, depends on: 0, job name: cutadapt, flag: cutadapt.$g  
    flag=1.0.cutadapt.$g
    flag=${flag//\//_}
    deps=""
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " cutadapt -f fastq -a $Adapter -O 3 -m 1 --error-rate=0.21  --length-tag 'length=' -o ${preout1} ${initialFilesDir}/${f}.fastq")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[1]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[0]="no"
        [ -z ${startNewLoop[1]} ] && jobIDs[1]="" && startNewLoop[1]="no" && echo starting new loop for 1! 
        jobID[1]=$id
        jobIDs[1]=${jobIDs[1]}.$id
        depes[1]=$flag.jid
        depesm[1]="${depesm[1]} $flag.jid"
        depes[1]="${depes[1]} ${depes[0]}"
        depesm[1]="${depesm[1]} ${depesm[0]}"
        sorted_unique_ids=$(echo "${depes[1]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    # trim_right 1 to remove the A added by RT superscipt polymerase (when the fragment is smaller than the read length)
    #@2,1,clean: clean the fastq

    echo step: 2, depends on: 1, job name: clean, flag: clean.$g  
    flag=2.1.clean.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[1]}"
    tempids=""
    tempids="$tempids ${depes[1]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " perl ${prinseq} -fastq ${preout1}  -out_good ${outDir}/${f}_cleaned -out_bad ${bad1}  -no_qual_header -min_len 7 -min_qual_mean 20 -trim_right 1  -trim_ns_right 1 -trim_qual_right 20 -trim_qual_type mean -trim_qual_window 3 -trim_qual_step 1")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[2]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[1]="no"
        [ -z ${startNewLoop[2]} ] && jobIDs[2]="" && startNewLoop[2]="no" && echo starting new loop for 2! 
        jobID[2]=$id
        jobIDs[2]=${jobIDs[2]}.$id
        depes[2]=$flag.jid
        depesm[2]="${depesm[2]} $flag.jid"
        depes[2]="${depes[2]} ${depes[1]}"
        depesm[2]="${depesm[2]} ${depesm[1]}"
        sorted_unique_ids=$(echo "${depes[2]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    #@3,2,barcodeXtract: extract molecular barcode

    echo step: 3, depends on: 2, job name: barcodeXtract, flag: barcodeXtract.$g  
    flag=3.2.barcodeXtract.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[2]}"
    tempids=""
    tempids="$tempids ${depes[2]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " python /groups/churchman/jd187/NETseq/script/extractMolecularBarcode.py  ${outDir}/${f}_cleaned.fastq  ${outDir}/${f}_cleaned_noBarcode.fastq  ${outDir}/${f}_barcodeDistribution.txt  ${outDir}/${f}_ligationDistribution.txt")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[3]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[2]="no"
        [ -z ${startNewLoop[3]} ] && jobIDs[3]="" && startNewLoop[3]="no" && echo starting new loop for 3! 
        jobID[3]=$id
        jobIDs[3]=${jobIDs[3]}.$id
        depes[3]=$flag.jid
        depesm[3]="${depesm[3]} $flag.jid"
        depes[3]="${depes[3]} ${depes[2]}"
        depesm[3]="${depesm[3]} ${depesm[2]}"
        sorted_unique_ids=$(echo "${depes[3]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    ### Align reads without barcode ###
    mkdir -p ${baseDir}/TopHat2/${f}/LogErr
    outDir=${baseDir}/TopHat2/${f}
    index=${indexDir}/${g#*-}
    reads1=${baseDir}/postCleaning/${f}/${f}_cleaned_noBarcode.fastq
    seg=20
    #@4,3,tophat_no_barcode: align without barcode

    echo step: 4, depends on: 3, job name: tophat_no_barcode, flag: tophat_no_barcode.$g  
    flag=4.3.tophat_no_barcode.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[3]}"
    tempids=""
    tempids="$tempids ${depes[3]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " tophat2 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 -o $outDir --min-anchor-length 8  --max-insertion-length 3  --max-deletion-length 3 --num-threads 4 --max-multihits 100  --library-type fr-firststrand --segment-mismatches 2 --no-coverage-search --segment-length ${seg}  --b2-sensitive  $index ${reads1}")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[4]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[3]="no"
        [ -z ${startNewLoop[4]} ] && jobIDs[4]="" && startNewLoop[4]="no" && echo starting new loop for 4! 
        jobID[4]=$id
        jobIDs[4]=${jobIDs[4]}.$id
        depes[4]=$flag.jid
        depesm[4]="${depesm[4]} $flag.jid"
        depes[4]="${depes[4]} ${depes[3]}"
        depesm[4]="${depesm[4]} ${depesm[3]}"
        sorted_unique_ids=$(echo "${depes[4]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    ### Remove PCR duplicates ###
    BAMdir=${baseDir}/TopHat2/${f}
    #@5,4,removePCRdups: remove PCR duplicates

    echo step: 5, depends on: 4, job name: removePCRdups, flag: removePCRdups.$g  
    flag=5.4.removePCRdups.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[4]}"
    tempids=""
    tempids="$tempids ${depes[4]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " python removePCRdupsFromBAM.py ${BAMdir}/accepted_hits.bam ${BAMdir}/accepted_hits_noPCRdup.bam")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[5]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[4]="no"
        [ -z ${startNewLoop[5]} ] && jobIDs[5]="" && startNewLoop[5]="no" && echo starting new loop for 5! 
        jobID[5]=$id
        jobIDs[5]=${jobIDs[5]}.$id
        depes[5]=$flag.jid
        depesm[5]="${depesm[5]} $flag.jid"
        depes[5]="${depes[5]} ${depes[4]}"
        depesm[5]="${depesm[5]} ${depesm[4]}"
        sorted_unique_ids=$(echo "${depes[5]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    #@6,5,sort_bam: sort bam file

    echo step: 6, depends on: 5, job name: sort_bam, flag: sort_bam.$g  
    flag=6.5.sort_bam.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[5]}"
    tempids=""
    tempids="$tempids ${depes[5]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " samtools sort ${BAMdir}/accepted_hits_noPCRdup.bam ${BAMdir}/accepted_hits_noPCRdup_sorted")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[6]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[5]="no"
        [ -z ${startNewLoop[6]} ] && jobIDs[6]="" && startNewLoop[6]="no" && echo starting new loop for 6! 
        jobID[6]=$id
        jobIDs[6]=${jobIDs[6]}.$id
        depes[6]=$flag.jid
        depesm[6]="${depesm[6]} $flag.jid"
        depes[6]="${depes[6]} ${depes[5]}"
        depesm[6]="${depesm[6]} ${depesm[5]}"
        sorted_unique_ids=$(echo "${depes[6]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
    ### Calculate coverage ###
    mkdir -p ${baseDir}/Coverage/${f}
    script="/groups/churchman/ue4/Scripts/customCoverage.py"
    #@7,6,coverage: calculate coverage

    echo step: 7, depends on: 6, job name: coverage, flag: coverage.$g  
    flag=7.6.coverage.$g
    flag=${flag//\//_}
    deps=""
    deps="${jobID[6]}"
    tempids=""
    tempids="$tempids ${depes[6]}"
    sorted_unique_ids=$(echo "$tempids" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
    for d in ""$sorted_unique_ids""; do 
        [[ -f $cwd/${d%.*}.redo && -f $cwd/$flag.success ]] && rm $cwd/$flag.success 
    done 
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(bsubRun $xsub -flag $deps $cwd $flag " samtools view -q 50 ${baseDir}/TopHat2/${f}/accepted_hits_noPCRdup_sorted.bam |  python $script ${baseDir}/Coverage/${f}/${f}_TH")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[7]=""
    else
        alljobs="$alljobs $id"
        nmap[$id]=$flag
        mmap[$id]=$deps
        startNewLoop[6]="no"
        [ -z ${startNewLoop[7]} ] && jobIDs[7]="" && startNewLoop[7]="no" && echo starting new loop for 7! 
        jobID[7]=$id
        jobIDs[7]=${jobIDs[7]}.$id
        depes[7]=$flag.jid
        depesm[7]="${depesm[7]} $flag.jid"
        depes[7]="${depes[7]} ${depes[6]}"
        depesm[7]="${depesm[7]} ${depesm[6]}"
        sorted_unique_ids=$(echo "${depes[7]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') 
        for d in ""$sorted_unique_ids""; do 
           echo $id ${mmap[$id]} ${nmap[$id]} >> $cwd/$d
        done 
    fi
#loopEnd
done
cd $cwd/..
echo all submitted jobs: 
printf "%-10s   %-20s   %-10s\n" job_id depend_on job_flag
echo ---------------------------------------------------------
for i in $alljobs; do
    printf "%-10s | %-20s | %-10s\n" $i ${mmap[$i]} ${nmap[$i]}
done
printf "%-10s   %-20s   %-10s\n" job_id depend_on job_flag > flag/alljobs.jid
echo ---------------------------------------------------------
for i in $alljobs; do
    printf "%-10s  %-20s  %-10s\n" $i ${mmap[$i]} ${nmap[$i]} >> flag/alljobs.jid
done
echo bjobs -w output:
bjobs -w
