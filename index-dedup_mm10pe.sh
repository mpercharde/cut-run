#!/bin/bash

###################################################
#       CnR analysis pipeline                     #
#       Michelle Percharde, PhD 2018              #
#                                                 #
#                v Mm10 **PE                      #
###################################################

#~~~~~~~~~~EDIT BELOW THIS LINE ~~~~~~~~~~~~~~~~~~#

# ONLY DEDUP AND INDEXING - REQUIRES SORTED BAM FILES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flagcheck=0

while getopts ':gchi:' flag; do
    case ${flag} in
      i) dir=$OPTARG
        flagcheck=1 ;;
        g) gz='true' ;;
        c) clontech='true' ;;
        h) echo ""
           echo "Usage: $0 [-h] [-g] [-c] [-i <path/to/files/>]"
           echo ""
           echo "    -h        Help mode, prints Usage"
           echo "    -g        Input FASTQ is .gz compressed"
           echo "    -c        Libraries are CLONTECH"
           echo "    -i        input file directory/. use "./" for current dir (not recommended)"
           echo ""
           flagcheck=1
           exit 1 ;;
        \?)
          echo ""
          echo "Invalid option, type -h for help"
          echo ""
          flagcheck=1
          exit 1 ;;
    esac
done

if [ "$flagcheck" == 0 ] ; then
  echo ""
  echo "Incorrect or no -i specified, please read help -h and try again!"
  echo ""
  exit 1
fi

if [ "$gz" == "true" ]; then
  echo ""
  echo "your file is compressed, trim and align will be run on .gz files"
fi

mkdir -p trimmed/fastqc/
mkdir -p sorted_bam/
mkdir -p alignment_summaries/
# echo "indexing the transcriptome data one time - prefix hg19" ###DO THIS FIRST TIME ###
# tophat -G /data/refs/hg19/genes.gtf --transcriptome-index=transcriptome_data/known_genes /data/refs/hg19/genome

# echo "do test alignment start"
# tophat -o test_aligned/ -p 8 -g 20 --no-coverage-search --library-type fr-unstranded \
# --no-novel-indels --transcriptome-index=transcriptome_data/known_genes /data/refs/hg19/genome \
# trimmed/164_fslg1_1_val_1.fq.gz trimmed/164_fslg1_2_val_2.fq.gz

for file in "$dir"*_1* ; do
    echo ""
    if [ "$gz" == "true" ]; then
      base=$(basename $file _1.fq.gz)
      name1=$(basename $file .fq.gz)
    	name2=$(basename $file _1.fq.gz)_2
    	file1=$file
    	file2="$dir"${name2}.fq.gz
      trimfile1=${name1}_val_1.fq.gz
      trimfile2=${name2}_val_2.fq.gz
    else
      base=$(basename $file _1.fq)
      name1=$(basename $file .fq)
    	name2=$(basename $file _1.fq)_2
    	file1=$file
    	file2="$dir"${name2}.fq
      trimfile1=${name1}_val_1.fq
      trimfile2=${name2}_val_2.fq #needs to be _val_1 or _val_2
    fi
    echo "analysing file: $name1, $name2"
    echo "trimmed will be $trimfile1, $trimfile2"
    echo "file1 is $file1, file2 is $file2"

    if [ "$clontech" == "true" ]; then
      echo ""
      echo "$name1 and $name2 are from a clontech library, adaptor trimming will include clipping 3bp"
      echo ""
      trim="trim_galore --clip_R1 3"
    else
      echo ""
      echo "$name1 and $name2 are an nebnext library"
      echo ""
      trim="trim_galore"
    fi

    echo ""
    echo "#######################################"
    echo "## Cut&Run-seq pipeline - mm10 PE    ##"
    echo "##   by Michelle Percharde, PhD      ##"
    echo "##   mpercharde@gmail.com            ##"
    echo "#######################################"
    # echo "1. trimming paired files for $base"
    echo ""
    # $trim --fastqc --fastqc_args " --outdir trimmed/fastqc/" --illumina --paired -o trimmed/ $file1 $file2

    # echo ""
    # echo "2. aligning $base to mm10 and sorting in one"
    # echo ""
    # mkdir -p ${base}_aligned/

    # (bowtie2 -p8 -x /data/refs/mm10/mm10 -1 trimmed/$trimfile1 -2 trimmed/$trimfile2 \
    # --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
    # | samtools view -Suo - - | samtools sort - -o sorted_bam/${base}.sorted.bam) 2> alignment_summaries/${base}_alignment.txt
    #
    # echo "2. aligning $base to yeast repeat-masked genome and sorting in one"
    # (bowtie2 -p8 -x /data/refs/scR64/SC_rm -1 trimmed/$trimfile1 -2 trimmed/$trimfile2 \
    # --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
    # | samtools view -Suo - - | samtools sort - -o sorted_bam/${base}_SC.sorted.bam) 2> alignment_summaries/${base}_SC_alignment.txt

    # tophat -o ${base}_aligned/ -p 8 -g 20 --no-coverage-search --library-type fr-unstranded \
    # --no-novel-indels --transcriptome-index=/data/refs/hg19/transcriptome_data/known_genes /data/refs/hg19/genome trimmed/$trimfile1 trimmed/$trimfile2
    # echo ""
    # echo "3. sorting $base bam file ready for DL"
    # echo ""
    # samtools sort -o sorted_bam/${base}.sorted.bam ${base}_aligned/accepted_hits.bam
    # echo "cleaning up: removing raw files and removing non-sorted file"
    # rm ${base}_aligned/accepted_hits.bam
    # rm $file1
    # rm $file2

   echo ""
   echo "3. deduplicating $base bam file"
   echo ""
   mkdir -p dedup_bam
   samtools rmdup sorted_bam/${base}.sorted.bam dedup_bam/${base}.sorted.dedup.bam

   echo ""
   echo "4. Generating bam index for $base"
   echo ""
   samtools index sorted_bam/${base}.sorted.bam sorted_bam/${base}.sorted.bai
   samtools index dedup_bam/${base}.sorted.dedup.bam dedup_bam/${base}.sorted.dedup.bai
    echo ""

done
