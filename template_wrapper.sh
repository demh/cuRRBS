#!/bin/bash

# SOMATYPUS: A PLATYPUS-BASED VARIANT CALLING PIPELINE FOR CANCER DATA
# Adrian Baez-Ortega, Transmissible Cancer Group, University of Cambridge
# 2016

# somatypus
# Core pipeline script

# INPUT
# -i: path to folder containing the input BAM files (accompanied by BAI indices)
# -g: path to genome FASTA (must have FAI index)
# -o: path to global output folder
# -r: path to file of regions in chr:start-end format (optional)
# -c: number of CPUs (processes) for Platypus (*SHOULD NOT EXCEED 8 DUE TO BUG*) (optional)
# -e: extra options for Platypus (within quotes, separated by spaces) (optional)



VERSION=1.3


####################################### FUNCTIONS #######################################

# AUXILIARY FUNCTIONS
# print_help()
# Prints a help guide if -h, or no arguments, are input
print_help() {
    echo
    echo
    echo "| SOMATYPUS"
    echo "| A Platypus-based variant calling pipeline for cancer data"
    echo "| Version $VERSION"
    echo "|"    
    echo "| Required input:"
    echo "|    -i  Absolute path to folder containing the input BAM files (accompanied by BAI indices)."
    echo "|    -g  Absolute path to reference genome FASTA file (accompanied by FAI index)."
    echo "|    -o  Absolute path to the output folder (it will be created if needed)."
    echo "|"
    echo "| Optional input:"
    echo "|    -r  Absolute path to file of regions to use, one per line in CHR:START-END format."
    echo "|    -c  Number of CPUs (processes) for Platypus *(should not exceed 8 due to a bug)*."
    echo "|    -p  Additional options for Platypus, within quotes and separated by spaces."
    echo "|"
    echo "| Options:"
    echo "|    -h  Print this usage information and exit."
    echo "|    -v  Print version and exit."
    echo "|"
    echo "| Usage:"
    echo "|    somatypus -i /path/to/bams_dir -o /path/to/out_dir -g /path/to/genome.fna -r /path/to/regions.txt -c <1-8> -p \"--option=VAL --option=VAL\""
    echo
    echo
}


# check_file()
# Checks if a file exists and is not empty. In that case, it displays an error message and exits
# Used for checking the output of each step
check_file() {

    if [ ! -s $1 ]; then
        echo -e "\nERROR: Output file $1 was not correctly generated. Please check the logs folder for more information.\n" >&2
        exit 1
    fi
    
}


# PIPELINE STEPS
# 1-2) individual_calling()
# Calls variants using Platypus from each sample in a sample set, individually
# INPUT: $1 - Settings to use for calling: default (0) or alternative (1)
individual_calling() {

    # If the user input a regions file: include the --regions argument
    REGIONSARG=""
    if [ "$REGIONS" != "no" ]; then
        REGIONSARG="--regions=$REGIONS"
    fi

    # Create directories for output files and logs
    mkdir -p $OUTDIR/1-2_individual_calls
    mkdir -p $OUTDIR/logs/1_individual_default
    mkdir -p $OUTDIR/logs/2_individual_alternative

    # Run Platypus with the default/alternative settings, for each BAM file
    # Default settings:
    if [ "$1" -eq 0 ]; then
        for FILE in `ls $BAMSDIR/*.bam`; do 
            NAME=`basename $FILE`
            echo -e "\nCalling on $NAME\n"
            Platypus.py callVariants \
            --logFileName=$OUTDIR/logs/1_individual_default/"${NAME%.*}"_default.log \
            --refFile=$REFERENCE \
            --bamFiles=$FILE \
            $REGIONSARG \
            --minPosterior=0 \
            --minReads=3 \
            --nCPU=$CPUS \
            $EXTRA \
            -o $OUTDIR/1-2_individual_calls/platypusVariants_"${NAME%.*}"_default.vcf
        done

    # Alternative (minFlank=0) settings:
    else
        for FILE in `ls $BAMSDIR/*.bam`; do 
            NAME=`basename $FILE`
            echo -e "\nCalling on $NAME\n"
            Platypus.py callVariants \
            --logFileName=$OUTDIR/logs/2_individual_alternative/"${NAME%.*}"_alternative.log \
            --refFile=$REFERENCE \
            --bamFiles=$FILE \
            $REGIONSARG \
            --minPosterior=0 \
            --minReads=3 \
            --minFlank=0 \
            --trimReadFlank=10 \
            --nCPU=$CPUS \
            $EXTRA \
            -o $OUTDIR/1-2_individual_calls/platypusVariants_"${NAME%.*}"_alternative.vcf
        done
    fi
    
}


# 3) split_calls()
# Splits multi-allelic calls and MNPs in individual Platypus VCFs into bi-allelic SNVs
split_calls() {

    # Create directory for output split files
    mkdir -p $OUTDIR/3_individual_split

    # For each individual VCF, split calls
    for FILE in `ls $OUTDIR/1-2_individual_calls/platypusVariants_*`; do
        Somatypus_SplitMA-MNVs.py $FILE >> $OUTDIR/logs/3_split.log
        mv "${FILE%.*}".split.vcf $OUTDIR/3_individual_split/
    done

}


# 4) filter_calls()
# Filters Platypus SNV calls showing flags badReads, MQ, strandBias, SC or QD  
filter_calls() {

    # Create directory for output filtered VCFs
    mkdir -p $OUTDIR/4_individual_filtered

    # For each split VCF, remove calls with flags badReads, MQ, strandBias, SC or QD
    for FILE in `ls $OUTDIR/3_individual_split/platypusVariants_*`; do
        NAME=`basename $FILE`
        awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $FILE > $OUTDIR/4_individual_filtered/"${NAME%.*}".filtered.vcf
    done

}


# 5) indel_flag()
# Identifies SNVs close to indels in any sample, from a set of Platypus VCFs
indel_flag() {

    # Create directory for output indel-flagged SNVs list
    mkdir -p $OUTDIR/5-7_merged

    # Create a list of all split VCF files
    rm -f $OUTDIR/3_individual_split/list.txt
    ls -1 $OUTDIR/3_individual_split/*.split.vcf > $OUTDIR/3_individual_split/list.txt

    # Create list of indel-flagged SNVs
    Somatypus_IndelFlag.py $OUTDIR/3_individual_split/list.txt $OUTDIR/5-7_merged/indel_flagged_SNVs.txt > $OUTDIR/logs/5_indel_flag.log

}


# 6) merge_calls()
# Merges filtered Platypus SNV calls obtained from individual calling
merge_calls() {

    # Create a list of filtered individual VCFs
    ls -1 $OUTDIR/4_individual_filtered/*.filtered.vcf > $OUTDIR/4_individual_filtered/list.txt

    # Merge SNVs from all filtered VCFs
    Somatypus_SNVmerge.py $OUTDIR/4_individual_filtered/list.txt $OUTDIR/5-7_merged/indel_flagged_SNVs.txt $OUTDIR/5-7_merged > $OUTDIR/logs/6_SNV_merge.log

}


# 7) extract_indels()
# Extracts indels from Platypus VCF files obtained from individual calling
extract_indels() {

    # Create list of original individual VCFs
    # (Indel merging uses unfiltered, unsplit calls)
    ls -1 $OUTDIR/1-2_individual_calls/*.vcf > $OUTDIR/1-2_individual_calls/list.txt

    # Merge and filter indels; only bi-allelic, PASS-flagged indels are selected
    Somatypus_IndelMerge.py $OUTDIR/1-2_individual_calls/list.txt $OUTDIR/5-7_merged > $OUTDIR/logs/7_indel_merge.log

}


# 8) prepare_genotyping()
# Generates adequate VCF and region files for SNV and indel genotyping
prepare_genotyping() {

    # Sort, compress and index VCF files for inputting them to Platypus
    vcf-sort $OUTDIR/5-7_merged/MergedIndels.vcf > $OUTDIR/5-7_merged/MergedIndels.sorted.vcf 2> /dev/null
    bgzip -f $OUTDIR/5-7_merged/MergedIndels.sorted.vcf
    tabix -f -p vcf $OUTDIR/5-7_merged/MergedIndels.sorted.vcf.gz
    for FILE in `ls $OUTDIR/5-7_merged/MergedSNVs_allele?.vcf`; do
        vcf-sort $FILE > "${FILE%.*}".sorted.vcf 2> /dev/null
        bgzip -f "${FILE%.*}".sorted.vcf
        tabix -f -p vcf "${FILE%.*}".sorted.vcf.gz
    done

    # Create directory for new region files (and genotyping output)
    mkdir -p $OUTDIR/8-18_genotyped

    # If the user did not input a regions file, define regions of +/-200 bp around every variant
    if [ "$REGIONS" == "no" ]; then
        awk '{if ($2 < 200) start=0; else start=$2-200; print $1 ":" start "-" $2+200}' $OUTDIR/5-7_merged/*.vcf > $OUTDIR/8-18_genotyped/variant_regions_200bp.txt
        Somatypus_MergeRegions.py $OUTDIR/8-18_genotyped/variant_regions_200bp.txt
        REGIONS="$OUTDIR/8-18_genotyped/variant_regions_200bp_merged.txt"
    fi
    
    # Extract only regions containing variants
    #if [ "$REGIONS" != "no" ]; then
    Somatypus_ExtractRegions.py $REGIONS $OUTDIR/5-7_merged/MergedSNVs_allele1.vcf $OUTDIR/5-7_merged/MergedSNVs_allele2.vcf $OUTDIR/5-7_merged/MergedSNVs_allele3.vcf $OUTDIR/5-7_merged/MergedIndels.vcf $OUTDIR/8-18_genotyped 0 > $OUTDIR/logs/8_extract_regions.log
    #fi

    # Create list of BAM files for Platypus
    ls -1 $BAMSDIR/*.bam > $OUTDIR/8-18_genotyped/bam_list.txt

}


# 9-12) genotyping()
# Runs Platypus to genotype SNVs and indels obtained after individual calling and filtering
# INPUT: $1 - Number of the allele to be genotyped (1, 2 or 3; 0 means indels)
genotyping() {

    IND="$1"
    
    # If IND==0: genotype indels
    if [ "$IND" -eq 0 ]; then
    
        # If the user input a regions file: include the --regions argument
        REGIONSARG=""
        if [ "$REGIONS" != "no" ]; then
            REGIONSARG="--regions=$OUTDIR/8-18_genotyped/regions_indels.txt"
        fi

        # Run Platypus to genotype indels
        Platypus.py callVariants \
        --logFileName=$OUTDIR/logs/12.1_genotype_indels_first.log \
        --refFile=$REFERENCE \
        --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
        $REGIONSARG \
        --minPosterior=0 \
        --nCPU=$CPUS \
        --minReads=3 \
        --source=$OUTDIR/5-7_merged/MergedIndels.sorted.vcf.gz \
        --getVariantsFromBAMs=0 \
        $EXTRA \
        -o $OUTDIR/8-18_genotyped/GenotypedIndels_first.vcf

        # (Some calls may not be genotyped due to the way Platypus builds haplotypes)
        # Extract missing calls by comparing merged and genotyped VCFs
        tail -n +49 $OUTDIR/8-18_genotyped/GenotypedIndels_first.vcf | cut -f1,2,4,5 > $OUTDIR/geno_pos.txt
        cut -f1,2,4,5 $OUTDIR/5-7_merged/MergedIndels.vcf > $OUTDIR/merged_pos.txt
        grep -vxFf $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt > $OUTDIR/coords.txt

        # Create a new regions file containing only the bases of the missing variants
        # The size of the region is the length of the SNV/indel
        awk '{if (length($3) >= length($4)) { print $1 ":" $2 "-" $2+length($3)-1 } else { print $1 ":" $2 "-" $2+length($4)-1 }}' $OUTDIR/coords.txt > $OUTDIR/8-18_genotyped/varRegions_indels.txt
        Somatypus_MergeRegions.py $OUTDIR/8-18_genotyped/varRegions_indels.txt > $OUTDIR/logs/12.2_merge_indel_regions.log
        rm $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt $OUTDIR/coords.txt

        # If there are missing calls: run Platypus to re-genotype them
        if [ -s $OUTDIR/8-18_genotyped/varRegions_indels_merged.txt ]; then
            echo -e "\nGenotyping missing calls in indels\n"
            Platypus.py callVariants \
            --logFileName=$OUTDIR/logs/12.3_genotype_indels_second.log \
            --refFile=$REFERENCE \
            --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
            --regions=$OUTDIR/8-18_genotyped/varRegions_indels_merged.txt \
            --minPosterior=0 \
            --nCPU=$CPUS \
            --minReads=3 \
            --source=$OUTDIR/5-7_merged/MergedIndels.sorted.vcf.gz \
            --getVariantsFromBAMs=0 \
            $EXTRA \
            -o $OUTDIR/8-18_genotyped/GenotypedIndels_second.vcf
        else
            echo -e "\nNo missing calls"
        fi
        
    # SNV genotyping (allele $IND) 
    else

        # If the user input a regions file: include the --regions argument
        REGIONSARG=""
        if [ "$REGIONS" != "no" ]; then
            REGIONSARG="--regions=$OUTDIR/8-18_genotyped/regions_allele${IND}.txt"
        fi

        # Run Platypus to genotype the specified allele
        Platypus.py callVariants \
        --logFileName=$OUTDIR/logs/$(( 8 + $IND )).1_genotype_allele${IND}_first.log \
        --refFile=$REFERENCE \
        --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
        $REGIONSARG \
        --minPosterior=0 \
        --nCPU=$CPUS \
        --minReads=3 \
        --source=$OUTDIR/5-7_merged/MergedSNVs_allele${IND}.sorted.vcf.gz \
        --getVariantsFromBAMs=0 \
        $EXTRA \
        -o $OUTDIR/8-18_genotyped/GenotypedSNVs_allele${IND}_first.vcf

        # (Some calls may not be genotyped due to the way Platypus builds haplotypes)
        # Extract missing calls by comparing merged and genotyped VCFs
        tail -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_allele${IND}_first.vcf | cut -f1,2,4,5 > $OUTDIR/geno_pos.txt
        cut -f1,2,4,5 $OUTDIR/5-7_merged/MergedSNVs_allele${IND}.vcf > $OUTDIR/merged_pos.txt
        grep -vxFf $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt > $OUTDIR/coords.txt

        # Create a new regions file containing only the bases of the missing variants
        # The size of the region is the length of the SNV/indel
        awk '{if (length($3) >= length($4)) { print $1 ":" $2 "-" $2+length($3)-1 } else { print $1 ":" $2 "-" $2+length($4)-1 }}' $OUTDIR/coords.txt > $OUTDIR/8-18_genotyped/varRegions_allele${IND}.txt
        rm $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt $OUTDIR/coords.txt

        # If there are missing calls: run Platypus to re-genotype them
        if [ -s $OUTDIR/8-18_genotyped/varRegions_allele${IND}.txt ]; then
            echo -e "\nGenotyping missing calls in allele $IND\n"
            Platypus.py callVariants \
            --logFileName=$OUTDIR/logs/$(( 8 + $IND )).2_genotype_allele${IND}_second.log \
            --refFile=$REFERENCE \
            --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
            --regions=$OUTDIR/8-18_genotyped/varRegions_allele${IND}.txt \
            --minPosterior=0 \
            --nCPU=$CPUS \
            --minReads=3 \
            --source=$OUTDIR/5-7_merged/MergedSNVs_allele${IND}.sorted.vcf.gz \
            --getVariantsFromBAMs=0 \
            $EXTRA \
            -o $OUTDIR/8-18_genotyped/GenotypedSNVs_allele${IND}_second.vcf
        else
            echo -e "\nNo missing calls"
        fi
    fi
}


# 13) prepare_genotyping_indelflagged()
# Generates adequate VCF and region files for genotyping of SNVs near indels
prepare_genotyping_indelflagged() {

    # Sort, compress and index VCF files for inputting them to Platypus
    for FILE in `ls $OUTDIR/5-7_merged/IndelExcludedSNVs_allele?.vcf`; do
        vcf-sort $FILE > "${FILE%.*}".sorted.vcf 2> /dev/null
        bgzip -f "${FILE%.*}".sorted.vcf
        tabix -f -p vcf "${FILE%.*}".sorted.vcf.gz
    done

    # If the user input a regions file:
    # Extract only regions containing variants
    if [ "$REGIONS" != "no" ]; then
        Somatypus_ExtractRegions.py $REGIONS $OUTDIR/5-7_merged/IndelExcludedSNVs_allele1.vcf $OUTDIR/5-7_merged/IndelExcludedSNVs_allele2.vcf $OUTDIR/5-7_merged/IndelExcludedSNVs_allele3.vcf none $OUTDIR/8-18_genotyped 1 > $OUTDIR/logs/13_extract_regions_excluded.log
    fi

    # Create list of BAM files for Platypus
    ls -1 $BAMSDIR/*.bam > $OUTDIR/8-18_genotyped/bam_list.txt

}


# 14-16) genotyping_indelflagged()
# Runs Platypus to genotype SNVs excluded for being close to indels
# INPUT: $1 - Number of the allele to be genotyped (1, 2 or 3)
genotyping_indelflagged() {

    IND="$1"

    # If the user input a regions file: include the --regions argument
    REGIONSARG=""
    if [ "$REGIONS" != "no" ]; then
        REGIONSARG="--regions=$OUTDIR/8-18_genotyped/regions_allele${IND}_indelExcluded.txt"
    fi

    # Run Platypus to genotype the specified allele
    Platypus.py callVariants \
    --logFileName=$OUTDIR/logs/$(( 13 + $IND )).1_genotype_indelExcluded_allele${IND}_first.log \
    --refFile=$REFERENCE \
    $REGIONSARG \
    --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
    --minPosterior=0 \
    --nCPU=$CPUS \
    --minReads=3 \
    --source=$OUTDIR/5-7_merged/IndelExcludedSNVs_allele${IND}.sorted.vcf.gz \
    --getVariantsFromBAMs=0 \
    $EXTRA \
    -o $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele${IND}_first.vcf

    # (Some calls may not be genotyped due to the way Platypus builds haplotypes)
    # Extract missing calls by comparing merged and genotyped VCFs
    tail -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele${IND}_first.vcf | cut -f1,2,4,5 > $OUTDIR/geno_pos.txt
    cut -f1,2,4,5 $OUTDIR/5-7_merged/IndelExcludedSNVs_allele${IND}.vcf > $OUTDIR/merged_pos.txt
    grep -vxFf $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt > $OUTDIR/coords.txt

    # Create a new regions file containing only the bases of the missing variants
    # The size of the region is the length of the SNV/indel
    awk '{if (length($3) >= length($4)) { print $1 ":" $2 "-" $2+length($3)-1 } else { print $1 ":" $2 "-" $2+length($4)-1 }}' $OUTDIR/coords.txt > $OUTDIR/8-18_genotyped/varRegions_allele${IND}_indelExcluded.txt
    rm $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt $OUTDIR/coords.txt

    # If there are missing calls: run Platypus to re-genotype them
    if [ -s $OUTDIR/8-18_genotyped/varRegions_allele${IND}_indelExcluded.txt ]; then
        echo -e "\nGenotyping missing calls in allele $IND\n"
        Platypus.py callVariants \
        --logFileName=$OUTDIR/logs/$(( 13 + $IND )).2_genotype_indelExcluded_allele${IND}_second.log \
        --refFile=$REFERENCE \
        --bamFiles=$OUTDIR/8-18_genotyped/bam_list.txt \
        --regions=$OUTDIR/8-18_genotyped/varRegions_allele${IND}_indelExcluded.txt \
        --minPosterior=0 \
        --nCPU=$CPUS \
        --minReads=3 \
        --source=$OUTDIR/5-7_merged/IndelExcludedSNVs_allele${IND}.sorted.vcf.gz \
        --getVariantsFromBAMs=0 \
        $EXTRA \
        -o $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele${IND}_second.vcf
    else
        echo -e "\nNo missing calls"
    fi

}


# 17) merge_filter_indelflagged()
# Merges and filters genotyped SNVs which are flagged for being close to indels
merge_filter_indelflagged() {

    # Collect all the genotyped indel-excluded SNVs in a single file
    # Header
    head -48 $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele1_first.vcf > $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.vcf
    # Content
    tail -q -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele?_first.vcf >> $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.vcf 2> /dev/null
    tail -q -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele?_second.vcf >> $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.vcf 2> /dev/null

    # Filter calls with flags badReads, MQ, strandBias, SC or QD
    awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.vcf > $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.filtered.vcf

    # Filter indel-flagged SNVs with median read coverage <20, median VAF <0.2 or median VAF >0.9
    Somatypus_IndelRescuedFilter.py $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.filtered.vcf > $OUTDIR/logs/17_indel_rescued_filter.log

}


# 18) merge_filter_all()
# Merges and filters all genotyped variants to create the output VCF file
merge_filter_all() {

    # Collect all the genotyped SNVs and indels in two files
    # Indels
    cp $OUTDIR/8-18_genotyped/GenotypedIndels_first.vcf $OUTDIR/8-18_genotyped/GenotypedIndels_merged.vcf
    tail -n +49 $OUTDIR/8-18_genotyped/GenotypedIndels_second.vcf >> $OUTDIR/8-18_genotyped/GenotypedIndels_merged.vcf 2> /dev/null
    # SNVs
    head -48 $OUTDIR/8-18_genotyped/GenotypedSNVs_allele1_first.vcf > $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.vcf
    tail -q -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_allele?_first.vcf >> $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.vcf 2> /dev/null
    tail -q -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_allele?_second.vcf >> $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.vcf 2> /dev/null
    tail -n +49 $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.filtered.VAFfilt.vcf >> $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.vcf

    # Filter SNV and indel calls with flags badReads, MQ, strandBias, SC or QD
    awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.vcf > $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.filtered.vcf
    awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $OUTDIR/8-18_genotyped/GenotypedIndels_merged.vcf > $OUTDIR/8-18_genotyped/GenotypedIndels_merged.filtered.vcf

    # Filter variants with a VAF >0.9 in all samples
    Somatypus_VAFfilter.py $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.filtered.vcf $OUTDIR/8-18_genotyped/GenotypedIndels_merged.filtered.vcf > $OUTDIR/logs/18_VAF_filter.log

    # Sort VCFs
    vcf-sort -c $OUTDIR/8-18_genotyped/GenotypedSNVs_merged.filtered.VAFfilt.vcf > $OUTDIR/Somatypus_SNVs_final.vcf 2> /dev/null
    vcf-sort -c $OUTDIR/8-18_genotyped/GenotypedIndels_merged.filtered.VAFfilt.vcf > $OUTDIR/Somatypus_Indels_final.vcf 2> /dev/null
        
}

################################### END OF FUNCTIONS ####################################


# Check that dependencies (Platypus, tabix, bgzip, vcf-sort, Somatypus scripts) are installed
hash Platypus.py 2>/dev/null || { echo -e "\nERROR: Platypus.py: command not found. Please install Platypus and add its directory to your PATH.\n" >&2; exit 1; }    
hash tabix 2>/dev/null || { echo -e "\nERROR: tabix: command not found. Please install the SAMtools tabix package and add its directory to your PATH.\n" >&2; exit 1; }    
hash bgzip 2>/dev/null || { echo -e "\nERROR: bgzip: command not found. Please install the SAMtools tabix package and add its directory to your PATH.\n" >&2; exit 1; }    
hash vcf-sort 2>/dev/null || { echo -e "\nERROR: vcf-sort: command not found. Please install VCFtools and add its directory to your PATH.\n" >&2; exit 1; }    
hash Somatypus_SplitMA-MNVs.py 2>/dev/null || { echo -e "\nERROR: Somatypus directory not included in the PATH. Please add the somatypus/src directory to your PATH.\n" >&2; exit 1; }    



# If no arguments (or -h): print help
if [ "$#" -eq 0 ]; then
    print_help
    exit 0
fi


# Parse input
BAMSDIR=""
REFERENCE=""
REGIONS="no"
OUTDIR=""
CPUS=1
EXTRA=""
while getopts ":i:g:r:o:c:p:hv?" OPT; do
  case $OPT in
    i)
      BAMSDIR=$OPTARG
      ;;
    g)
      REFERENCE=$OPTARG
      ;;
    r)
      REGIONS=$OPTARG
      ;;
    o)
      OUTDIR=$OPTARG
      ;;
    c)
      CPUS=$OPTARG
      ;;
    p)
      EXTRA=$OPTARG
      ;;
    h)
      print_help
      exit 0
      ;;
    v)
      echo "Somatypus $VERSION"
      exit 0
      ;;
    \?)
      print_help
      echo -e "Invalid option: -$OPTARG\n" >&2
      exit 1
      ;;
  esac
done


# Check that all mandatory inputs are present
if [ -z "$BAMSDIR" ]; then
   print_help
   echo -e "Input BAM files folder (-i) is required\n" >&2
   exit 1
fi

if [ -z "$REFERENCE" ]; then
   print_help
   echo -e "Reference genome FASTA file (-g) is required\n" >&2
   exit 1
fi

if [ -z "$OUTDIR" ]; then
   print_help
   echo -e "Path to the output folder (-o) is required\n" >&2
   exit 1
fi


# Sanity checks: 
# Check that: Input files exist; FASTA is indexed; there are input BAMs; all BAMs are indexed; CPUS >0;
# extra Platypus options do not contain options used within the pipeline
if [ ! -s $REFERENCE ]; then
    echo -e "\nERROR: Reference FASTA file not found or empty. Please check the path.\n" >&2
    exit 1
fi

if [ ! -s ${REFERENCE}.fai ]; then
    echo -e "\nERROR: Reference FASTA file has no accompanying .fai index file. Please index the FASTA file using 'samtools faidx' or equivalent.\n" >&2
    exit 1
fi

if  [ "$REGIONS" != "no" ] && [ ! -s $REGIONS ]; then
    echo -e "\nERROR: Regions file not found or empty. Please check the path.\n" >&2
    exit 1
fi

if ! [[ $CPUS =~ ^[1-9]+[0-9]*$ ]]; then
    echo -e "\nERROR: Number of CPUs must be greater than 0\n" >&2
    exit 1
fi

if echo "$EXTRA" | grep -q -E "\-\-logFileName|\-\-refFile|\-\-bamFiles|\-\-regions|\-\-minPosterior|\-\-minReads|\-\-minFlank|\-\-trimReadFlank|\-\-source|\-\-getVariantsFromBAMs|\-\-nCPU|\-\-output=|\-o " ; then 
    echo -e "\nERROR: Additional Platypus options cannot include --logFileName, --refFile, --bamFiles, --regions, --minPosterior, --minReads, --minFlank, --trimReadFlank, --source, --getVariantsFromBAMs, --nCPU, --output, or -o.\n" >&2
    exit 1
fi

if [ ! -d $BAMSDIR ]; then
    echo -e "\nERROR: $BAMSDIR directory not found. Please check the path.\n" >&2
    exit 1
fi

FILES=`ls -1 $BAMSDIR/*.bam 2> /dev/null | wc -l`
if [ "$FILES" -eq 0 ]; then
    echo -e "\nERROR: $BAMSDIR contains no files with .bam extension. Please check the path.\n" >&2
    exit 1
fi

for FILE in `ls $BAMSDIR/*.bam`; do
    if [ ! -s ${FILE}.bai ]; then
        echo -e "\nERROR: The file $FILE has no accompanying .bai index file. Please index the file using 'samtools sort' and 'samtools index' or equivalent.\n" >&2
        exit 1
    fi
done



# START RUNNING
# Copy all standard out and standard error to log file
mkdir -p $OUTDIR/logs
exec &> >(tee -ia $OUTDIR/logs/SOMATYPUS_`date +"%y%m%d%H%M"`.log)

echo -e "\nThis is Somatypus $VERSION\n"

echo "Input BAMs directory:   $BAMSDIR"
echo "Input reference genome: $REFERENCE"
echo "Input regions file:     $REGIONS"
echo "Output directory:       $OUTDIR"
echo "Number of CPUs to use:  $CPUS"
if [ -z "$EXTRA" ]; then
    echo "Extra Platypus options: no"
else
    echo "Extra Platypus options: $EXTRA"
fi


# Check if there is a checkpoint file from a previous run in the output folder
STEP=0
if [ -s $OUTDIR/logs/CHECKPOINT ]; then
    CHK=`tail -1 $OUTDIR/logs/CHECKPOINT`
    STEP=`echo $CHK | cut -f1 -d" "`
    STEPNAME=`echo $CHK | cut -f2 -d" "`
    echo -e "\n*CHECKPOINT FILE FOUND*"
    echo -e "Resuming execution after last completed step: $STEPNAME"
fi


echo -e "\nExecution started on `date`"


# Each step is performed only if its index is higher than STEP (last finished step index)
# 1. RUN PLATYPUS (DEFAULT) INDIVIDUALLY ON EVERY SAMPLE
if [ "$STEP" -lt 1 ]; then

    echo -e "\n(1) RUNNING PLATYPUS (DEFAULT SETTINGS) INDIVIDUALLY ON EVERY SAMPLE"
    individual_calling 0
    
    # Check successful execution
    for FILE in `ls $BAMSDIR/*.bam`; do 
        NAME=`basename $FILE`
        check_file $OUTDIR/1-2_individual_calls/platypusVariants_"${NAME%.*}"_default.vcf
    done
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "1 individual_calling_default" >> $OUTDIR/logs/CHECKPOINT

fi


# 2. RUN PLATYPUS (ALTERNATIVE) INDIVIDUALLY ON EVERY SAMPLE
if [ "$STEP" -lt 2 ]; then

    echo -e "\n(2) RUNNING PLATYPUS (ALTERNATIVE SETTINGS) INDIVIDUALLY ON EVERY SAMPLE"
    individual_calling 1
    
    # Check successful execution
    for FILE in `ls $BAMSDIR/*.bam`; do 
        NAME=`basename $FILE`
        check_file $OUTDIR/1-2_individual_calls/platypusVariants_"${NAME%.*}"_alternative.vcf
    done
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "2 individual_calling_alternative" >> $OUTDIR/logs/CHECKPOINT

fi


# 3. SPLIT INDIVIDUAL CALLS
if [ "$STEP" -lt 3 ]; then

    echo -e "\n(3) SPLITTING MULTI-ALLELIC AND MNP CALLS"
    split_calls
    
    # Check successful execution
    for FILE in `ls $OUTDIR/1-2_individual_calls/platypusVariants_*`; do
        NAME=`basename $FILE`
        check_file $OUTDIR/3_individual_split/"${NAME%.*}".split.vcf
    done

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "3 split_calls" >> $OUTDIR/logs/CHECKPOINT

fi


# 4. FILTER INDIVIDUAL CALLS
if [ "$STEP" -lt 4 ]; then

    echo -e "\n(4) FILTERING INDIVIDUAL CALLS"
    filter_calls

    # Check successful execution
    for FILE in `ls $OUTDIR/3_individual_split/platypusVariants_*`; do
        NAME=`basename $FILE`
        check_file $OUTDIR/4_individual_filtered/"${NAME%.*}".filtered.vcf
    done

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "4 filter_calls" >> $OUTDIR/logs/CHECKPOINT

fi


# 5. IDENTIFY SNVS CLOSE TO INDELS IN ANY SAMPLE
if [ "$STEP" -lt 5 ]; then

    echo -e "\n(5) FLAGGING SNVS CLOSE TO INDELS IN ANY SAMPLE (USING UNFILTERED DATA)"
    indel_flag

    # Check successful execution
    check_file $OUTDIR/5-7_merged/indel_flagged_SNVs.txt

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "5 indel_flag" >> $OUTDIR/logs/CHECKPOINT

fi


# 6. MERGE THE VCFS ACCORDING TO CHR,POS,REF,ALT VALUES
if [ "$STEP" -lt 6 ]; then

    echo -e "\n(6) MERGING FILTERED SNVS"
    merge_calls

    # Check successful execution
    check_file $OUTDIR/5-7_merged/MergedSNVs_allele1.vcf
    check_file $OUTDIR/5-7_merged/IndelExcludedSNVs_allele1.vcf
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "6 merge_calls" >> $OUTDIR/logs/CHECKPOINT

fi


# 7. EXTRACT AND FILTER INDELS
if [ "$STEP" -lt 7 ]; then

    echo -e "\n(7) EXTRACTING AND FILTERING INDELS"
    extract_indels
    
    # Check successful execution
    check_file $OUTDIR/5-7_merged/MergedIndels.vcf

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "7 extract_indels" >> $OUTDIR/logs/CHECKPOINT

fi


# 8. PREPARE DATA FOR GENOTYPING
if [ "$STEP" -lt 8 ]; then
    
    echo -e "\n(8) PREPARING DATA FOR VARIANT GENOTYPING"
    prepare_genotyping
    
    # Check successful execution
    for FILE in `ls $OUTDIR/5-7_merged/MergedSNVs_allele?.vcf`; do
        check_file "${FILE%.*}".sorted.vcf.gz
        check_file "${FILE%.*}".sorted.vcf.gz.tbi
    done
    check_file $OUTDIR/5-7_merged/MergedIndels.sorted.vcf.gz
    check_file $OUTDIR/5-7_merged/MergedIndels.sorted.vcf.gz.tbi
    if [ "$REGIONS" != "no" ]; then
        check_file $OUTDIR/8-18_genotyped/regions_allele1.txt
        check_file $OUTDIR/8-18_genotyped/regions_indels.txt
    fi
    check_file $OUTDIR/8-18_genotyped/bam_list.txt
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "8 prepare_genotyping" >> $OUTDIR/logs/CHECKPOINT

fi


# 9-11. GENOTYPE ALLELE 1/2/3 SNVS
for ALL in `seq 1 3`; do
    if [ "$STEP" -lt $(( 8 + $ALL )) ]; then
        
        echo -e "\n($(( 8 + $ALL ))) GENOTYPING 'ALLELE ${ALL}' SNVS\n"
        genotyping $ALL
        
        # Check successful execution
        check_file $OUTDIR/8-18_genotyped/GenotypedSNVs_allele${ALL}_first.vcf
        
        # Update checkpoint file
        echo -e "\nSuccess"
        echo "$(( 8 + $ALL )) genotyping_allele$ALL" >> $OUTDIR/logs/CHECKPOINT

    fi
done



# 12. GENOTYPE INDELS
if [ "$STEP" -lt 12 ]; then
        
        echo -e "\n(12) GENOTYPING INDELS\n"
        genotyping 0
        
        # Check successful execution
        check_file $OUTDIR/8-18_genotyped/GenotypedIndels_first.vcf
        
        # Update checkpoint file
        echo -e "\nSuccess"
        echo "12 genotyping_indels" >> $OUTDIR/logs/CHECKPOINT

fi



# 13. PREPARE DATA FOR GENOTYPING OF INDEL-EXCLUDED SNVS
if [ "$STEP" -lt 13 ]; then
    
    echo -e "\n(13) PREPARING DATA FOR GENOTYPING OF INDEL-EXCLUDED SNVS"
    prepare_genotyping_indelflagged
    
    # Check successful execution
    for FILE in `ls $OUTDIR/5-7_merged/IndelExcludedSNVs_allele?.vcf`; do
        check_file "${FILE%.*}".sorted.vcf.gz
        check_file "${FILE%.*}".sorted.vcf.gz.tbi
    done
    if [ "$REGIONS" != "no" ]; then
        check_file $OUTDIR/8-18_genotyped/regions_allele1_indelExcluded.txt
    fi
    check_file $OUTDIR/8-18_genotyped/bam_list.txt
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "13 prepare_genotyping_indelflagged" >> $OUTDIR/logs/CHECKPOINT

fi


# 14-16. GENOTYPE ALLELE 1/2/3 INDEL-EXCLUDED SNVS
for ALL in `seq 1 3`; do
    if [ "$STEP" -lt $(( 13 + $ALL )) ]; then

        echo -e "\n($(( 13 + $ALL ))) GENOTYPING 'ALLELE ${ALL}' INDEL-EXCLUDED SNVS\n"
        genotyping_indelflagged $ALL

        # Check successful execution
        check_file $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_allele${ALL}_first.vcf

        # Update checkpoint file
        echo -e "\nSuccess"
        echo "$(( 13 + $ALL )) genotyping_indelflagged_allele$ALL" >> $OUTDIR/logs/CHECKPOINT

    fi
done


# 17. MERGE AND FILTER INDEL-EXCLUDED SNVS
if [ "$STEP" -lt 17 ]; then
    
    echo -e "\n(17) MERGING AND FILTERING GENOTYPED INDEL-EXCLUDED SNVS"
    merge_filter_indelflagged
    
    # Check successful execution
    check_file $OUTDIR/8-18_genotyped/GenotypedSNVs_indelExcluded_merged.filtered.VAFfilt.vcf
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "17 merge_filter_indelflagged" >> $OUTDIR/logs/CHECKPOINT

fi


# 18. MERGE AND FILTER ALL VARIANTS
if [ "$STEP" -lt 18 ]; then
    
    echo -e "\n(18) MERGING AND FILTERING ALL GENOTYPED VARIANTS"
    merge_filter_all
    
    # Check successful execution
    check_file $OUTDIR/Somatypus_SNVs_final.vcf
    check_file $OUTDIR/Somatypus_Indels_final.vcf
    
    # Update checkpoint file
    echo -e "\nSuccess"
    echo "18 merge_filter_all" >> $OUTDIR/logs/CHECKPOINT

fi


echo -e "\nExecution finished on `date`"

echo -e "\nALL DONE!"
echo -e "Output is in: $OUTDIR/Somatypus_SNVs_final.vcf\n              $OUTDIR/Somatypus_Indels_final.vcf\n"


