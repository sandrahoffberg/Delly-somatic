#!/usr/bin/env bash

if [ $# -eq 0 ]; then
    echo "No arguments supplied"
else
    echo "args:"
    for i in $*; do 
        echo $i 
    done
    echo ""
fi


if [ -z ${1} ]; then
    samplesheet=$(find -L ../data/ -name "samplesheet*.tsv")
else
    samplesheet=${1}
fi

# Check that there is only 1 reference
number_samplesheet=$(find -L ../data/ -name "*samplesheet*" | wc -l )
if [ ${number_samplesheet} -eq 0 ]; then
    echo "Check your input data: there is no samplesheet."
    exit 1
elif [ ${number_samplesheet} -gt 1 ]; then
    echo "Check your input data: there are multiple samplesheets."
    exit 1
fi

# make the path absolute
samplesheet=$PWD/${samplesheet}


# URL of the s3 bucket with bams for analysis
if [ -z ${2} ]; then
    num_bam_files=$(find -L ../data -name "*.bam" | wc -l)
    if [ ${num_bam_files} -gt 0 ]; then
        echo "Using bams in the /data folder"
        data_dir="../data"
    else
        echo "There are no sample files available for this capsule.  Please either add some to the /data folder or provide a URL to an S3 bucket."
    fi
else
    s3_url=${2}
    
    # location to download the bams
    if [ -z ${3} ]; then
        data_dir=../results/data
    else
        data_dir=${3}
    fi
    
    #get BAM file(s) that have been trimmed, aligned to reference, sorted, read groups adjusted if necessary, and indexed
    echo "Downloading data from the S3 bucket at the provided URL."
    aws s3 sync --no-sign-request --only-show-errors ${s3_url} ${data_dir}


    # Organize the controls and cases into separate folders
    control_bams=$(while IFS= read -r line; do
    if [[ ${line} =~ "control" ]]; then
        control_name=$(echo ${line} | cut -d' ' -f1)
        echo $(find -L ../data ../results -name ${control_name}*.bam)
    fi
    done < ${samplesheet})
    echo ${control_bams}

    mkdir -p ../results/data/controls
    mv $(dirname ${control_bams}) ../results/data/controls


    # find the tumors
    tumor_bams=$(while IFS= read -r line; do
        if [[ ${line} =~ "tumor" ]]; then
            tumor_name=$(echo ${line} | cut -d' ' -f1)
            echo $(find -L ../data ../results -name ${tumor_name}*.bam)
        fi
    done < ${samplesheet})
    echo ${tumor_bams}

    mkdir -p ../results/data/cases
    mv $(dirname ${tumor_bams}) ../results/data/cases
fi


if [ -z $AWS_BATCH_JOB_ID ]; then
    directory_with_cases=../results/data/cases/*
    directory_with_controls=../results/data/controls/*
else
    directory_with_cases=${data_dir}/cases/*
    directory_with_controls=${data_dir}/controls/*
fi


if [ -z ${4} ]; then
    comparesheet=$(find -L ../data/ -name "comparesheet*.csv")
else
    comparesheet=${4}
fi

# Check that there is only 1 reference
number_comparesheet=$(find -L ../data/ -name "*compare*" | wc -l )
if [ ${number_comparesheet} -eq 0 ]; then
    echo "Check your input data: there is no comparesheet."
    exit 1
elif [ ${number_comparesheet} -gt 1 ]; then
    echo "Check your input data: there are multiple comparesheets."
    exit 1
fi

# make the path absolute
comparesheet=$PWD/${comparesheet}





# If a reference file has not been specified in the app panel, find it
if [ -z ${5} ]; then
    reference=$(find -L ../data/ -name "*.fa")
else
    reference=${5}
fi

# Check that there is only 1 reference
number_references=$(find -L ../data/ -name "*.fa" | wc -l )
if [ ${number_references} -eq 0 ]; then
    echo "Check your input data: there is no reference .fa file."
    exit 1
elif [ ${number_references} -gt 1 ]; then
    echo "Check your input data: there are multiple reference .fa files."
    exit 1
fi

# make the path absolute
reference=$PWD/${reference}


# If a bed file of regions to exclude has not been specified in the app panel, find it
if [ -z ${6} ]; then
  SVbed=$(find -L ../data/ -name "*.bed")
else
  SVbed=${6}
fi

# Check that there is only 1 bed file
number_SVbed=$(find -L ../data/ -name "*.bed" | wc -l )
if [ ${number_SVbed} -gt 1 ]; then
  echo "Check your input data: there are multiple bed files with regions to exclude."
  exit 1
fi

# make the path absolute
SVbed=$PWD/${SVbed}


# If an annotation file has not been specified in the app panel, find it
if [ -z ${7} ]; then
    CNVmap=$(find -L ../data/ -name "*blacklist.gz")
else
    CNVmap=${7}
fi

number_CNVmap=$(find -L ../data/ -name "*blacklist.gz" | wc -l )
if [ ${number_CNVmap} -gt 1 ]; then
    echo "Check your input data: there are multiple blacklist.gz files."
    exit 1
fi

# make the path absolute
CNVmap=$PWD/${CNVmap}



#Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples.
cores=$(get_cpu_count.py)

if [ -z ${8} ]; then
    num_threads=$((cores * 2))
else
    num_threads=${8}
fi

num_attached_bams=$(find -L ${data_dir} -name "*.bam" | wc -l)

# Check to make sure all of the necessary files are attached. 
if [ ${num_attached_bams} == 0 ]; then
    echo "Check your input data: there are no (indexed) bam files."
    exit 1
fi

min=$(( ${num_attached_bams} < ${num_threads} ? ${num_attached_bams} : ${num_threads} ))

export OMP_NUM_THREADS=${min}
