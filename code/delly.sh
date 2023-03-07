#!/usr/bin/bash
set -ex

source ./config.sh


mkdir -p ../results/dellySV/1calls
mkdir -p ../results/dellySV/2prefilter
mkdir -p ../results/dellySV/3genotypes
mkdir -p ../results/dellySV/4filter

mkdir -p ../results/dellyCNV/1tumor
mkdir -p ../results/dellyCNV/2control
mkdir -p ../results/dellyCNV/3merge
mkdir -p ../results/dellyCNV/4filter
mkdir -p ../results/dellyCNV/5plot



# 1. SV discovery: Use at least 1 tumor and matched control
for line in $(cat ${comparesheet}); do
    control=$(echo $line | awk -F, '{print $1}')
    case=$(echo $line | awk -F, '{print $2}')

    control_file=$(find -L ${data_dir} -name "${control}"*.bam)
    case_file=$(find -L ${data_dir} -name "${case}"*.bam)
    
    /delly/src/delly call -g ${reference} -x ${SVbed} -q ${map_quality} -r ${qual_translocation} -c ${min_clip} -z ${min_clique_size} -m ${min_split_read_dist} -n ${max_split_read_dist} -o ../results/dellySV/1calls/sample_pair_${case}_${control}.bcf ${case_file} ${control_file}   
done



# 2. Somatic prefiltering
for pairedbcfs in ../results/dellySV/1calls/sample_pair*.bcf
do
    case=$(basename ${pairedbcfs} | sed 's/_[^_]*$//g' | uniq)
    /delly/src/delly filter -f somatic -o ../results/dellySV/2prefilter/${case}_prefilter.bcf -s ${samplesheet} ${pairedbcfs}
done 



# 3. Genotype pre-filtered somatic sites across a larger panel of control samples to efficiently filter false postives and germline SVs. (1 tumor and lots of controls)

# # Re-specify these variables
# control_bams=$(while IFS= read -r line; do
#     if [[ ${line} =~ "control" ]]; then
#         control_name=$(echo ${line} | cut -d' ' -f1)
#         echo $(find -L /data /results/data/controls -name ${control_name}*.bam)
#     fi
# done < ${samplesheet})
# echo $controls


# # find the tumors
# tumors=$(while IFS= read -r line; do
#     if [[ ${line} =~ "tumor" ]]; then
#         echo ${line} | cut -d' ' -f1
#     fi
# done < ${samplesheet})
# echo ${tumors}

for tumorsamp in ${tumors}; do
    #tumorsamp=$(basename ${file} .bam)
    mkdir -p ../results/dellySV/3genotypes/${tumorsamp}
    tumor_bcf=$(find -L ../results/ -name "sample_pair_${tumorsamp}_*prefilter.bcf")
    /delly/src/delly call -g ${reference} -v ${tumor_bcf} -x ${SVbed} -q ${map_quality} -r ${qual_translocation} -c ${min_clip} -z ${min_clique_size} -m ${min_split_read_dist} -n ${max_split_read_dist} -o ../results/dellySV/3genotypes/${tumorsamp}/${tumorsamp}_geno.bcf $(find -L /data /results -name ${tumorsamp}*.bam) ${control_bams} #RG

    # 4. Post-filter for somatic SVs using all control samples.
    /delly/src/delly filter -f somatic -o ../results/dellySV/4filter/${tumorsamp}/${tumorsamp}_somatic.bcf -s ${samplesheet} ../results/dellySV/3genotypes/${tumorsamp}/${tumorsamp}_geno.bcf
done



############################################
# CNV


# 1. Somatic copy-number alterations detection (-u is required). Depending on the coverage, tumor purity and heterogeneity you can adapt parameters -z, -t and -x which control the sensitivity of SCNA detection.

for tumorsamp in ${tumors}; do
    /delly/src/delly cnv -u -o ../results/dellyCNV/1tumor/tumor_${tumorsamp}.bcf -c ../results/dellyCNV/1tumor/tumor_${tumorsamp}_out.cov.gz -g ${reference} -m ${CNVmap} $(find -L ../data ../results -name ${tumorsamp}*.bam) #RG
done 


# 2. Then these tumor SCNAs are genotyped in the control sample (-u is required).

for line in $(cat ${comparesheet}); do
    control=$(echo $line | awk -F, '{print $1}')
    case=$(echo $line | awk -F, '{print $2}')

    control_file=$(find -L ${data_dir} -name "${control}"*.bam)
    
    /delly/src/delly cnv -u -v ../results/dellyCNV/1tumor/tumor_${case}.bcf -o ../results/dellyCNV/2control/control_${control}.bcf -g ${reference} -m ${CNVmap} ${control_file}
done



# 3. The VCF IDs are matched between tumor and control. Thus, you can merge both files using bcftools.

for line in $(cat ${comparesheet}); do
    control=$(echo $line | awk -F, '{print $1}')
    case=$(echo $line | awk -F, '{print $2}')

    bcftools merge -m id --threads ${num_threads} -O b -o ../results/dellyCNV/3merge/${case}_${control}.bcf ../results/dellyCNV/1tumor/tumor_${case}.bcf ../results/dellyCNV/2control/control_${control}.bcf
    bcftools index -c --threads ${num_threads} ../results/dellyCNV/3merge/${case}_${control}.bcf

done



# 4. Somatic filtering requires a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

for line in $(cat ${comparesheet}); do
    control=$(echo $line | awk -F, '{print $1}')
    case=$(echo $line | awk -F, '{print $2}')

    /delly/src/delly classify -f somatic -o ../results/dellyCNV/4filter/somatic_${case}_${control}.bcf -s ${samplesheet} ../results/dellyCNV/3merge/${case}_${control}.bcf
done



# 5. Plot the SCNAs using bcftools and R.

for line in $(cat ${comparesheet}); do
    echo $line 
    control=$(echo $line | awk -F, '{print $1}')
    case=$(echo $line | awk -F, '{print $2}')

    # A. Plot copy-number distribution for large number of samples (>>100)
    # First convert formats, then cd because Rscript will print to pwd
    if [ ${num_attached_bcfs} -gt 100 ]; then
        bcftools query -f "%ID[\t%RDCN]\n" -o ../results/plot.tsv ../data/somatic_${case}_${control}.bcf 
        cd ../results
        R CMD BATCH "--args ./plot.tsv" /delly/R/cnv.R
        cd ../code
    fi

    # B. Plot segmented read depth profiles
    # First convert formats, then cd because Rscript will print to pwd.  Must change back to /code/ so that relative paths are correct. 
    bcftools query -s ${case} -f "%CHROM\t%POS\t%INFO/END\t%ID\t[%RDCN]\n" -o ../results/dellyCNV/5plot/${case}_segmentation.bed ../results/dellyCNV/4filter/somatic_${case}_${control}.bcf
    mkdir -p ../results/dellyCNV/5plot/${case}_plots_seg
    cd ../results/dellyCNV/5plot/${case}_plots_seg
    R CMD BATCH "--args ../../1tumor/tumor_${case}_out.cov.gz ../${case}_segmentation.bed" /delly/R/rd.R 
    cd ../../../../code
done
