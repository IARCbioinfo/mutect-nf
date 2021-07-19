#! /usr/bin/env nextflow

// Copyright (C) 2010 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.mem           = 8
params.cpu           = 4
params.output_folder = "mutect_pon_results"
params.mutect_args   = ""
params.nsplit        = 1
params.java          = "java"
params.known_snp     = "NO_SNP_FILE"
params.region        = null
params.bed           = null
params.gatk_version  = "4"
params.normal_file       = null
params.ext                      = "cram"
//#known_snp: /data/gcs/mesomics/files/databases/BundleMutect2/af-only-gnomad.hg38.vcf.gz
params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  mutect-nf 2.2.0: Mutect pipeline for PON creation with Nextflow "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {
log.info '-------------------------------------------------------------'
    log.info ' USAGE  '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'nextflow run mutect_build_pon.nf --ref ref.fasta --normal_file normals.tab [OPTIONS] '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --normal_file            FILE                input tabulation-separated values file with columns sample (sample name),'
    log.info '                                                 and normal (full path to matched normal cram);'
    log.info '    --ref                FILE (with indexes)     Reference fasta file.'
    log.info 'Optional arguments:'
    log.info '    --nsplit             INTEGER                 Split the region for calling in nsplit pieces and run in parallel (default: 1).'
    log.info '    --known_snp          FILE                    VCF file with known variants and frequency (e.g., from gnomad).'
    log.info '    --mutect_args        STRING                  Arguments you want to pass to mutect.'
    log.info '    --gatk_version 	   INTEGER                 gatk version, used to call Mutect and add the appropriate options (default: 4).'
    log.info '    --cpu                INTEGER                 Number of cpu used (default: 4).'
    log.info '    --mem                INTEGER                 Java memory passed to mutect in GB (default: 8).'
    log.info '    --output_folder      FOLDER                  Output folder (default: mutect_pon_results).'
    log.info '    --java               PATH                    Name of the JAVA command  (default: java).'
    log.info '    --ext                STRING                  Type of alignment file [Def:cram]'
    log.info ''
    exit 0
}else{
    /* Software information */
    log.info "mem                    = ${params.mem}"
    log.info "cpu                    = ${params.cpu}"
    log.info "output_folder          = ${params.output_folder}"
    log.info "mutect_args            = ${params.mutect_args}"
    log.info "nsplit                 = ${params.nsplit}"
    log.info "known_snp              = ${params.known_snp}"
    log.info "gatk_version           = ${params.gatk_version}"
    log.info "normal_file            = ${params.normal_file}"
    log.info "ref                    = ${params.ref}"
    log.info "ext                   = ${params.ext}"
}

//load reference
fasta_ref      = file( params.ref )
fasta_ref_fai  = file( params.ref+'.fai' )
fasta_ref_gzi  = file( params.ref+'.gzi' )
fasta_ref_dict = file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )


ext_ind = ".crai"
if(params.ext=="bam"){ ext_ind=".bai"}

//load know snps
known_snp     = file(params.known_snp)
known_snp_tbi = file(params.known_snp+".tbi")



//load input files
if (params.normal_file) {
    // FOR INPUT AS A TAB DELIMITED FILE
    pairs = Channel.fromPath(params.normal_file).splitCsv(header: true, sep: '\t', strip: true)
                       .map{ row -> [ row.sample , file(row.normal), file(row.normal+ext_ind) ] }
    tn_bambai = pairs.groupTuple(by: 0)
                              .map { row -> tuple(row[0] , row[1] , row[2]) }

} else {
    exit 0;
}



/* manage input positions to call (bed or region or whole-genome) */
input_region = 'whole_genome'

if (params.region) {
    input_region = 'region'
} else if (params.bed) {
    input_region = 'bed'
    bed = file(params.bed)
} else {
    input_region = 'whole_genome'
}


/* process to create a bed file from region or from faidx if whole-genome, otherwise return the input bed file */
process bed {
      output:
      file "temp.bed" into outbed, all_regions

      shell:
      if (input_region == 'region')
      '''
      echo !{params.region} | sed -e 's/[:|-]/	/g' > temp.bed
      '''
      else if (input_region == 'whole_genome')
      '''
      cat !{fasta_ref_fai} | awk '{print $1"	"0"	"$2 }' | grep -v -P "alt|random|Un|chrEBV|HLA" > temp.bed
      '''
  }

/* split bed file into nsplit regions */
process split_bed {

      input:
      file bed from outbed

      output:
      file '*_regions.bed' into split_bed, count_split_bed mode flatten

      shell:
      '''
      grep -v '^track' !{bed} | sort -T $PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1" "$2" "$3}' | cut_into_small_beds.r !{params.nsplit}
      '''
}

( split_bed1 , split_bed2) = split_bed.into(2)

process mutect {
    memory params.mem+'GB'
    cpus params.cpu

    tag { printed_tag }

    input:
    set val(sample), file(bamN), file(baiN), file(bed)  from tn_bambai.combine(split_bed1)
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_dict
    file known_snp
    file known_snp_tbi

    output:
    set val(sample), file("${printed_tag}_*.vcf") into mutect_output1

    shell:
    bed_tag0 = bed.baseName
    bed_tag = bed_tag0.replaceAll("[^a-zA-Z0-9 _-]+","")
    printed_tag = "${sample}_" + bed_tag
    if("${params.gatk_version}" == "4"){
    '''
    gatk Mutect2 --java-options "-Xmx!{params.mem}G" -R !{fasta_ref} -L !{bed} !{params.mutect_args} -I !{bamN} --max-mnp-distance 0 -O !{printed_tag}_pon.vcf
    '''
    }else{
      '''
      echo "Mutect version not supported\n"
      exit 1;
      '''
    }
}

beds_length = count_split_bed.count().val

process mergeMuTectOutputs {

    tag { normal_tag }

    publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
                 if (filename.indexOf(".vcf.gz") > 0)     "intermediate_calls/raw_calls/$filename"
    }

    input:
    set val(normal_tag), file(vcf_files) from mutect_output1.groupTuple(size: beds_length)

    output:
    file("${normal_tag}_pon.vcf.gz") into res_merged
    file("${normal_tag}_pon.vcf.gz") into res_merged_vcfs
    file("${normal_tag}_pon.vcf.gz.tbi") into res_merged_vcfs_indexs

    shell:
    '''
    #rm -f !{normal_tag}_pon.vcf
    #rm -f !{normal_tag}_pon.vcf.gz
    #rm -f !{normal_tag}_pon.vcf.gz.tbi
    # MERGE VCF FILES
    sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt
    # Check if sort command allows sorting in natural order (chr1 chr2 chr10 instead of chr1 chr10 chr2)
    if [ `sort --help | grep -c 'version-sort' ` == 0 ]
    then
        sort_ops="-k1,1d"
    else
        sort_ops="-k1,1V"
    fi
    # Add all VCF contents and sort
    grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -T $PWD -t '	' $sort_ops -k2,2n >> header.txt
    mv header.txt !{normal_tag}_pon.vcf
    #we index the created PON of normals
    bcftools view  -O z !{normal_tag}_pon.vcf > !{normal_tag}_pon.vcf.gz
    bcftools index -t !{normal_tag}_pon.vcf.gz

    '''
}

process create_Genomics_DB{
  tag { "Build_PON" }
  publishDir params.output_folder+'/PON/', mode: 'copy'
  input:
  file(vcfs) from res_merged_vcfs.collect()
  file(vcfs_index) from res_merged_vcfs_indexs.collect()
  file(intervals) from all_regions
  file fasta_ref
  file fasta_ref_fai
  file fasta_ref_dict
  file known_snp
  file known_snp_tbi
  output:
    file("pon.vcf.gz") into mutect_pon_result
    file("pon.vcf.gz.tbi") into mutect_pon_result_index
  shell:
  input_vcf=""
  for( vcf in vcfs ){
      input_vcf=input_vcf+" -V ${vcf}"
  }
    '''
    # We create the GenomicsDBinport database
    gatk GenomicsDBImport -R !{fasta_ref} -L !{intervals} \
        --genomicsdb-workspace-path pon_db !{input_vcf}
    # We create the panel of normals
    gatk CreateSomaticPanelOfNormals -R !{fasta_ref} \
     --germline-resource !{known_snp} \
     -V gendb://pon_db \
     -O pon.vcf.gz
    '''
}
