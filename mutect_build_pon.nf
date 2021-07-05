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

params.suffix_tumor  = "_T"
params.suffix_normal = "_N"
params.mem           = 8
params.cpu           = 4
params.output_folder = "mutect_results"
params.mutect_args   = ""
params.nsplit        = 1
params.region        = null
params.bed           = null
params.java          = "java"
params.known_snp     = "NO_SNP_FILE"
params.snp_contam    = "NO_FILE"
params.cosmic        = "NO_COSMIC_FILE"
params.mutect_jar    = null
params.mutect2_jar   = null
params.gatk_version  = "4"
params.tn_file       = null
params.tumor_bam_folder  = null
params.normal_bam_folder = null
params.PON           = null
params.estimate_contamination = null
params.filter_readorientation = null
params.genotype      = null
params.ref_RNA       = "NO_REF_RNA_FILE"
params.ext                      = "cram"

params.help = null

log.info "" 
log.info "--------------------------------------------------------"
log.info "  mutect-nf 2.2.0: Mutect pipeline for somatic variant calling with Nextflow "
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
    log.info 'nextflow run mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --ref ref.fasta [OPTIONS] '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --tumor_bam_folder   FOLDER                  Folder containing tumor BAM files to be called.'
    log.info '    --normal_bam_folder  FOLDER                  Folder containing matched normal BAM files.'
    log.info '    OR'
    log.info '    --tn_file            FILE                    input tabulation-separated values file with columns sample (sample name),'
    log.info '                                                 tumor (full path to tumor bam), normal (full path to matched normal bam);'
    log.info '                                                 optionally (for --genotype mode), columns preproc (is the bam RNAseq needing'
    log.info '                                                 preprocessing: yes or no) and vcf (full path to vcf file containing alleles to genotype)'
    log.info '                                                 where each line contains the path of two matched BAM files.'    
    log.info '    --ref                FILE (with indexes)     Reference fasta file.'
    log.info '    --mutect_jar         FILE                    mutect*.jar explicit path.'
    log.info '    OR'
    log.info '    --mutect2_jar        FILE                    gatk*.jar explicit path to run mutect2 (integrated to GATK).'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --region             REGION                  A region defining the calling, in the format CHR:START-END.'
    log.info '    NOTE: if neither --bed or --region, will perform the calling on whole genome, based on the faidx file.'
    log.info '    --nsplit             INTEGER                 Split the region for calling in nsplit pieces and run in parallel (default: 1).'
    log.info '    --known_snp          FILE                    VCF file with known variants and frequency (e.g., from gnomad).'
    log.info '    --snp_contam 		   FILE                    VCF file with known germline variants to genotype for contamination estimation'
    log.info '                                                 (requires --estimate_contamination)'
    log.info '    --PON 	           FILE                    VCF file of GATK panel of normals used to filter calls.'
    log.info '    --cosmic             FILE                    Cosmic VCF file required by mutect (CAUTION: not a symbolic link); not in gatk4.'
    log.info '    --mutect_args        STRING                  Arguments you want to pass to mutect.'
    log.info '                                                 WARNING: form is " --force_alleles " with spaces between quotes.'
    log.info '    --gatk_version 	   INTEGER                 gatk version, used to call Mutect and add the appropriate options (default: 4).'
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '    --ref_RNA 		   PATH                    fasta reference for preprocessing RNA (required when preproc column contains yes'
    log.info '                                                 in input tn_file).'
    log.info '    --cpu                INTEGER                 Number of cpu used (default: 4).'
    log.info '    --mem                INTEGER                 Java memory passed to mutect in GB (default: 8).'
    log.info '    --output_folder      FOLDER                  Output folder (default: mutect_results).'
    log.info '    --java               PATH                    Name of the JAVA command  (default: java).'
    log.info ''
    log.info 'Flags:'
    log.info '    --estimate_contamination 	                   Run extra step of estimating contamination and use results to filter calls; only for gatk4'
    log.info '    --filter_readorientation                     Run extra step learning read orientation model and using it to filter reads.'
    log.info '    --genotype                                   Use genotyping from vcf mode instead of usual variant calling;'
    log.info '                                                 requires tn_file with vcf column and gatk4, and if RNA-seq included, requires preproc column'
    log.info ''
    exit 0
}else{
    /* Software information */
    log.info "suffix_tumor           = ${params.suffix_tumor}"
    log.info "suffix_normal          = ${params.suffix_normal}"
    log.info "mem                    = ${params.mem}"
    log.info "cpu                    = ${params.cpu}"
    log.info "output_folder          = ${params.output_folder}"
    log.info "mutect_args            = ${params.mutect_args}"
    log.info "nsplit                 = ${params.nsplit}"
    log.info "region                 = ${params.region}"
    log.info "bed                    = ${params.bed}"
    log.info "java                   = ${params.java}"
    log.info "known_snp              = ${params.known_snp}"
    log.info "snp_contam             = ${params.snp_contam}"
    log.info "cosmic                 = ${params.cosmic}"
    log.info "mutect_jar             = ${params.mutect_jar}"
    log.info "mutect2_jar            = ${params.mutect2_jar}"
    log.info "gatk_version           = ${params.gatk_version}"
    log.info "tn_file                = ${params.tn_file}"
    log.info "tumor_bam_folder       = ${params.tumor_bam_folder}"
    log.info "normal_bam_folder      = ${params.normal_bam_folder}"
    log.info "PON                    = ${params.PON}"
    log.info "estimate_contamination = ${params.estimate_contamination}"
    log.info "filter_readorientation = ${params.filter_readorientation}"
    log.info "genotype               = ${params.genotype}"
    log.info "ref                    = ${params.ref}"
    log.info "ref_RNA                = ${params.ref_RNA}"
    log.info "ext                   = ${params.ext}"	
}

//load reference
fasta_ref      = file( params.ref )
fasta_ref_fai  = file( params.ref+'.fai' )
fasta_ref_gzi  = file( params.ref+'.gzi' )
fasta_ref_dict = file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )


ext_ind = ".crai"
if(params.ext=="bam"){ ext_ind=".bai"}


if(params.genotype){
    if(params.ref_RNA == "NO_REF_RNA_FILE"){
        fasta_ref_RNA      = file( params.ref )
        fasta_ref_RNA_fai  = file( params.ref+'.fai' )
        fasta_ref_RNA_dict = file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )
    }else{
        fasta_ref_RNA      = file( params.ref_RNA )
        fasta_ref_RNA_fai  = file( params.ref_RNA+'.fai' )
        fasta_ref_RNA_dict = file( params.ref_RNA.replace(".fasta",".dict").replace(".fa",".dict") )
    }
}

//load jar files for gatk<=3
jar  = file('NO_MUTECT_JAR_FILE')
if(params.mutect_jar){
    jar = file(params.mutect_jar)
}
jar2 = file('NO_MUTECT2_JAR_FILE')
if(params.mutect2_jar){
    jar2 = file(params.mutect2_jar)
}
mutect_version = params.mutect_jar ? 1 : 2

//load VCFs
snp_contam = file(params.snp_contam)
snp_contam_tbi = file(params.snp_contam+'.tbi')

cosmic = file(params.cosmic)
if (params.cosmic == "NO_COSMIC_FILE") { cosmic_option = "" } else { cosmic_option = "--cosmic ${cosmic.getName()}" }

known_snp     = file(params.known_snp)
known_snp_tbi = file(params.known_snp+".tbi")
if (params.known_snp == "NO_SNP_FILE") { known_snp_option = "" } else {
        if(params.gatk_version == "4"){
                known_snp_option = "--germline-resource ${known_snp.getName()}"
        }else{
                known_snp_option = "--dbsnp ${known_snp.getName()}"
        }
}
if (params.PON) { 
	PON = file(params.PON)
	PON_tbi = file(params.PON+'.tbi')
} else { 
	PON = file('NO_FILE')
	PON_tbi = file('NO_FILE2')
}

//load input files
if (params.tn_file) {
    // FOR INPUT AS A TAB DELIMITED FILE
    pairs = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
                       .map{ row -> [ row.sample , file(row.tumor), file(row.tumor+ext_ind), file(row.normal), file(row.normal+ext_ind) ] }

	pairs2 = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
                       .map{ row -> [ row.sample , file(row.tumor), file(row.tumor+ext_ind), file(row.normal), file(row.normal+ext_ind) ] }

    tn_bambai2 = pairs2.groupTuple(by: 0)
                              .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0]  ) }

    tn_bambai = pairs.groupTuple(by: 0)
                              .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0]  ) }

    if(params.estimate_contamination){
	    pairsT4cont = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
                       .map{ row -> [ row.sample , 'T' , file(row.tumor), file(row.tumor+ext_ind) ] }
	    pairsN4cont = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
                         .map{ row -> [ row.sample , 'N', file(row.normal), file(row.normal+ext_ind) ] }
		       .unique()
	    pairs4cont  = pairsT4cont.concat( pairsN4cont )
    }
} else {
    // FOR INPUT AS TWO FOLDER
    // recovering of bam files
	ext_ind
    tumor_bams = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.'+params.ext )
              .ifEmpty { error "Cannot find any bam/cram file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}."+params.ext,""), path ] }

    // recovering of bai files
    tumor_bais = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.'+params.ext+ext_ind )
              .ifEmpty { error "Cannot find any bai/crai file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}."+params.ext+ext_ind,""), path ] }

    // building bam-bai pairs
    tumor_bam_bai = tumor_bams
    	      .phase(tumor_bais)
    	      .map { tumor_bam, tumor_bai -> [ tumor_bam[0], tumor_bam[1], tumor_bai[1] ] }

    // FOR NORMAL
    // recovering of bam files
    normal_bams = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.'+params.ext )
              .ifEmpty { error "Cannot find any bam/cram file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}."+params.ext,""), path ] }

    // recovering of bai files
    normal_bais = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.'+params.ext+ext_ind )
              .ifEmpty { error "Cannot find any bai/crai file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}."+params.ext+ext_ind,""), path ] }

    // building bam-bai pairs
    normal_bam_bai = normal_bams
    	      .phase(normal_bais)
    	      .map { normal_bam, normal_bai -> [ normal_bam[0], normal_bam[1], normal_bai[1] ] }

    // building 4-uplets corresponding to {tumor_bam, tumor_bai, normal_bam, normal_bai}
    tn_bambai = tumor_bam_bai
	      .phase(normal_bam_bai)
	      .map {tumor_bb, normal_bb -> [ tumor_bb[0], tumor_bb[1], tumor_bb[2], normal_bb[1], normal_bb[2] ] }
    // here each element X of tn_bambai channel is a 4-uplet. X[0] is the tumor bam, X[1] the tumor bai, X[2] the normal bam and X[3] the normal bai.
    if(params.estimate_contamination){
	    pairsT4cont = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.'+params.ext )
                             .map {  path -> [ path.name.replace("${params.suffix_tumor}."+params.ext,""), 'T',
                             file(path), file(path +ext_ind) ] }
        pairsN4cont = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.'+params.ext )
                             .map {  path -> [ path.name.replace("${params.suffix_normal}."+params.ext,""), 'N',
                             file(path), file(path +ext_ind ) ] }
		                    .unique()
	    pairs4cont  = pairsT4cont.concat( pairsN4cont )
    }
}


/* manage input positions to call (bed or region or whole-genome) */
  if (params.region) {
      input_region = 'region'
  } else if (params.bed) {
      input_region = 'bed'
      bed = file(params.bed)
  } else {
      input_region = 'whole_genome'
}

//genotyping mode
if(params.genotype){
    pairs2 = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
                       .map{ row -> [ row.sample , row.preproc, file(row.tumor), 
                       file(row.tumor+ext_ind), file(row.normal), 
                       file(row.normal+ext_ind), file(row.vcf) ] }

    pairs2.branch{
                    bam2preproc: it[1]=="yes"
                    bam2nopreproc: it[1]!="yes"
                 }
        .set{ bams }

//pre-processing of RNAseq data
process RNAseq_preproc_fixMCNDN_fixMQ{
    memory params.mem+'GB'
    cpus params.cpu
    tag { sample }

    input:
    set val(sample), val(preproc), file(bam), file(bai), file(bamN), file(baiN), file(vcf) from bams.bam2preproc

    output:
    set val(sample), file("*_MCNDNfixed_rehead.bam"), file("*_MCNDNfixed_rehead.bai"), file(bamN), file(baiN), file(vcf) into bampreproc_mcndn

    shell:
    '''
    if [ -L "None" ]; then unlink None; unlink None.bai; touch None;touch None.bai; fi
    if [ -L "none" ]; then unlink none; unlink none.bai; touch none;touch none.bai; fi
    SM=`samtools view -H !{bam} | grep "^@RG" | head -1 | awk '{print $NF}' | cut -c 4-`
    python !{baseDir}/bin/correctNDN.py !{bam} !{sample}_$SM"_MCNDNfixed.bam"
    samtools view -H !{sample}_$SM"_MCNDNfixed.bam" | sed -e "s/SM:"$SM"/SM:"$SM"_MCNDNfixed/" | samtools reheader - !{sample}_$SM"_MCNDNfixed.bam" > !{sample}_$SM"_MCNDNfixed_rehead.bam"
    samtools index !{sample}_$SM"_MCNDNfixed_rehead.bam" !{sample}_$SM"_MCNDNfixed_rehead.bai"
    '''
}

process RNAseq_preproc_split{
    memory params.mem+'GB'
    cpus '2'
    tag { sample }

    input:
    set val(sample), file(bam), file(bai), file(bamN), file(baiN), file(vcf) from bampreproc_mcndn
    file fasta_ref_RNA
    file fasta_ref_RNA_fai
    file fasta_ref_RNA_dict

    output:
    set val(sample), file("*_split*.bam"), file("*_split*.bai"), file(bamN), file(baiN), file(vcf) into bampreproc

    shell:
    new_tag = sample+"_MCNDNfixed_split"
    '''
    SM=`samtools view -H !{bam} | grep "^@RG" | head -1 | awk '{print $NF}' | cut -c 4-`
    gatk SplitNCigarReads --java-options "-Xmx!{params.mem}G -Djava.io.tmpdir=$PWD" --add-output-sam-program-record  -fixNDN true -R !{fasta_ref_RNA} -I !{bam} -O !{new_tag}_$SM.bam
    '''
}

bam2nopreproc2 = bams.bam2nopreproc.map{ row -> tuple(row[0], row[2], row[3], row[4], row[5], row[6]) }
bams2 = bam2nopreproc2.concat(bampreproc)

tn_bambaivcf0 = bams2.groupTuple(by: 0)
                            .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0] , row[5][0] ) }

tn_bambaivcf0.into{tn_bambaivcf; tn_bambaivcf4print}

process genotype{
    memory params.mem+'GB'
    cpus params.cpu
    tag { sample }

    input:
    set val(sample), file(bamT), file(baiT), file(bamN), file(baiN), file(vcf)  from tn_bambaivcf
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_dict
    file PON
    file PON_tbi
    file known_snp
    file known_snp_tbi

    output:
    set val(sample), file(vcf) , file("${printed_tag}*.vcf") into mutect_geno
    set val(sample), file("${printed_tag}*stats*") into mutect_output2

    publishDir "${params.output_folder}/stats", mode: 'copy', pattern: '{*stats*}' 
  
    shell:
    printed_tag = "${sample}"
    if("${params.gatk_version}" == "4"){
    input_t=""
    for( bamTi in bamT ){
	input_t=input_t+" -I ${bamTi}"
    }
    if(bamN.baseName == 'None' ) {
	input_n=" "
    }else{
        input_n="-I ${bamN} -normal \$normal_name"
    }
    if (params.PON) {
        PON_option = "--panel-of-normals ${PON}"
    } else {
        PON_option = ""
    }
    '''
    !{baseDir}/bin/prep_vcf_bed.sh !{known_snp} !{PON}
    normal_name=`samtools view -H !{bamN} | grep "^@RG" | head -1 | awk '{print $NF}' | cut -c 4-`
    gatk IndexFeatureFile -I !{vcf}
    gatk Mutect2 --java-options "-Xmx!{params.mem}G" -R !{fasta_ref} !{known_snp_option} !{PON_option} !{input_t} !{input_n} \
    -O !{printed_tag}_genotyped.vcf !{params.mutect_args} --alleles !{vcf} -L regions.bed --disable-read-filter NonChimericOriginalAlignmentReadFilter --disable-read-filter NotDuplicateReadFilter \
    --disable-read-filter ReadLengthReadFilter --disable-read-filter WellformedReadFilter \
    --force-call-filtered-alleles --genotype-filtered-alleles --genotype-germline-sites --genotype-pon-sites --active-probability-threshold 0.000 --min-base-quality-score 0 --initial-tumor-lod -100000000000  --tumor-lod-to-emit \
    -100000000000 --force-active --max-reads-per-alignment-start 0
    '''
    }
}

process CompressAndIndex {
    tag { tumor_normal_tag }

    input:
    set val(tumor_normal_tag), file(vcf) , file(vcf_geno) from mutect_geno

    output:
    set val(tumor_normal_tag), file("*.vcf.gz"), file("*.vcf.gz.tbi") into res_filtered_PASS

    publishDir params.output_folder, mode: 'copy'

    shell:
    vcf_name = vcf_geno[0].baseName
    '''
    bcftools view -O z !{vcf_geno[0]} > !{vcf_name}.vcf.gz
    bcftools index -t !{vcf_name}.vcf.gz
    '''
}

}else{
/* process to create a bed file from region or from faidx if whole-genome, otherwise return the input bed file */
process bed {
      output:
      file "temp.bed" into outbed

      shell:
      if (input_region == 'region')
      '''
      echo !{params.region} | sed -e 's/[:|-]/	/g' > temp.bed
      '''

      else if (input_region == 'bed')
      '''
      ln -s !{bed} temp.bed
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
    set val(sample), file(bamT), file(baiT), file(bamN), file(baiN), file(bed)  from tn_bambai.combine(split_bed1)
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_dict
    file PON
    file PON_tbi
    file jar2
    file jar
    file known_snp
    file known_snp_tbi
    file cosmic

    output:
    set val(sample), file("${printed_tag}_*.vcf") into mutect_output1
    set val(sample), file("${printed_tag}*stats*") into mutect_output2
    set val(sample), file("*_f1r2.tar.gz") optional true into f1r2

    shell:
    bed_tag0 = bed.baseName 
    bed_tag = bed_tag0.replaceAll("[^a-zA-Z0-9 _-]+","")
    printed_tag = "${sample}_" + bed_tag
    if("${params.gatk_version}" == "4"){
    input_t=""
    for( bamTi in bamT ){
	input_t=input_t+" -I ${bamTi}"
    }
    if(bamN.baseName == 'None' )  input_n=" "
    else input_n="-I ${bamN} -normal \$normal_name"
    if (params.PON) {
        PON_option = "--panel-of-normals ${PON}"
    } else {
        PON_option = ""
    }
    '''
    normal_name=`samtools view -H !{bamN} | grep "^@RG" | head -1 | awk '{print $NF}' | cut -c 4-`
    gatk Mutect2 --java-options "-Xmx!{params.mem}G" -R !{fasta_ref} !{known_snp_option} !{PON_option} \
    !{input_t} !{input_n} -O !{printed_tag}_calls.vcf -L !{bed} !{params.mutect_args} --f1r2-tar-gz !{printed_tag}_f1r2.tar.gz
    '''
    }else{
    '''
    if [ "!{mutect_version}" == "2" ]
		        then
		            !{params.java} -Xmx!{params.mem}g -jar !{jar2} -T MuTect2  -R !{fasta_ref} !{known_snp_option} \
                    !{cosmic_option} -I:normal !{bamN} -I:tumor !{bamT} \
                    -o "!{sample}_!{bed_tag}_calls.vcf" -L !{bed} !{params.mutect_args}
			    touch !{sample}_!{bed_tag}_calls_stats.txt 
	        	else
	        	    !{params.java} -Xmx!{params.mem}g -jar !{jar} --analysis_type MuTect --reference_sequence !{fasta_ref} \
                    !{known_snp_option} !{cosmic_option} --intervals !{bed} \
                    --input_file:tumor !{bamT} --input_file:normal !{bamN} --out "!{sample}_!{bed_tag}_calls_stats.txt" \
                    --vcf "!{sample}_!{bed_tag}_calls.vcf" !{params.mutect_args}
    fi
    '''
    }
}

beds_length = count_split_bed.count().val

process mergeMuTectOutputs {

    tag { tumor_normal_tag }

    publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
                 if (filename.indexOf(".stats") > 0)   "stats/$filename"
            else if (filename.indexOf(".vcf") > 0)     "intermediate_calls/raw_calls/$filename"
    }
    
    input:
    set val(tumor_normal_tag), file(vcf_files) from mutect_output1.groupTuple(size: beds_length)
    set val(tumor_normal_tag2), file(txt_files) from mutect_output2.groupTuple(size: beds_length)

    output:
    set val(tumor_normal_tag), file("${tumor_normal_tag}_calls.vcf"), file("${tumor_normal_tag}_calls.vcf.stats") into res_merged

    shell:
    input_stats=""
    for( txt in txt_files ){
        input_stats=input_stats+" -stats ${txt}"
    }

    '''
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
    mv header.txt !{tumor_normal_tag}_calls.vcf

    # MERGE TXT FILES if mutect1
    if [ "!{mutect_version}" == "1" ]
        then
          head -n2 `ls -1 *.txt | head -1` > header.txt
          sed -i '1,2d' *calls_stats.txt
          cat *calls_stats.txt >> header.txt
          mv header.txt !{tumor_normal_tag}_calls.vcf.stats
        else
            if [ "!{params.gatk_version}" == "4" ] 
            then
	            gatk MergeMutectStats !{input_stats} -O !{tumor_normal_tag}_calls.vcf.stats
            else  
                touch !{tumor_normal_tag}_calls.vcf.stats
            fi
    fi  
    '''
}

if(params.gatk_version=="4"){
    if(params.filter_readorientation){
	//println("Filtering output")
	process ReadOrientationLearn {
            tag { tumor_normal_tag }

            publishDir params.output_folder+'/stats', mode: 'copy'

            input:
            set val(tumor_normal_tag), file(f1r2_files) from f1r2.groupTuple()

            output:
            set val(tumor_normal_tag), file("*model.tar.gz") into ROmodel

            shell:
            input_f1r2=""
            for( f1r2 in f1r2_files ){
                input_f1r2=input_f1r2+" -I ${f1r2}"
            }
            '''
            gatk LearnReadOrientationModel !{input_f1r2} -O !{tumor_normal_tag}_read-orientation-model.tar.gz
            '''
        }
        res_merged_RO = res_merged.join(ROmodel)
    }else{
        res_merged_RO = res_merged.map{row -> [row[0], row[1], row[2] , null]}
    }

if(params.estimate_contamination){
	process ContaminationEstimationPileup {
	    memory params.mem+'GB'
        cpus '2'
	    tag { tumor_normal_tag }

	    input:
	    set val(tumor_normal_tag), val(TN), file(bam), file(bai) from pairs4cont
	    file fasta_ref
	    file fasta_ref_fai
	    file fasta_ref_dict
	    file snp_contam
	    file snp_contam_tbi

	    output:
	    set val(tumor_normal_tag), val(TN) , file("*.table") into pileups
	    

	    shell:
	    basename=bam.baseName
	
	    if(basename != 'None'){
	    '''
	    gatk --java-options "-Xmx!{params.mem}G" GetPileupSummaries -R !{fasta_ref} -I !{bam} -V !{snp_contam} -L !{snp_contam} -O !{basename}_pileups.table
	    '''
	    }else{
	    '''
	    touch empty.table
	    '''
	}
	}
	pileupsN0   = Channel.create()
	pileupsT0   = Channel.create()
	pileups.choice( pileupsN0,pileupsT0 ) { a -> a[1] == "N" ? 0 : 1 }

	pileupsN   = Channel.create()
	pileupsT   = Channel.create()
	pileupsN4pr   = Channel.create()
        pileupsT4pr   = Channel.create()
	pileupsN0.into(pileupsN,pileupsN4pr)
	pileupsT0.into(pileupsT,pileupsT4pr)
	pileups4cont = pileupsN.cross(pileupsT)
			       .map { row -> tuple(row[0][0] , row[0][2] , row[1][2]  ) }
				//pileups.groupTuple(by: 0)
			     // .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0]  ) }


	
	pileups4contpr = pileupsN4pr.cross(pileupsT4pr)
				    .map { row -> tuple(row[0][0] , row[0][2] , row[1][2]  ) }
				    .subscribe{ row -> println "${row}" }

	process ContaminationEstimation {
   	    memory params.mem+'GB'
	    cpus '2'
	    tag { tumor_normal_tag }
	    publishDir params.output_folder+'/stats', mode: 'copy'

	    input:
	    set val(tumor_normal_tag), file(pileupN) , file(pileupT) from pileups4cont

	    output:
	    set val(tumor_normal_tag), file("*contamination.table") into contam

	    shell:
	    basename=pileupT.baseName
        if(pileupN.baseName == 'empty' )  input_n=" "
        else input_n="-matched $pileupN"
	    '''
	    gatk --java-options "-Xmx!{params.mem}G" CalculateContamination -I !{pileupT} !{input_n} -O !{basename}_calculatecontamination.table
	    '''
	}
   res_merged_RO_contam = res_merged_RO.join(contam)
}else{
   res_merged_RO_contam = res_merged_RO.map{row -> [row[0], row[1], row[2] , row[3], null]}
}

process FilterMuTectOutputs {
    tag { tumor_normal_tag }

    input:
    set val(tumor_normal_tag), file(vcf), file(stats), file(ROmodel), file(contam_tables) from res_merged_RO_contam
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_gzi
    file fasta_ref_dict

    output:
    set val(tumor_normal_tag), file("*filtered.vcf*") into res_filtered

    publishDir params.output_folder+'/intermediate_calls/filtered', mode: 'copy'

    shell:
    RO=""
    if(params.filter_readorientation){
	for(model in ROmodel){
		RO=RO+"--ob-priors ${model} "
	}
    }
    contam=""
    if(params.estimate_contamination){
	for(contmp in contam_tables){
		contam=contam+"--contamination-table ${contmp} "
	}
    }
    '''
    gatk FilterMutectCalls -R !{fasta_ref} -V !{vcf} !{contam} !{RO} -O !{tumor_normal_tag}_filtered.vcf
    '''
}

process FilterMuTectOutputsOnPass {
    tag { tumor_normal_tag }

    input:
    set val(tumor_normal_tag), file(vcf_filtered) from res_filtered

    output:
    file("*_PASS.vcf*") into res_filtered_PASS

    publishDir params.output_folder, mode: 'copy'

    shell:
    vcf_name = vcf_filtered[0].baseName
    '''
    bcftools view -f PASS -O z !{vcf_filtered[0]} > !{vcf_name}_PASS.vcf.gz
    bcftools index -t !{vcf_name}_PASS.vcf.gz
    '''
}

}

}
