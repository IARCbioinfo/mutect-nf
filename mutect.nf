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
params.mem           = 16
params.cpu           = 4
params.output_folder = "mutect_results"
params.mutect_args   = ""
params.nsplit        = 1
params.bed           = null
params.java          = "java"
params.snp_vcf       = null
params.snp_contam    = null
params.tn_file       = null
params.PON           = null
params.filter_readorientation = null
params.genotype      = null
params.ref_RNA       = "NO_REF_RNA_FILE"

params.help = null

log.info "" 
log.info "--------------------------------------------------------"
log.info "  mutect-nf 2.2.0: Mutect pipeline for somatic variant calling with Nextflow DSL2"
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
    log.info '    --tn_file            FILE                    input tabulation-separated values file with columns sample (sample name),'
    log.info '                                                 tumor (full path to tumor bam), normal (full path to matched normal bam);'
    log.info '                                                 optionally (for --genotype mode), columns preproc (is the bam RNAseq needing'
    log.info '                                                 preprocessing: yes or no) and vcf (full path to vcf file containing alleles to genotype)'
    log.info '                                                 where each line contains the path of two matched BAM files.'    
    log.info '    --ref                FILE (with indexes)     Reference fasta file.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --bed                FILE                    Bed or region file containing intervals.'
    log.info '    NOTE: if neither bed or region file is provided, will perform the calling on whole genome, based on the faidx file.'
    log.info '    --nsplit             INTEGER                 Split the region for calling in nsplit pieces and run in parallel (default: 1).'
    log.info '    --snp_vcf            FILE                    VCF file with known variants and frequency (e.g., from gnomad).'
    log.info '    --snp_contam 	       FILE                    VCF file with known germline variants to genotype for contamination estimation'
    log.info '    --PON                FILE                    VCF file of GATK panel of normals used to filter calls.'
    log.info '    --mutect_args        STRING                  Arguments you want to pass to mutect.'
    log.info '                                                 WARNING: form is " --force_alleles " with spaces between quotes.'
    log.info '    --ref_RNA 	       PATH                    fasta reference for preprocessing RNA (required when preproc column contains yes'
    log.info '                                                 in input tn_file).'
    log.info '    --cpu                INTEGER                 Number of cpu used (default: 4).'
    log.info '    --mem                INTEGER                 Java memory passed to mutect in GB (default: 8).'
    log.info '    --output_folder      FOLDER                  Output folder (default: mutect_results).'
    log.info '    --java               PATH                    Name of the JAVA command  (default: java).'
    log.info ''
    log.info 'Flags:'
    log.info '    --filter_readorientation                     Run extra step learning read orientation model and using it to filter reads.'
    log.info '    --genotype                                   Use genotyping from vcf mode instead of usual variant calling;'
    log.info '                                                 requires tn_file with vcf column and gatk4, and if RNA-seq included, requires preproc column'
    log.info ''
    exit 0
}else{
    /* Software information */
    log.info "mem                    = ${params.mem}"
    log.info "cpu                    = ${params.cpu}"
    log.info "output_folder          = ${params.output_folder}"
    log.info "mutect_args            = ${params.mutect_args}"
    log.info "nsplit                 = ${params.nsplit}"
    log.info "bed                    = ${params.bed}"
    log.info "java                   = ${params.java}"
    log.info "snp_vcf                = ${params.snp_vcf}"
    log.info "snp_contam             = ${params.snp_contam}"
    log.info "tn_file                = ${params.tn_file}"
    log.info "PON                    = ${params.PON}"
    log.info "filter_readorientation = ${params.filter_readorientation}"
    log.info "genotype               = ${params.genotype}"
    log.info "ref                    = ${params.ref}"
    log.info "ref_RNA                = ${params.ref_RNA}"	
}


/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/


// reference file and its indexes
ref = tuple file(params.ref), file(params.ref+'.fai'), file(params.ref+'.gzi'), 
    file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )

fasta_ref_fai = file( params.ref+'.fai' )

// reference file and its indexes for genotyping
ref_RNA = tuple file( params.ref_RNA ), 
    file( params.ref_RNA+'.fai' ),
    file( params.ref_RNA.replace(".fasta",".dict").replace(".fa",".dict") )

ref_RNA = (params.ref_RNA == "NO_REF_RNA_FILE") ? ref : ref_RNA

// know snp VCFs
snp_vcf = params.snp_vcf ? (tuple file(params.snp_vcf), file(params.snp_vcf+".tbi")) : (tuple file("NO_SNP"), file("NO_SNP.tbi"))
snp_vcf_option = params.snp_vcf ? "" : "--germline-resource ${snp_vcf.get(0)}"

// contamination VCFs
snp_contam = ((params.snp_contam) && (!params.snp_contam=="null")) ? (tuple file(params.snp_contam), file(params.snp_contam+'.tbi')) : null
estimate_contamination = snp_contam ? true : false
log.info "snp_contam = ${snp_contam}"


// Pannel Of Normal
PON = ((params.PON) && (!params.PON=="null")) ? (tuple file(params.PON), file(params.PON +'_TBI')) : (tuple file("NO_FILE"), file("NO_FILE_TBI"))
PON_option = ((params.PON) && (!params.PON=="null")) ? "--panel-of-normals ${PON.get(0)}" : ""

// manage input positions to call (bed or region or whole-genome)
intervals = params.bed ? params.bed : 'whole_genome'


/***************************************************************************************/
/************************  Process : RNAseq_preproc_fixMCNDN_fixMQ *********************/
/***************************************************************************************/

process RNAseq_preproc_fixMCNDN_fixMQ{
   
    memory params.mem+'GB'
    cpus params.cpu

    input:
        tuple val(sample), val(preproc), path(bam), path(bai), path(bamN), path(baiN), path(vcf)

    output:
        tuple val(sample), path("*_MCNDNfixed_rehead.bam"), path("*_MCNDNfixed_rehead.bai"), path(bamN), path(baiN), path(vcf)

    shell:
        """
        if [ -L "None" ]; then unlink None; unlink None.bai; touch None;touch None.bai; fi
        if [ -L "none" ]; then unlink none; unlink none.bai; touch none;touch none.bai; fi
        SM=`samtools view -H $bam | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        python $baseDir/bin/correctNDN.py $bam ${sample}_\$SM"_MCNDNfixed.bam"
        samtools view -H ${sample}_\$SM"_MCNDNfixed.bam" | sed -e "s/SM:"\$SM"/SM:"\$SM"_MCNDNfixed/" | samtools reheader - ${sample}_\$SM"_MCNDNfixed.bam" > ${sample}_$SM"_MCNDNfixed_rehead.bam"
        samtools index ${sample}_\$SM"_MCNDNfixed_rehead.bam" ${sample}_\$SM"_MCNDNfixed_rehead.bai"
        """
}


/***************************************************************************************/
/************************  Process : RNAseq_preproc_split ******************************/
/***************************************************************************************/

process RNAseq_preproc_split{

    memory params.mem+'GB'
    cpus '2'

    input:
        tuple val(sample), path(bam), path(bai), path(bamN), path(baiN), path(vcf)
        tuple path(fasta_ref_RNA), path(fasta_ref_RNA_fai), path(fasta_ref_RNA_dict)

    output:
        tuple val(sample), path("*_split*.bam"), path("*_split*.bai"), path(bamN), path(baiN), path(vcf)

    shell:
        new_tag = sample+"_MCNDNfixed_split"
        """
        SM=`samtools view -H ${bam} | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        gatk SplitNCigarReads --java-options "-Xmx${params.mem}G -Djava.io.tmpdir=$PWD" --add-output-sam-program-record  -fixNDN true -R $fasta_ref_RNA -I $bam -O ${new_tag}_\$SM.bam
        """

    stub:
        new_tag = sample+"_MCNDNfixed_split"
        """
        SM=`samtools view -H ${bam} | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        touch ${new_tag}_\$SM.bam ${new_tag}_\$SM.bam.bai
        """
}

/***************************************************************************************/
/************************  Process : genotype ******************************************/
/***************************************************************************************/

process genotype{

    memory params.mem+'GB'
    cpus params.cpu

    publishDir "${params.output_folder}/stats", mode: 'copy', pattern: '{*stats*}' 

    input:
        tuple val(sample), path(bamT), path(baiT), path(bamN), path(baiN), path(vcf)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_gzi), path(fasta_ref_dict)
        tuple path(PON), path(PON_tbi)
        tuple path(snp_vcf), path(snp_vcf_tbi)

    output:
        tuple val(sample), path(vcf) , path("${sample}*.vcf"), emit : vcfs
        tuple val(sample), path("${sample}*stats*")

    shell:
        input_t = "-I " + bamT.join(" -I ")
        input_n = (bamN.baseName == 'None') ? "" : "-I ${bamN} -normal \$normal_name"
        //PON_option = params.PON ? "--panel-of-normals ${PON}" : ""
        """
        ${baseDir}/bin/prep_vcf_bed.sh $snp_vcf $PON
        normal_name=`samtools view -H $bamN | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        gatk IndexFeatureFile -I $vcf
        gatk Mutect2 --java-options "-Xmx${params.mem}G" -R $fasta_ref $snp_vcf_option $PON_option $input_t $input_n \
        -O ${sample}_genotyped.vcf $params.mutect_args --alleles $vcf -L regions.bed --disable-read-filter NonChimericOriginalAlignmentReadFilter --disable-read-filter NotDuplicateReadFilter \
        --disable-read-filter ReadLengthReadFilter --disable-read-filter WellformedReadFilter \
        --force-call-filtered-alleles --genotype-filtered-alleles --genotype-germline-sites --genotype-pon-sites --active-probability-threshold 0.000 --min-base-quality-score 0 --initial-tumor-lod -100000000000  --tumor-lod-to-emit \
        -100000000000 --force-active --max-reads-per-alignment-start 0
        """

    stub:
        """
        touch ${sample}_genotyped.vcf ${sample}_genotyped.vcf.stats
        """
}


/***************************************************************************************/
/************************  Process : CompressAndIndex **********************************/
/***************************************************************************************/

process CompressAndIndex {

    publishDir params.output_folder, mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf) , path(vcf_geno)

    output:
        tuple val(tumor_normal_tag), path("*.vcf.gz"), path("*.vcf.gz.tbi")

    shell:
        vcf_name = vcf_geno[0].baseName
        """
        bcftools view -O z $vcf_geno[0] > ${vcf_name}.vcf.gz
        bcftools index -t ${vcf_name}.vcf.gz
        """
    
    stub:
        vcf_name = vcf_geno[0].baseName
        """
        touch ${vcf_name}.vcf.gz ${vcf_name}.vcf.gz.tbi
        """
}

/***************************************************************************************/
/************************  Process : make_bed ******************************************/
/***************************************************************************************/

/* process to create a bed file from region or from faidx if whole-genome, otherwise return the input bed file */
/* FORMALLY CALLED bed !!!!! */

process make_bed {

    input:
        path fasta_ref_fai

    output:
        file "temp.bed"

    shell:
        if (intervals =~ /^[0-9a-zA-Z]+:[0-9]+-[0-9]+/){
        """
        echo $intervals | sed -e 's/[:|-]/	/g' > temp.bed
        """
        
        }else if( intervals =~ /^[0-9a-zA-Z]+\t[0-9]+\t[0-9]+/){
        """
        ln -s $intervals temp.bed
        """

        }else{
        """
        cat $fasta_ref_fai | awk '{print \$1"	"0"	"\$2 }' | grep -v -P "alt|random|Un|chrEBV|HLA" > temp.bed
        """
        }

    stub:
        """
        touch temp.bed
        """
  }


/***************************************************************************************/
/************************  Process : split_bed *****************************************/
/***************************************************************************************/

process split_bed {

    input:
        path bed

    output:
        path '*_regions.bed'

    shell:
        """
        grep -v '^track' $bed | sort -T \$PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print \$1" "\$2" "\$3}' | cut_into_small_beds.r $params.nsplit
        """
    
    stub:
        """
        grep -v '^track' $bed | sort -T \$PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print \$1" "\$2" "\$3}' | cut_into_small_beds.r $params.nsplit
        """
}


/***************************************************************************************/
/************************  Process : mutect ********************************************/
/***************************************************************************************/


process mutect {
    
    memory params.mem+'GB' 
    cpus params.cpu

    input:
        tuple val(sample), path(bamT), path(baiT), path(bamN), path(baiN), path(bed)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_gzi), path(fasta_ref_dict)
        tuple path(PON), path(PON_tbi)
        tuple path(snp_vcf), path(snp_vcf_tbi)

    output:
        tuple val(sample), path("${printed_tag}_*.vcf"), path("${printed_tag}*stats*"), emit : vcfs
        tuple val(sample), path("*_f1r2.tar.gz"), emit : f1r2

    shell:
        bed_tag = bed.baseName.replaceAll("[^a-zA-Z0-9 _-]+","")
        printed_tag = "${sample}_" + bed_tag
        input_t = "-I " + bamT.join(" -I ")
        input_n = (bamN.baseName == 'None') ? "" : "-I ${bamN} -normal \$normal_name"
        """
        normal_name=`samtools view -H $bamN | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        echo \$normal_name
        gatk Mutect2 --java-options "-Xmx${params.mem}G" -R $fasta_ref $snp_vcf_option $PON_option \
        $input_t $input_n -O ${printed_tag}_calls.vcf -L $bed $params.mutect_args --f1r2-tar-gz ${printed_tag}_f1r2.tar.gz
        """

    stub:
        bed_tag = bed.baseName.replaceAll("[^a-zA-Z0-9 _-]+","")
        printed_tag = "${sample}_" + bed_tag
        """
        touch ${printed_tag}_calls.vcf ${printed_tag}_calls.vcf.stats ${printed_tag}_calls_f1r2.tar.gz
        """
}

/***************************************************************************************/
/************************  Process : mergeMuTectOutputs ********************************/
/***************************************************************************************/

process mergeMuTectOutputs {

    publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
        if (filename.indexOf(".stats") > 0) "stats/$filename"
        else if (filename.indexOf(".vcf") > 0) "intermediate_calls/raw_calls/$filename"
    }
    
    input:
        tuple val(tumor_normal_tag), path(vcf_files), path(txt_files)

    output:
        tuple val(tumor_normal_tag), path("${tumor_normal_tag}_calls.vcf"), path("${tumor_normal_tag}_calls.vcf.stats")

    shell:
        input_stats = "-stats " + txt_files.join(" -stats ")
        """
        # MERGE VCF FILES
        sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt

        # Determine sort options based on the availability of version-sort in sort command
        sort_ops=\$(sort --help | grep -q 'version-sort' && echo "-k1,1V" || echo "-k1,1d")

        # Concatenate VCF contents and sort
        grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -T \$PWD -t '	' \$sort_ops -k2,2n >> header.txt

        # Rename the merged file
        mv header.txt ${tumor_normal_tag}_calls.vcf

        # MERGE STAT FILES
        gatk MergeMutectStats $input_stats -O ${tumor_normal_tag}_calls.vcf.stats
        """

    stub:
        """
        touch ${tumor_normal_tag}_calls.vcf ${tumor_normal_tag}_calls.vcf.stats
        """
}

/***************************************************************************************/
/************************  Process : ReadOrientationLearn ******************************/
/***************************************************************************************/

process ReadOrientationLearn {
            
    publishDir params.output_folder+'/stats', mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(f1r2_files)

    output:
        tuple val(tumor_normal_tag), path("*model.tar.gz")

    shell:
        input_f1r2 = "-I " + f1r2_files.join(" -I ")
        """
        gatk LearnReadOrientationModel ${input_f1r2} -O ${tumor_normal_tag}_read-orientation-model.tar.gz
        """
    
    stub:
        """
        touch ${tumor_normal_tag}_read-orientation-model.tar.gz
        """
}

/***************************************************************************************/
/************************  Process : ContaminationEstimationPileup *********************/
/***************************************************************************************/

process ContaminationEstimationPileup {

    cpus '16'
    memory '164 GB'

	input:
        tuple val(tumor_normal_tag), val(TN), path(bam), path(bai)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_gzi), path(fasta_ref_dict)
        tuple path(snp_contam), path(snp_contam_tbi)

	output:
	    tuple val(tumor_normal_tag), val(TN) , path("*.table")
	    
	shell:
        """
        gatk --java-options "-Xmx${params.mem}G" GetPileupSummaries -R ${fasta_ref} -I $bam -V $snp_contam -L $snp_contam -O ${bam.baseName}_pileups.table
        """
    
    stub:
        """
        touch ${bam.baseName}_pileups.table
        """
}


/***************************************************************************************/
/************************  Process : process ContaminationEstimation *******************/
/***************************************************************************************/

process ContaminationEstimation {
   	
    memory params.mem+'GB'
	cpus '2'
	    
    publishDir params.output_folder+'/stats', mode: 'copy'

	input:
	    tuple val(tumor_normal_tag), path(pileupN), path(pileupT)

	output:
	    tuple val(tumor_normal_tag), path("*contamination.table")
	
    shell:
        input_n = (pileupN.baseName == 'empty' ) ? " " : "-matched $pileupN"
	    """
	    gatk --java-options "-Xmx${params.mem}G" CalculateContamination -I $pileupT $input_n -O ${pileupT.baseName}_calculatecontamination.table
	    """
    
    stub:
        """
        touch ${pileupT.baseName}_calculatecontamination.table
        """
}


/***************************************************************************************/
/************************  Process : process FilterMuTectOutputs ***********************/
/***************************************************************************************/

process FilterMuTectOutputs {

    publishDir params.output_folder+'/intermediate_calls/filtered', mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf), path(stats), path(ROmodel), path(contam_tables)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_gzi), path(fasta_ref_dict)

    output:
        tuple val(tumor_normal_tag), path("*filtered.vcf*")

    shell:
        RO = (ROmodel.baseName=="NO_ROmodel") ? "": "--ob-priors " + ROmodel.join(" --ob-priors ")
        contam = (contam_tables.baseName == "NO_contam") ? "" : "--contamination-table " + contam_tables.join(" --contamination-table ")
        """
        gatk FilterMutectCalls -R $fasta_ref -V $vcf $contam $RO -O ${tumor_normal_tag}_filtered.vcf
        """
    
    stub:
        """
        touch ${tumor_normal_tag}_filtered.vcf
        """
}

/***************************************************************************************/
/************************  Process : process FilterMuTectOutputsOnPass *****************/
/***************************************************************************************/

process FilterMuTectOutputsOnPass {

    publishDir params.output_folder, mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf_filtered)

    output:
        path("*_PASS.vcf*")

    shell:
        vcf_name = vcf_filtered[0].baseName
        """
        bcftools view -f PASS -O z ${vcf_filtered[0]} > ${vcf_name}_PASS.vcf.gz
        bcftools index -t ${vcf_name}_PASS.vcf.gz
        """

    stub:
        vcf_name = vcf_filtered[0].baseName
        """
        touch ${vcf_name}_PASS.vcf.gz
        """
}


/***************************************************************************************/
/************************  Workflow : main *********************************************/
/***************************************************************************************/


workflow {

    //load input files
    assert(params.tn_file) : "Error: please provide tn_file"
    pairs = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
        .map{ row -> 
        assert(row.sample) : "Error: sample tag missing check your input_file"
        assert(row.tumor) : "Error: tumor bam is missing check your input_file"
        assert(row.normal) : "Error: normal bam is missing check your input_file"
        tuple(
            row.sample,
            row.preproc ? row.preproc : "no", 
            file(row.tumor), file("${row.tumor}.{bai,crai}")[0],
            file(row.normal), file("${row.normal}.{bai,crai}")[0],
            row.vcf ? file(row.vcf) : file("NO_VCF"),
            row.vcf ? file(row.vcf+".tbi")[0] : file("NO_VCF_TBI")
        )}
    

    // contamination
    if(estimate_contamination){
        pairsT4cont = pairs.map{ row -> tuple( row[0] , 'T' , row[2], row[3] )}
        pairsN4cont = pairs.map{ row -> tuple( row[0] , 'N' , row[4], row[5] )}.unique()
        pairs4cont  = pairsT4cont.concat( pairsN4cont )
    }

    if(params.genotype){

        bams = pairs.branch{
            preproc: it[1]=="yes"
            nopreproc: it[1]!="yes"
        }
        bampreproc = RNAseq_preproc_fixMCNDN_fixMQ(bams.preproc) | RNAseq_preproc_split
        bamnopreproc = bams.nopreproc.map{ row -> tuple(row[0], row[2], row[3], row[4], row[5], row[6]) }
        
        bams = bamnopreproc.concat(bampreproc).groupTuple(by: 0).map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0] , row[5][0] ) }
        genotype(bams,ref,PON,snp_vcf)
        pass = CompressAndIndex(genotype.out.vcfs)

    } else{

        // mutect
        bams = pairs.groupTuple(by: 0).map { row -> tuple(row[0] , row[2], row[3] , row[4][0] , row[5][0]  ) }
        bed = make_bed(fasta_ref_fai) | split_bed | flatten
        mutect(bams.combine(bed),ref,PON,snp_vcf)

        
        mutectOutput = mergeMuTectOutputs(mutect.out.vcfs.groupTuple(size: params.nsplit ))

        // read orientation
        if(params.filter_readorientation){
            ReadOrientationLearn(mutect.out.f1r2.groupTuple(size: params.nsplit ))
            mutectOutput = mutectOutput.join( ReadOrientationLearn.out )
        }else{
            mutectOutput = mutectOutput.map{row -> [row[0], row[1], row[2] , file("NO_ROmodel") ]}
        }

        // contamination
        if(estimate_contamination){
            pileups4cont  = ContaminationEstimationPileup(pairs4cont,ref,snp_contam)
                .groupTuple(by:[0,1]).groupTuple()
                .map{ group -> tuple( group[0], group[2][0], group[2][1] )}
                .transpose(by:1)
            contamination = ContaminationEstimation(pileups4cont)
            mutectOutput = mutectOutput.join(contamination)
        }else{
            mutectOutput = mutectOutput.map{row -> [row[0], row[1], row[2] , row[3], file("NO_contam") ]}
        }

        // filter
        FilterMuTectOutputs(mutectOutput, ref) | FilterMuTectOutputsOnPass
        
    }


}
