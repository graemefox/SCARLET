#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

nextflow_version="v.0.1"

// default behaviour is to NOT:
// run rapidCNS2 classifier
// run sturgeon classifier 
// run nanoDx classifier
// run NanoPlot

// default the various optional parameters to false
params.rapidcns2 = false
params.sturgeon = false
params.nanoplot = false
params.nanodx = false

log.info """\

        ================================================================
        SCARLET - Nextflow P I P E L I N E ${nextflow_version}
        ================================================================

        INPUTS
        ================================================================
        sample                  : ${params.sample}
        input_bam               : ${params.bam}
        outdir                  : ${params.outdir}
        reference               : ${params.reference}
        annotations             : ${params.annotations}
        threads                 : ${params.threads}
        bam_min_coverage        : ${params.bam_min_coverage}
        min_mgmt_coverage       : ${params.minimum_mgmt_cov}
        rapidcns2               : ${params.rapidcns2}
        sturgeon                : ${params.sturgeon}
        nanoplot                : ${params.nanoplot}
        nanodx                  : ${params.nanodx}

        ================================================================
        To run with SLURM, add -process.executor='slurm' to your nextflow command.
        ================================================================
        """
        .stripIndent()

process get_versions {
    input:

    output:
        path "software_versions1.txt", emit: versions_file1

    script:
        """
        echo -e "Package and version" > software_versions1.txt
        samtools --version | head -n 1 >> software_versions1.txt
        /modkit --version >> software_versions1.txt
        NanoPlot -v >> software_versions1.txt
        /mosdepth --version -v >> software_versions1.txt
        python3 /methylartist/methylartist -v >> software_versions1.txt
        bcftools -v | head -n 1 >> software_versions1.txt
        echo "clairS_to v0.1.0" >> software_versions1.txt
        echo "igv_reports v1.12.0" >> software_versions1.txt
        echo "cnvpytor v1.3.1" >> software_versions1.txt
        vcftools | head -n 2 | tail -n 1 >> software_versions1.txt
        pigz -V >> software_versions1.txt
        gzip -V | head -n 1 >> software_versions1.txt
        /annovar/convert2annovar.pl | grep Version | sed 's/^/annovar: /; s/,//g' >> "software_versions1.txt"
        R --version | head -n 1 >> software_versions1.txt
        source /sturgeon-0.4.4/venv/bin/activate
        sturgeon -v >> software_versions1.txt
        deactivate
        echo "nanoDx - TODO" >> software_versions1.txt
        """
}

process get_versions_outside_docker {
    input:

    output:
        path "software_versions2.txt", emit: versions_file2

    script:
        """
        docker -v | sed 's/,//g' > "software_versions2.txt"
        nextflow -v >> software_versions2.txt
        nextflow run epi2me-labs/wf-human-variation | grep '^Launching' | awk '{for (i=1; i<=NF; i++) if (\$i ~ /^revision:/) print "epi2me-labs/wf-human-variation " \$i " " \$(i+1)}' >> software_versions2.txt
        """
}

process merge_versions_files {
    input:
        path(versions_file1)
        path(versions_file2)

    output:  
        path "all_software_versions.txt", emit: versions_file

    script:
        """
        cat \
        ${versions_file1} \
        ${versions_file2} \
        > all_software_versions.txt
        """
}

process check_bam_has_meth_data {
    input:
        path(input_bam)
        val(threads)
    
    output:
        stdout emit: meth_check

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${input_bam} \
        | grep -m 1 MM:Z
        """
}

process modkit_adjust_mods {
    input:
        path(input_bam)
        val(sample)
        val(threads)

    output:
        path "*_modkit_merge.bam", emit: modkit_merged_bam

    script:
        """
        /modkit \
        adjust-mods \
        --convert h m \
        ${input_bam} \
        ${sample}_modkit_merge.bam \
        --threads ${threads}
        """
}

process index_input_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: indexed_bam
        tuple path(input_bam), path("*.bai"), emit: indexed_bam_tuple

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process index_merged_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        tuple path(input_bam), path("*.bai"), emit: indexed_bam

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process nanoplot {
    input:
        path(input_bam)
        val(threads)
        val(sample)
        path(indexed_bam)

    publishDir("${params.outdir}", mode: 'copy')

    output:
        path "*NanoPlot-report.html", emit: nanoplot

    script:
        """
        NanoPlot \
        -t ${threads} \
        --bam ${input_bam} \
        -p ${sample}
        """
}

process check_mgmt_coverage {
    input:
        path(input_bam)
	path(mgmt_bed)
	val(minimum_mgmt_coverage)
        path(indexed_bam)
        val(threads)

    output:
        stdout emit: mgmt_avg_cov

    script:
        """
        /mosdepth \
        -t ${threads} \
        -n \
        --by ${mgmt_bed} \
        mgmt_cov ${input_bam}        
        cov="\$(grep "^chr10_region" mgmt_cov.mosdepth.summary.txt | awk '{ print \$4 }')"
        echo \${cov}
	"""
}

process draw_mgmt_methylartist {
    maxRetries 3
    errorStrategy { (task.attempt <= maxRetries) ? 'retry' : 'ignore' }

    input:
        tuple path(bam), path(bai)
        path(reference)
        val(params.outdir)
        val(mgmt_avg_cov)
        val(minimum_mgmt_coverage)

    output:
        path "*.png", emit: mgmt_plot optional true
    
    script:        
          if ( mgmt_avg_cov.toFloat() > minimum_mgmt_coverage.toFloat() )
            """
            python3 /methylartist/methylartist \
            locus \
            -i chr10:129466536-129467536 \
            -b ${bam} \
            --ref ${reference} \
            --motif CG \
            --mods m
            """
            else

            """
            echo "dummy_file" > dummy_plot.png
            """
            
}

process mosdepth {
    input:
        val(threads)
        path(targets)
        path(input_bam)
        val(sample)
        path(indexed_bam)

    output:
        path "*.mosdepth.summary.txt", emit: mosdepth_out

    script:
        """
        /mosdepth \
        -t ${threads} \
        -n \
        --by ${targets} \
        --fast-mode \
        ${sample} \
        ${input_bam}
        """
}

process index_subsetted_bam {
    input:
        path(input_bam)
        val(threads)
    
    output:
        tuple path(input_bam), path("*.bai"), emit: subsetted_bam_index
        val(true)

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process human_variation_mods {
    input:
        tuple path(merged_bam), path(bai)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        val(threads)
    
    output:
        path("*/*wf_mods.bedmethyl.gz"), emit: bedmethyl_gz
    
    script:
        """
        nextflow run epi2me-labs/wf-human-variation -r master \
        -profile standard \
        --ref ${reference} \
        --mod \
        --bam ${merged_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --threads ${threads} \
        """ 
}

process human_variation_cnv {

    // must use --use_qdnaseq for the CNV module as SPECTRE not suitable
    // for data generated by adaptive sampling
    // see https://github.com/epi2me-labs/wf-human-variation/issues/198

    // low depth of coverage will cause this module to fail - not handled gracefully by the sub-workflow
    // other than providing the CNV report, the outputs are not used and so errors can be ignored
    errorStrategy 'ignore'

    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        val(threads)

    publishDir("${params.outdir}", mode: 'copy')

    output:
        path("*/*wf-human-cnv-report.html"), emit: human_variation_cnv_report

    script:
        """
        nextflow run epi2me-labs/wf-human-variation -r master \
        -profile standard \
        --ref ${reference} \
        --cnv \
        --use_qdnaseq \
        --bam ${input_bam} \
        --out_dir wf-human-variation_reports \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --threads ${threads}
        """
}

process human_variation_sv {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        path(indexed_bam)
        val(threads)

    publishDir("${params.outdir}", mode: 'copy')

    output:
        path("*/*wf-human-sv-report.html"), emit: human_variation_sv_report

    script:
        """
        nextflow run epi2me-labs/wf-human-variation -r master \
        -profile standard \
        --ref ${reference} \
        --sv \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --out_dir wf-human-variation_reports \
        --bam_min_coverage ${bam_min_coverage} \
        --threads ${threads} \
        --sniffles_args="--non-germline"
        """
}

process subset_bam_by_bed {
    input:
        path(input_bam)
        path(input_bed)
        val(sample)
        val(threads)

    output:
        path "*subset.bam", emit: subsetted_bam

    script:
        """
        samtools view -@{threads} \
        -b \
        -h \
        -L ${input_bed} \
        ${input_bam} \
        > ${sample}_subset.bam
        """
}

process get_wanted_read_ids {
    input:
        path(subset_bam)
        val(threads)

    output:
        path "temp_ids.txt", emit: wanted_ids

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${subset_bam} \
        | cut -f 1 \
        > temp_ids.txt
        """

}

process add_suppl_alignments_to_subset {
    input:
        path(subset_bam)
        path(input_bam)
        val(threads)
        val(sample)
        path(list_of_ids)

    output:
        path "*_subset.suppl.bam", emit: suppl_subset_bam

    script:
        """
        samtools \
        view \
        -@${threads} \
        -h \
        -f 2048 \
        -N ${list_of_ids} \
        ${input_bam} \
        | samtools \
        merge \
        -@${threads} \
        -o ${sample}_subset.suppl.bam \
        ${subset_bam} \
        -
        """
}

process index_suppl_subset_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: suppl_subset_indexed_bam

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}   
        """ 
}

process human_variation_snp {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        val(threads)

    publishDir("${params.outdir}", mode: 'copy', pattern: 'wf-human-variation_reports/*wf-human-snp-report.html')

    output:
        path("*/*wf-human-snp-report.html"), emit: human_variation_snp_report
        path("*/*wf_snp.vcf.gz"), emit: human_variation_snp_vars

    script:
        """
        nextflow run epi2me-labs/wf-human-variation -r master \
        -profile standard \
        --ref ${reference} \
        --snp \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir wf-human-variation_reports \
        --threads ${threads}
        """
}

process index_reference {
    input:
        path(reference)

    output:
        tuple path(reference), path("*.fai"), emit: reference_index

    script:
        """
        samtools faidx ${reference}
        """
}

process clairS_to_variants {
    input:
        tuple path(subsetted_bam), path(bai)
        tuple path(reference), path(fai)
        val(params.outdir)
        val(threads)

    output:
        path("*/snv.vcf.gz"), emit: clairS_snv
        path("*/indel.vcf.gz"), emit: clairS_indel

    script:
        """
        /opt/bin/run_clairs_to \
        --tumor_bam ${subsetted_bam} \
        --ref_fn ${reference} \
        --threads ${threads} \
        --platform ont_r10_guppy_hac_5khz \
        --output_dir clairS_to_output
        """
}
 
process filter_to_somatic_variants_only {
    input:
        path(clairS_snv)
        path(clairS_indel)

    output:
        path("snv.somatic.vcf.gz"), emit: clairS_somatic_snv
        path("indel.somatic.vcf.gz"), emit: clairS_somatic_indel

    script:
        """
        bcftools view \
        -f PASS \
        ${clairS_snv} \
        | gzip > snv.somatic.vcf.gz
        bcftools view \
        -f PASS \
        ${clairS_indel} \
        | gzip > indel.somatic.vcf.gz
        """
}

process igv_reports {
    input:
        path(clair3_report)  // filter-report output
        path(somatic_clair3_report)  // filter-report somatic output
        val(sample)
        path(reference)
        path(input_bam)
        path(indexed_bam)
        path(annotations)

    output:
        path("*_igv-report.html"), emit: igv_report
 
    script:
        """
        echo -e "Chr\tStart\tEnd\tFunc\tGene\tExonicFunc\tAAChange\tcytoBand\t1000g_EUR\tCOSMIC\tMUTATION_CLASS" > ${sample}_clair3_report.merged.csv
        if [ -f ${somatic_clair3_report} ]; then
            sed -e 's/,/\t/g' -e 's/\"//g' ${sample}_somatic_clair3_report.csv > ${sample}_somatic_clair3_report.fmt.csv 
            awk 'BEGIN {FS=OFS="\\t"} NR==1 {\$0=\$0 OFS "MUTATION_CLASS"} NR>1 {\$(NF+1)="SOMATIC"}1' ${sample}_somatic_clair3_report.fmt.csv > ${sample}_somatic_clair3_report.new.fmt.csv
            sed 1,1d ${sample}_somatic_clair3_report.new.fmt.csv >> ${sample}_clair3_report.merged.csv
        fi
        if [ -f ${clair3_report} ]; then
            sed -e 's/,/\t/g' -e 's/\"//g' ${sample}_clair3_report.csv > ${sample}_clair3_report.fmt.csv 
            awk 'BEGIN {FS=OFS="\\t"} NR==1 {\$0=\$0 OFS "MUTATION_CLASS"} NR>1 {\$(NF+1)="ALL"}1' ${sample}_clair3_report.fmt.csv > ${sample}_clair3_report.new.fmt.csv
            sed 1,1d ${sample}_clair3_report.new.fmt.csv >> ${sample}_clair3_report.merged.csv
        fi
        if [ -f ${sample}_clair3_report.merged.csv ]; then
            create_report ${sample}_clair3_report.merged.csv \
            --fasta ${reference} \
            --sequence 1 \
            --begin 2 \
            --end 3 \
            --flanking 1000 \
            --info-columns Chr Start End Func Gene ExonicFunc AAChange cytoBand 1000g_EUR COSMIC MUTATION_CLASS \
            --output ${sample}_igv-report.html \
            --standalone \
            --tracks ${input_bam} ${annotations}
        fi
        """
}

process cnvpytor {
    input:
        val(sample)
        path(input_bam)
        val(threads)

    output:
        path("*_cnvpytor_100k.global.0000.png"), emit: cnv_plot
    
    script:
        """
        cnvpytor -root ${sample}_CNV.pytor -rd ${input_bam.toRealPath()} -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -his 1000 10000 100000 -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -partition 1000 10000 100000 -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -call 1000 -j ${threads} > ${sample}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 10000 -j ${threads} > ${sample}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 100000 -j ${threads} > ${sample}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${sample}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${sample}_cnvpytor_100k.png
        """
}

process vcftools {
    maxRetries 3
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }

    input:
        path(human_variation_vars) // human-variation SNP module
        val(sample)

    output:
        path "*.vcf", emit: variants_out

    script:
        """
        vcftools \
        --gzvcf ${human_variation_vars} \
        --remove-filtered-all \
        --recode \
        --out ${sample}
        """
}

process vcftools_somatic {
    maxRetries 3
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }

    input:
        path(somatic_snv) // filtered snvs
        path(somatic_indel) // filtered indels
        val(sample)

    output:
        path "*snv.recode.vcf", emit: snv_somatic_variants_out
        path "*indel.recode.vcf", emit: indel_somatic_variants_out

    script:
        """
        vcftools \
        --gzvcf ${somatic_snv} \
        --remove-filtered-all \
        --recode \
        --out ${sample}_snv
        vcftools \
        --gzvcf ${somatic_indel} \
        --remove-filtered-all \
        --recode \
        --out ${sample}_indel
        """
}

process gzip {
    input:
        path(input_file)

    output:
        path "*.gz", emit: compressed_out

    script:
        """
        pigz \
        -1 \
        -c \
        ${input_file} \
        > ${input_file}.gz
        """
}

process gzip_somatic {
    input:
        path(input_file_snv)
        path(input_file_indel)

    output:
        path "*snv.recode.vcf.gz", emit: snv_compressed_out
        path "*indel.recode.vcf.gz", emit: indel_compressed_out

    script:
        """
        pigz \
        -1 \
        -c \
        ${input_file_snv} \
        > ${input_file_snv}.gz
        cat ${input_file_indel} \
        | grep -v ^# \
        | gzip \
        > ${input_file_indel}.gz
        """
}

process merge_somatic_snvs_and_indels {
    input:
        path(somatic_snvs)
        path(somatic_indels)
        val(sample)
    
    output:
        path "*all.somatic.recode.vcf.gz", emit: all_compressed_out

    script:
        """
        cat ${somatic_snvs} ${somatic_indels} > ${sample}.all.somatic.recode.vcf.gz
        """
}

process bedtools_intersect {
    input:
        path(input1)
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools \
        intersect \
        -a ${input1} \
        -b ${input2} > ${output_file}${ext}
        """
}

process bedtools_intersect_somatic {
    input:
        path(somatic_vars) // merged somatic snvs and indels
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools intersect -a ${somatic_vars} -b ${input2} > ${output_file}${ext} 
        """
}

process convert2annovar {
    input:
        path(input)
        val(output_file)
        val(ext)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        /annovar/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo \
        > ${output_file}${ext}
        """
}

process convert2annovar_somatic {
    input:
        path(input)
        val(output_file)
        val(ext)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        /annovar/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo \
        > ${output_file}${ext}
        """
}

process table_annovar {
    input:
        path(input)
        val(annovar_ver)
        val(output_file)
        val(ext)
        val(threads)
    
    output:
        path "*_multianno.csv", emit: clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} \
        /annovar/humandb/ \
        -buildver ${annovar_ver} \
        -out ${output_file}${ext} \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo \
        -thread ${threads}
        """
}

process table_annovar_somatic {
    input:
        path(input)
        val(annovar_ver)
        val(output_file)
        val(ext)
        val(threads)

    output:
        path "*_somatic_clair3_panel.hg38_multianno.csv", emit: somatic_clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} \
        /annovar/humandb/ \
        -buildver ${annovar_ver} \
        -out ${output_file}_somatic_${ext} \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo \
        -thread ${threads}
        """
}

process bedtools_intersect2 {
    input:
        path(bed_methyl) // filtered bedmethyl file from filter to just 5mC
        path(input2)
        val(output_file)
        val(ext)
        val(sample)

    output:   
        path "*.bed", emit: intersect_bed
        
    script:
        """
        bedtools \
        intersect \
        -a ${bed_methyl} \
        -b ${input2} > ${output_file}.${ext}
        """
}

process mgmt_pred {
    input:
        path(mgmt_pred)
        file(intersect_bed)
        path(probes)
        path(model)
        val(sample)
        val(threads)
        val(mgmt_avg_cov)
        val(minimum_mgmt_coverage)

    output:
        path "*mgmt_status.csv", emit: mgmt_status

    script:
        if ( mgmt_avg_cov.toFloat() > minimum_mgmt_coverage.toFloat() )
            """
            Rscript ${mgmt_pred} \
            --input ${intersect_bed} \
            --probes ${probes} \
            --model ${model} \
            --sample ${sample} \
            --threads ${threads}
            """
        else
            """
            echo "NA" > low_cov_mgmt_status.csv
            """
}

process rapid_cns2_meth_classification {
    input:
        path(meth_class)
        val(sample)
        path(topprobes)
        path(trainingdata)
        path(arrayfile)
        val(threads)
        path(bed_methyl) // output from human-variation mods

    output:
        path "*rf_details.tsv", emit: rf_details
        path "*votes.tsv", emit: votes

    script:
        """
        Rscript ${meth_class} \
        --sample ${sample} \
        --in_file ${bed_methyl} \
        --probes ${topprobes} \
        --training_data ${trainingdata} \
        --array_file ${arrayfile} \
        --threads ${threads}
        """
}

process filter_report {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
        val(params.outdir)
    
    output:
        path "*_clair3_report.csv", emit: clair3_report_csv
    
    script:
        """
        Rscript ${filterreport} \
        --input ${clair3_multianno} \
        --sample ${sample} > ${sample}_clair3_report.csv \
        --output ${sample}_clair3_report.csv 
        """        
}

process filter_report_somatic {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
        val(params.outdir)
    
    output:
        path "*_clair3_report.csv", emit: somatic_clair3_report_csv

    script:
        """
        Rscript ${filterreport} \
        --input ${clair3_multianno} \
        --sample ${sample} > ${sample}_somatic_clair3_report.csv \
        --output ${sample}_somatic_clair3_report.csv 
        """        
}

process find_seq_platform {
    input:
        path(find_seq_platform_py)
        path(input_bam)

    output:
        stdout emit: seq_platform

    script:
        """
        platform=\$(python3 ${find_seq_platform_py} -i ${input_bam})
        echo -n "\$platform"
        """
}

process make_report {
    // cache false forces it to regenerate the report each time and not use the cache
    // necessary for the report as you add/remove different classifiers
    cache false

    input:
        path(makereport)
        path(cnv_plot) // cnv pytor
        path(mgmt_status) // mgmt_pred
        val(rf_details) // rapidcns2 rf_details
        val(votes) // rapidcns2 votes
        path(clair3_report) // filter_report
        path(somatic_clair3_report) // filter report somatic
        val(sample)
        path(report_UKHD)
        path(mosdepth_plot_data) // mosdepth
        val(mgmt_cov)
        val(mgmt_minimum_cov)
        val(sturgeon_class) // sturgeon classification
        path(igv_report) //igv_reports
        val(nextflow_version)
        path(input_bam)
        val(seq)
        val(nanodx_votes) // nanodx classifier output
        path(versions_file) // file of software version info
        path(unique_genes_cov) // unique genes cov
        path(mgmt_plot)
        val(avg_gene_cov)
        path(user_params)
    
    output:
	val(true)

    script:
        """
        Rscript ${makereport} \
        --prefix ${sample} \
        --mutations ${clair3_report} \
        --somatic_mutations ${somatic_clair3_report} \
        --cnv_plot ${cnv_plot} \
        --rf_details ${rf_details} \
        --votes ${votes} \
        --output_dir ${PWD}/${params.outdir} \
        --coverage ${mosdepth_plot_data} \
        --sample ${sample} \
        --report_UKHD ${report_UKHD} \
        --methylartist ${mgmt_plot} \
        --sturgeon_csv ${sturgeon_class} \
        --igv_report ${igv_report} \
        --nextflow_ver ${nextflow_version} \
        --seq ${seq} \
        --nanodx_votes ${nanodx_votes} \
        --software_versions ${versions_file} \
        --mgmt ${mgmt_status} \
        --avg_gene_cov ${avg_gene_cov} \
        --uniq_genes_cov ${unique_genes_cov} \
        --user_params ${user_params} \
        --mgmt_minimum_cov ${mgmt_minimum_cov} \
        --promoter_mgmt_coverage ${mgmt_cov} 
        """
}

process STURGEON_modkit_extract {
    input:
        path(mod_merged_bam)
        val(sample)
        val(threads)

    output:
        path "*_modkit_output.txt", emit: modkit_extract_output

    script:
        """
        /modkit \
        extract \
        ${mod_merged_bam} \
        ${sample}_modkit_output.txt \
        --threads ${threads}
        """
}

process STURGEON_inputtobed {
    input:
        path(modkit_out)

    output:
        path "*/merged_probes_methyl_calls.bed", emit: sturgeon_bed_convert

    script:
        """
        source /sturgeon-0.4.4/venv/bin/activate
        sturgeon \
        inputtobed \
        --reference-genome hg38 \
        -i ${modkit_out} \
        -s modkit \
        -o output 
        deactivate
        """
}

process STURGEON_predict {
    input:
        path(input_bed) // sturgeon inputtobed output

    output:
        path "*/merged_probes_methyl_calls_CAPPER_MODEL.csv", emit: sturgeon_prediction

    script:
        """
        source /sturgeon-0.4.4/venv/bin/activate
        sturgeon \
        predict \
        -i ${input_bed} \
        -o output \
        --model-files /sturgeon-0.4.2/sturgeon/include/models/CAPPER_MODEL.zip \
        --plot-results        
        """
}

process nanoDx_modkit {
    input:
        tuple path(input_bam), path(bam_index)
        path(reference)
        path(mapping)
        val(sample)
        val(threads)

    output:
        path "*CpG.bed", emit: nanoDx_modkit_out

    script:
        """
        /modkit pileup \
        ${input_bam} - \
        --ref ${reference} \
        --include-bed ${mapping} \
        --preset traditional \
        --only-tabs \
        --edge-filter 0,27 \
        --threads ${threads} \
        > ${sample}.CpG.bed
        """
}

process gzip_nanoDx_modkit {
    input:
        path(nanoDx_modkit_out)
        val(sample)

    output:
        path "*CpG.bed.gz", emit: nanoDx_modkit_out_gz

    script:
        """
        gzip -c ${nanoDx_modkit_out} > ${sample}.CpG.bed.gz
        """
}

process nanoDx_annotate_bedMethyl {
    input:
        path(nanodx_pileup_out)
        path(mapping)
        val(sample)
        path(gzip_nanoDx) // gzipped nanoDx modkit out

    output:
        path "*CpG.450K.bed", emit: nanoDx_annotate_out

    script:
        """
        zcat ${gzip_nanoDx} \
        | cut -f1-11 \
        | bedtools intersect \
        -a - \
        -b ${mapping} \
        -wa \
        -wb \
        | awk -v OFS='\t' '\$4=\$15' \
        | cut -f1-11 \
        > ${sample}.CpG.450K.bed
        """
}

process nanoDx_NN_classify {
    input:
        path(classify_NN_bedMethyl_py)
        path(model)
        path(input_bed)
        val(sample)

    output:
        path "*nanoDx_votes.txt", emit: nanoDx_votes

    script:
        """
        python3 ${classify_NN_bedMethyl_py} \
        -m ${model} \
        -i ${input_bed} \\
        -v ${sample}_nanoDx_votes.txt \
        -o ${sample}_nanoDx_output.txt
      """
}

process unique_genes_cov {
    input:
        path(unique_genes_bed)
        path(input_bam)
        val(threads)
        path(bam_index)

    output:
        path "*bed.gz", emit: unique_genes_coverage
        
    script:
        """
        /mosdepth \
        -n \
        --by ${unique_genes_bed} \
        uniq_genes_cov_ \
        ${input_bam} \
        -t ${threads}
        """
}

process format_genes_cov {
    input:
        path(mosdepth_out)

    output:
        path "*csv", emit: format_unique_genes_coverage
 
    script:
        """
        echo "Chr,Start,End,Gene_ID,Coverage" > unique_genes_coverage.csv
        zcat ${mosdepth_out} | awk -F'\t' '{ gsub(/,/, "-", \$4); print \$1 "," \$2 "," \$3 "," \$4 "," \$5 }' >> unique_genes_coverage.csv
        """
}

process filter_gene_cov {
    input:
        path(filter_script) // R script to filter gene cov table
        path(raw_gene_table)

    output:
        path "filtered_gene_cov.csv", emit: filtered_unique_genes_coverage
        path "avg_gene_cov.csv", emit: avg_gene_cov

    script:
        """
        Rscript ${filter_script} -i ${raw_gene_table} -o filtered_gene_cov.csv
        """
}

process write_user_params {
    input:
        val(sample)
        path(input_bam)
        path(reference)
        val(outdir)
        val(threads)
        val(minimum_mgmt_cov)
        path(annotations)
        val(rapidcns2)
        val(sturgeon)
        val(nanodx)
        val(nanoplot)
 
    output:
        path "user_parameters.csv", emit: user_params

    script:
        """
        echo -e "SAMPLE,${sample}" > user_parameters.csv
        echo -e "INPUT_BAM,${input_bam}" >> user_parameters.csv
        echo -e "REFERENCE,${reference}" >> user_parameters.csv
        echo -e "OUTDIR,${outdir}" >> user_parameters.csv
        echo -e "THREADS,${threads}" >> user_parameters.csv
        echo -e "MINIMUM_MGMT_COV,${minimum_mgmt_cov}" >> user_parameters.csv
        echo -e "ANNOTATIONS,${annotations}" >> user_parameters.csv
        echo -e "RAPIDCNS2,${rapidcns2}" >> user_parameters.csv
        echo -e "STURGEON,${sturgeon}" >> user_parameters.csv
        echo -e "NANODX,${nanodx}" >> user_parameters.csv
        echo -e "NANOPLOT,${nanoplot}" >> user_parameters.csv
        """
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {

    // Read all input parameters

    // User defined parameters

    Channel.fromPath(params.bam, checkIfExists: true)
    .set {input_bam}

    Channel.fromPath(params.reference, checkIfExists: true)
    .set {reference}

    Channel.from(params.sample)
    .set {sample}

    Channel.from(params.outdir)
    .set {outdir}

    Channel.from(params.threads)
    .set {threads}

    Channel.from(params.minimum_mgmt_cov)
    .set {minimum_mgmt_cov}

    Channel.fromPath(params.annotations, checkIfExists: true)
    .set {annotations}

    // Collect variables and scripts from bin

    Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed", checkIfExists: true)
    .set {targets}

    Channel.fromPath("${projectDir}/bin/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmt_bed}

    Channel.fromPath("${projectDir}/bin/mgmt_probes.Rdata", checkIfExists: true)
    .set {probes}

    Channel.fromPath("${projectDir}/bin/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {model}

    Channel.fromPath("${projectDir}/bin/mgmt_pred_v0.4.R", checkIfExists: true)
    .set {mgmt_pred}

    Channel.fromPath("${projectDir}/bin/methylation_classification_nanodx_v0.3.R", checkIfExists: true)
    .set {meth_class}

    Channel.fromPath("${projectDir}/bin/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topprobes}
    
    Channel.fromPath("${projectDir}/bin/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingdata}

    Channel.fromPath("${projectDir}/bin/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayfile}

    Channel.fromPath("${projectDir}/bin/filter_report_v0.3.R", checkIfExists: true)
    .set {filterreport}

    Channel.fromPath("${projectDir}/bin/make_report_v0.7.R", checkIfExists: true)
    .set {makereport}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_v0.8.Rmd", checkIfExists: true)
    .set {report_UKHD}

    Channel.fromPath("${projectDir}/bin/find_seq_platform.py", checkIfExists: true)
    .set {find_seq_platform_py}

    Channel.fromPath("${projectDir}/bin/hglft_genome_260e9_91a970_clean.bed", checkIfExists: true)
    .set {nanoDx_mapping}

    Channel.fromPath("${projectDir}/bin/classify_NN_bedMethyl.py", checkIfExists: true)
    .set {classify_NN_bedMethyl_py}

    Channel.fromPath("${projectDir}/bin/Capper_et_al_NN.pkl", checkIfExists: true)
    .set {nanoDx_model}

    Channel.fromPath("${projectDir}/bin/unique_genes.bed", checkIfExists: true)
    .set {unique_genes}

    Channel.fromPath("${projectDir}/bin/filter_gene_cov.R", checkIfExists: true)
    .set {filter_gene_cov}

    //////////////
    // WORKFLOW //
    //////////////

    // get software versions
    get_versions_ch = get_versions()
    get_versions_outside_docker_ch = get_versions_outside_docker()
    merge_versions_files_ch = merge_versions_files(get_versions_ch.versions_file1, get_versions_outside_docker_ch.versions_file2)

    // check input bam file for methylation tags
    check_ch = check_bam_has_meth_data(input_bam, threads)

    // use modkit to merge 'm' and 'h' mods into a single score
    modkit_adjust_ch = modkit_adjust_mods(input_bam, sample, threads)

    // index the input bam
    index_ch = index_input_bam(input_bam, threads)

    // index the merged bam file 
    merged_index_ch = index_merged_bam(modkit_adjust_ch.modkit_merged_bam, threads)

    if ( params.nanoplot != false) {
        // generate NanoPlot QC Report
        nanoplot(input_bam, threads, sample, index_ch.indexed_bam)
    }

    // check the coverage of the mgmt region in the sequence data
    mgmt_coverage_ch = check_mgmt_coverage(input_bam, mgmt_bed, minimum_mgmt_cov, index_ch.indexed_bam, threads)

    // methylartist mgmt plot
    // methylartist only runs if the coverage is detected to be high enough (> minimum_mgmt_cov)    
    methyl_artist_ch = draw_mgmt_methylartist(merged_index_ch.indexed_bam, reference, outdir, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov)

    // mosdepth coverage plots
    mosdepth_ch = mosdepth(threads, targets, input_bam, sample, index_ch.indexed_bam)

    // subset the input bam file by the bed file - to speed up the SNP portion of human-variation
    subsetted_bam_ch = subset_bam_by_bed(modkit_adjust_ch.modkit_merged_bam, targets, sample, threads)

    // index the subsetted bam generated above
    index_subsetted_bam_ch = index_subsetted_bam(subsetted_bam_ch.subsetted_bam, threads)

    // get the read IDs out of the subsetted bam, find all alignmnets associated with IDs and merge them back in
    wanted_ids_ch = get_wanted_read_ids(subsetted_bam_ch.subsetted_bam, threads)

    // using the wanted IDs, pull out the supplementary alignments from the input bam and merge them into the subsetted bam
    subset_suppl_bam_ch = add_suppl_alignments_to_subset(subsetted_bam_ch.subsetted_bam, input_bam, threads, sample, wanted_ids_ch.wanted_ids)

    // index the subsetted bam with the suppl reads added back in
    index_suppl_subset_bam_ch=index_suppl_subset_bam(subset_suppl_bam_ch.suppl_subset_bam, threads)

    // call and run the epi2me-labs/wf-human-variation : mods
    human_variation_mods_ch = human_variation_mods(merged_index_ch.indexed_bam, targets, reference, sample, outdir, 1, threads)
 
    // call and run the epi2me-labs/wf-human-variation : sv
    human_variation_sv(subset_suppl_bam_ch.suppl_subset_bam, targets, reference, sample, outdir, 1, index_suppl_subset_bam_ch.suppl_subset_indexed_bam, threads)

    // call and run the epi2me-labs/wf-human-variation : cnv
    human_variation_cnv(index_ch.indexed_bam_tuple, targets, reference, sample, outdir, 1, threads)

    // run the SNP human variation workflow on the subsetted bam
    human_variation_snp_ch = human_variation_snp(index_subsetted_bam_ch.subsetted_bam_index, targets, reference, sample, outdir, 1, threads)

    // index the supplied reference fasta
    index_reference_ch = index_reference(reference)

    // run clairS_to for SNP calling
    clairS_to_ch = clairS_to_variants(index_subsetted_bam_ch.subsetted_bam_index, index_reference_ch.reference_index, outdir, threads)

    // filter clairS generated vcf
    filter_to_somatic_ch = filter_to_somatic_variants_only(clairS_to_ch.clairS_snv, clairS_to_ch.clairS_indel)

    // run cnvPYTOR
    cnvpytor_ch = cnvpytor(sample, input_bam, threads)

    // run vcftools
    variants_ch = vcftools(human_variation_snp_ch.human_variation_snp_vars, sample)
    variants_somatic_ch = vcftools_somatic(filter_to_somatic_ch.clairS_somatic_snv, filter_to_somatic_ch.clairS_somatic_indel, sample)

    // compress output from vcftools
    compressed_variants_ch = gzip(variants_ch.variants_out)
    compressed_variants_somatic_ch = gzip_somatic(variants_somatic_ch.snv_somatic_variants_out, variants_somatic_ch.indel_somatic_variants_out)

    // merge SNVs and INDELS into a single var file
    merged_somatic_vars_ch = merge_somatic_snvs_and_indels(compressed_variants_somatic_ch.snv_compressed_out, compressed_variants_somatic_ch.indel_compressed_out, sample)

    // run bedtools intersect on output of vcftools
    vcf_intersect_ch = bedtools_intersect(compressed_variants_ch.compressed_out, targets, sample, '_clair3_panel.vcf')
    vcf_intersect_somatic_ch = bedtools_intersect_somatic(merged_somatic_vars_ch.all_compressed_out, targets, sample, '_somatic_clair3_panel.vcf')
	
    // convert vcf to annovar input
    converted_annovar_ch = convert2annovar(vcf_intersect_ch.intersect_vcf, sample, '_clair3_panel.avinput')
    converted_annovar_somatic_ch = convert2annovar_somatic(vcf_intersect_somatic_ch.intersect_vcf, sample, '_somatic_clair3_panel.avinput')

    // run table_annovar on the converted annovar input
    clair3_annovar_ch = table_annovar(converted_annovar_ch.annovar_input, 'hg38', sample, '_clair3_panel', threads)
    clair3_annovar_somatic_ch = table_annovar_somatic(converted_annovar_somatic_ch.annovar_input, 'hg38', sample, 'clair3_panel', threads)

    // run bedtools intersect
    intersect_bed_ch = bedtools_intersect2(human_variation_mods_ch.bedmethyl_gz, mgmt_bed, 'mgmt_5mC.hg38', 'bed', sample)

    // run the mgmt_pred script
    mgmt_pred_ch = mgmt_pred(mgmt_pred, intersect_bed_ch.intersect_bed, probes, model, sample, threads, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov)

    if (params.rapidcns2 != false) {
        // run the meth classification script
        rapid_cns2_ch = rapid_cns2_meth_classification(meth_class, sample, topprobes, trainingdata, arrayfile, threads, human_variation_mods_ch.bedmethyl_gz)
    }
    
    // collect report data and generate report
    filter_report_ch = filter_report(filterreport, clair3_annovar_ch.clair3_output, sample, params.outdir)  
    filter_report_somatic_ch = filter_report_somatic(filterreport, clair3_annovar_somatic_ch.somatic_clair3_output, sample, params.outdir)  

    //  produce igv_report for each SNP in the clair3 report
    igv_reports_ch = igv_reports(filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, reference, input_bam, index_ch.indexed_bam, annotations)

    // find sequence platform for report
    find_seq_platform_ch = find_seq_platform(find_seq_platform_py, input_bam)

    // coverage of list of target genes
    unique_genes_cov_ch = unique_genes_cov(unique_genes, input_bam, threads, index_ch.indexed_bam)

    format_genes_cov_ch = format_genes_cov(unique_genes_cov_ch.unique_genes_coverage)

    filter_gene_cov_ch = filter_gene_cov(filter_gene_cov, format_genes_cov_ch.format_unique_genes_coverage)

    write_user_params_ch = write_user_params(sample, input_bam, reference, outdir, threads, minimum_mgmt_cov, annotations, params.rapidcns2, params.sturgeon, params.nanodx, params.nanoplot)

    ////////////////////////////////////////
    // Decode various classifier options //
    ///////////////////////////////////////

    if (params.rapidcns2 != false ) {
        if ( params.sturgeon != false ) {
            // modkit extract values
            modkit_extract_ch = STURGEON_modkit_extract(modkit_adjust_ch.modkit_merged_bam, sample, threads)

            // sturgeon convert input to bed file
            sturgeon_inputtobed_ch = STURGEON_inputtobed(modkit_extract_ch.modkit_extract_output)

            // sturgeon predict
            sturgeon_predict_ch = STURGEON_predict(sturgeon_inputtobed_ch.sturgeon_bed_convert)

            // nanoDx section
            if ( params.nanodx != false) {
                nanoDx_modkit_ch = nanoDx_modkit(index_ch.indexed_bam_tuple, reference, nanoDx_mapping, sample, threads)

                gzip_nanoDx_modkit_ch = gzip_nanoDx_modkit(nanoDx_modkit_ch.nanoDx_modkit_out, sample)

                nanoDx_annotate_bed_ch = nanoDx_annotate_bedMethyl(nanoDx_modkit_ch.nanoDx_modkit_out, nanoDx_mapping, sample, gzip_nanoDx_modkit_ch.nanoDx_modkit_out_gz)

                nanoDx_classify_ch = nanoDx_NN_classify(classify_NN_bedMethyl_py, nanoDx_model, nanoDx_annotate_bed_ch.nanoDx_annotate_out, sample)

                // generate report inc. rapid_CNS2, inc. sturgeon and inc. nanodx 
                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, rapid_cns2_ch.rf_details, rapid_cns2_ch.votes, filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, sturgeon_predict_ch.sturgeon_prediction, igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, nanoDx_classify_ch.nanoDx_votes, merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)
            } else {
                // generate report inc. rapid_CNS2, inc. sturgeon without nanoDx
                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, rapid_cns2_ch.rf_details, rapid_cns2_ch.votes, filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, sturgeon_predict_ch.sturgeon_prediction, igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, "NULL", merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)
            }
        } else { 
            // nanoDx section
            if ( params.nanodx != false) {
                // inc. rapidCNS2, nanoDx, but no sturgeon

                nanoDx_modkit_ch = nanoDx_modkit(index_ch.indexed_bam_tuple, reference, nanoDx_mapping, sample, threads)

                gzip_nanoDx_modkit_ch = gzip_nanoDx_modkit(nanoDx_modkit_ch.nanoDx_modkit_out, sample)

                nanoDx_annotate_bed_ch = nanoDx_annotate_bedMethyl(nanoDx_modkit_ch.nanoDx_modkit_out, nanoDx_mapping, sample, gzip_nanoDx_modkit_ch.nanoDx_modkit_out_gz)

                nanoDx_classify_ch = nanoDx_NN_classify(classify_NN_bedMethyl_py, nanoDx_model, nanoDx_annotate_bed_ch.nanoDx_annotate_out, sample)

                // generate report inc. nanoDx without sturgeon
                
                  make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, rapid_cns2_ch.rf_details, rapid_cns2_ch.votes, filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, "NULL", igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, nanoDx_classify_ch.nanoDx_votes, merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)
            } else {
                // collect data and generate final report inc. rapidCNS2,  no sturgeon and no nanoDx
                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, rapid_cns2_ch.rf_details, rapid_cns2_ch.votes, filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, "NULL", igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, "NULL", merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)

            // close else nanoDx == false
            }
        }
    } else {
        if ( params.sturgeon != false ) {
            // modkit extract values
            modkit_extract_ch = STURGEON_modkit_extract(modkit_adjust_ch.modkit_merged_bam, sample, threads)

            // sturgeon convert input to bed file
            sturgeon_inputtobed_ch = STURGEON_inputtobed(modkit_extract_ch.modkit_extract_output)

            // sturgeon predict
            sturgeon_predict_ch = STURGEON_predict(sturgeon_inputtobed_ch.sturgeon_bed_convert)

            if ( params.nanodx != false) {
                nanoDx_modkit_ch = nanoDx_modkit(index_ch.indexed_bam_tuple, reference, nanoDx_mapping, sample, threads)

                gzip_nanoDx_modkit_ch = gzip_nanoDx_modkit(nanoDx_modkit_ch.nanoDx_modkit_out, sample)

                nanoDx_annotate_bed_ch = nanoDx_annotate_bedMethyl(nanoDx_modkit_ch.nanoDx_modkit_out, nanoDx_mapping, sample, gzip_nanoDx_modkit_ch.nanoDx_modkit_out_gz)

                nanoDx_classify_ch = nanoDx_NN_classify(classify_NN_bedMethyl_py, nanoDx_model, nanoDx_annotate_bed_ch.nanoDx_annotate_out, sample)

                // generate report no rapidCNS2, inc. nanoDx inc. sturgeon

                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, "NULL", "NULL", filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, sturgeon_predict_ch.sturgeon_prediction, igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, nanoDx_classify_ch.nanoDx_votes, merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)
            } else {
            // generate report no rapidCNS2, inc. sturgeon without nanoDx
            make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, "NULL", "NULL", filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, STURGEON_predict.out, igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, "NULL", merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)
        }
        } else { 
            if ( params.nanodx != false) {
                nanoDx_modkit_ch = nanoDx_modkit(index_ch.indexed_bam_tuple, reference, nanoDx_mapping, sample, threads)

                gzip_nanoDx_modkit_ch = gzip_nanoDx_modkit(nanoDx_modkit_ch.nanoDx_modkit_out, sample)

                nanoDx_annotate_bed_ch = nanoDx_annotate_bedMethyl(nanoDx_modkit_ch.nanoDx_modkit_out, nanoDx_mapping, sample, gzip_nanoDx_modkit_ch.nanoDx_modkit_out_gz)

                nanoDx_classify_ch = nanoDx_NN_classify(classify_NN_bedMethyl_py, nanoDx_model, nanoDx_annotate_bed_ch.nanoDx_annotate_out, sample)

                // generate report no rapidCNS2, inc. nanoDx without sturgeon

                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, "NULL", "NULL", filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, "NULL", igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, nanoDx_classify_ch.nanoDx_votes, merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)

            } else {
            // get here when all classifier are false
                make_report(makereport, cnvpytor_ch.cnv_plot, mgmt_pred_ch.mgmt_status, "NULL", "NULL", filter_report_ch.clair3_report_csv, filter_report_somatic_ch.somatic_clair3_report_csv, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, minimum_mgmt_cov, "NULL", igv_reports_ch.igv_report, nextflow_version, input_bam, find_seq_platform_ch.seq_platform, "NULL", merge_versions_files_ch.versions_file, filter_gene_cov_ch.filtered_unique_genes_coverage, methyl_artist_ch.mgmt_plot, filter_gene_cov_ch.avg_gene_cov, write_user_params_ch.user_params)

            }
        }
    }
}
