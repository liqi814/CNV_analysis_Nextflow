#!/usr/bin/env nextflow
/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 *
 * This nextflow pipeline is for comprehensive CNV analysis.
 * Author: Qi Li (ql2387)
 */

nextflow.enable.dsl=2

params.sampleList = '/nfs/external/az-ipf-garcia/CNVanalysis_shortTLwoQV/WGS_BAM_File/*.bam'
params.reference = "/nfs/seqscratch09/AZ-IPF/reference/hs37d5.fa" 
params.Interval300bp = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_300bp/Interval_list/hs37d5.preprocessed_300bp.interval_list"
params.Interval_1kb = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_1kb/hs37d5.preprocessed.1000bp.primary_contigs.interval_list"
params.outdir = "/nfs/external/az-ipf-garcia/CNVanalysis_shortTLwoQV"
params.WindowSizes = [300]

WindowSizes = Channel.fromList(params.WindowSizes)

log.info """\
 C N V - N F   P I P E L I N E
 ===================================
 samplelist     : ${sampleBAM_ch}
 reference      : ${params.reference}
 outdir         : ${params.outdir}
 """

Channel.fromPath(params.sampleList)
	.map { it -> [it.getSimpleName, it] }
	.set{sampleBAM_ch}

process CaseCollectRCs {
	publishDir "${sample_file.baseName}/case_RCs", mode: 'copy'
	tag "${sample_file.baseName}"
	executor 'sge'
	cpus 2
	penv 'threaded'
	clusterOptions = '-S /bin/bash' 

	input:
	file(sample_file)
	each mode
	val params.Interval300bp
	val params.Interval_1kb
	val params.reference

	output:
	file ('*clean_counts.tsv')

	script:
	if( mode == 300)
	"""
	GATK="/nfs/goldstein/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
	java -jar \${GATK} CollectReadCounts -I $sample_file  -L $params.Interval300bp -R $params.reference --format TSV -imr OVERLAPPING_ONLY -O ${sample_file.baseName}.300bp_clean_counts.tsv
	"""
	else
	"""
	GATK="/nfs/goldstein/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
    java -jar \${GATK} CollectReadCounts -I $sample_file  -L $params.Interval_1kb -R $params.reference --format TSV -imr OVERLAPPING_ONLY -O ${sample_file.baseName}.1kb_clean_counts.tsv
	"""

}

workflow{
	CaseCollectRCs(sampleBAM_ch, WindowSizes, params.Interval300bp, params.Interval_1kb, params.reference)
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
