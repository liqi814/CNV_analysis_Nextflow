#!/usr/bin/env nextflow
/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 *
 * This nextflow pipeline is for comprehensive CNV analysis.
 * Author: Qi Li (ql2387)
 */

params.sampleList = '/nfs/external/az-ipf-garcia/CNVanalysis_shortTLwoQV/WGS_BAM_File/*.bam'
params.reference = "/nfs/seqscratch09/AZ-IPF/reference/hs37d5.fa" 
params.Interval300bp = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_300bp/Interval_list/hs37d5.preprocessed_300bp.interval_list"
params.Interval_1kb = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_1kb/hs37d5.preprocessed.1000bp.primary_contigs.interval_list"
params.outdir = "/nfs/external/az-ipf-garcia/CNVanalysis_shortTLwoQV"
params.GATK="/nfs/goldstein/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
params.WindowSizes = [1000]
params.FemalePoN_1kb = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_1kb/PoN/cnv_FemaleHFM_1kb.pon.hdf5"
params.FemalePoN_300bp = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_300bp/PoN/cnv_FemaleHFM.pon.hdf5"


log.info """\
 C N V - N F   P I P E L I N E
 ===================================
 samplelist     : ${Channel.fromPath(params.sampleList)}
 reference      : ${params.reference}
 outdir         : ${params.outdir}
 GATK			: ${params.GATK}
 """

process CaseCollectRCs {
	publishDir "${sample_file.baseName}/case_RCs", mode: 'copy'
	tag "${sample_file.baseName}"

	input:
	file(sample_file) from Channel.fromPath(params.sampleList)
	each mode from Channel.fromList(params.WindowSizes)

	output:
	file ('*clean_counts.tsv') into case_RCs_ch

	script:
	if(mode == 300)
	"""
	java -jar ${params.GATK} CollectReadCounts -I $sample_file \
	-L $params.Interval300bp -R $params.reference --format TSV -imr OVERLAPPING_ONLY \
	-O ${sample_file.baseName}.${mode}_clean_counts.tsv
	"""
	else
	"""
	java -jar ${params.GATK} CollectReadCounts -I $sample_file \
	-L $params.Interval_1kb -R $params.reference --format TSV -imr OVERLAPPING_ONLY \
	-O ${sample_file.baseName}.${mode}_clean_counts.tsv
	"""
}

process DenoiseReadCounts {
	publishDir "${sample_file.baseName}/Denoise", mode: 'copy'
	tag "${sample_file.baseName}"

	input:
	file(sample_file) from Channel.fromPath(params.sampleList)
	file(readCount_file) from case_RCs_ch
	each mode from Channel.fromList(params.WindowSizes)

	output:
	file ('*_clean.denoisedCR.tsv') into denoised_RCs_ch
	file ('*_clean.standardizedCR.tsv')

	script:
	if( mode == 300)
	"""
	java -jar ${params.GATK} DenoiseReadCounts -I $readCount_file \
	--count-panel-of-normals ${params.FemalePoN_300bp} \
	--standardized-copy-ratios ${sample_file.baseName}.${mode}_clean.standardizedCR.tsv \
	--denoised-copy-ratios ${sample_file.baseName}.${mode}_clean.denoisedCR.tsv
	"""
	else
	"""
	java -jar ${params.GATK} DenoiseReadCounts -I $readCount_file \
	--count-panel-of-normals ${params.FemalePoN_1kb} \
	--standardized-copy-ratios ${sample_file.baseName}.${mode}_clean.standardizedCR.tsv \
	--denoised-copy-ratios ${sample_file.baseName}.${mode}_clean.denoisedCR.tsv
	"""
}

process CollectAllelicCounts {
	publishDir "${sample_file.baseName}/AllelicCounts", mode: 'copy'
	tag "${sample_file.baseName}"

	input:
	file(sample_file) from Channel.fromPath(params.sampleList)
	each mode from Channel.fromList(params.WindowSizes)

	output:
	file ('*_allelicCounts.tsv') into alleleCount_ch

	script:
	if( mode == 300)
	"""
	java -jar ${params.GATK} CollectAllelicCounts -I $sample_file \
	-L $params.Interval300bp -R $params.reference \
	-O ${sample_file.baseName}.${mode}_allelicCounts.tsv
	"""
	else
	"""
	java -jar ${params.GATK} CollectAllelicCounts -I $sample_file \
	-L $params.Interval_1kb -R $params.reference \
	-O ${sample_file.baseName}.${mode}_allelicCounts.tsv
	"""
}

process ModelSegments {
	tag "${sample_file.baseName}"

	input:
	file(sample_file) from Channel.fromPath(params.sampleList)
	file(denoised_RCs) from denoised_RCs_ch
	file(alleleCounts) from alleleCount_ch

	output:
	file ('*.cr.seg') into segments_ch

	script:
	"""
	java -jar ${params.GATK} ModelSegments \
	--denoised-copy-ratios $denoised_RCs --allelic-counts $alleleCounts \
	--output ${params.outdir}/${sample_file.baseName}/ModelSegments \
	--output-prefix ${sample_file.baseName}
	"""
}

process CallCopyRatioSegments {
	publishDir "${sample_file.baseName}/CopyRatioSegments", mode: 'copy'
	tag "${sample_file.baseName}"

	input:
	file(segment_file) from segments_ch
	file(sample_file) from Channel.fromPath(params.sampleList)

	output:
	file('*.called.seg')

	scripts:
	"""
	java -jar ${params.GATK} CallCopyRatioSegments \
	--input $segment_file
	--output ${sample_file.baseName}.called.seg
	"""

}
workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}

