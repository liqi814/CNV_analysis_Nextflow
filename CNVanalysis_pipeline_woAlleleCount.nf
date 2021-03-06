#!/usr/bin/env nextflow
/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 *
 * This nextflow pipeline is for comprehensive CNV analysis.
 * Author: Qi Li (ql2387)
 */

params.gender = "F"
if(params.gender == "F") {
        params.PoN_1kb = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_1kb_old/PoN/cnv_FemaleHFM_1kb.pon.hdf5"
        params.PoN_300bp = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_300bp/PoN/cnv_FemaleHFM.pon.hdf5"
}
else {
        params.PoN_1kb = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_1kb_old/PoN/cnv_MaleHFM_1kb.pon.hdf5"
        params.PoN_300bp = "/nfs/projects/CNV_WGS/GATK_gCNV/RCs/schizo_HFM_300bp/PoN/cnv_MaleHFM.pon.hdf5"
}

Channel
	.fromPath(params.sampleList)
	.splitCsv(sep:'')
	.into{bamFileLoc1; bamFileLoc2}


log.info """\
 C N V - N F   P I P E L I N E
 ===================================
 samplelist     : ${Channel.fromPath(params.sampleList)}
 reference      : ${params.reference}
 GATK           : ${params.GATK}
 """

process CaseCollectRCs {
	publishDir "${sample_file.baseName}/case_RCs", mode: 'copy'
	tag "${sample_file.baseName}"

	input:
	path(sample_file) from bamFileLoc1
	each mode from Channel.fromList(params.WindowSizes)

	output:
	set val(sample_file.baseName), val(mode), file('*clean_counts.tsv') into case_RCs_ch

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
	publishDir "${sample_ID}/Denoise", mode: 'copy'
	tag "${sample_ID}"

	input:
	set val(sample_ID), val(mode), file(readCount_file) from case_RCs_ch

	output:
	set val(sample_ID), val(mode), file('*_clean.denoisedCR.tsv') into denoised_RCs_ch
	file ('*_clean.standardizedCR.tsv')

	script:
	if( mode == 300)
	"""
	java -jar ${params.GATK} DenoiseReadCounts -I $readCount_file \
	--annotated-intervals ${params.Annotated_300bp} \
	--count-panel-of-normals ${params.PoN_300bp} \
	--standardized-copy-ratios ${sample_ID}.${mode}_clean.standardizedCR.tsv \
	--denoised-copy-ratios ${sample_ID}.${mode}_clean.denoisedCR.tsv
	"""
	else
	"""
	java -jar ${params.GATK} DenoiseReadCounts -I $readCount_file \
	--annotated-intervals ${params.Annotated_1kb} \
	--count-panel-of-normals ${params.PoN_1kb} \
	--standardized-copy-ratios ${sample_ID}.${mode}_clean.standardizedCR.tsv \
	--denoised-copy-ratios ${sample_ID}.${mode}_clean.denoisedCR.tsv
	"""
}

process ModelSegments {
	tag "${sample_ID}"
	publishDir "$sample_ID/ModelSegments", mode: 'copy'

	input:
	set val(sample_ID), val(mode), val(count_files) from denoised_RCs_ch

	output:
	set val(sample_ID), val(mode), file ('*.cr.seg') into segments_ch
	file "*"

	script:
	"""
	java -jar ${params.GATK} ModelSegments \
	--denoised-copy-ratios ${count_files} \
	--output . \
	--output-prefix ${sample_ID}.${mode}
	"""
}

process CallCopyRatioSegments {
	publishDir "$sample_ID/CopyRatioSegments", mode: 'copy'
	tag "$sample_ID"

	input:
	set val(sample_ID), val(mode), file(segment_file) from segments_ch

	output:
	set val(sample_ID), file('*.called.seg') into calledSegments_ch
	file "*"

	script:
	"""
	java -jar ${params.GATK} CallCopyRatioSegments \
	--input $segment_file \
	--output ${sample_ID}.${mode}.called.seg
	"""
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}

