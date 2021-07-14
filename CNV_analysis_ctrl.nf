#!/usr/bin/env nextflow
/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 *
 * This nextflow pipeline is for comprehensive CNV analysis.
 * Author: Qi Li (ql2387)
 */

Channel
	.fromPath(params.sampleListF)
	.splitCsv(sep:'')
	.set{ bamFileLocF }

Channel
	.fromPath(params.sampleListM)
	.splitCsv(sep:'')
	.set{ bamFileLocM }

process PreprocessIntervals {
	publishDir "Interval_List", mode: 'copy'

	input:
	each mode from Channel.fromList(params.WindowSizes)

	output:
	set val(mode), file('*bp.interval_list') into Preprocess_ch

	script:
	if(mode == 300)
	"""
	java -jar ${params.GATK}  PreprocessIntervals -R $params.reference \
	-XL  $params.BlackList --sequence-dictionary $params.DICT \
	--bin-length ${mode}  -imr OVERLAPPING_ONLY \
	-O hs37d5.preprocessed.${mode}bp.interval_list
	"""
	else
	"""
	java -jar ${params.GATK}  PreprocessIntervals -R $params.reference \
	-XL  $params.BlackList --sequence-dictionary $params.DICT \
	--bin-length ${mode}  -imr OVERLAPPING_ONLY \
	-O hs37d5.preprocessed.${mode}bp.interval_list
	"""
}

process AnnotateIntervals {
	publishDir "Interval_List", mode: 'copy'

	input:
	set val(mode), file(PreprocessedInterval) from Preprocess_ch

	output:
	set val(mode), file('*bp.interval_list') into Annotate_ch

	script:
	"""
	java -jar ${params.GATK}  AnnotateIntervals -R $params.reference \
	-L $PreprocessedInterval --interval-merging-rule OVERLAPPING_ONLY \
	-O hs37d5.annotated.${mode}bp.interval_list
	"""
}

process CtrlCollectRCsF {
	publishDir "ctrl_RCs/Female", mode: 'copy'
	tag "${sample_fileF.baseName}"

	input:
	path(sample_fileF) from bamFileLocF
	each mode from Channel.fromList(params.WindowSizes)

	output:
	file('*_F_clean_counts.hdf5')

	script:
	"""
	java -jar ${params.GATK} CollectReadCounts -I $sample_fileF \
	-L ${params.OUTDIR}/Interval_List/hs37d5.preprocessed.${mode}bp.interval_list \
	-R $params.reference --format HDF5 -imr OVERLAPPING_ONLY \
	-O ${sample_fileF.baseName}.${mode}_F_clean_counts.hdf5
	"""
}

process CtrlCollectRCsM {
	publishDir "ctrl_RCs/Male", mode: 'copy'
	tag "${sample_fileM.baseName}"

	input:
	path(sample_fileM) from bamFileLocM
	each mode from Channel.fromList(params.WindowSizes)

	output:
	file('*_M_clean_counts.hdf5')

	script:
	"""
	java -jar ${params.GATK} CollectReadCounts -I $sample_fileM \
	-L ${params.OUTDIR}/Interval_List/hs37d5.preprocessed.${mode}bp.interval_list \
	-R $params.reference --format HDF5 -imr OVERLAPPING_ONLY \
	-O ${sample_fileM.baseName}.${mode}_M_clean_counts.hdf5
	"""
}
