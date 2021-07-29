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
        set val(mode), file('*bp.interval_list') into(Preprocess_ch1, Preprocess_ch2, Preprocess_ch3)

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
        set val(mode), file(PreprocessedInterval) from Preprocess_ch1

        output:
        set val(mode), file('*bp.interval_list') into(Annotate_ch1, Annotate_ch2)

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
        each(sample_fileF) from bamFileLocF
        set val(mode), file(PreprocessedInterval) from Preprocess_ch2

        output:
        file('*_F_clean_counts.hdf5') into RCsFileF_ch

        script:
        """
        java -jar ${params.GATK} CollectReadCounts -I $sample_fileF \
        -L $PreprocessedInterval \
        -R $params.reference --format HDF5 -imr OVERLAPPING_ONLY \
        -O ${sample_fileF.baseName}.${mode}_F_clean_counts.hdf5
        """
}

process CtrlCollectRCsM {
        publishDir "ctrl_RCs/Male", mode: 'copy'
        tag "${sample_fileM.baseName}"

        input:
        each(sample_fileM) from bamFileLocM
        set val(mode), file(PreprocessedInterval) from Preprocess_ch3

        output:
        file('*_M_clean_counts.hdf5') into RCsFileM_ch

        script:
        """
        java -jar ${params.GATK} CollectReadCounts -I $sample_fileM \
        -L $PreprocessedInterval \
        -R $params.reference --format HDF5 -imr OVERLAPPING_ONLY \
        -O ${sample_fileM.baseName}.${mode}_M_clean_counts.hdf5
        """
}

/*
Channel
        .fromPath( 'ctrl_RCs/Male/*_M_clean_counts.hdf5' )
        .collect()
        .set{ RCsFileM_ch }

Channel
        .fromPath( 'ctrl_RCs/Female/*_F_clean_counts.hdf5' )
        .collect()
        .set{ RCsFileF_ch }
*/

process PoN_M {
        publishDir "PoNs", mode: 'copy'

        input:
        val(MaleRCsList) from RCsFileM_ch.collect()
        set val(mode), file(AnnotatedInterval) from Annotate_ch1

        output:
        file('*')

        script:
        """
        echo $MaleRCsList | tr -d '[] '  | tr ',' '\n' >> Male_ReadCount.list
        java -jar ${params.GATK} CreateReadCountPanelOfNormals \
        --minimum-interval-median-percentile 5.0 \
        --output cnv_MaleHFM_${mode}bp.pon.hdf5 \
        --annotated-intervals $AnnotatedInterval \
        --input Male_ReadCount.list
        """
}

process PoN_F {
        publishDir "PoNs", mode: 'copy'

        input:
        val(FemaleRCsList) from RCsFileF_ch.collect()
        set val(mode), file(AnnotatedInterval) from Annotate_ch2

        output:
        file('*')

        script:
        """
        echo $FemaleRCsList | tr -d '[] '  | tr ',' '\n' >> Female_ReadCount.list
        java -jar ${params.GATK} CreateReadCountPanelOfNormals \
        --minimum-interval-median-percentile 5.0 \
        --output cnv_FemaleHFM_${mode}bp.pon.hdf5 \
        --annotated-intervals $AnnotatedInterval \
        --input Female_ReadCount.list
        """
}

