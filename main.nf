#!/usr/bin/env nextflow
/*
 * Author       :
 *      Tyler K. Chafin
 *      tkchafin@uark.edu
 *
 *  Written as a part of funded research
 *  by the Arkansas Game and Fish Commission
 *  building predictive models to understand
 *  and forecast the spread of chronis wasting
 *  disease in white-tailed deer
 *  2019
 *
 * Collaborators:
 *      Zach D. Zbinden
 *      Bradley T. Martin
 *      Michael E. Douglas
 *      Marlis R. Douglas
 *
 * Description  : Nextflow pipeline for ResistanceGA
 *
 */

rasters = Channel.fromPath(params.rasters)


process ss_optimization{

    publishDir "${workflow.workDir}/${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    file r from rasters

    output:
    file("*RESULTS_*.SS.rds") into ss_results
    file("*_K.SS.rds") into k_list
    file("*_CD.SS.rds") into cd_list
    file("*_MLPE.SS.rds") into mlpe_list
    file("*_distMat.csv") into distmat_list
    file("*.asc") into ss_asc_files

    script:
    """
    Rscript --vanilla ${params.code_dir}/bin/pipeline_1_SSoptimization.R $workflow.workDir/$params.publish_dir $params.cores_per_job_ss $params.R_source_dir $params.gen_dist $params.coords $params.max_iter $r
    """
}

process ss_gather{

    publishDir "${workflow.workDir}/${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    file ss from ss_results.collect()

    output:
    file("single_surface.RESULTS_ALL.rds") into ss_results_gathered

    script:
    """
    Rscript --vanilla ${params.code_dir}/bin/pipeline_2_SSgather.R $workflow.workDir/$params.publish_dir $params.R_source_dir
    """
}

Channel
     .fromPath( params.rasters )
     .set { raster_set }
raster_set.into  { raster_files_input }

process ms_make_batches{

    publishDir "${workflow.workDir}/${params.publish_dir}", mode: 'copy', overwrite: true

    input:
	file input_files from raster_files_input.toList()

    output:
    file("ms.batch") into ms_batch
	file("gdist-inputs.MS.rds") into gdist_inputs
	file("GA-inputs.MS.rds") into GA_inputs

    script:
    """
    Rscript --vanilla ${params.code_dir}/bin/pipeline_3_MS-get-batches.R $workflow.workDir/$params.publish_dir $params.R_source_dir $params.max_combination $params.max_iter $params.cores_per_job_ms $params.gen_dist $params.coords ${input_files}
    """
}

ms_batch.splitText(by: params.ms_jobs_per_node, file: true).set{chunks}

gdist_inputs.into {gd_opt; gd_gather}
GA_inputs.into {ga_opt; ga_gather}

process ms_optimization{

    publishDir "${workflow.workDir}/${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    file (chunk) from chunks
    file gd from gd_opt
    file ga from ga_opt

    output:
    file("*.MS.rds") into ms_results
	file("*.asc") into ms_asc

    script:
    """
    cat $chunk | parallel -j ${params.ms_jobs_parallel} "Rscript --vanilla ${params.code_dir}/bin/pipeline_4_MSoptimization.R $workflow.workDir/$params.publish_dir $params.R_source_dir $params.max_combination $gd $ga {}"
    """
}

process ms_gather{

    publishDir "${workflow.workDir}/${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    file ms from ms_results.collect()
    file ss from ss_results_gathered
    file gd from gd_gather
    file ga from ga_gather

    output:
    file("multi_surface.RESULTS_ALL.rds") into ms_results_gathered
	file("All_Combinations_Summary.csv") into ms_summary
	file("Bootstrap_Results.csv") into boot_summary

    script:
    """
    Rscript --vanilla ${params.code_dir}/bin/pipeline_5_MSgather.R $workflow.workDir/$params.publish_dir $params.R_source_dir $ss $gd $ga $params.boots
    """

}
