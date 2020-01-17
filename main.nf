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

    input:
    file r from rasters

    output:
    file("*RESULTS_*.SS.rds") into ss_results
    file("*_K.SS.rds") into k_list
    file("*_CD.SS.rds") into cd_list
    file("*_MLPE.SS.rds") into mlpe_results

    script:
    """
    Rscript --vanilla ${params.code_dir}/bin/pipeline_1_SSoptimization.R $workflow.workDir 4 $params.R_source_dir $params.gen_dist $params.coords $params.max_iter $r
    """
}

process ss_gather{

    input:
    file ss from ss_results.collect()

    output:
    file("single_surface.RESULTS_ALL.rds") into ss_results_gathered

    script:
    """

    """

}
