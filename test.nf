#!/usr/bin/env nextflow

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
    Rscript --vanilla ${params.code_dir}/bin/pipeline_5_MSgather.R $workflow.workDir/$params.publish_dir $params.R_source_dir $gd $ga
    """

}
