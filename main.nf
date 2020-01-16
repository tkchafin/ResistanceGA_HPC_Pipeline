#!/usr/bin/env nextflow
/*
 * Author       :
 *      Tyler K. Chafin
 *      tkchafin@uark.edu
 #
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

process prep_dirs{

    input:
    val x from 1..params.replicates

    output:
    file("${params.all_comb}/rep_*") into rep_dirs

    script:
    """
    mkdir -p $params.all_comb"/rep_"$x
    """
}

process ss_optimization{

    input:
    file rep from rep_dirs

    output:
    file("${rep}/single.surface") into ss_dirs

    script:
    """
    echo $rep
    mkdir $rep"/single.surface"
    """

}
