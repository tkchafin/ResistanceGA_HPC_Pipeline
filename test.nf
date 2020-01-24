#!/usr/bin/env nextflow

Channel
	.fromPath("/Users/tkchafin/ResistanceGA_HPC_Pipeline/ms.batch")
	.splitText(by: 1, file: true)
	.set{test}

test.view()
