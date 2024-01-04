#!/usr/bin/env nextflow

params.infile = "$baseDir/data/sample.fa"
params.matrix = "$baseDir/data/GSE72857_umitab.txt"
params.outSeurat = "tmp/tmp_seurat.rds" // within nextflow output folder


//  set variables
matrix_file = file(params.matrix)
 
/*
 * Run data download and setup script to convert gene count matrix to Seurat object
 */
process data_setup {
		tag { setup }

    input:
	    path "input"

    output:
    	path("tmp_seurat.rds")

    """
    Rscript --vanilla '$baseDir/R/00_setup.R' $input
    """
}


process run_slingshot {
		input:
				path("input")
    output:
        path("test")
		
    """
    Rscript --vanilla '$baseDir/R/01_slingshot.R' tmp_seurat.rds
    """
}

 
/*
 * Define the workflow
 */
workflow {
    data_setup(params.matrix) |
    run_slingshot | view
}