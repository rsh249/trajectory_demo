#!/usr/bin/env nextflow

params.infile = "$baseDir/data/sample.fa"
params.matrix = "$baseDir/data/GSE72857_umitab.txt"


//  set variables
matrix_file = file(params.matrix)
 
/*
 * Run data download and setup script to convert gene count matrix to Seurat object
 */
process data_setup {
    input:
    		file(infile)
    output:
        file('*.rds')
		
    """
    Rscript --vanilla '${workflow.projectDir}/R/00_setup.R' '${params.matrix}'
    """
}


process run_slingshot {
    input:
        file('infile')

    output:
        file('*.rds')
		
    """
    Rscript --vanilla '${workflow.projectDir}/R/01_slingshot.R'
    """
}

 
/*
 * Define the workflow
 */
workflow {
    data_setup(params.matrix) | run_slingshot | view
}