
nextflow.enable.dsl=2

/*
    Parameters
*/

date = new Date().format('yyyyMMdd')
log.info("Source: ${params.source}")

if (params.debug) {
    println """
    *** Using debug mode ***
    """
    params.output = "assembly-${date}-debug"
    params.pbdata = "${workflow.projectDir}/test_data"
} else {
    if (params.source == 'default') {
        println """
            Running on default source.
            """
        if (params.sample_sheet == null) {
            println """
            Please specify a sample sheet with option --sample_sheet.
            """
            exit 1
        }
        if (params.outdir == null) {
            println """
            Please specify an output directory name with option --outdir.
            """
            exit 1
        }
        params.output = "${params.outdir}-assembly"
        // params.pbdata = "${params.data_path}"
    } else if (params.source == 'umd') {
        println """
            Running on UMD source.
            """
        params.output = "${params.raw_dir}-assembly"
        params.pbdata = "${params.data_path}/${params.raw_dir}"
    } else {
	log.info("Missing source, exiting.")
        exit 1
    }
}

def log_summary() {
    // Corrected log summary function to print information instead of recursive call
    log.info("Workflow summary: \n" + 
             "Debug mode: ${params.debug}\n" + 
             "Output directory: ${params.output}\n" + 
             "PB Data path: ${params.pbdata}")
    
    // Show help if requested
    if (params.help) {
        log.info("Help requested, exiting.")
        exit 1
    }
}


/*
    Workflow
*/
workflow {
    
    // add a BRAKER3 process - make a different nextflow script that runs BRAKER3, creates proteomes, runs BUSCO on these proteomes, and creates proteomes that are the longest isoform for OrthoFinder
    braker_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                    .splitCsv(sep: "\t", header: false)
                    .map { row -> [row.species, row.strain, row.filt_asm] } 
    
    braker3(braker_ch)

    busco_p_ch = braker3.out.geneAnno()
            
    buscoProt()
    
}


process braker3 {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'braker'

    input:
    tuple val(species), val(strain), path(filt_asm)

    output:
    tuple val(species), val(strain), path(filt_asm), path(".gff3"), emit: geneAnno

    script:
    """
    
    """
}


process busco_prot {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'busco'

    input:
    tuple val(strain), path(filt_asm), val (species), path(.gff3), path(.txt)

    output:
    path(".txt")

    script:
    """
    run busco proteome and append to the filtered stats TXT
    """
}
