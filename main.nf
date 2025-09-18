
nextflow.enable.dsl=2

/*
    Parameters
*/

date = new Date().format('yyyyMMdd')
log.info("Source: ${params.source}.")

println '''
--samplesheet expects the filtered assembly stats table output from assembly-nf by defualt.
To run ONLY BRAKER3, use "--source only-braker" and provide a TSV to --sample_sheet that contains columns: species, strain, asm_path
'''

/*
                                                                    __ 
                          /\                                       / _|
  __ _  ___ _ __   ___   /  \   _ __  _ __   ___             _ __ | |_ 
 / _` |/ _ \ '_ \ / _ \ / /\ \ | '_ \| '_ \ / _ \   ______  | '_ \|  _|
| (_| |  __/ | | |  __// ____ \| | | | | | | (_) | |______| | | | | |  
 \__, |\___|_| |_|\___/_/    \_\_| |_|_| |_|\___/           |_| |_|_|  
  __/ |                                                                
 |___/                                                                 

*/

if (params.debug) {
    println """
    *** Using debug mode ***
    """
    params.output = "annotation-${date}-debug"
    params.pbdata = "${workflow.projectDir}/test_data"
} else {
    if (params.source == 'default') {
        println """
            Running on default source. BRAKER3 and downstream analyses will be executed.
            """
        if (params.sample_sheet == null) {
            println """
            Please specify a sample sheet with parameter --sample_sheet.
            """
            exit 1
        }
        if (params.outdir == null) {
            println """
            Please specify an output directory name with option --outdir.
            """
            exit 1
        }
        params.output = "${params.outdir}-braker_annotation"

    } else if (params.source == 'only-braker') {
            println """
                Running on "only-braker" source. Only BRAKER3 will be executed.
                """
            if (params.sample_sheet == null) {
                println """
                Please specify a sample sheet with parameter --sample_sheet.
                """
                exit 1
            }
            if (params.outdir == null) {
                println """
                Please specify an output directory name with option --outdir.
                """
                exit 1
            }
            params.output = "${params.outdir}-braker_annotation"
        
    } else {
        log.info("Missing source, exiting. Please provide '--source default' or '--source only-braker'.")
            exit 1
    }
}


def log_summary() {
    // Corrected log summary function to print information instead of recursive call
    log.info("Workflow summary: \n" + 
             "Debug mode: ${params.debug}\n" + 
             "Output directory: ${params.output}\n")
    
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

    if (params.source == "only-braker") {
        braker_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        .ifEmpty { exit 1, "Please provide a TSV sample sheet that contains: species, strain, asm_path" }
                        .splitCsv(sep: "\t", header: true)
                        .map { row -> [row.species, row.strain, row.asm_path] } 
        braker3(braker_ch)

    } else if (params.source == "default") {
        // Takes in the entire table of filtered assembly stats!
        braker_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        .ifEmpty { exit 1, "Please provide the filtered assemlby stats sheet output from assembly-nf" }
                        .splitCsv(sep: "\t", header: true)
                        .map { row -> [row.species, row.strain, row.asm_path] } 
        
        braker3(braker_ch)

        busco_p_ch = braker3.out.geneAnno()
        
        agat_ch = braker3.out.geneAnno
                    .map { species, strain, asm_path, gff3 -> tuple(species, strain, gff3) }
        
        longestIso(agat_ch)

        agat_output_ch = (agat_ch.geneAnno)
        
        proteome(agat_output_ch, GENOME)

        buscoProt(busco_p_ch)

        gatherAllStats()
    }
    
}


process braker3 {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
        pattern: "**/*.gff3",
    )
    
    label 'braker'
    container "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/loconn13999-braker3_20250724.sif"
    beforeScript 'module load singularity'

    input:
    tuple val(species), val(strain), path(asm_path)

    output:
    tuple val(species), val(strain), path(asm_path), path("${species}/${strain}/braker/output/${strain}.braker.gff3"), emit: geneAnno

    script:
    def prot_fa = params.prot_path
    """
    mkdir -p ${species}/${strain}/augustus_config 
    mkdir -p ${species}/${strain}/braker/output
    
    # Copy Augustus config (inside the container... nextflow handles this)
    cp -r /opt/Augustus/config/* ${species}/${strain}/augustus_config/
    
    # Adjust permissions
    chmod -R u+w ${species}/${strain}/augustus_config

    # Set environment variable so braker knows where to find the augustus config folder
    export AUGUSTUS_CONFIG_PATH=\$PWD/${species}/${strain}/augustus_config

    # Run BRAKER3
    braker.pl \
        --genome ${asm_path} \
        --species ${species}_${strain} \
        --prot_seq ${prot_fa} \
        --threads ${task.cpus} \
        --busco_lineage=nematoda_odb10 \
        --gff3 \
        --workingdir ${species}/${strain}/braker/output

    # Renaming "braker.gff3" that is produced:
    mv ${species}/${strain}/braker/output/braker.gff3 ${species}/${strain}/braker/output/${strain}.braker.gff3

    """
}

/*
process longestIso {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'agat'

    input:
    tuple val(species), val(strain), path(gff3)

    output:
    tuple val(species), val(strain), path(asm_path), path(".gff3"), emit: geneAnno

    script:
    """
    agat_sp_keep_longest_isoform.pl  -f $file -o ${file%.*}.longest.gff3

    """
}

process proteome {

    input:
    tuple val(species), val(strain), path(asm_path), path(".gff3")


}

process busco_prot {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'busco'

    input:
    tuple val(strain), path(asm_path), val (species), path(.gff3)

    output:
    path(".txt")

    script:
    """
    run busco proteome and append to the filtered stats TXT
    """
}

process gatherAllStats {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'gatherAllStats'

    input:
    tuple path(.txt), (busco_p)

    output:
    path(".txt")

    script:
    """
    # Append busco proteome scores to final, filtered asm_stats table output from assembly-nf




    """
}
*/