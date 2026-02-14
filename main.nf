
nextflow.enable.dsl=2

/*
    Parameters
*/

date = new Date().format('yyyyMMdd')
log.info("Source: ${params.source}.")

println '''
--samplesheet expects the filtered assembly stats table output from assembly-nf by default.
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
    geneAnno-nf does not have a debug mode yet.....
    """
    params.output = "braker.annotation-${date}-debug"
    params.test_data = "${workflow.projectDir}/test_data"
} else {
    if (params.source == 'default') {
        println """
            Running on default source. BRAKER3 and downstream analyses will be executed.
            """
        if (params.sample_sheet == null) {
            println """
            Please specify a sample sheet with parameter --sample_sheet. The sample sheet should be the filtered asm stats sheet that is output from assemlby-nf.
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
                Please specify a sample sheet with parameter --sample_sheet. The sample sheet should be a TSV with header: species, strain, asm_path
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


// Which proteome BRAKER3 should use
def proteome_map = [
    "c_elegans": "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/libraries/N2.WBonly.WS283.PConly.prot.fa",
    "c_tropicalis": "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/libraries/N2.WBonly.WS283.PConly.prot.fa",
    "c_briggsae": "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/libraries/CBCE_mixed_custom_library.prot.fa",
    "c_nigoni": "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/libraries/Eukaryota.fa" ]

def braker_proteome = proteome_map[params.species]

if (!braker_proteome) {
    error """
    Please provide a valid species to --species parameters. Either 'c_elegans', 'c_tropicalis', 'c_briggsae', or 'c_nigoni'.
    """ }

// def rna_seq_map = ["c_briggsae": ""]

// def braker_rna_seq = rna_seq_map[params.species]

// if (!braker_rna_seq) {
//     error """
//     As of right now, only species "c_briggsae" is compatible with this version of main.nf which incorporates CGC2 SR RNA-seq to help guide gene model prediction
//     """}

// Which fasta to use for masking
def mask_map= [
    "c_elegans": "/vast/eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Ce_repeats_bed_file/final_bed/Ce.clust.class.noProt.fa",
    "c_tropicalis": "/vast/eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Ct_repeats_bed_file/final_bed/Ct.clust.class.noProt.fa",
    "c_briggsae": "/vast/eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Cb_repeats_bed_file/final_bed/Cb.clust.class.noProt.fa"]

def mask_file = mask_map[params.species]


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

        busco_download() 
        busco_db_lineage = busco_download.out

        braker_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        .ifEmpty { exit 1, "Please provide a TSV sample sheet that contains: species, strain, asm_path" }
                        .splitCsv(sep: "\t", header: true)
                        .map { row -> [row.species, row.strain, row.asm_path] } 
        braker3(braker_ch, busco_db_lineage)

    } else if (params.source == "default") {
                                                                                                                                                    //////////////////////////////// NEED TO ADD IFELSE STATMENT FOR NIGONI - WHICH DOES NOT HAVE A MASK FILE.... execute unmasked braker3 version and the "successful_annotations" channel
        busco_download() 
        busco_db_lineage = busco_download.out

        mask_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        .ifEmpty { exit 1, "Please provide the filtered assemlby stats sheet output from assembly-nf" }
                        .splitCsv(sep: "\t", header: true)
                        .map { row -> [row.species, row.strain, row.asm_path] }

        softMask(mask_ch)
	


	// add conditional on if row.rna_seq1 and row.rna_seq2 are present
	//	then perform braker with RNA seq
	//	    process( STAR RNA aln ) 
		//	also need a mode param to specify only RNA, RNA + prot, or only prot..... so maybe not column conditional



        braker_ch = softMask.out.masked.map { species, strain, asm_masked -> tuple(species, strain, asm_masked) }
        
        braker3(braker_ch, busco_db_lineage)

        // successful_annotations = braker3.out.geneAnno
        //     .filter { species, strain, asm_path, gff3_file ->
        //         def success = gff3_file.size() > 100
        //         if (!success) {
        //             log.info "Filtering out failed sample ${species}_${strain} because braker.gff3 could not be made (most likely >100K protein prediction)"
        //         }
        //         return success
        // }
    
        // agat_ch = successful_annotations
        //             .map { species, strain, asm_path, gff3 -> tuple(species, strain, asm_path, gff3) }
        
        agat_ch = braker3.out.geneAnno
                    .map { species, strain, asm_masked, gff3 -> tuple(species, strain, asm_masked, gff3) }
        
        // agat_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true) /// REMOVE - THIS WAS JUST FOR NIC
        //                 .ifEmpty { exit 1, "Please provide a TSV sample sheet that contains: species, strain, asm_path" }
        //                 .splitCsv(sep: "\t", header: true)
        //                 .map { row -> [row.species, row.strain, row.asm_path, row.gff3] } 
        
        longestIso(agat_ch)

        agat_output_ch = longestIso.out.longest 
                            .map { species, strain, asm_masked, gff3 -> tuple(species, strain, asm_masked, gff3) }

        
        proteome(agat_output_ch)

        busco_p_ch = proteome.out.prot
                        .map { species, strain, asm_masked, prot_path -> tuple(species, strain, prot_path) }

        busco_prot(busco_p_ch)

        asm_filt_table_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true).splitCsv(sep: "\t", header: true)
        
        // gff_path_ch = agat_ch.map { species, strain, asm_path, gff3 -> tuple(strain, gff3) } /// REMOVE - THIS WAS JUST FOR NIC

        // gff_path_ch = successful_annotations.map { species, strain, asm_path, gff3 -> 
        //                 def full_path = "${workflow.launchDir}/${params.output}/${species}/${strain}/braker/output/${gff3.name}" // name removes the full path and keeps only the file name - gff3 is stored in a cached directory in /scratch4 so we want to remove this 
        //                 tuple(strain, full_path) } 
        gff_path_ch = braker3.out.geneAnno.map { species, strain, asm_masked, gff3 -> 
                        def full_path = "${workflow.launchDir}/${params.output}/${species}/${strain}/braker/output/${gff3.name}" // name removes the full path and keeps only the file name - gff3 is stored in a cached directory in /scratch4 so we want to remove this 
                        tuple(strain, full_path) } 

        busco_stats_ch = busco_prot.out.buscoStat.map { species, strain, stats_file, _ -> file(stats_file) }.collectFile(name: "all_busco_scores.tsv", keepHeader: true, skip: 1).splitCsv(sep: "\t", header: true).map { row -> tuple(row.strain, row.busco_completeness_protein, row.proteome_path) }

        combined_ch = asm_filt_table_ch.map { row -> tuple(row.strain, row) }  // the second value in the tuple, row, contains the entire row: [strain [entire_row]]
            .join(gff_path_ch)
            .join(busco_stats_ch)
            .map { strain, sample_row, gff_path, busco_prot, prot_path ->
                // Add the new columns to the sample row
                sample_row.gff3_path = gff_path
                sample_row.proteome_path = prot_path
                sample_row.protein_busco = busco_prot
                return sample_row
            }


        final_table = combined_ch
            .collect()  // This ensures all data is gathered before processing
            .map { rows ->
                if (rows.isEmpty()) return ""
                
                // Create header from first row keys
                def header = rows[0].keySet().join('\t')
                
                // Create data rows
                def data = rows.collect { row -> 
                    row.values().join('\t') 
                }.join('\n')
                
                return header + '\n' + data + '\n'
            }
            .collectFile(
                name: "${params.outdir}_all_stats.tsv",
                storeDir: "${workflow.launchDir}/${params.output}")
            .view { "Final table created: $it" }    
    }
    
}



process busco_download {
    
    label 'braker'
    container "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/loconn13999-braker3_20250724.sif"
    beforeScript 'module load singularity'

    output:
    path("mb_downloads")

    script:
    """
    /opt/compleasm_kit/compleasm.py download nematoda_odb10

    """
}

process softMask {
    publishDir(
        path: "${params.output}",
        mode: 'copy',
        pattern: "**/*.fa",
    )

    label 'masking'

    input:
    tuple val(species), val(strain), path (asm_path)

    output:
    tuple val(species), val(strain), path("${species}/${strain}/${asm_path.baseName}_softMasked.fa"), emit: masked

    script:
    """
    mkdir -p ${species}/${strain}

    RepeatMasker -s -xsmall -lib $mask_file -pa ${task.cpus} $asm_path 

    mv ${asm_path.name}.masked ${species}/${strain}/${asm_path.baseName}_softMasked.fa

    """

}


process rna_aln {
	wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"
source activate star

cd $wkdir/${GENOME%%.*}/alignments/
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 600000000000 \
--genomeDir . \
--genomeFastaFiles $wkdir/$GENOME \
--genomeSAindexNbases 12 \
--alignIntronMax 10000
STAR \
--runThreadN 24 \
--genomeDir . \
--outSAMtype BAM Unsorted SortedByCoordinate \
--twopassMode Basic \
--readFilesCommand zcat \
--alignIntronMax 10000 \
--readFilesIn $wkdir/shortreads/${GENOME%%.*}/${GENOME%%.*}.reads.f.fq.gz $wkdir/shortreads/${GENOME%%.*}/${GENOME%%.*}.reads.r.fq.gz


--outSAMstrandField intronMotif # needed for StringTie compatibility with BRAKER3!
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
    tuple val(species), val(strain), path(asm_masked)
    path(busco_db_lineage)

    output:
    tuple val(species), val(strain), path(asm_masked), path("${species}/${strain}/braker/output/${strain}.softMasked.braker.gff3"), emit: geneAnno

    script:
    """    
    mkdir -p ${species}/${strain}/augustus_config 
    mkdir -p ${species}/${strain}/braker/output

    # -L flag to follow symlinks and copy actual content
    cp -rL ${busco_db_lineage} ${species}/${strain}/braker/output/
    #ls -la ${species}/${strain}/braker/output/mb_downloads/

    # Copy Augustus config (inside the container... nextflow handles this)
    cp -r /opt/Augustus/config/* ${species}/${strain}/augustus_config/
    
    # Adjust permissions
    chmod -R u+w ${species}/${strain}/augustus_config

    # Set environment variable so braker knows where to find the augustus config folder
    export AUGUSTUS_CONFIG_PATH=\$PWD/${species}/${strain}/augustus_config

    # Run BRAKER3
    braker.pl \
        --genome ${asm_masked} \
        --species ${species}_${strain} \
        --prot_seq ${braker_proteome} \
        --threads ${task.cpus} \
        --busco_lineage=nematoda_odb10 \
        --gff3 \
        --workingdir ${species}/${strain}/braker/output 

    # Renaming "braker.gff3" that is produced:
    mv ${species}/${strain}/braker/output/braker.gff3 ${species}/${strain}/braker/output/${strain}.softMasked.braker.gff3

    """
}


// process braker3 { // rename to "braker3_unmasked"

//     publishDir(
//         path: "${params.output}",
//         mode: 'copy',
//         pattern: "**/*.gff3",
//     )
    
//     label 'braker'
//     container "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/loconn13999-braker3_20250724.sif"
//     beforeScript 'module load singularity'
//     errorStrategy 'ignore'  // Continue pipeline even if some samples fail

//     input:
//     tuple val(species), val(strain), path(asm_path)
//     path(busco_db_lineage)

//     output:
//     tuple val(species), val(strain), path(asm_path), path("${species}/${strain}/braker/output/${strain}.braker*.gff3"), emit: geneAnno, optional: true

//     script:
//     """    
//     mkdir -p ${species}/${strain}/augustus_config 
//     mkdir -p ${species}/${strain}/braker/output

//     # -L flag to follow symlinks and copy actual content
//     cp -rL ${busco_db_lineage} ${species}/${strain}/braker/output/
//     #ls -la ${species}/${strain}/braker/output/mb_downloads/

//     # Copy Augustus config (inside the container... nextflow handles this)
//     cp -r /opt/Augustus/config/* ${species}/${strain}/augustus_config/
    
//     # Adjust permissions
//     chmod -R u+w ${species}/${strain}/augustus_config

//     # Set environment variable so braker knows where to find the augustus config folder
//     export AUGUSTUS_CONFIG_PATH=\$PWD/${species}/${strain}/augustus_config

//     # Run BRAKER3
//     braker.pl \
//         --genome ${asm_path} \
//         --species ${species}_${strain} \
//         --prot_seq ${braker_proteome} \
//         --threads ${task.cpus} \
//         --busco_lineage=nematoda_odb10 \
//         --gff3 \
//         --workingdir ${species}/${strain}/braker/output && BRAKER_SUCCESS=true || BRAKER_SUCCESS=false

//     # Handle output based on success/failure
//     if [ "\$BRAKER_SUCCESS" = "true" ] && [ -f "${species}/${strain}/braker/output/braker.gff3" ] && [ -s "${species}/${strain}/braker/output/braker.gff3" ]; then
//         mv ${species}/${strain}/braker/output/braker.gff3 ${species}/${strain}/braker/output/${strain}.braker.gff3
//     else
//         echo "BRAKER failed for ${species}_${strain} (likely >100K protein issue)"
//         touch ${species}/${strain}/braker/output/${strain}.braker.EMPTY.gff3
//     fi

//     # Always exit successfully to prevent Nextflow task failure
//     exit 0

//     # Renaming "braker.gff3" that is produced:
//     #mv ${species}/${strain}/braker/output/braker.gff3 ${species}/${strain}/braker/output/${strain}.braker.gff3

//     """
// }


process longestIso {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
        pattern: "**/*.gff3",
    )
    
    label 'agat'

    input:
    tuple val(species), val(strain), path(asm_masked), path(gff3)

    output:
    tuple val(species), val(strain), path(asm_masked), path("${species}/${strain}/${gff3.baseName}.longest.gff3"), emit: longest

    script:
    """
    mkdir -p ${species}/${strain}
    agat_sp_keep_longest_isoform.pl  -f $gff3 -o ${species}/${strain}/${gff3.baseName}.longest.gff3

    """
}


process proteome {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
        pattern: "**/*.protein.fa",
    )
    
    label 'prot'

    input:
    tuple val(species), val(strain), path(asm_masked), path(gff3)

    output:
    tuple val(species), val(strain), path(asm_masked), path("${species}/${strain}/protein/${gff3.baseName}.protein.fa"), emit: prot

    script:
    """
    mkdir -p ${species}/${strain}/protein

    gffread -S $gff3 -g $asm_masked -y ${species}/${strain}/protein/${gff3.baseName}.protein.fa

    """ 
}


process busco_prot {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'busco'

    input:
    tuple val(species), val(strain), path(prot_path)

    output:
    tuple val(species), val(strain), path("${species}/${strain}/busco/${prot_path.baseName}.busco/${prot_path.baseName}.busco.stat.tsv"), path("${species}/${strain}/busco/${prot_path.baseName}.busco/short_summary.specific.nematoda_odb10.${prot_path.baseName}.busco.txt"), emit: buscoStat

    script:
    """
    mkdir -p ${species}/${strain}/busco

    busco -i $prot_path -c 12 -m prot -l /vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/busco_downloads/lineages/nematoda_odb10/ -o ${species}/${strain}/busco/${prot_path.baseName}.busco -c ${task.cpus} --offline

    echo -e "strain\tbusco_completeness_protein\tproteome_path" > header.tsv
    grep "C:" ${species}/${strain}/busco/${prot_path.baseName}.busco/short_summary.specific.nematoda_odb10.${prot_path.baseName}.busco.txt > ${species}/${strain}/busco/${prot_path.baseName}.busco/tmp.tsv
    awk '{ match(\$0, /C:([0-9.]+)%/, a); print a[1] }' ${species}/${strain}/busco/${prot_path.baseName}.busco/tmp.tsv > ${species}/${strain}/busco/${prot_path.baseName}.busco/tmp2.tsv 
    paste -d '\t' <(echo "$strain") ${species}/${strain}/busco/${prot_path.baseName}.busco/tmp2.tsv <(echo "${workflow.launchDir}/${params.output}/${species}/${strain}/protein/${prot_path.baseName}.fa") > strain_busco.tsv
    
    cat header.tsv strain_busco.tsv > ${species}/${strain}/busco/${prot_path.baseName}.busco/${prot_path.baseName}.busco.stat.tsv
    """
}
