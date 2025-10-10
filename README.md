This is a nextflow process that predicts gene models in *Caenorhabditis de novo* assembled genomes using BRAKER3 guided with reference/Eukaryota proteomes. 

```
                                                                    __ 
                          /\                                       / _|
  __ _  ___ _ __   ___   /  \   _ __  _ __   ___             _ __ | |_ 
 / _` |/ _ \ '_ \ / _ \ / /\ \ | '_ \| '_ \ / _ \   ______  | '_ \|  _|
| (_| |  __/ | | |  __// ____ \| | | | | | | (_) | |______| | | | | |  
 \__, |\___|_| |_|\___/_/    \_\_| |_|_| |_|\___/           |_| |_|_|  
  __/ |                                                                
 |___/                                                                 
```
This Nextflow process is meant to be run after [assembly-nf](https://github.com/AndersenLab/assembly-nf), and it will add full paths to files produced and busco protein stats directly to the stats sheet output from assembly-nf.

## --source
The parameter ```--source``` can be equal to "only-braker" or "default". If ```--source``` is set to "only-braker", then only BRAKER3 will run. Otherwise, if ```--source``` is equal to "default", then the entire workflow will run, which is inclusive of BRAKER3, agat for extracting the longest isoform of each gene model prediction, gffread for creating a proteome from the longest isoform GFF and the supplied genome, BUSCO on the created proteome, and then gathering all paths and stats and adding to the supplied ```--sample_sheet```.
   
## --sample_sheet
The formatting of the supplied ```--sample_sheet``` is dependent on what ```--source``` is equal to. If ```--source``` is equal to "only-braker", then the provided ```--sample_sheet``` should be a TSV with headers: species, strain, asm_path. If ```--source``` is equal to "default", then the provided ```--sample_sheet``` can either be the stats TSV output from assembly-nf, or a custom TSV with any number of headers, but it must contain the required columns of: species, strain, asm_path. The absolute path to the ```--sample_sheet``` needs to be provided.
  
## --species
The parameter ```--species``` can be equal to "c_elegans", "c_tropicalis", "c_briggsae", or "c_nigoni". The ```--species``` specifies which proteome/protein database to use for guided gene model predictions. Additionally, the specified species will determine which file is used for soft-masking an assembly prior to running BRAKER3.

## --outdir & -profile
The parameter ```-profile``` specifies the configuration of the pipeline, inclusive of default parameters, resource allocation, and process-specific environments/containers. The only option for ```-profile``` right now is "rockfish". The user can created their own config file in /geneAnno-nf/conf/ if they wish to run on their HPC.


## Use cases
### Running on newly assembled genomes output from assembly-nf:
```nextflow run main.nf --source default --sample_sheet $PWD/nigoni_asm-nf_filteredStats.tsv --outdir nigoni_geneModels --species c_nigoni -profile rockfish ```
