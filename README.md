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
The parameters ```--source``` can be equal to "only-braker", or "default". If ```--source``` is set to "only-braker", then only BRAKER3 will run. Otherwise, if ```--source``` is equal to "default", then the entire workflow will run, which is inclusive of BRAKER3, agat for extracting the longest isoform of each gene model prediction, gffread for creating a proteome from the longest isoform GFF and the supplied genome, BUSCO on the created proteome, and then gathering all paths and stats and adding to the supplied ```--sample_sheet```.
   
## --sample_sheet

## --species

## --outdir & -profile
