## To start on the HPC, run:
`nextflow run funghub/PFAS_project --input PRJNA604830 -profile spartan_hpc -latest -resume`
- -profile spartan_hpc sets the profile to run on specific nodes on the HPC
- **code_acc_bowtie.nf** (a pipeline) runs through Nextflow by pulling from GitHub that had **nextflow.config** configured already.
**NOTE: please download the TempO-Seq Manifest (Human S1500+ Surrogate v2.0 or another version) into a folder "manifest" in the directory you are running nextflow.**

## About command to start
### Usage to run command:

`nextflow run funghub/PFAS_project --input [SRA accession number] --output_dir [folder for results] -profile spartan_hpc -latest -resume`

### Options:
- `--input`   Input SRA accession number (no default value)
- `--output_dir`   Directory to save results (default: results)

### Input option
The input for the command is the SRA accession number in the SRA database on NCBI that you want to retrieve (ex: PRJNA1137368).

## Simplicity of pipeline
- enter only SRA accession number you find on NCBI database into --input option
- pipeline will save important files in your current directory

## Important Step (depends on single read or paired end read): To Verify Before Trimming The FASTQ Files
**this pipeline assumes single read in the FASTP process, you may change the process block in the pipeline to fit your data*

In the process if you would like to verify before trimming the fastq files, you can comment out processes in the workflow block in main: and publish: and output block after MULTIQC_pretrim.
- you can proceed to viewing the Multiqc files of pretrim runs
- you may edit the fastp process after if you have any further specifications for options
- proceed by uncommenting and resume pipeline run