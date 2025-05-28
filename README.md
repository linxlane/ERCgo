ERCgo is an experimental, proof of concept python program which analyzes the gene ontology of ERCnet/interactome networks.
```
python3 ERCgo -i network_file -g tair.gaf
```
or
```
python3 ERCgo -d my_input_files
```

Please note that either the -i and -g flags or the -d flag must be included in the command for ERCgo to properly run. Do not use all three.

| Short flag | Long flag         | Description | Required? |
|------------|-------------------|-------------|-----------|
| -j         | --job_name        | Job name for this run of ERCgo. If a directory with this job name already exists at the output path, it will be erased and rewritten. Avoid including spaces or special characters ("_" is ok) | Yes |
| -i         | --interactome     | Path to interactome that will be used for the GO analysis | No |
| -g         | --gaf             | Path to gene association file that will be used for the GO analysis (https://geneontology.org/docs/download-go-annotations/)| No |
| -d         | --directory_input | Path to directory which contains the interactome and gene association files that will be used in the GO analysis. GAF file must have .gaf extension. Interactome file must have keyword "interactome" in file name. | No |
| -o         | --output          | Path where new directory for ERCgo output will be created with the job name. If not included, it will be written in the ERCgo OUTPUT directory. If this path already exists, it will be deleted and a new directory will be created at the output_directory/job_name path. | No |
