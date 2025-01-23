ERCgo is an experimental, proof of concept python program to analyze the gene ontology of ERCnet output.

python3 ERCgo -i my_erc_net_ouput -g tair.gaf -o analysis_OUT

| Short flag | Long flag         | Description | Required? |
|------------|-------------------|-------------|-----------|
| -i         | --input           | Path to ERCnet output files. This must contain one edge and one vertices file to analyze. | Yes |
| -g         | --gaf             | Path to gaf file. https://geneontology.org/docs/download-go-annotations/| Yes |
| -o         | --output          | Path to output folder where files generated during analysis will be written. This directory must be created prior to running analysis. | Yes |
