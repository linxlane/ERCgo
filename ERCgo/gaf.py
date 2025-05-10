import pandas
from goatools.anno.gaf_reader import GafReader

def goatoolsReadGaf(gafPath):
  ogaf = GafReader(gafPath)
  ns2assc = ogaf.get_ns2assc()
  geneAndGoList = []

  for namespace, associations in ns2assc.items():
    for protein_id, go_ids in sorted(associations.items()):
      geneGoPair = []
      geneGoPair.append(protein_id)
      geneGoPair.append(sorted(go_ids))
      geneAndGoList.append(geneGoPair)

  geneAndGoListDf = pandas.DataFrame(geneAndGoList, columns=['Gene_ID', 'GO_Terms'])
  return geneAndGoListDf


def generateAssocDict(df):
  assoc_dict = df.set_index('Gene_ID')['GO_Terms'].to_dict()
  return assoc_dict


def processGaf(gafPath, masterOut):
  #Get GO terms for each gene in gene association file, returns dataframe
  print('> Utilize GOATOOLS package to extract GO terms for each gene from GAF', flush=True)
  goAssocDF = goatoolsReadGaf(gafPath)
  print(' > DONE')


  #Save go terms for each gene in tsv
  print('> Write [COMP_ID, GO_Terms] table to tsv: ID_GO_TERMS_TABLE.tsv', flush=True)
  compGoPath = masterOut + '/[COMP_ID, GO_Terms]_TABLE.tsv'
  goAssocDF.to_csv(compGoPath, sep='\t', index=False)
  print(' > DONE')


  #Turn dataframe into dictionary and return for later gene-go look up
  print('> Generate dictionary of {Gene_ID:GO_Terms} to find GO terms associated with genes in ERCnet identified gene pairs', flush=True)
  assoc_dict = generateAssocDict(goAssocDF)
  print(' > DONE')

  return assoc_dict
