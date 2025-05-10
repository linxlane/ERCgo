import pandas
from goatools.anno.gaf_reader import GafReader
import itertools

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


def generateAssocDict(filePath):
  goTermsDict = {}
  with open(filePath, 'r') as goTerms:
    next(goTerms)
    for line in goTerms:
      lineData = line.strip().split('\t')
      key = lineData[0]
      value = eval(lineData[1])

      if key in goTermsDict:
        currentVals = goTermsDict[key]
        #print('currentVals:' + str(currentVals))
        for item in value:
          #print('Item:' + str(item))
          currentVals.append(str(item))
          #print('Update:' + str(currentVals))
      else:
        goTermsDict[key] = value
  
  return goTermsDict


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
  assoc_dict = generateAssocDict(compGoPath)
  print(' > DONE')

  return assoc_dict
