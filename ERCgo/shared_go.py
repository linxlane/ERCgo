import pandas

def generateAssocDict(df):
  assoc_dict = df.set_index('Gene_ID')['GO_Terms'].to_dict()
  return assoc_dict

def genePairGO(genePairsPath, goTermsDict, outputPath):
  geneGOList = []

  with open(genePairsPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      matchingGoList = []
      genePair = line.strip().split('\t')
      matchingGoList.append(genePair[0])
      matchingGoList.append(genePair[1])

      goTermsA = goTermsDict.get(genePair[0])
      goTermsB = goTermsDict.get(genePair[1])

      matchingGoList.append(goTermsA)
      matchingGoList.append(goTermsB)
      geneGOList.append(matchingGoList)

  geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'GO_Terms_A', 'GO_Terms_B'])
  geneGoDF.to_csv(outputPath + '/COMP_GO_TABLE.tsv', sep='\t', index=False)
  return geneGoDF
  