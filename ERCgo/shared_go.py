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

      if goTermsA == None:
        goTermsA = []

      if goTermsB == None:
        goTermsB = []

      matchingGoList.append(goTermsA)
      matchingGoList.append(goTermsB)
      geneGOList.append(matchingGoList)

  geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'GO_Terms_A', 'GO_Terms_B'])
  print('Write [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B] table to tsv: COMP_GO_TABLE.tsv')
  geneGoDF.to_csv(outputPath + '/COMP_GO_TABLE.tsv', sep='\t', index=False)
  return geneGoDF
  
def compareGoTerms(geneGoPath, geneGoDF, outputPath, edgeFileName):
  goTermIntersectionList = []
  sharedGoLen = []
  maxPossibleShared = []
  propSharedList = []
  with open(geneGoPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      lineData = line.strip().split('\t')
      #print(lineData)

      goListA = eval(lineData[2])
      goListB = eval(lineData[3])

      goTermIntersection = set(goListA) & set(goListB)
      goTermIntersectionList.append(goTermIntersection)
      sharedGoLen.append(len(goTermIntersection))

      shortestList = min(len(goListA), len(goListB))
      maxPossibleShared.append(shortestList)

      if shortestList == 0:
        propSharedGO = None
      else:  
        propSharedGO = len(goTermIntersection)/shortestList
      propSharedList.append(propSharedGO)

  geneGoDF['Shared_GO'] = goTermIntersectionList
  geneGoDF['Number_of_Shared_GO'] = sharedGoLen
  geneGoDF['Max_Shared_GO'] = maxPossibleShared
  geneGoDF['Observed/Max_Shared_GO'] = propSharedList
  print('Write [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B, Shared_GO, Number_of_Shared_GO, Max_Shared_GO, Observed/Max_Shared_GO] table to tsv: SHARED_GO_TABLE')
  geneGoDF.to_csv(outputPath + '/SHARED_GO_TABLE_' + edgeFileName + '.tsv', sep='\t', index=False, na_rep='N/A')