import pandas
from collections import Counter

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
  print('> Write [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B] table to tsv: [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B].tsv')
  geneGoDF.to_csv(outputPath, sep='\t', index=False)
  return geneGoDF

def generatePopulation(uniqueIds, assoc_dict):
  goTermsPopulation = []
  keyErrorCounter = 0
  for gene in uniqueIds:
    if gene in assoc_dict.keys():
      goTermsPopulation.extend(assoc_dict[gene])
    else:
      keyErrorCounter += 1
  #print(goTermsPopulation)
  print(f'Number of keys not found: {keyErrorCounter}')

  populationCounts = Counter(goTermsPopulation)
  return populationCounts

def overlapScore(setA, setB, population):
  weightedIntersection = 0
  
  totalSetLengths = len(setA) + len(setB)

  if totalSetLengths == 0:
    overlapScore = 0.0
  else:
    for term in setA & setB:
      weightedIntersection += 1/population[term]
    
    overlapScore = (2 * weightedIntersection / totalSetLengths) * (1 - (1 / totalSetLengths))

  return overlapScore

def sanityCheck(baseDF, writePath, geneGoPath, counts):
  sanityDF = baseDF.copy(deep=True)

  lengthGoAList = []
  lengthGoBList = []
  goTermIntersectionList = []
  sharedGoLen = []
  countsA = []
  countsA_List = []
  countsB = []
  countsB_List = []
  weightedIntersectionList = []
  weightedIntersection = 0.0
  totalSetLengthsList = []
  overlapScore = 0.0
  overlapScoreList = []
  formula = ''
  formulaList = []
  termntersectionCount = []
  intersectionCountsList = []

  with open(geneGoPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      lineData = line.strip().split('\t')
      #print(lineData)

      goListA = eval(lineData[2])
      goListB = eval(lineData[3])

      goSetA = set(goListA)
      goSetB = set(goListB)

      lengthGoAList.append(len(goListA))
      lengthGoBList.append(len(goListB))

      goTermIntersection = goSetA & goSetB
      goTermIntersectionList.append(goTermIntersection)
      sharedGoLen.append(len(goTermIntersection))

      for termA in goListA:
        #print(term, counts[term])
        countsA.append(counts[termA])
      #print(countsA)
      countsA_List.append(countsA)
      countsA = []

      for termB in goListB:
        countsB.append(counts[termB])
      countsB_List.append(countsB)
      countsB = []

      totalSetLengths = len(goSetA) + len(goSetB)
      totalSetLengthsList.append(totalSetLengths)

      if totalSetLengths == 0:
        overlapScore = 0.0
        weightedIntersection = 0.0
      else:
        for term in goSetA.intersection(goSetB):
          termCount = counts[term]
          termntersectionCount.append(termCount)
          weightedIntersection += 1/termCount
          overlapScore = (2 * weightedIntersection / totalSetLengths) * (1 - (1 / totalSetLengths))
      weightedIntersectionList.append(weightedIntersection)
      weightedIntersection = 0.0
      intersectionCountsList.append(termntersectionCount)
      termntersectionCount = []
      overlapScoreList.append(overlapScore)
      overlapScore = 0.0

      formula = f'(2 * {weightedIntersection} / {totalSetLengths}) * (1 - (1 / {totalSetLengths}))'
      formulaList.append(formula)

  sanityDF['Length_Go_A'] = lengthGoAList
  sanityDF['Length_Go_B'] = lengthGoBList
  sanityDF['Shared_GO'] = goTermIntersectionList
  sanityDF['Number_of_Shared_GO'] = sharedGoLen  
  sanityDF['Population_Counts_A'] = countsA_List
  sanityDF['Population_Counts_B'] = countsB_List
  sanityDF['Intersection_Counts'] = intersectionCountsList
  sanityDF['Weighted_Intersection'] = weightedIntersectionList
  sanityDF['Length_SetA_and_SetB'] = totalSetLengthsList
  sanityDF['Formula'] = formulaList
  sanityDF['Overlap_Score'] = overlapScoreList

  sanityDF.to_csv(writePath, sep='\t', index=False, na_rep='N/A')


def compareGoTerms(geneGoPath, geneGoDF, outputPath, edgeFileName, population):
  goTermIntersectionList = []
  sharedGoLen = []
  maxPossibleShared = []
  propSharedList = []
  overlapScoreList = []

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

      overlapScoreList.append(overlapScore(set(goListA), set(goListB), population))

  geneGoDF['Shared_GO'] = goTermIntersectionList
  geneGoDF['Number_of_Shared_GO'] = sharedGoLen
  geneGoDF['Max_Shared_GO'] = maxPossibleShared
  geneGoDF['Observed/Max_Shared_GO'] = propSharedList
  geneGoDF['Overlap_Score'] = propSharedList
  geneGoDF['label'] = edgeFileName
  print('> Write [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B, Shared_GO, Number_of_Shared_GO, Max_Shared_GO, Observed/Max_Shared_GO, Overlap_Score, label] table to tsv: SHARED_GO_TABLE')
  geneGoDF.to_csv(outputPath, sep='\t', index=False, na_rep='N/A')
