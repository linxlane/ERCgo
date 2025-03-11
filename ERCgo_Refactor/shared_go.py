import os
import hog_comp_ids
import population
import pandas

def generateSharedGOTable(masterOutPath, edgeFilePath, hogCompDict, geneGoDict, edgeFileName):
  #Make directory to write conversion files for later reference
  intermediateFilesPath = masterOutPath + '/Intermediate_Files_' + edgeFileName
  os.makedirs(intermediateFilesPath)

  ##Convert HOG and comp ids
  print('1. Convert and prepare gene pairs for analysis')
  print('  > Generate table of [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]')
  hogCompPath = intermediateFilesPath + '/[HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]_FULL_' + edgeFileName + '.tsv'
  #Create and write tsv that contains equivalent HOG-COMP id values
  geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict, hogCompPath)

  print('  > Drop rows that contain None values, ie there was no COMP_ID match for a given HOG_ID')
  hogCompNaPath = intermediateFilesPath + '/[HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]_DROP_NA_' + edgeFileName + '.tsv'
  compPairsNaPath = intermediateFilesPath + '/[COMP_GENE_A, COMP_GENE_B]_DROP_NA_' + edgeFileName + '.tsv'
  geneLookupDF_DROP = hog_comp_ids.dropNaRows(geneLookupDF_FULL, hogCompNaPath, compPairsNaPath)

  ##Summarize GO term population for all GO terms associated with genes in the current edge file
  print('2. Gather GO term population of all genes in edge file')
  #Get unique genes from column A and B, returns set of gene comp ids
  print('  > Get unique genes from column A and B')
  uniqueIds = population.populationUniqueGenes(compPairsNaPath)

  #Use counter container to get frequencies of each go term in the population
  print('  > Create counter for GO term frequencies')
  goTermFreq = population.countGoTermsFrequency(uniqueIds, geneGoDict)
  print('  > Write frequency counter to file')
  populationWritePath = masterOutPath + '/Population_Frequencies_' + edgeFileName + '.tsv'
  population.writePopulationFrequencies(populationWritePath, goTermFreq)

  ##Get GO terms for each gene in gene pairs in edge file, write table to tsv
  print('3. Collect GO terms for each gene in a pair and construct table')
  compPairsGoWritePath = intermediateFilesPath + '/[COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B]_TABLE_' + edgeFileName + '.tsv'
  genePairWithGoDF = genePairGO(compPairsNaPath, geneGoDict, compPairsGoWritePath)

  return genePairWithGoDF, compPairsGoWritePath, goTermFreq


def genePairGO(genePairsPath, goTermsDict, outputPath):
  ##For each gene pair, get GO terms for each gene and write all associations to tsv file

  #Create new list to store each row of gene pairs and GO terms
  geneGOList = []

  #Iterate through each gene pair and get GO terms
  with open(genePairsPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      #New list to represent new row in table
      matchingGoList = []
  
      #Get each gene in gene pair and append each to row representation list
      genePair = line.strip().split('\t')
      matchingGoList.append(genePair[0])
      matchingGoList.append(genePair[1])

      #Get GO terms for each gene
      goTermsA = goTermsDict.get(genePair[0])
      goTermsB = goTermsDict.get(genePair[1])

      #Check for no GO terms to assign empty list
      if goTermsA == None:
        goTermsA = []

      if goTermsB == None:
        goTermsB = []

      #Append GO terms for each gene to row representation list
      matchingGoList.append(goTermsA)
      matchingGoList.append(goTermsB)

      #Append 'new row' to table
      geneGOList.append(matchingGoList)

  ##Write table to file
  geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'GO_Terms_A', 'GO_Terms_B'])
  print('  > Write [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B] table to tsv: [COMP_GENE_A, COMP_GENE_B, GO_TERMS_A, GO_TERMS_B].tsv')
  geneGoDF.to_csv(outputPath, sep='\t', index=False)
  return geneGoDF


def analyzeSharedGo(baseDF, masterOutPath, geneGoPath, frequencies, edgeFileName):
  print('4. Analyze shared GO terms intersection for each pair and construct table')

  ##Create copy of gene pairs with GO terms dataframe to add on analysis data
  sharedStatsDF = baseDF.copy(deep=True)

  ##Initialize new rows and variables to conduct calculations/analysis
  lengthGoAList = []
  lengthGoBList = []
  goTermIntersectionList = []
  sharedGoLen = []
  frequenciesA = []
  frequenciesA_List = []
  frequenciesB = []
  frequenciesB_List = []
  weightedIntersectionList = []
  weightedIntersection = 0.0
  totalSetLengthsList = []
  overlapScore = 0.0
  overlapScoreList = []
  formula = ''
  formulaList = []
  termIntersectionCount = []
  intersectionFrequenciesList = []

  ##Iterate through gene pairs and GO terms and conduct analysis
  with open(geneGoPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      #Separate 'row' data
      lineData = line.strip().split('\t')

      #Get GO terms for gene A and gene B
      goListA = eval(lineData[2])
      goListB = eval(lineData[3])

      #Convert GO term data to sets
      goSetA = set(goListA)
      goSetB = set(goListB)

      #Get the length of each GO term set and append to appropriate column 
      lengthGoAList.append(len(goListA))
      lengthGoBList.append(len(goListB))

      #Determine intersection, if any, of GO terms for each gene in the gene pair
      goTermIntersection = goSetA & goSetB
      #Append intersecting GO terms to column list
      goTermIntersectionList.append(goTermIntersection)
      #Get and append number of intersecting GO terms to column list
      sharedGoLen.append(len(goTermIntersection))

      ##Get frequency counts for all GO terms for each gene in gene pair
      #Frequency counts for each GO term associated with gene A
      for termA in goListA:
        #Get each frequency and append to list
        frequenciesA.append(frequencies[termA])
      #Append frequency list to column
      frequenciesA_List.append(frequenciesA)
      #Reset for next row
      frequenciesA = []

      #Frequency counts for each GO term associated with gene B
      for termB in goListB:
        #Get each frequency and append to list
        frequenciesB.append(frequencies[termB])
      #Append frequency list to column
      frequenciesB_List.append(frequenciesB)
      #Reset for next row
      frequenciesB = []

      ##Calculate 
      totalSetLengths = len(goSetA) + len(goSetB)
      totalSetLengthsList.append(totalSetLengths)

      uniqueGoTermsLenth = len(goSetA.union(goSetB))

      if uniqueGoTermsLenth == 0:
        overlapScore = 0.0
        weightedIntersection = 0.0
      else:
        for term in goSetA.intersection(goSetB):
          termCount = frequencies[term]
          termIntersectionCount.append(termCount)
          weightedIntersection += 1/termCount
          overlapScore = (weightedIntersection / uniqueGoTermsLenth)
      formula = f'({weightedIntersection} / {uniqueGoTermsLenth})'
      formulaList.append(formula)
      weightedIntersectionList.append(weightedIntersection)
      weightedIntersection = 0.0
      intersectionFrequenciesList.append(termIntersectionCount)
      termIntersectionCount = []
      overlapScoreList.append(overlapScore)
      overlapScore = 0.0

  sharedStatsDF['Length_Go_A'] = lengthGoAList
  sharedStatsDF['Length_Go_B'] = lengthGoBList
  sharedStatsDF['Shared_GO'] = goTermIntersectionList
  sharedStatsDF['Number_of_Shared_GO'] = sharedGoLen  
  sharedStatsDF['Population_Frequencies_A'] = frequenciesA_List
  sharedStatsDF['Population_Frequencies_B'] = frequenciesB_List
  sharedStatsDF['Intersection_Frequencies'] = intersectionFrequenciesList
  sharedStatsDF['Weighted_Intersection'] = weightedIntersectionList
  sharedStatsDF['Length_SetA_and_SetB'] = totalSetLengthsList
  sharedStatsDF['Formula'] = formulaList
  sharedStatsDF['Overlap_Score'] = overlapScoreList
  sharedStatsDF['label'] = edgeFileName

  print('  > Write analysis table to tsv: GO_ANALYSIS_.tsv')
  analysisWritePath = masterOutPath + '/GO_ANALYSIS_' + edgeFileName + '.tsv'
  sharedStatsDF.to_csv(analysisWritePath, sep='\t', index=False, na_rep='N/A')
