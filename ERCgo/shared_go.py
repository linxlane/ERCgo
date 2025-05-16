import pandas

def collectGoTerms(genePairsPath, goTermsDict, intermediateFilesPath, argsDict):
  ##For each gene pair, get GO terms for each gene and write all associations to tsv file
  print('2. Collect GO terms for each gene in a pair and construct table', flush=True)

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

      #Get GO terms for each gene
      goTermsA = goTermsDict.get(genePair[0])
      goTermsB = goTermsDict.get(genePair[1])

      if goTermsA != None and goTermsB != None:
        matchingGoList.append(genePair[0])
        matchingGoList.append(genePair[1])
        #Append R2 values and p values
        matchingGoList.append(genePair[2])
        matchingGoList.append(genePair[3])
        matchingGoList.append(genePair[4])
        matchingGoList.append(genePair[5])
        if argsDict['analysis'] == 'full':
          matchingGoList.append(genePair[6])

        #Append GO terms for each gene to row representation list
        matchingGoList.append(goTermsA)
        matchingGoList.append(goTermsB)

        #Append 'new row' to table
        geneGOList.append(matchingGoList)

  ##Write table to file
  if argsDict['analysis'] == 'hits':
    geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'P_R2', 'P_Pval', 'S_R2', 'S_Pval','GO_Terms_A', 'GO_Terms_B'])
    print('  > Write [COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval, GO_TERMS_A, GO_TERMS_B] table to tsv', flush=True)
  
  if argsDict['analysis'] == 'full':
    geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'Slope', 'P_R2', 'P_Pval', 'S_R2', 'S_Pval', 'GO_Terms_A', 'GO_Terms_B'])
    print('  > Write [COMP_GENE_A, COMP_GENE_B, Slope, P_R2, P_Pval, S_R2, S_Pval, GO_TERMS_A, GO_TERMS_B] table to tsv', flush=True)
  writePath = intermediateFilesPath + '/gene_pairs_w_GO_terms_TABLE_' + argsDict['job_name'] + '.tsv'
  geneGoDF.to_csv(writePath, sep='\t', index=False)
  print(' > DONE', flush=True)
  return geneGoDF, writePath

def colorCode(geneA, geneB, interestGenes, clpGenes):
  colorStr = ''
  if geneA in clpGenes and geneB in clpGenes:
    colorStr = 'Clp-Clp interaction'
  elif geneA in clpGenes and geneB in interestGenes:
    colorStr = 'Clp-Interest interaction'
  elif geneA in interestGenes and geneB in clpGenes:
    colorStr = 'Interest-Clp interaction'
  elif geneA in interestGenes and geneB in interestGenes:
    colorStr = 'Interest-Interest interaction'
  else:
    colorStr = 'Not of interest'
      
  return colorStr

def analyzeSharedGo(baseDF, masterOutPath, geneGoPath, frequencies, argsDict):
  print('3. Analyze shared GO terms intersection for each pair and construct table', flush=True)

  ##Create copy of gene pairs with GO terms dataframe to add on analysis data
  sharedStatsDF = baseDF.copy(deep=True)

  #color code genes
   #Clp genes and other genes of interest for plot coloring
  interestGenes = ['AT1G06950', 'AT5G53080', 'AT2G04270', 'AT1G26460', 'AT3G12380', 
                'AT3G60830', 'AT1G10510', 'AT3G19800', 'AT1G36320']
  
  clpGenes = ['ATCG00670', 'AT1G49970', 'AT1G12410', 'AT1G09130', 
                   'AT4G17040', 'AT1G66670', 'AT4G12060']

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
  color = []

  ##Iterate through gene pairs and GO terms and conduct analysis
  with open(geneGoPath, 'r') as genePairs:
    #Skip first line with column titles
    next(genePairs)
    for line in genePairs:
      #Separate 'row' data
      lineData = line.strip().split('\t')

      geneA = lineData[0]
      geneB = lineData[1]

      color.append(colorCode(geneA, geneB, interestGenes, clpGenes))

      #Get GO terms for gene A and gene B
      if argsDict['analysis'] == 'hits':
        goListA = eval(lineData[6])
        goListB = eval(lineData[7])

      if argsDict['analysis'] == 'full':
        goListA = eval(lineData[7])
        goListB = eval(lineData[8])

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
  sharedStatsDF['Color'] = color
  sharedStatsDF['label'] = argsDict['job_name']

  print('  > Write analysis table to tsv', flush=True)
  analysisWritePath = masterOutPath + '/GO_ANALYSIS_' + argsDict['job_name'] + '.tsv'
  sharedStatsDF.to_csv(analysisWritePath, sep='\t', index=False, na_rep='N/A')
  print('   > DONE', flush=True)
  print(' > DONE', flush=True)