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
    geneGoDF = pandas.DataFrame(geneGOList, columns=['COMP_GENE_A', 'COMP_GENE_B', 'P_R2', 'P_Pval', 'S_R2', 'S_Pval', 'Slope', 'GO_Terms_A', 'GO_Terms_B'])
    print('  > Write [COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval, Slope, GO_TERMS_A, GO_TERMS_B] table to tsv', flush=True)
  writePath = intermediateFilesPath + '/gene_pairs_w_GO_terms_TABLE_' + argsDict['job_name'] + '.tsv'
  geneGoDF.to_csv(writePath, sep='\t', index=False)
  print(' > DONE', flush=True)
  return geneGoDF, writePath

def colorCode(geneA, geneB, alphaList, betaList, rpnList, rptList, allInterestGenes):
  colorStr = ''
  if geneA in allInterestGenes or geneB in allInterestGenes:
    if geneA in alphaList:
      colorStr += 'Alpha-'
    elif geneA in betaList:
      colorStr += 'Beta-'
    elif geneA in rpnList:
      colorStr += 'RPN-'
    elif geneA in rptList:
      colorStr += 'RPT-'
    else:
      colorStr += 'Other-'
    
    if geneB in alphaList:
      colorStr += 'Alpha'
    elif geneB in betaList:
      colorStr += 'Beta'
    elif geneB in rpnList:
      colorStr += 'RPN'
    elif geneB in rptList:
      colorStr += 'RPT'
    else:
      colorStr += 'Other'

  else:
    colorStr = 'Not of interest'
      
  return colorStr

def analyzeSharedGo(baseDF, masterOutPath, geneGoPath, frequencies, argsDict):
  print('3. Analyze shared GO terms intersection for each pair and construct table', flush=True)

  ##Create copy of gene pairs with GO terms dataframe to add on analysis data
  sharedStatsDF = baseDF.copy(deep=True)

  #color code genes
  alphaList = ['AT2G05840',
                'AT5G35590',
                'AT1G16470',
                'AT1G79210',
                'AT3G22110',
                'AT4G15165',
                'AT3G51260',
                'AT5G66140',
                'AT1G53850',
                'AT3G14290',
                'AT1G47250',
                'AT5G42790',
                'AT2G27020']

  betaList = ['AT4G31300',
              'AT3G27430',
              'AT5G40580',
              'AT1G21720',
              'AT1G77440',
              'AT3G22630',
              'AT4G14800',
              'AT1G13060',
              'AT3G26340',
              'AT3G60820',
              'AT1G56450']

  rpnList = ['AT4G38630',
              'AT5G23540',
              'AT1G64520',
              'AT5G42040',
              'AT2G20580',
              'AT4G28470',
              'AT1G04810',
              'AT2G32730',
              'AT1G20200',
              'AT1G75990',
              'AT5G09900',
              'AT5G64760',
              'AT1G29150',
              'AT4G24820',
              'AT3G11270',
              'AT5G05780',
              'AT4G19006',
              'AT5G45620']

  rptList = ['AT1G53750',
              'AT1G53780',
              'AT2G20140',
              'AT4G29040',
              'AT5G58290',
              'AT1G45000',
              'AT5G43010',
              'AT1G09100',
              'AT3G05530',
              'AT5G19990',
              'AT5G20000']
  
  allInterestGenes = alphaList + betaList + rpnList + rptList

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

      color.append(colorCode(geneA, geneB, alphaList, betaList, rpnList, rptList, allInterestGenes))

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