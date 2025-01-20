import pandas

def generate_lookup_dict(path):
  df = pandas.read_csv(path, sep='\t', usecols=['HOG_ID', 'Comprehensive_ID'])
  hog_dict = df.set_index('HOG_ID')['Comprehensive_ID'].to_dict()
  return hog_dict

def lookup(hogGene, hogCompDict):
  compGene = hogCompDict.get(hogGene)
  compGene = formatCompID(compGene)
  return compGene

def formatCompID(compGene):
  if compGene.startswith('HOG'):
    compGene = None
  
  else:
    compGene = compGene[6:]
    if len(compGene) > 10:
      compGene = compGene[:9]

  return compGene