import pandas

def generate_lookup_dict(path):
  df = pandas.read_csv(path, sep='\t', usecols=['HOG_ID', 'Comprehensive_ID'])
  hog_dict = df.set_index('HOG_ID')['Comprehensive_ID'].to_dict()
  return hog_dict

def lookup(hogGene, hogCompDict):
  compGene = hogCompDict.get(hogGene)
  if compGene.startswith('HOG'):
    compGene = 'N/A'
  
  else:
    compGene = compGene[6:]
    if len(compGene) > 10:
      compGene = compGene[:9]

  return compGene
