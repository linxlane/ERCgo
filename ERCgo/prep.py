import pandas

def generate_lookup_dict(path):
  df = pandas.read_csv(path, sep='\t', usecols=['HOG_ID', 'Comprehensive_ID'])
  hog_dict = df.set_index('HOG_ID')['Comprehensive_ID'].to_dict()
  return hog_dict

