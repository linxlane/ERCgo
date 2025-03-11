import argparse

def runParser():
  parser = argparse.ArgumentParser()

  parser.add_argument('-i', '--input', required=True, metavar='dir_path',
    help='''Path to ERCnet output files and GAF file which will be used in the GO analysis''')
  
  parser.add_argument('-o', '--output', required=True, metavar='dir_path',
    help='''Path to output folder where ERCnet GO analysis information will be written''')
  
  parser.add_argument('-r', '--random', required=True, type = int, metavar='int',
    help='''Number of random replicates to generate and analyze.''')
  
  parser.add_argument('-m', '--rand_method', required=False, default='B', choices=['AB', 'B'],
    help='''Randomize both A and B columns or just B column in randomized edge files''')
  
  return parser.parse_args()
