import argparse

def runParser():
  parser = argparse.ArgumentParser()

  parser.add_argument('-j', '--job_name', required=True,
    help='''Job name for this run of ERCgo. If a directory with this job name already exists at the output path, it will be erased and rewritten.
            Avoid including spaces or special characters ("_" is ok)''')

  parser.add_argument('-i', '--input', required=True, metavar='dir_path',
    help='''Path to ERCnet output files and GAF file which will be used in the GO analysis''')
  
  parser.add_argument('-a', '--analysis', required=True, choices=['full', 'hits'],
    help='''Specifies type of ERCgo analysis. 'full' will analyze entire ERCnet analysis and requires the output ERC_results tsv from an ERCnet run.
            'hits' will analyze only the hits of an ERCnet run. This mode requires the edges and vertices files. 'both' performs a GO analysis on the full ERCnet output
            and the hit data.''')
  
  parser.add_argument('-o', '--output', required=False, metavar='dir_path', default=None,
    help='''Path where new directory for ERCgo output will be created with the job name. 
            If not included, it will be written in the ERCgo directory.
            If this path already exists, it will be deleted and a new directory will be created at the output_directory/job_name path.
            ''')
  
  return parser.parse_args()
