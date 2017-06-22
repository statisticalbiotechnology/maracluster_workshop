import sys
import csv
import os
import getopt
import collections

def main(argv):
  helpMessage = 'extract_unidentified_spectra.py --spectra <mgf_file> --identifications <cluster_identifications> --outfile <output_file> --minsize <min_cluster_size>'
  minClusterSize = 70
  try:
    opts, args = getopt.getopt(argv, "s:i:o:m:", ["spectra=","identifications=","outfile=","minsize="])
  except getopt.GetoptError:
    print(helpMessage)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(helpMessage)
      sys.exit()
    elif opt in ("-s", "--spectra"):
      mgfFile = arg
    elif opt in ("-i", "--identifications"):
      clusterIdentificationsFile = arg
    elif opt in ("-o", "--outfile"):
      outputFile = arg
    elif opt in ("-m", "--minsize"):
      minClusterSize = int(arg)
  
  if not os.path.isfile(mgfFile):
    print("Could not find mgf spectra file:", mgfFile)
    sys.exit(2)
  
  if not os.path.isfile(clusterIdentificationsFile):
    print("Could not find cluster identifications file:", clusterIdentificationsFile)
    sys.exit(2)
  
  if len(outputFile) == 0:
    print("No output file specified")
    print(helpMessage)
    sys.exit(2)
     
  unidentifiedScanNrs = list()
  for clusterId in parseclusterIdentificationsFile(clusterIdentificationsFile):
    if clusterId.qvalue > 0.01 and clusterId.num_spectra_in_cluster >= minClusterSize:
      unidentifiedScanNrs.append(clusterId.consensus_scannr)
  unidentifiedScanNrs = set(unidentifiedScanNrs)
  
  with open(outputFile, 'w') as f:
    for scannr, spectrum in parseMgf(mgfFile):
      if scannr in unidentifiedScanNrs:
        for line in spectrum:
          f.write(line)
          if line.startswith("TITLE"):
            f.write("SCANS=%d\n" % scannr)

ClusterIdentificationRowHeaders = ["consensus_scannr", "num_spectra_in_cluster", "qvalue", "peptide", "proteins"]
ClusterIdentificationRow = collections.namedtuple("ClusterIdentificationRow", ClusterIdentificationRowHeaders)

def parseclusterIdentificationsFile(clusterIdentificationsFile):
  reader = csv.reader(open(clusterIdentificationsFile, 'r'), delimiter = '\t')
  next(reader) # skip header
  
  for row in reader:
    yield ClusterIdentificationRow(int(row[0]), int(row[1]), float(row[2]), row[3], row[4])

def parseMgf(inputFile):
  scannr = -1
  spectrum = list()
  with open(inputFile, 'r') as f:
    for line in f:
      if line.startswith('BEGIN IONS'):
        if len(spectrum) > 0:
          yield scannr, spectrum
        spectrum = list()
      if line.startswith('TITLE'):
        scannr = int(int(line.split('=')[2]) / 100)
      spectrum.append(line)
      
if __name__ == "__main__":
   main(sys.argv[1:])  
