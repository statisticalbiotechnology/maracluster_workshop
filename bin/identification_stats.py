import sys
import csv
import os
import getopt
import collections

ClusterIdentificationRowHeaders = ["consensus_scannr", "num_spectra_in_cluster", "qvalue", "peptide", "proteins"]
ClusterIdentificationRow = collections.namedtuple("ClusterIdentificationRow", ClusterIdentificationRowHeaders)

def main(argv):
  helpMessage = 'identification_stats.py --clusters <cluster_file> --identifications <percolator_psms> --outfile <output_file> --minsize <min_cluster_size>'
  minClusterSize = 6
  percolatorPsmsFile, clusterFile, outputFile = "", "", ""
  try:
    opts, args = getopt.getopt(argv, "i:c:o:", ["identifications=","clusters=","outfile=","minsize="])
  except getopt.GetoptError:
    print(helpMessage)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(helpMessage)
      sys.exit()
    elif opt in ("-i", "--identifications"):
      percolatorPsmsFile = arg
    elif opt in ("-c", "--clusters"):
      clusterFile = arg
    elif opt in ("-o", "--outfile"):
      outputFile = arg
    elif opt in ("-m", "--minsize"):
      minClusterSize = int(arg)
  
  if not os.path.isfile(percolatorPsmsFile):
    print("Could not find percolator identifications file:", percolatorPsmsFile)
    sys.exit(2)
  
  if not os.path.isfile(clusterFile):
    print("Could not find cluster file:", clusterFile)
    sys.exit(2)
  
  if len(outputFile) == 0:
    print("No output file specified")
    print(helpMessage)
    sys.exit(2)
  
  psms = dict()
  for psm in parsePsmsPout(percolatorPsmsFile, qThresh = 0.01, msgf = True):
    if not psm.scannr in psms:
      psms[psm.scannr] = psm
  
  outputRows = list()
  for i, cluster in enumerate(parseClusterFile(clusterFile)):
    clusterIdx = i + 1
    if len(cluster) >= minClusterSize:
      psm = psms.get(clusterIdx, getDefaultPoutPsms())
      outputRows.append(ClusterIdentificationRow(clusterIdx, len(cluster), psm.qvalue, psm.peptide, ",".join(psm.proteins)))
  
  outputRows = sorted(outputRows, key = lambda x : x[1], reverse = True)
  
  writer = csv.writer(open(outputFile, 'w'), delimiter = '\t')
  writer.writerow(ClusterIdentificationRowHeaders)
  for row in outputRows:
    writer.writerow(row)

ClusterSpectrumEntry = collections.namedtuple("ClusterSpectrumEntry", "filePath scannr")

def parseClusterFile(clusterFile):
  reader = csv.reader(open(clusterFile, 'r'), delimiter='\t')
  rows = list()
  for row in reader:
    if len(row) >= 2:
      rows.append(ClusterSpectrumEntry(row[0], int(row[1])))
    else:
      if len(rows) > 0:
        yield rows
      rows = list()

PercolatorPoutPsmsBase = collections.namedtuple("PercolatorPoutPsms", "id filename scannr charge svm_score qvalue PEP peptide proteins")

class PercolatorPoutPsms(PercolatorPoutPsmsBase):
  def toList(self):
    return [self.id, self.svm_score, self.qvalue, self.PEP, self.peptide] + self.proteins
  
  def toString(self):
    return "\t".join(map(str, self.toList()))
    
def getDefaultPoutPsms():
  return PercolatorPoutPsms("", "", 0, 0, 0, 1.0, 1.0, "NA", ["NA"])
  
# works for peptide and psms pouts
def parsePsmsPout(poutFile, qThresh = 1.0, msgf = True):
  reader = csv.reader(open(poutFile, 'r'), delimiter='\t')
  headers = next(reader) # save the header
  
  for row in reader:
    if float(row[2]) <= qThresh:
      proteins = row[5:]
      yield PercolatorPoutPsms(row[0], getFileName(row[0], msgf), getId(row[0], msgf), getCharge(row[0]), float(row[1]), float(row[2]), float(row[3]), row[4], proteins)
    else:
      break

def getId(PSMid, msgf = False):
  if msgf:
    return int(int((PSMid.split('_'))[-3]) / 100)
  else:
    return int((PSMid.split('_'))[-3])

def getCharge(PSMid):
  return int((PSMid.split('_'))[-2])

def getFileName(PSMid, msgf = False):
  if msgf:
    return '_'.join(PSMid.split('_')[:-6])
  else:
    return '_'.join(PSMid.split('_')[:-3])
    
if __name__ == "__main__":
   main(sys.argv[1:])  
