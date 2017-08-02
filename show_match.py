from pyteomics import mgf, mass
import os, sys, getopt
from urllib import urlretrieve
import pylab
import re

# Code based on the tutorial in
# https://pythonhosted.org/pyteomics/examples/example_msms.html

def main(argv):
  helpMessage = 'show_match.py --moda <moda_file> --mgf <mgf_file> --mass <theoretical_mass> --scan <scan nr>'
  modaFile, mgfFile, wantedScannr = "consensus/Latosinska_LFQ_consensus.moda.psms.tsv", "consensus/Latosinska_LFQ_consensus.filtered.mgf", 1
  try:
    opts, args = getopt.getopt(argv, "a:s:g:m:", ["moda=","mgf=","scan="])
  except getopt.GetoptError:
    print(helpMessage)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(helpMessage)
      sys.exit()
    elif opt in ("-a", "--moda"):
      modaFile = arg
    elif opt in ("-g", "--mgf"):
      mgfFile = arg
    elif opt in ("-s", "--scan"):
      wanted_scan = int(arg)
    elif opt in ("-m", "--mass"):
      wanted_theoretical_mass = arg

  if not os.path.isfile(modaFile):
    print("Could not find ModA identifications file:", modaFile)
    sys.exit(2)

  if not os.path.isfile(mgfFile):
    print("Could not find mgf file:", mgfFile)
    sys.exit(2)

  peptide = ""
  mf = moda_read(modaFile)
  for scan,charge,peptide,mass in mf:
      if mass == wanted_theoretical_mass:
          break
  peptide = peptide[2:-2]
  print(peptide,scan,charge)

  with mgf.read(mgfFile) as spectra:
    for spectrum in spectra:
        if int(spectrum['params']['scans']) == scan:
            print("found the spectrum")
            pylab.figure()
            pylab.title('Theoretical and experimental spectra for '
                + peptide)
            pylab.xlabel('m/z, Th')
            pylab.ylabel('Intensity, rel. units')

            pylab.bar(spectrum['m/z array'], spectrum['intensity array'], width=0.1, linewidth=2,
                edgecolor='black')

            theor_spectrum = list(fragments(peptide, maxcharge=charge))
            pylab.bar(theor_spectrum,
                [spectrum['intensity array'].max()]*len(theor_spectrum),
                width=0.1, edgecolor='red', alpha=0.7)
            pylab.show()
            break


def moda_read(fname):
    record = None
    scan = None
    with open(fname,"r") as mfile:
        for line in mfile:
            if line[0] == ' ' or line[0] == '\t' or len(line)<=1:
                scan = None
                continue
            if scan:
                words = line.split('\t')
                mass = words[0]
                peptide = words[4]
                record = (scan,charge,peptide,mass)
                yield record
            if line[0:2] == ">>":
                words = line.split('\t')
                scan = int(words[4])
                charge = int(words[3])

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    varmods = []
    varmodexp = re.compile("[+-][0-9]+")
    mod = varmodexp.search(peptide)
    while mod:
        varmods.append((mod.start()-1+1,int(mod.group())))
        peptide = varmodexp.sub("", peptide,1)
        mod = varmodexp.search(peptide)
    for i in xrange(1, len(peptide)-1):
        pre,post = 0,0
        for pos,diff in varmods:
            if pos<=i:
                pre += diff
            else:
                post += diff
        for ion_type in types:
            for charge in xrange(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)+pre/charge
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)+post/charge

if __name__ == "__main__":
   main(sys.argv[1:])
