import os
import sys
import gzip
import time

def ExtractDSSP(f):
  flagISDSSP = False
  flagFOUNDRES = False
  flagFOUNDSS = False
  flagFOUNDCH = False

  residueIndex = 0
  ssIndex = 0
  chIndex = 0
  index = 0

  chains = {}
  currSeq = ""
  currSS = ""

  for line in f.readlines():
    line = line.decode('utf-8')
    #Check if in a loop
    if line.find("_dssp_struct_summary") != -1 and flagISDSSP == False:
      flagISDSSP = True
      index = 0

    elif line.strip() == "#":
      flagISDSSP = False
      flagFOUNDRES = False
      flagFOUNDSS = False
      flagFOUNDCH = False

    if flagISDSSP:
      # print(line, flagFOUNDRES, flagFOUNDSS)
      point = line.find(".")
      if line[point:].strip() == ".label_comp_id":
        flagFOUNDRES = True
        residueIndex = index
        # print(f"res {index}")
      if line[point:].strip() == ".secondary_structure":
        flagFOUNDSS = True
        ssIndex = index
        # print(f"ss {index}")
      if line[point:].strip() == ".label_asym_id":
        flagFOUNDCH = True
        chIndex = index
      index += 1

      if line.find("_dssp_struct_summary") == -1 and flagFOUNDRES and flagFOUNDSS:
        spl = line.split()
        if flagFOUNDCH:
          chain = spl[chIndex]
        else:
          chain = '$'

        if chain not in chains.keys():
          chains[chain] = ["", ""]

        chains[chain][0] += spl[residueIndex] + " "
        chains[chain][1] += spl[ssIndex]
  return chains



AAs = {
    "ARG" : "R",
    "HIS" : "H",
    "LYS" : "K",
    "ASP" : "D",
    "GLU" : "E",
    "SER" : "S",
    "THR" : "T",
    "ASN" : "N",
    "GLN" : "Q",
    "CYS" : "C",
    "GLY" : "G",
    "PRO" : "P",
    "ALA" : "A",
    "VAL" : "V",
    "ILE" : "I",
    "LEU" : "L",
    "MET" : "M",
    "PHE" : "F",
    "TYR" : "Y",
    "TRP" : "W",
}

def Three2One(chains):
  for key in chains.keys():
    newStr = ""
    spl = chains[key][0].split()
    for aa in spl:
      if aa in AAs.keys():
        newStr += AAs[aa]
      else:
        newStr += "X"

    chains[key][0] = newStr

def WriteChains(chains, cif, todir):
  if len(chains) == 0:
    return

  with open(f"{todir}/{cif}.fbs", "w") as fp:
    for key in chains.keys():
      fp.write(f">{cif}|{key}|sequence\n")
      fp.write(f"{chains[key][0]}\n")
      fp.write(f">{cif}|{key}|secondary_structure\n")
      fp.write(f"{chains[key][1]}\n")




path = sys.argv[1]
todir = sys.argv[2]

print("Preparing directory")
if not os.path.exists(todir) and not os.path.isdir(todir):
  os.mkdir(todir)
print("Preparing to extract files")
print("Extracting...")
print("Getting data")
total = 0
t1 = time.time()
for ext in os.listdir(path=path):
  extended_path = path + f"/{ext}"
  total += len(os.listdir(path=extended_path))
print(f"Completed initial check in {time.time() - t1} seconds")

print("Beginning DSSP extraction")
    

iterations = 0
for ext in os.listdir(path=path):
  extended_path = path + f"/{ext}"
  for cif in os.listdir(path=extended_path):
    iterations += 1
    full_path = extended_path + f"/{cif}"
    print(f"({iterations}/{total}) | Total time: {time.time() - t1 : .4f} |Reading {full_path}", end="\r")
    if os.path.exists(f"{todir}/{cif[:4]}.fbs"):
        continue
    
    with gzip.open(full_path, 'rb') as f:
      chains = ExtractDSSP(f)
      Three2One(chains)
      WriteChains(chains, cif[:4], todir)
