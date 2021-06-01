import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = "Get reference sequence filename")
  parser.add_argument("fpath", metavar="fpath", type = str, help = "Path to full human reference sequence")
  args = parser.parse_args()
  fpath = args.fpath
  nucout = "nuc/nuc.fna"
  mitout ="mito/mito.fna"

  nuc = open(nucout,"w")
  mit = open(mitout,"w")

  towrite = nuc
 
  with open(fpath) as fp:
    line = fp.readline()
    while line:  
      if ">" in line:
        if "mitochondrion" in line:
          towrite = mit
          print("Mitochondrial sequence")
        else:
          towrite = nuc
          print("Nuclear sequence")
        print(line.rstrip())
      towrite.write(line)
      line = fp.readline()

  nuc.close()
  mit.close()