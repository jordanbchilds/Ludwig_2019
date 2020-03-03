mito = ">NC_012920.1 Homo sapiens mitochondrion, complete genome"

fpath = "GCA_000001405.28_GRCh38.p13_genomic.fna"
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