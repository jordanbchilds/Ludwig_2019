import pysam
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def makesumm(bamfname):
  bamfile = pysam.AlignmentFile(bamfname, "rb")
  Nbases = bamfile.lengths[0]

  nucs = ["A","C","G","T"]
  summ = {"A": [0]*Nbases,"C":[0]*Nbases,"G": [0]*Nbases,"T": [0]*Nbases}
  
  i = 0
  for i,pileupcolumn in enumerate(bamfile.pileup()):
      pileups = [p for p in pileupcolumn.pileups if ((not p.is_del) and (not p.is_refskip))]
      reads = [pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pileups]
      for nuc in nucs:
        nnuc = len([r for r in reads if r==nuc])
        summ[nuc][i] = nnuc
  bamfile.close()

  summ["Coverage"] = [summ["A"][i] + summ["C"][i] + summ["G"][i] + summ["T"][i] for i,j in enumerate(summ["A"])]
  return(summ)

mtDNA_file = "mtDNA.fa"
acc_file = "ExamplePath.txt"

if not os.path.exists("reports"):
  os.mkdir("reports")

with open(mtDNA_file) as f:
   lines = f.readlines()
del lines[0]
mtseq = "".join([l.rstrip() for l in lines])

with open(acc_file) as f:
  rootstrs = f.readlines()
rootstrs = [r.rstrip() for r in rootstrs]

roots = [r.split("\t")[0] for r in rootstrs]
labs = [r.split("\t")[1] for r in rootstrs]
res = {}

for root,lab in zip(roots,labs):
  bamfname = os.path.join("bam",root+"_header.bam")
  print(bamfname)
  summ = makesumm(bamfname)
  nucs = ["A","C","G","T"]
  summ["MutFrac"] = [sum([summ[m][i] for m in [n for n in nucs if n != mtseq[i]]])/summ["Coverage"][i] if summ["Coverage"][i]>0 else float('NaN') for i,j in enumerate(summ["Coverage"])]
  res[root] = summ
  
colormap = plt.cm.rainbow #gist_ncar 
colors = [colormap(i) for i in np.linspace(0, 1, len(roots))]
xposition = [1495, 2110, 8003, 15089]

for i,(root,lab) in enumerate(zip(roots,labs)):
  summ = res[root]
  fig1 = plt.figure(figsize = (16,8))
  ax = fig1.add_subplot(111)
  plt.subplots_adjust(hspace=0.4)
  p1 = plt.subplot(2,1,1)
  l1 = plt.plot(summ["Coverage"],color=colors[i], label=lab)
  x1 = plt.xlabel('Chromosome coordinate (bases)')
  y1 = plt.ylabel('Coverage')
  for xc in xposition:
    plt.axvline(x=xc, color='k', linestyle='--',alpha=0.35)
  p2 = plt.subplot(2,1,2)
  l2 = plt.plot(summ["MutFrac"],color=colors[i], label=lab)
  x2 = plt.xlabel('Chromosome coordinate (bases)')
  y2 = plt.ylabel('Point mutation load (all point mutations)')
  for xc in xposition:
    plt.axvline(x=xc, color='k', linestyle='--',alpha=0.35)
  sttl = plt.suptitle(lab+" "+root)
  ax.legend(loc="center left", bbox_to_anchor=(1,0.5))
  plt.savefig(os.path.join("reports",root+'_Report.pdf'))

fig1 = plt.figure(figsize = (16,8))
ax = fig1.add_subplot(111)
x2 = plt.xlabel('Chromosome coordinate (bases)')
y2 = plt.ylabel('Point mutation load (all point mutations)')

for i,(root,lab) in enumerate(zip(roots,labs)):
  summ = res[root]
  ax.plot(summ["MutFrac"],alpha=0.5, label=lab)

for i,j in enumerate(ax.lines):
    j.set_color(colors[i])
	
for xc in xposition:
    plt.axvline(x=xc, color='k', linestyle='--',alpha=0.35)

ax.set_yscale('linear')
ax.set_ylim(bottom=0.0, top=1.0)
ax.legend(loc="center left", bbox_to_anchor=(1,0.5))

plt.savefig(os.path.join("reports",'Linear_Report.pdf'))

fig1 = plt.figure(figsize = (16,8))
ax = fig1.add_subplot(111)
x2 = plt.xlabel('Chromosome coordinate (bases)')
y2 = plt.ylabel('Point mutation load (all point mutations)')

for i,(root,lab) in enumerate(zip(roots,labs)):
  summ = res[root]
  ax.plot(summ["MutFrac"],alpha=0.5, label=lab)

for i,j in enumerate(ax.lines):
    j.set_color(colors[i])

for xc in xposition:
    plt.axvline(x=xc, color='k', linestyle='--',alpha=0.35)

ax.set_yscale('log')
ax.set_ylim(top=1.0)
ax.legend(loc="center left", bbox_to_anchor=(1,0.5))

plt.savefig(os.path.join("reports",'Log_Report.pdf'))

fig1 = plt.figure(figsize = (16,8))
ax = fig1.add_subplot(111)
x2 = plt.xlabel('Chromosome coordinate (bases)')
y2 = plt.ylabel('Point mutation load (all point mutations)')

for i,(root,lab) in enumerate(zip(roots,labs)):
  summ = res[root]
  ax.plot(summ["MutFrac"],alpha=0.5, label=lab)

for i,j in enumerate(ax.lines):
    j.set_color(colors[i])

for xc in xposition:
    plt.axvline(x=xc, color='k', linestyle='--',alpha=0.35)

ax.set_yscale('linear')
ax.set_ylim(bottom=0.0, top=0.3)
ax.legend(loc="center left", bbox_to_anchor=(1,0.5))

plt.savefig(os.path.join("reports",'Trimmed_Report.pdf'))

quants = {}
for root in roots:
  covs = res[root]["MutFrac"]
  covsinv = [1.0-c for c in covs]
  covmin = [min(a,b) for a,b in zip(covs,covsinv)]
  upper = np.percentile(covmin,95)
  quants[root]=upper

uppers = [quants[root] for root in roots]

#for xpos in range(0,16569):
for xpos in [821,1492,1557,5006,5261,6226,6961,6964,7162,7344,7345,7346,8001,8699,10329,10397,10399,10404,11718,12348,12852,13204,15048,15290,16128,16155,16171]:
  covs = [res[root]["MutFrac"][xpos] for root in roots]
  covsinv = [1.0-c for c in covs]
  fig, ax = plt.subplots()
  fig.canvas.draw()
  ax.set_xticks(range(0,10))
  ax.set_xticklabels(labs)
  ax.set_ylim(bottom=0.0, top=0.3)
  y2 = plt.ylabel('Point mutation load')
  plt.axhline(y=np.mean(uppers), color='k', linestyle='--',alpha=0.35)
  plt.plot(covs)
  plt.plot(covsinv)
  cc = str(xpos).zfill(5)
  plt.suptitle("Chromosome coordinate of point mutation: "+cc)
  plt.savefig(os.path.join("frames_examine",cc+".png"), dpi=600)

# Surprised to see so many point mutations with very high mutation loads in one sample
# Print out mutation loads and check in pileup file
#for i,j in enumerate(summ["MutFrac"][0:200]):
#    print(i+1,j)



