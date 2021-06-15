from matplotlib.pyplot import bar, title
from classes import CaloParticle, LayerCluster, RecHit, Trackster
from readTree import *
from tqdm import tqdm
from plotsUtils import *
path ="./plotsAnalysis/"
dirs = [path+"PCAWeightedNewPlot0Rechits", path+"PCAWeightedNewPlot1Rechits", path+"PCAWeightedNewPlot2Rechits", path+"PCAWeightedNewPlot3Rechits"]
files_LC = []
files_RH = []
label = []
label_RH = []
j = 1
for d in dirs[1:]:
    files_LC.append(d+"/cp_fraction005.npy")
    files_RH.append(d+"/cp_fraction005_rh.npy")
    label.append("PCA-LCs Rechits >="+str(j))
    label_RH.append("PCA-RHs Rechits >="+str(j))
    j+=1

outputPath = path+ "/Results/"
mkdir_p(outputPath)
mkdir_p(outputPath+"pickle/")

l_LC = []
l_RH = []

for f in files_LC:
    n = np.load(f)
    l_LC.append(n)

for f in files_RH:
    n = np.load(f)
    l_RH.append(n)

print(l_LC)

ranges = (0.9996, 1)
plotMultiHisto(
    l_LC,
    "Angle between PCA direction and CaloParticle direction",
    x_title="Cosine",
    y_title="Entries",
    saveName="Multi_Fraction005_zoom",
    labels=label,
    bins=100,
    ranges=ranges,
    outputdir=outputPath,
)

plotMultiHisto(
    l_RH,
    "Angle between PCA direction and CaloParticle direction",
    x_title="Cosine",
    y_title="Entries",
    saveName="Multi_Fraction005_RH",
    labels=label_RH,
    bins=100,
    ranges=ranges,
    outputdir=outputPath,
)