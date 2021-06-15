from matplotlib.pyplot import bar, title
from numpy.lib.financial import ipmt
from classes import CaloParticle, LayerCluster, RecHit, Trackster
from readTree import *
from tqdm import tqdm
from plotsUtils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--plot', type = bool, help='make plots')
parser.add_argument('--debug', type = bool, help='debugging mode')
args = parser.parse_args()
plotflag = bool(args.plot)
debug = bool(args.debug)

plt.style.use([hep.style.ROOT, hep.style.firamath])

hep.rcParams.label.data = True
hep.rcParams.label.paper = False


# Loading Tree
print("Loading Tree \n")
dataDir = "./data/"
rootFileName = "output_singlephoton_1k_vertexing_cpNew_eta1p7.root"
# rootFileName = "output_SingleGamma_BS.root"
fileRoot = ROOT.TFile(dataDir + rootFileName)
tree = fileRoot.Get("trackstersAnalyzer/tree")
outputDir = "./plotsAnalysis/"
# outputSubDir = "PCAWeightedDZVertex/"
outputSubDir = "PCAWeightedVertexingEta1p73RhitsCutFraction0p05/"
outputPath = outputDir + outputSubDir
outputZoomSubDir = outputPath + "Zoom/"
outputPathExample = outputPath + "Example/"
mkdir_p(outputPath)
mkdir_p(outputPath + "pickle")
mkdir_p(outputZoomSubDir)
mkdir_p(outputZoomSubDir + "pickle")
mkdir_p(outputPathExample)

getBranches(tree)

if(debug):
    Nevts = 10
else:
    Nevts = tree.GetEntries()

vertex_positions = []
CaloParticleCollection = []
TrackstersCollection = []
beam_spots = []
vertices = []
lx = []
ly = []
lz = []
ind = 0
for evt in range(Nevts):
    tree.GetEntry(evt)

    # print(f"\nEvent number: {evt}")

    # Loop over number of trackster
    for i in range(N_trackster.at(0)):  # loop over tracksters
        # print(i)

        t = Trackster(evt, ind, [])
        ind += 1
        tmp_hit_idx = 0
        if N_lcs.size() > 0:
            for j in range(N_lcs.at(0)):  # loop over layerclusters
                # print(N_lcs.at(0))
                x = lc_x.at(j)
                y = lc_y.at(j)
                z = lc_z.at(j)
                E = lc_E.at(j)
                lc = LayerCluster(x, y, z, E, i, [])
                if N_hit_per_lc.at(j) >= 3:
                    for k in range(N_hit_per_lc.at(j)):  # loop over hit
                        # be carefull with the index, is a single array containing the rechits for all the layer clusters
                        new_k = tmp_hit_idx + k
                        hit_x = rh_x.at(new_k)
                        hit_y = rh_y.at(new_k)
                        hit_z = rh_z.at(new_k)
                        hit_E = rh_E.at(new_k)
                        hit_fraction = rh_fraction.at(new_k)
                        hit = RecHit(hit_x, hit_y, hit_z, hit_E, hit_fraction)
                        lc.addRechit(hit)
                    lx.append(lc.x)
                    ly.append(lc.y)
                    lz.append(abs(lc.z))
                    t.addLayerCluster(lc)
                tmp_hit_idx += N_hit_per_lc.at(j)
                
            cp = CaloParticle(
                cp_px.at(i),
                cp_py.at(i),
                cp_pz.at(i),
                cp_pt.at(i),
                cp_E.at(i),
                cp_eta.at(i),
                cp_phi.at(i),
                cp_x.at(i),
                cp_y.at(i),
                cp_z.at(i),
                evt,
                i,
            )
            # cp = CaloParticle(
            #     cp_px.at(i),
            #     cp_py.at(i),
            #     cp_pz.at(i),
            #     cp_pt.at(i),
            #     cp_E.at(i),
            #     cp_eta.at(i),
            #     cp_phi.at(i),
            #     0,
            #     0,
            #     0,
            #     evt,
            #     i,
            # )
            # print(cp_eta.at(i))
            CaloParticleCollection.append(cp)
            TrackstersCollection.append(t)
            vertex_pos = [v_x.at(0), v_y.at(0), v_z.at(0)]
            vertex_positions.append(vertex_pos)
            vtxs = []
            '''
            for k in range(vtxs_x.size()):
                vtx =  [vtxs_x.at(k),vtxs_y.at(k),vtxs_z.at(k)]
                vtxs.append(vtx)
            vertices.append(vtxs)
            # for k in range(bs_x.size()):
            bs = [bs_x.at(i),bs_y.at(i),bs_z.at(i)]
            beam_spots.append(bs)
            '''

lcsdict = {"lcx":lx, "lcy":ly, "lcz":lz}
lcsdf = pd.DataFrame(lcsdict)
lcsgroupz = lcsdf.groupby(['lcz'])
list_of_layerz = []
for key, item in lcsgroupz:
    list_of_layerz.append(float(key))



print(len(CaloParticleCollection), len(TrackstersCollection), len(beam_spots), len(vertices))

list_of_layerz_sorted = sorted(list_of_layerz)
layers_number = [i+1 for i in range(len(list_of_layerz_sorted))]
list_of_layer_list = [[] for _ in range(len(list_of_layerz_sorted))]
list_of_layer_list_dx = [[] for _ in range(len(list_of_layerz_sorted))]
list_of_layer_list_dy = [[] for _ in range(len(list_of_layerz_sorted))]
# Process all the data and make all the plots
cp_vs_tracksterRH_filter = []  # list of list
cp_vs_tracksterLC_filter = []
# fractions = np.arange(0.001, 0.005, step=0.001)
# fractions_lc = [0.05, 0.1, 0.2, 0.3]
fractions_lc = [0]
fractions = [0]
real_fraction_used = []


weighting = True
scaling = False
nonlinearweighting = False
bestLCS = True
filter_energy = 0.1


angle_yz_planeLC = []
angle_xy_planeLC = []
angle_xz_planeLC = []

angle_yz_planeRH = []
angle_xy_planeRH = []
angle_xz_planeRH = []

eigen_value2 = []
eigen_value3 = []
eigen_value2_bestLC = []
eigen_value3_bestLC = []
eigen_value2_filterLC = []
eigen_value3_filterLC = []

profile_lc_pca_cp_pt = ROOT.TProfile("Angles1", "Angles1", 30, 0, 150, 0.5, 1.5)
profile_bestlc_pca_cp_pt = ROOT.TProfile("Angles2", "Angles2", 30, 0, 150, 0.5, 1.5)
profile_bestlc_pca_filter_cp_pt = ROOT.TProfile(
    "Angles3", "Angles3", 30, 0, 150, 0.5, 1.5
)

profile_lc_pca_cp_eta = ROOT.TProfile("Angles11", "Angles11", 30, 1.5, 2.8, 0.5, 1.5)
profile_bestlc_pca_cp_eta = ROOT.TProfile("Angles22", "Angles22", 30, 1.5, 2.8, 0.5, 1.5)
profile_bestlc_pca_filter_cp_eta = ROOT.TProfile(
    "Angles33", "Angles33", 30, 1.5, 2.8, 0.5, 1.5
)

profile_rh_pca_cp_pt = ROOT.TProfile("Angles12", "Angles1", 30, 0, 150, 0.5, 1.5)
profile_bestrh_pca_cp_pt = ROOT.TProfile("Angles21", "Angles2", 30, 0, 150, 0.5, 1.5)
profile_bestrh_pca_filter_cp_pt = ROOT.TProfile(
    "Angles31", "Angles3", 30, 0, 150, 0.5, 1.5
)

lcs_layers_dxy = ROOT.TProfile(
    "Angles31", "Angles3", 1000, 310, 380, -10, 10
)

cp_vs_tracksterLC = []  # CP vs WPCA-LC standard
cp_vs_tracksterRH = []  # CP vs WPCA-RH standard
cp_vs_tracksterBary = []  # CP vs Barycenter-Vtx direction
cp_vs_tracksterLC_BestLCs = (
    []
)  # CP vs WPCA-LC standard, with best LCs for each HGCAL layer
cp_vs_tracksterRH_BestLCs = (
    []
)  # CP vs WPCA-RH standard, with best LCs for each HGCAL layer
cp_vs_tracksterLC_Filter_0p1 = []  # CP vs WPCA-LC standard, on LCs filtered by energy
cp_vs_tracksterLC_Filter_0p05 = []  # CP vs WPCA-LC standard, on LCs filtered by energy
caloparticles = makeCaloParticlesDF(CaloParticleCollection, save=False)
p = []
cp_eta = []
eff = 0.0
eff_best = 0.0
eff_filter = 0.0
tot = 0.0
dz_verticesLC = []
dz_verticesLC_Best = []
dz_verticesLC_Filter = []


dx_cp_LC_Best = []
dy_cp_LC_Best = []
dxy_cp_LC = []
dxy_cp_LC_Fraction = []
dxy_cp_LC_Best = []
dxLayers = []
dyLayers = []
dxyLayers = []
k = 0
LCs_E = []
LCs_x = []
LCs_y = []
LCs_z = []
INT_x = []
INT_y = []
INT_z = []
E_bad_LC = []
for evt in tqdm(range(len(TrackstersCollection))):
    tot += 1
    t = TrackstersCollection[evt]
    # vtxs = vertices[evt]
    # bs = beam_spots[evt]
    vtxs = []
    bs = []
    
    evt_id = t.event_id
    cp = caloparticles[caloparticles.loc[:, "evt_id"] == evt_id]
    if(cp.cp_pt.values[0] > 0):
        cp_eta.append(cp.cp_eta.values[0])
        # print(len(t.layerclusters[1].rechits))
        # Assign PCA on LayerClusters
        
        t.assignPCAToTracksterLC(
            3,
            weighting,
            filter_percentage=0,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=False,
        )  # weighted
        vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
        if t.valid == True:
            eff += 1
            dot = dot_product(vec, t.eigenvectors[0])
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterLC.append(dot)

            angle_yz_planeLC.append(findYZAngle(vec, t.eigenvectors[0]))
            angle_xy_planeLC.append(findXYAngle(vec, t.eigenvectors[0]))
            angle_xz_planeLC.append(findXZAngle(vec, t.eigenvectors[0]))

            profile_lc_pca_cp_pt.Fill(cp.cp_pt.values[0], dot)
            profile_lc_pca_cp_eta.Fill(cp.cp_eta.values[0], dot)
            eigen_value2.append(t.eigenvalues[1] / t.eigenvalues[0])
            eigen_value3.append(t.eigenvalues[2] / t.eigenvalues[0])
            # lpd = line_to_point_distance(t.eigenvectors[0], t.barycenter, bs)
            # dzs = []
            # for j in range(len(vtxs)):
            #     dz_vtx = abs(lpd[-1] - vtxs[j][-1])
            #     dzs.append(dz_vtx)
            # dz_verticesLC.append(min(dzs))
            #############################################################

        # Assign PCA on RecHits
        t.assignPCAToTracksterRH(
            3,
            weighting,
            filter_percentage=0,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=False,
        )
        if t.valid == True:
            dot = dot_product(vec, t.eigenvectors[0])
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterRH.append(dot)
            profile_rh_pca_cp_pt.Fill(cp.cp_pt.values[0], dot)
            # angle_yz_planeRH.append(findYZAngle(t.eigenvectors[0], vec))
            # angle_xy_planeRH.append(findXYAngle(t.eigenvectors[0], vec))
            # angle_xz_planeRH.append(findXZAngle(t.eigenvectors[0], vec))
            #############################################################
        
        # Assign PCA on LayerClusters
        t.assignPCAToTracksterLC(
            3,
            weighting,
            filter_percentage=0,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=True,
        )  # weighted
        vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
        lcs_x = []
        lcs_y = []
        lcs_z = []
        lcs_E = []
        int_x = []
        int_y = []
        int_z = []
        if t.valid == True:
            eff_best += 1
            dot = dot_product(vec, t.eigenvectors[0])
            # print(f"Eigenvector {t.eigenvectors[0]}")
            # print(f"Barycenter {t.barycenter}")
            # print(f"CaloParticle {cp.cp_x.values[0],cp.cp_y.values[0],cp.cp_z.values[0]}")
            layers_best = t.bestLCS

            for l in layers_best:
                # print(f"LC Energy {l.E}")
                index = list_of_layerz_sorted.index(abs(l.z))
                xyzIntersection = intersectionPlane(t.eigenvectors[0], t.barycenter, l.z) #intersection @ lc z plane
                # dxy_layer = pointDistance(xyzIntersection[:2], [l.x, l.y])
                dx_layer =  xyzIntersection[0] - l.x
                dy_layer =  xyzIntersection[1] - l.y
                dxy_layer = m.sqrt(dx_layer**2 + dy_layer**2)
                if(dxy_layer > 2):
                    E_bad_LC.append(l.E)

                lcs_x.append(l.x)
                lcs_y.append(l.y)
                lcs_z.append(l.z)
                lcs_E.append(l.E)
                int_x.append(xyzIntersection[0])
                int_y.append(xyzIntersection[1])
                int_z.append(xyzIntersection[2])

                dxLayers.append(dx_layer)
                dyLayers.append(dy_layer)
                dxyLayers.append(dxy_layer)

                list_of_layer_list[index].append(dxy_layer)
                list_of_layer_list_dx[index].append(dx_layer)
                list_of_layer_list_dy[index].append(dy_layer)

            LCs_x.append(lcs_x)
            LCs_y.append(lcs_y)
            LCs_z.append(lcs_z)
            LCs_E.append(lcs_E)
            INT_x.append(int_x)
            INT_y.append(int_y)
            INT_z.append(int_z)

            xyzIntersectionCP = intersectionPlane(t.eigenvectors[0], t.barycenter, cp.cp_z.values[0])
            dxCP = xyzIntersectionCP[0] - cp.cp_x.values[0]
            dyCP = xyzIntersectionCP[1] - cp.cp_y.values[0]
            # dxyCP = pointDistance(xyzIntersectionCP[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
            dxyCP = m.sqrt(dxCP**2 + dyCP**2)
            dx_cp_LC_Best.append(dxCP)
            dy_cp_LC_Best.append(dyCP)
            dxy_cp_LC_Best.append(dxyCP)
            # print(f"intersectionCP {xyzIntersectionCP}")
                # lcs_layers_dxy.Fill(abs(l.z), dxy_layer)
            # xyzIntersection = intersectionPlane(t.eigenvectors[0], t.barycenter, cp.cp_z.values[0])
            # dXYCP = pointDistance(xyzIntersection[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
            # dxy_cp_LC_Best.append(dXYCP)
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterLC_BestLCs.append(dot)
            profile_bestlc_pca_cp_pt.Fill(cp.cp_pt.values[0], dot)
            profile_bestlc_pca_cp_eta.Fill(cp.cp_eta.values[0], dot)
            eigen_value2_bestLC.append(t.eigenvalues[1])
            eigen_value3_bestLC.append(t.eigenvalues[2])


            #DXY for each layer cluster position
            # lcs_best = t.bestLCS
            # lcs_z = []
            # for l in lcs_best:
            #     lcs_z.append(l.z)
            # print(lcs_z)
            # lpd = line_to_point_distance(t.eigenvectors[0], t.barycenter, bs)
            # dzs = []
            # for j in range(len(vtxs)):
            #     dz_vtx = abs(lpd[-1] - vtxs[j][-1])
            #     dzs.append(dz_vtx)
            # dz_verticesLC_Best.append(min(dzs))
            #############################################################
        
        # Assign PCA on RecHits
        t.assignPCAToTracksterRH(
            3,
            weighting,
            filter_percentage=0,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=True,
        )
        if t.valid == True:
            eff_filter += 1
            dot = dot_product(vec, t.eigenvectors[0])
            xyzIntersection = intersectionPlane(t.eigenvectors[0], t.barycenter, cp.cp_z.values[0])
            dXYCP = pointDistance(xyzIntersection[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
            dxy_cp_LC_Best.append(dXYCP)
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterRH_BestLCs.append(dot)
            profile_bestrh_pca_cp_pt.Fill(cp.cp_pt.values[0], dot)
            #############################################################

        t.assignPCAToTracksterLC(
            3,
            weighting,
            filter_percentage=0.1,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=False,
        )  # weighted
        vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
        if t.valid == True:
            dot = dot_product(vec, t.eigenvectors[0])
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterLC_Filter_0p1.append(dot)

            #############################################################

        t.assignPCAToTracksterLC(
            3,
            weighting,
            filter_percentage=0.05,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=False,
        )  # weighted
        vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
        if t.valid == True:
            dot = dot_product(vec, t.eigenvectors[0])
            xyzIntersection = intersectionPlane(t.eigenvectors[0], t.barycenter, cp.cp_z.values[0])
            dXYCP = pointDistance(xyzIntersection[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
            dxy_cp_LC_Fraction.append(dXYCP)
            # print(xyzIntersection,(cp.cp_x.values[0], cp.cp_y.values[0], cp.cp_z.values[0]))
            if(dXYCP > 50):
                print(f"dXY {dXYCP}, intersection {xyzIntersection}, cp_pos {(cp.cp_x.values[0], cp.cp_y.values[0], cp.cp_z.values[0])}")
                print(t.eigenvectors[0])
            # print(xyzIntersection)
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterLC_Filter_0p05.append(dot)
            profile_bestlc_pca_filter_cp_pt.Fill(cp.cp_pt.values[0], dot)
            profile_bestlc_pca_filter_cp_eta.Fill(cp.cp_eta.values[0], dot)
            eigen_value2_filterLC.append(t.eigenvalues[1])
            eigen_value3_filterLC.append(t.eigenvalues[2])

            # lpd = line_to_point_distance(t.eigenvectors[0], t.barycenter, bs)
            # dzs = []
            # for j in range(len(vtxs)):
            #     dz_vtx = abs(lpd[-1] - vtxs[j][-1])
            #     dzs.append(dz_vtx)
            # dz_verticesLC_Filter.append(min(dzs))
            #############################################################

        # Assign PCA on RecHits
        t.assignPCAToTracksterRH(
            3,
            weighting,
            filter_percentage=0.05,
            scaling=scaling,
            nonlinearweight=nonlinearweighting,
            bestLCS=False,
        )
        # print(t.bestLCSFilter)
        # print("NEW")
        # for l in (t.bestLCSFilter):
        #     print(l.E, t.dataframe.tot_energy.values[0])
        if t.valid == True:
            eff_filter += 1
            dot = dot_product(vec, t.eigenvectors[0])
            if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                dot = -dot
            else:
                dot = dot
            cp_vs_tracksterRH_filter.append(dot)
            profile_bestrh_pca_filter_cp_pt.Fill(cp.cp_pt.values[0], dot)
            #############################################################

            # Direction Barycenter-Vtx WRT CaloParticle
            barycenter_wrt_vertex = t.barycenter - vertex_positions[evt]
            dot_cp_tracksterBary = dot_product(barycenter_wrt_vertex, vec)
            cp_vs_tracksterBary.append(m.acos(dot_cp_tracksterBary))
            
        #############################################################

cp_vs_tracksterLC_Filter_0p05_np = np.array(cp_vs_tracksterLC_Filter_0p05)
cp_vs_tracksterRH_np = np.array(cp_vs_tracksterRH)
np.save(outputPath + "cp_fraction005.npy", cp_vs_tracksterLC_Filter_0p05_np )
np.save(outputPath + "cp_fraction005_rh.npy", cp_vs_tracksterRH_np )


mean_dxy_layer = []
std_dxy_layer = []

mean_dx_layer = []
std_dx_layer = []

mean_dy_layer = []
std_dy_layer = []

for d in list_of_layer_list:
    mean_dxy_layer.append(np.array(d).mean())
    std_dxy_layer.append(np.array(d).std())

for d in list_of_layer_list_dx:
    mean_dx_layer.append(np.array(d).mean())
    std_dx_layer.append(np.array(d).std())

for d in list_of_layer_list_dy:
    mean_dy_layer.append(np.array(d).mean())
    std_dy_layer.append(np.array(d).std())

# fig = plt.figure(figsize = (30,25))
# # fig  = plt.figure()
# ax = fig.add_subplot(1,1,1, projection = '3d')    
# cm = plt.cm.get_cmap('hot')
# noise_idx_l = -1
# ax.scatter(provax, provay, provaz, marker='o')
# cm = plt.cm.get_cmap('hot')
# ax.scatter(provaxLC, provayLC, provaZLC, c = provaELC, cmap = cm, marker='o')
# plt.show()


for i in range(10):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection  = '3d')
    cm = plt.cm.get_cmap('winter')
    ax.scatter(INT_x[i],INT_y[i],INT_z[i], c = 'red')
    ax.scatter(LCs_x[i],LCs_y[i],LCs_z[i], c = LCs_E[i], cmap = cm)
    plt.savefig(outputPathExample + f"example_trackster{i}.png")

plt.figure(figsize=(15,15))
plt.hist(E_bad_LC, bins = 100)
plt.xlabel("Energy [GeV]")
plt.ylabel("Entries")
plt.savefig(outputPath + "E_bad_lc.png")
# plt.show()

plt.figure(figsize = (15,15))
plt.errorbar(layers_number, mean_dxy_layer, yerr=std_dxy_layer,fmt='o')
plt.xlabel("Layer number")
plt.ylabel("dxy [cm]")
plt.savefig(outputPath + "layers_dxy_fromPCA.png")
# plt.savefig( outputPath + "pickle/" + "layers_dxy_fromPCA.png")
# plt.show()
plt.figure(figsize = (15,15))
plt.errorbar(layers_number, mean_dx_layer, yerr=std_dx_layer,fmt='o')
plt.xlabel("Layer number")
plt.ylabel("$(x_{PCA} - x_{LCs})$")
plt.savefig(outputPath + "layers_dx_fromPCA.png")
# plt.savefig( outputPath + "pickle/" + "layers_dx_fromPCA.pickle")
# plt.show()
plt.figure(figsize = (15,15))
plt.errorbar(layers_number, mean_dy_layer, yerr=std_dy_layer,fmt='o')
plt.xlabel("Layer number")
plt.ylabel("$(y_{PCA} - y_{LCs})$")
plt.savefig(outputPath + "layers_dy_fromPCA.png")
# plt.savefig( outputPath + "pickle/" + "layers_dy_fromPCA.picle")

plt.figure(figsize = (15,15))
plt.scatter(dxLayers, dyLayers)
plt.xlabel("dx [cm]")
plt.ylabel("dy [cm]")
plt.xlim(-max(dxLayers),max(dxLayers))
plt.ylim(-max(dxLayers),max(dxLayers))
# plt.show()
plt.savefig(outputPath + "scatter_layers_dx_dy_fromPCA.png")
# plt.savefig( outputPath + "pickle/" + "scatter_layers_dx_dy_fromPCA.pickle")

plt.figure(figsize = (15,15))
plt.scatter(dxLayers, dxyLayers)
plt.xlabel("dx [cm]")
plt.ylabel("dxy [cm]")
# plt.show()
plt.savefig(outputPath + "scatter_layers_dx_dxy_fromPCA.png")
# plt.savefig( outputPath + "pickle/" + "scatter_layers_dx_dxy_fromPCA.pickle")

plt.figure(figsize = (15,15))
plt.scatter(dyLayers, dxyLayers)
plt.xlabel("dy [cm]")
plt.ylabel("dxy [cm]")
# plt.show()
plt.savefig(outputPath + "scatter_layers_dy_dxy_fromPCA.png")
# plt.savefig(outputPath + "pickle/" + "scatter_layers_dy_dxy_fromPCA.pickle")

plotHisto(
        dxy_cp_LC_Best,
        "dxy CaloParticle and PCA-Direction",
        "dxy [cm]",
        "Entries",
        "PCA-dxy_CP-PCA_Best",
        range=(0,5),
        outputdir=outputPath,
    )

fitHisto(
        dx_cp_LC_Best,
        "dx CaloParticle and PCA-Direction",
        "dx [cm]",
        "Entries",
        "PCA-dx_CP-PCA_Best",
        bins = 150,
        range=(-2,2),
        outputdir=outputPath,
    )

fitHisto(
        dy_cp_LC_Best,
        "dy CaloParticle and PCA-Direction",
        "dy [cm]",
        "Entries",
        "PCA-dy_CP-PCA_Best",
        bins = 150,
        range=(-2,2),
        outputdir=outputPath,
    )


# c2 = ROOT.TCanvas("c", "c")
# leg = ROOT.TLegend(0.7, 0.7, 1, 0.9)
# lcs_layers_dxy.SetTitle(
#     "dxy layercluster"
# )
# # profile_bestlc_pca_filter_cp_pt.SetLineColor(ROOT.kAzure)
# # profile_bestlc_pca_cp_pt.SetLineColor(ROOT.kOrange)
# # profile_lc_pca_cp_pt.SetLineColor(ROOT.kGreen)
# lcs_layers_dxy.SetStats(0000)
# lcs_layers_dxy.SetLineWidth(3)
# lcs_layers_dxy.Draw()
# c2.SaveAs(outputPath + "dxy_layers.png", "png")
# a = input()


# print(f"Efficiency Normal {eff/tot} , BestLCS  {eff_best/tot}, Filter {eff_filter/tot}")
# # plt.hist(dz_verticesLC, label = 'standard',histtype='step')
# # plt.hist(dz_verticesLC_Best, label = 'best',histtype='step')
# # plt.hist(dz_verticesLC_Filter, label = 'filter', histtype='step')
# # # plt.hist(cp_eta, bins = 50)
# # plt.legend()
# # plt.show()

# plt.hist(dxy_cp_LC_Fraction, bins = 100, range = (0, 5))
# plt.xlabel("dxy")
# plt.ylabel("Entries")
# plt.title("dxy distance from CaloParticle")
# plt.show()

if(plotflag == True):
    ranges = (0.99, 1)
    c = ROOT.TCanvas("c", "c", 1000, 1000)
    profile_lc_pca_cp_pt.GetYaxis().SetRangeUser(0.995,1.005)
    profile_bestlc_pca_cp_pt.GetYaxis().SetRangeUser(0.995,1.005)
    profile_bestlc_pca_filter_cp_pt.GetYaxis().SetRangeUser(0.995,1.005)
    # Plot vs CP pt
    singlePlot(
        c,
        profile_lc_pca_cp_pt,
        "Angle between PCA-LC direction wrt CP direction VS CP momentum",
        "pT [GeV]",
        "Cosine",
        "Profile_PCA-CP_pT",
        output_dir=outputPath,
    )
    singlePlot(
        c,
        profile_bestlc_pca_cp_pt,
        "Angle between PCA-LC direction wrt CP direction VS CP momentum",
        "pT [GeV]",
        "Cosine",
        "Profile_BestPCA-CP_pT",
        output_dir=outputPath,
    )
    singlePlot(
        c,
        profile_bestlc_pca_filter_cp_pt,
        f"Angle between PCA-Filter Fraction {filter_energy} direction wrt CP direction VS CP momentum",
        "pT [GeV]",
        "Cosine",
        "Profile_PCA-CP_Filter_pT",
        output_dir=outputPath,
    )
    c2 = ROOT.TCanvas("c", "c")
    leg = ROOT.TLegend(0.7, 0.7, 1, 0.9)
    profile_bestlc_pca_filter_cp_pt.SetTitle(
        "Angle between PCA-LC direction wrt CP direction vs CP momentum"
    )
    profile_bestlc_pca_filter_cp_pt.SetLineColor(ROOT.kAzure)
    profile_bestlc_pca_cp_pt.SetLineColor(ROOT.kOrange)
    profile_lc_pca_cp_pt.SetLineColor(ROOT.kGreen)
    profile_bestlc_pca_filter_cp_pt.SetStats(0000)
    profile_lc_pca_cp_pt.SetStats(0000)
    profile_bestlc_pca_cp_pt.SetStats(0000)
    profile_bestlc_pca_filter_cp_pt.GetYaxis().SetRangeUser(0.9999,1)
    profile_bestlc_pca_filter_cp_pt.SetLineWidth(3)
    profile_lc_pca_cp_pt.SetLineWidth(3)
    profile_bestlc_pca_cp_pt.SetLineWidth(3)
    leg.AddEntry(profile_bestlc_pca_filter_cp_pt, "PCA-LC-Fraction 0.05")
    leg.AddEntry(profile_bestlc_pca_cp_pt, "PCA-LC-Best")
    leg.AddEntry(profile_lc_pca_cp_pt, "PCA-LC")
    profile_bestlc_pca_filter_cp_pt.Draw()
    profile_bestlc_pca_cp_pt.Draw("SAME")
    profile_lc_pca_cp_pt.Draw("SAME")
    leg.Draw("SAME")
    c2.SaveAs(outputPath + "Profile_multi_pT.png", "png")
    a = input()
    c22 = ROOT.TCanvas("c", "c")
    leg2 = ROOT.TLegend(0.7, 0.7, 1, 0.9)
    profile_bestlc_pca_filter_cp_pt.SetLineColor(ROOT.kViolet)
    profile_bestlc_pca_cp_pt.SetLineColor(ROOT.kRed)
    profile_lc_pca_cp_pt.SetLineColor(ROOT.kGreen)
    profile_bestrh_pca_cp_pt.SetLineColor(ROOT.kOrange)
    profile_rh_pca_cp_pt.SetLineColor(ROOT.kAzure)
    # profile_bestrh_pca_cp_pt.SetLineColor(ROOT.kBlack)
    profile_bestlc_pca_filter_cp_pt.SetStats(0000)
    profile_lc_pca_cp_pt.SetStats(0000)
    profile_bestlc_pca_cp_pt.SetStats(0000)
    profile_rh_pca_cp_pt.SetStats(0000)
    profile_bestrh_pca_cp_pt.SetStats(0000)
    profile_bestlc_pca_filter_cp_pt.SetLineWidth(3)
    profile_lc_pca_cp_pt.SetLineWidth(3)
    profile_bestlc_pca_cp_pt.SetLineWidth(3)
    profile_rh_pca_cp_pt.SetLineWidth(3)
    profile_bestrh_pca_cp_pt.SetLineWidth(3)
    leg2.AddEntry(profile_bestlc_pca_filter_cp_pt, "PCA-LC-Fraction 0.05")
    leg2.AddEntry(profile_bestlc_pca_cp_pt, "PCA-LC-Best")
    leg2.AddEntry(profile_lc_pca_cp_pt, "PCA-LC")
    leg2.AddEntry(profile_bestrh_pca_cp_pt, "PCA-RH-Best")
    leg2.AddEntry(profile_rh_pca_cp_pt, "PCA-RH")
    profile_bestlc_pca_filter_cp_pt.SetTitle("Angle between PCA direction wrt CP direction vs CP momentum")
    profile_bestlc_pca_filter_cp_pt.GetYaxis().SetRangeUser(0.9999,1)
    profile_bestlc_pca_filter_cp_pt.Draw()
    profile_bestlc_pca_cp_pt.Draw("SAME")
    profile_lc_pca_cp_pt.Draw("SAME")
    profile_bestrh_pca_cp_pt.Draw("SAME")
    profile_rh_pca_cp_pt.Draw("SAME")
    leg2.Draw("SAME")
    c22.SaveAs(outputPath + "Profile_multi_pT_withRecHits.png", "png")
    a = input()

    c3 = ROOT.TCanvas("c", "c")
    leg23 = ROOT.TLegend(0.7, 0.7, 1, 0.9)
    profile_bestlc_pca_filter_cp_eta.SetLineColor(ROOT.kBlue)
    profile_bestlc_pca_cp_eta.SetLineColor(ROOT.kGreen)
    profile_lc_pca_cp_eta.SetLineColor(ROOT.kRed)
    profile_bestlc_pca_filter_cp_eta.SetStats(0000)
    profile_lc_pca_cp_eta.SetStats(0000)
    profile_bestlc_pca_cp_eta.SetStats(0000)
    profile_bestlc_pca_filter_cp_eta.SetLineWidth(3)
    profile_lc_pca_cp_eta.SetLineWidth(3)
    profile_bestlc_pca_cp_eta.SetLineWidth(3)
    leg23.AddEntry(profile_bestlc_pca_filter_cp_eta, "PCA-LC-Fraction 0.05")
    leg23.AddEntry(profile_bestlc_pca_cp_eta, "PCA-LC-Best")
    leg23.AddEntry(profile_lc_pca_cp_eta, "PCA-LC")
    profile_bestlc_pca_filter_cp_eta.SetTitle("Angle between PCA-LC direction wrt CP direction vs CP #eta")
    profile_bestlc_pca_filter_cp_eta.GetYaxis().SetRangeUser(0.9999,1)
    profile_bestlc_pca_filter_cp_eta.GetXaxis().SetTitle("#eta")
    profile_bestlc_pca_filter_cp_eta.Draw()
    profile_bestlc_pca_cp_eta.Draw("SAME")
    profile_lc_pca_cp_eta.Draw("SAME")
    leg23.Draw("SAME")
    c3.SaveAs(outputPath + "Profile_multi_eta.png", "png")
    a = input()
    i = input("a")

    # Plot PCA-LC and RH, standard
    plotHisto(
        cp_vs_tracksterLC,
        "Angle between PCA-LC direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-LC_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotHisto(
        dxy_cp_LC_Fraction,
        "x-y distance from CaloParticle",
        "dxy [cm]",
        "Entries",
        "dxy_CaloParticle_Fraction005",
        range=(0,5),
        outputdir=outputPath,
    )


    plotHisto(
        dxy_cp_LC_Best,
        "x-y distance from CaloParticle",
        "dxy [cm]",
        "Entries",
        "dxy_CaloParticle_Best",
        range=(0,5),
        outputdir=outputPath,
    )

    plotHisto(
        cp_vs_tracksterRH,
        "Angle between PCA-RH direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-RH_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-RH_CaloParticle",
        labels=["PCA-LC", "PCA-RH"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
    )
    ################################################3


    # Plot eigenvalues importance
    plotMultiHisto(
        [eigen_value2, eigen_value2_bestLC, eigen_value2_filterLC],
        "Importance of 2nd Eigenvalue",
        x_title="Importance",
        y_title="Entries",
        saveName="2ndEigenvalue_PCA_PCA-Best_PCA-Fraction",
        labels=["PCA-LC", "PCA-BestLCs", "PCA-Fraction > 0.05"],
        bins=100,
        ranges=(0,1),
        outputdir=outputPath,
        density=False,
    )
    plotMultiHisto(
        [eigen_value3, eigen_value3_bestLC, eigen_value3_filterLC],
        "Importance of 3nd Eigenvalue",
        x_title="Importance",
        y_title="Entries",
        saveName="3rdEigenvalue_PCA_PCA-Best_PCA-Filter",
        labels=["PCA-LC", "PCA-BestLCs", "PCA-Fraction > 0.05"],
        bins=100,
        ranges=(0,1),
        density=False,
        outputdir=outputPath,
    )

    # Plot PCA-LC and PCA-RH standard, Best LayerCluster on each hgcal layer
    plotHisto(
        cp_vs_tracksterLC_BestLCs,
        "Angle between PCA-LC direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-LC_BestLCs_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotHisto(
        cp_vs_tracksterRH_BestLCs,
        "Angle between PCA-RH direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-RH_BestLCs_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_BestLCs, cp_vs_tracksterRH_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-RH-BestLCs_CaloParticle",
        labels=["PCA-LC-Best", "PCA-RH-Best"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )

    # Combine PCA-LC and PCA-LC-Best
    plotMultiHisto(
        [cp_vs_tracksterLC, cp_vs_tracksterLC_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-LC_BestLCs_CaloParticle",
        labels=["PCA-LC", "PCA-LC-Best"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )
    # Combine PCA-RH and PCA-RH-Best
    plotMultiHisto(
        [cp_vs_tracksterRH, cp_vs_tracksterRH_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_CaloParticle",
        labels=["PCA-RH", "PCA-RH-Best"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )


    plotMultiHisto(
        [
            cp_vs_tracksterRH,
            cp_vs_tracksterRH_BestLCs,
            cp_vs_tracksterLC,
            cp_vs_tracksterLC_BestLCs,
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_-LC-LC-BestLCs_CaloParticle",
        labels=["PCA-RH", "PCA-RH-Best", "PCA-LC", "PCA-LC-Best"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )

    plotMultiHisto(
        [
            cp_vs_tracksterRH,
            cp_vs_tracksterRH_BestLCs,
            cp_vs_tracksterLC,
            cp_vs_tracksterLC_BestLCs,
            cp_vs_tracksterLC_Filter_0p05
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_-LC-LC-BestLCs-LC_Fraction0p05_CaloParticle",
        labels=["PCA-RH", "PCA-RH-Best", "PCA-LC", "PCA-LC-Best", "PCA-LC-Fraction 0.05"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )

    # PCA-LC filtered by energy
    number_dec1 = str(filter_energy - int(0.1))[2:]
    number_dec2 = str(filter_energy - int(0.05))[2:]
    plotHisto(
        cp_vs_tracksterLC_Filter_0p1,
        f"Angle between PCA-LC direction and CaloParticle direction - Energy fraction > {0.1}",
        "Cosine",
        "Entries",
        f"PCA-LC_Filter0p{number_dec1}_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_Filter_0p1, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilter0p{number_dec1}-LC_BestLCs-LC_CaloParticle",
        labels=[f"PCA-LC-Filter Fraction > {0.1}", "PCA-LC-Best", "PCA-LC"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
    )

    # PCA-LC filtered by energy
    # number_dec = str(filter_energy-int(0.05))[2:]
    plotHisto(
        cp_vs_tracksterLC_Filter_0p05,
        f"Angle between PCA-LC direction and CaloParticle direction - Energy fraction > {0.05}",
        "Cosine",
        "Entries",
        f"PCA-LC_Filter0p{number_dec2}_CaloParticle",
        range=ranges,
        outputdir=outputPath,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_Filter_0p05, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilter0p{number_dec2}-LC_BestLCs-LC_CaloParticle",
        labels=[f"PCA-LC-Fraction > {0.05}", "PCA-LC-Best", "PCA-LC"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )

    plotMultiHisto(
        [cp_vs_tracksterRH_filter, cp_vs_tracksterRH_BestLCs, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-RHFilter0p{number_dec2}-RH_BestLCs-RH_CaloParticle",
        labels=[f"PCA-RH-Fraction > {0.05}", "PCA-RH-Best", "PCA-RH"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )

    plotMultiHistoDifferent(
        [cp_vs_tracksterLC_Filter_0p05, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC, cp_vs_tracksterRH_filter, cp_vs_tracksterRH_BestLCs, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        diff_h = 3,
        saveName=f"MultiHisto_Angle_PCA-LCFilter0p{number_dec2}-LC_BestLCs-LC_CaloParticle_With_Rechit",
        labels=[f"PCA-LC-Fraction > {0.05}", "PCA-LC-Best", "PCA-LC",f"PCA-RH-Fraction > {0.05}", "PCA-RH-Best", "PCA-RH"],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
        density=False,
    )
    plotMultiHisto(
        [
            cp_vs_tracksterLC_Filter_0p05,
            cp_vs_tracksterLC_Filter_0p1,
            cp_vs_tracksterLC_BestLCs,
            cp_vs_tracksterLC,
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilters-LC_BestLCs-LC_CaloParticle",
        labels=[
            f"PCA-LC-Fraction > {0.05}",
            f"PCA-LC-Fraction > {0.1}",
            "PCA-LC-Best",
            "PCA-LC",
        ],
        bins=100,
        ranges=ranges,
        outputdir=outputPath,
    )


    ### ZOOM

    ranges = (0.999, 1)

    # Plot PCA-LC and RH, standard
    plotHisto(
        cp_vs_tracksterLC,
        "Angle between PCA-LC direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-LC_CaloParticle_zoom",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotHisto(
        cp_vs_tracksterRH,
        "Angle between PCA-RH direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-RH_CaloParticle_zoom",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-RH_CaloParticle_zoom",
        labels=["PCA-LC", "PCA-RH"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
    )
    ################################################3

    # Plot PCA-LC and PCA-RH standard, Best LayerCluster on each hgcal layer
    plotHisto(
        cp_vs_tracksterLC_BestLCs,
        "Angle between PCA-LC direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-LC_BestLCs_CaloParticle_zoom",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotHisto(
        cp_vs_tracksterRH_BestLCs,
        "Angle between PCA-RH direction and CaloParticle direction",
        "Cosine",
        "Entries",
        "PCA-RH_BestLCs_CaloParticle_zoom",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_BestLCs, cp_vs_tracksterRH_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-RH-BestLCs_CaloParticle_zoom",
        labels=["PCA-LC-Best", "PCA-RH-Best"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    # Combine PCA-LC and PCA-LC-Best
    plotMultiHisto(
        [cp_vs_tracksterLC, cp_vs_tracksterLC_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LC-LC_BestLCs_CaloParticle_zoom",
        labels=["PCA-LC", "PCA-LC-Best"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )
    # Combine PCA-RH and PCA-RH-Best
    plotMultiHisto(
        [cp_vs_tracksterRH, cp_vs_tracksterRH_BestLCs],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_CaloParticle_zoom",
        labels=["PCA-RH", "PCA-RH-Best"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    # PCA-LC filtered by energy
    number_dec = str(filter_energy - int(0.1))[2:]
    plotHisto(
        cp_vs_tracksterLC_Filter_0p1,
        f"Angle between PCA-LC direction and CaloParticle direction - Energy fraction > {0.1}",
        "Cosine",
        "Entries",
        f"PCA-LC_Filter0p{number_dec}_CaloParticle",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_Filter_0p1, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilter0p{number_dec1}-LC_BestLCs-LC_CaloParticle",
        labels=[f"PCA-LC-Fraction > {0.01}", "PCA-LC-Best", "PCA-LC"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
    )

    # PCA-LC filtered by energy

    plotHisto(
        cp_vs_tracksterLC_Filter_0p05,
        f"Angle between PCA-LC direction and CaloParticle direction - Energy fraction > {0.05}",
        "Cosine",
        "Entries",
        f"PCA-LC_Filter0p{number_dec2}_CaloParticle",
        range=ranges,
        outputdir=outputZoomSubDir,
    )

    plotMultiHisto(
        [cp_vs_tracksterLC_Filter_0p05, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilter{number_dec2}-LC_BestLCs-LC_CaloParticle",
        labels=[f"PCA-LC-Fraction > {0.05}", "PCA-LC-Best", "PCA-LC"],
        bins=60,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    plotMultiHisto(
        [
            cp_vs_tracksterLC_Filter_0p05,
            cp_vs_tracksterLC_Filter_0p1,
            cp_vs_tracksterLC_BestLCs,
            cp_vs_tracksterLC,
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-LCFilter-LC_BestLCs-LC_CaloParticle_Filters_Zoom",
        labels=[
            f"PCA-LC-Fraction > {0.05}",
            f"PCA-LC-Fraction > {0.1}",
            "PCA-LC-Best",
            "PCA-LC",
        ],
        bins=100,
        ranges=ranges,
        outputdir=outputZoomSubDir,
    )


    plotMultiHisto(
        [
            cp_vs_tracksterRH,
            cp_vs_tracksterRH_BestLCs,
            cp_vs_tracksterLC,
            cp_vs_tracksterLC_BestLCs,
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_-LC-LC-BestLCs_CaloParticle_zoom",
        labels=["PCA-RH", "PCA-RH-Best", "PCA-LC", "PCA-LC-Best"],
        bins=100,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    plotMultiHisto(
        [
            cp_vs_tracksterRH,
            cp_vs_tracksterRH_BestLCs,
            cp_vs_tracksterLC,
            cp_vs_tracksterLC_BestLCs,
            cp_vs_tracksterLC_Filter_0p05
        ],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName="MultiHisto_Angle_PCA-RH-RH_BestLCs_-LC-LC-BestLCs-LC_Fraction0p05_CaloParticle_zoom",
        labels=["PCA-RH", "PCA-RH-Best", "PCA-LC", "PCA-LC-Best", "PCA-LC-Fraction 0.05"],
        bins=100,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    plotMultiHistoDifferent(
        [cp_vs_tracksterLC_Filter_0p05, cp_vs_tracksterLC_BestLCs, cp_vs_tracksterLC, cp_vs_tracksterRH_filter, cp_vs_tracksterRH_BestLCs, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-LCFilter0p{number_dec2}-LC_BestLCs-LC_CaloParticle_With_Rechit_Zoom",
        diff_h = 3,
        labels=[f"PCA-LC-Fraction > {0.05}", "PCA-LC-Best", "PCA-LC",f"PCA-RH-Fraction > {0.05}", "PCA-RH-Best", "PCA-RH"],
        bins=100,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    plotMultiHisto(
        [cp_vs_tracksterRH_filter, cp_vs_tracksterRH_BestLCs, cp_vs_tracksterRH],
        "Angle between PCA direction and CaloParticle direction",
        x_title="Cosine",
        y_title="Entries",
        saveName=f"MultiHisto_Angle_PCA-RHFilter0p{number_dec2}-RH_BestLCs-RH_CaloParticle_Zoom",
        labels=[f"PCA-RH-Fraction > {0.05}", "PCA-RH-Best", "PCA-RH"],
        bins=100,
        ranges=ranges,
        outputdir=outputZoomSubDir,
        density=False,
    )

    plotHisto(
        dxy_cp_LC_Fraction,
        "x-y distance from CaloParticle",
        "dxy [cm]",
        "Entries",
        "dxy_CaloParticle_Fraction005",
        range=(0,3),
        outputdir=outputZoomSubDir,
    )

    plotHisto(
        dxy_cp_LC_Best,
        "x-y distance from CaloParticle",
        "dxy [cm]",
        "Entries",
        "dxy_CaloParticle_Best",
        range=(0,3),
        outputdir=outputPath,
    )
