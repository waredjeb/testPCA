from matplotlib.pyplot import bar, title
from numpy.lib.function_base import angle
from classes import CaloParticle, LayerCluster, RecHit, Trackster
from readTree import *
from tqdm import tqdm
from plotsUtils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--plot", type=bool, help="make plots")
parser.add_argument("--debug", type=bool, help="debugging mode")
args = parser.parse_args()
plotflag = bool(args.plot)
debug = bool(args.debug)

plt.style.use([hep.style.ROOT, hep.style.firamath])

hep.rcParams.label.data = True
hep.rcParams.label.paper = False


# Loading Tree


dataDir = "./data/"
# rootFileName = "output_singlephoton_1k_vertexing_cpNew.root"

file_names = ["output_singlephoton_1k_vertexing_cpNew.root"]#,"output_singlephoton_1k_vertexing_cpNew_eta1p7.root","output_singlephoton_1k_vertexing_cpNew_eta2p7.root","output_singlephoton_PU140_300.root"]
# file_names = ["output_singlephoton_PU140_300.root"]
etas = [0]
pileups = [0]
# etas = [0, 1.7, 2.7, 0]
# pileups = [0,0,0,140]

for i in range(len(file_names)):
    fileRoot = ROOT.TFile(dataDir + file_names[i])
    print(f"Loading Tree {dataDir + file_names[i]}\n")
    tree = fileRoot.Get("trackstersAnalyzer/tree")
    pileup = pileups[i]
    eta_val = etas[i]
    outputDir = "./plotsAnalysisSingleNew/"
    baseDir = "PCAWeightedVertexing_DEBUGVertexSmallAngles"
    eta_s = str(eta_val).replace(".", "p")
    pu_s = str(pileup)
    baseDir += eta_s + "_"
    baseDir += pu_s + "_"
    baseDir += "3HitsCuFractionCut0p005/"
    outputSubDir = baseDir
    # outputSubDir = "PCAWeightedDZVertex/"
    # outputSubDir = "PCAWeightedVertexing_3RhitsCutFractionCut2GeV/"
    outputPath = outputDir + outputSubDir
    outputZoomSubDir = outputPath + "Zoom/"
    outputPathExample = outputPath + "Example/"
    mkdir_p(outputPath)
    mkdir_p(outputPath + "pickle")
    mkdir_p(outputZoomSubDir)
    mkdir_p(outputZoomSubDir + "pickle")
    mkdir_p(outputPathExample)

    getBranches(tree)

    if debug:
        Nevts = 200
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

            t = Trackster(evt, ind, [])
            ind+=1
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
                
                for k in range(vtxs_x.size()):
                    vtx =  [vtxs_x.at(k),vtxs_y.at(k),vtxs_z.at(k)]
                    vtxs.append(vtx)
                vertices.append(vtxs)
                # for k in range(bs_x.size()):
                bs = [bs_x.at(i),bs_y.at(i),bs_z.at(i)]
                beam_spots.append(bs)
                

    lcsdict = {"lcx": lx, "lcy": ly, "lcz": lz}
    lcsdf = pd.DataFrame(lcsdict)
    lcsgroupz = lcsdf.groupby(["lcz"])
    list_of_layerz = []
    for key, item in lcsgroupz:
        list_of_layerz.append(float(key))


    print(
        len(CaloParticleCollection),
        len(TrackstersCollection),
        len(beam_spots),
        len(vertices),
    )

    list_of_layerz_sorted = sorted(list_of_layerz)
    print(f"LAYERZ {list_of_layerz_sorted[10:25]}")
    layers_number = [i + 1 for i in range(len(list_of_layerz_sorted))]
    list_of_layer_list = [[] for _ in range(len(list_of_layerz_sorted))]
    list_of_layer_list_dx = [[] for _ in range(len(list_of_layerz_sorted))]
    list_of_layer_list_dy = [[] for _ in range(len(list_of_layerz_sorted))]
    # Process all the data and make all the plots
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

    eigen_value2_bestLC = []
    eigen_value3_bestLC = []


    profile_bestlc_pca_cp_pt = ROOT.TProfile("Angles2", "Angles2", 30, 0, 150, 0.5, 1.5)


    profile_bestlc_pca_cp_eta = ROOT.TProfile(
        "Angles22", "Angles22", 30, 1.5, 2.8, 0.5, 1.5
    )

    profile_dz_vertex_eta = ROOT.TProfile(
        "Angles222", "Angles222", 10, 1.7, 2.7, -30, 30
    )


    lcs_layers_dxy = ROOT.TProfile("Angles31", "Angles3", 1000, 310, 380, -10, 10)
    cp_vs_tracksterLC_BestLCs = (
        []
    )  # CP vs WPCA-LC standard, with best LCs for each HGCAL layer
    cp_vs_tracksterLC_BestLCs_super = (
        []
    )  # CP vs WPCA-LC standard, with best LCs for each HGCAL layer

    caloparticles = makeCaloParticlesDF(CaloParticleCollection, save=False)
    p = []
    cp_ETA = []
    eff = 0.0
    eff_best = 0.0
    eff_filter = 0.0
    tot = 0.0
    dz_verticesLC_Best = []


    dx_cp_LC_Best = []
    dy_cp_LC_Best = []
    dxy_cp_LC_Best = []
    E_map_dxy_dz = []
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
    CP_x = []
    CP_y = []
    CP_z = []
    V_x = []
    V_y = []
    V_z = []

    TO_CP_x = []
    TO_CP_y = []
    TO_CP_z = []
    FROM_CP_x = []
    FROM_CP_y = []
    FROM_CP_z = []
    E_bad_LC = []
    dz_vertex = []
    lcs_energy = []
    tot_energies = []

    dz_wanted = []
    dxy_wanted = []
    angle_wanted = []
    dxy_caloparticle_w0 = []
    dxy_caloparticle_w0_std = []
    # w0s = np.arange(4, 6, 0.5)
    w0s = [5]
    bad_dz_vert = []
    bad_dot_vert = []
    print(w0s)
    for w0 in w0s:
        for evt in tqdm(range(len(TrackstersCollection))):
            tot += 1
            t = TrackstersCollection[evt]
            vtxs = vertices[evt]
            bs = beam_spots[evt]
            v = vtxs[0]
            # vtxs = []
            # bs = []
            # print(t)
            evt_id = t.event_id
            cp = caloparticles[caloparticles.loc[:, "evt_id"] == evt_id]
            if cp.cp_pt.values[0] > 0:
                cp_ETA.append(cp.cp_eta.values[0])
                t.assignPCAToTracksterLC(
                    3,
                    True,
                    filter_percentage=0,
                    scaling=scaling,
                    nonlinearweight=True,
                    bestLCS=True,
                    startingLC = 0,
                    endLC = -1,
                    w0 = w0
                )  # weighted
                vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
                vec_pos = cp.loc[:, ["cp_x", "cp_y", "cp_z"]].values[0]

                if t.valid == True:
                    # print(direction)
                    # print(t.barycenter)
                    df = t.dataframe
                    layers_best = t.bestLCS
                    layers_bestdf = t.bestLCSdf
                    max_index = layers_bestdf.loc[:,'lc_E'].idxmax()
                    highELC = layers_bestdf.loc[max_index, :].values[1:4]

                    # t.eigenvectors[0] = vec
                    # t.barycenter = vec_pos
                    # t.barycenter = highELC
                    # t.barycenter = v
                    direction = t.eigenvectors[0]
                    print(f"Direction {direction}")
                    print(f"Eigenvectors {t.eigenvectors[0]}")
                    print(f"Barycenter {t.barycenter}")

                    tot_energy = df.loc[0,'tot_energy']                    
                    eff_best += 1
                    dot = dot_product(vec, direction)               
                    if(abs(dot) > 0.999999):     
                        lcs_x = []
                        lcs_y = []
                        lcs_z = []
                        lcs_E = []
                        int_x = []
                        int_y = []
                        int_z = []
                        k_vert = 0
                        for l in layers_best:
                            # print(f"LC Energy {l.E}")
                            lcs_energy.append(l.E)
                            tot_energies.append(tot_energy)
                            index = list_of_layerz_sorted.index(abs(l.z))
                            xyzIntersection = intersectionPlane(
                                direction, t.barycenter, l.z
                            )  # intersection @ lc z plane
                            # dxy_layer = pointDistance(xyzIntersection[:2], [l.x, l.y])
                            dx_layer = xyzIntersection[0] - l.x
                            dy_layer = xyzIntersection[1] - l.y
                            dxy_layer = m.sqrt(dx_layer ** 2 + dy_layer ** 2)
                            if dxy_layer > 0.4:
                                E_bad_LC.append(l.E/tot_energy)

                            dxLayers.append(dx_layer)
                            dyLayers.append(dy_layer)
                            dxyLayers.append(dxy_layer)

                            list_of_layer_list[index].append(dxy_layer)
                            list_of_layer_list_dx[index].append(dx_layer)
                            list_of_layer_list_dy[index].append(dy_layer)

                            xyzIntersectionCP_bad = intersectionPlane(
                                direction, t.barycenter, cp.cp_z.values[0]
                            )
                            dxCP_bad = xyzIntersectionCP_bad[0] - cp.cp_x.values[0]
                            dyCP_bad = xyzIntersectionCP_bad[1] - cp.cp_y.values[0]
                            # dxyCP = pointDistance(xyzIntersectionCP[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
                            dxyCP_bad = m.sqrt(dxCP_bad ** 2 + dyCP_bad ** 2)

                            cpa_bad = findCPA(direction, t.barycenter, bs)
                            cpa_z_bad= cpa_bad[0][-1]
                            dz_v_bad = abs(cpa_z_bad - v[-1])
                            if(dz_v_bad > 0.0001):
                                lcs_x.append(l.x)
                                lcs_y.append(l.y)
                                lcs_z.append(l.z)
                                lcs_E.append(l.E)
                                int_x.append(xyzIntersection[0])
                                int_y.append(xyzIntersection[1])
                                int_z.append(xyzIntersection[2])
                                xyzIntersection = intersectionPlane(
                                    direction, t.barycenter, l.z
                                )  # intersection @ lc z plane
                                # dxy_layer = pointDistance(xyzIntersection[:2], [l.x, l.y])
                                dx_layer = xyzIntersection[0] - l.x
                                dy_layer = xyzIntersection[1] - l.y
                                dxy_layer = m.sqrt(dx_layer ** 2 + dy_layer ** 2)
                                # if(k_vert == 0):
                                #     print(f"Dz {dz_v_bad}")
                                # print(dxy_layer)        
                                

                                if(k_vert == 0):
                                    bad_dz_vert.append(dz_v_bad)
                                    bad_dot_vert.append(abs(dot))
                                k_vert+=1
                                    
                                if(l.z != xyzIntersection[2]):
                                    print(f"LC_Z {l.z} int_z {xyzIntersection[2]}")
                        to_CP_x = []
                        to_CP_y = []
                        to_CP_z = []
                        from_CP_x = []
                        from_CP_y = []
                        from_CP_z = []
                        if(len(lcs_x) > 0):
                            LCs_x.append(lcs_x)
                            LCs_y.append(lcs_y)
                            LCs_z.append(lcs_z)
                            LCs_E.append(lcs_E)
                            INT_x.append(int_x)
                            INT_y.append(int_y)
                            INT_z.append(int_z)
                            CP_x.append(cp.cp_x.values[0])
                            CP_y.append(cp.cp_y.values[0])
                            CP_z.append(cp.cp_z.values[0])
                            V_x.append(v[0])
                            V_y.append(v[1])
                            V_z.append(v[2])
                            if(lcs_z[0] < 0):
                                to_cp = np.linspace(lcs_z[0], 20, 100)
                            else:
                                to_cp = np.linspace(-20, lcs_z[0], 100)
                            print(to_cp)
                            for to_z in to_cp:
                                xyzIntersection_toCP = intersectionPlane(
                                    direction, t.barycenter, to_z
                                )
                                xyzIntersection_fromCP = intersectionPlane(
                                    vec, vec_pos, to_z
                                )
                                to_CP_x.append(xyzIntersection_toCP[0])
                                to_CP_y.append(xyzIntersection_toCP[1])
                                to_CP_z.append(xyzIntersection_toCP[2])
                                from_CP_x.append(xyzIntersection_fromCP[0])
                                from_CP_y.append(xyzIntersection_fromCP[1])
                                from_CP_z.append(xyzIntersection_fromCP[2])
                            TO_CP_x.append(to_CP_x)
                            TO_CP_y.append(to_CP_y)
                            TO_CP_z.append(to_CP_z)
                            FROM_CP_x.append(from_CP_x)
                            FROM_CP_y.append(from_CP_y)
                            FROM_CP_z.append(from_CP_z)

                        

                        xyzIntersectionCP = intersectionPlane(
                            direction, t.barycenter, cp.cp_z.values[0]
                        )
                        dxCP = xyzIntersectionCP[0] - cp.cp_x.values[0]
                        dyCP = xyzIntersectionCP[1] - cp.cp_y.values[0]
                        # dxyCP = pointDistance(xyzIntersectionCP[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
                        dxyCP = m.sqrt(dxCP ** 2 + dyCP ** 2)
                        dx_cp_LC_Best.append(dxCP)
                        dy_cp_LC_Best.append(dyCP)
                        dxy_cp_LC_Best.append(dxyCP)
                        E_map_dxy_dz.append(cp.cp_E.values[0])
                        # print(f"intersectionCP {xyzIntersectionCP}")
                        # lcs_layers_dxy.Fill(abs(l.z), dxy_layer)
                        # xyzIntersection = intersectionPlane(direction, t.barycenter, cp.cp_z.values[0])
                        # dXYCP = pointDistance(xyzIntersection[:2], [cp.cp_x.values[0], cp.cp_y.values[0]])
                        # dxy_cp_LC_Best.append(dXYCP)
                        if (t.barycenter[2] * t.eigenvectors[2][2]) < 0:
                            dot = -dot
                        else:
                            dot = dot
                        cp_vs_tracksterLC_BestLCs.append(dot)
                        
                        cp_vs_tracksterLC_BestLCs_super.append(dot)
                        profile_bestlc_pca_cp_pt.Fill(cp.cp_pt.values[0], dot)
                        profile_bestlc_pca_cp_eta.Fill(cp.cp_eta.values[0], dot)
                        # eigen_value2_bestLC.append(t.eigenvalues[1])
                        # # eigen_value3_bestLC.append(t.eigenvalues[2])
                        # print(direction)
                        # print(t.barycenter)
                    
                        cpa = findCPA(direction, t.barycenter, bs)
                        # print(f"CPA {cpa}")
                        # print(f"cp {[cp.cp_x.values[0], cp.cp_y.values[0], cp.cp_z.values[0]]}")
                        cpa_z = cpa[0][-1]
                        dz_v = abs(cpa_z - v[-1])
                        dz_vertex.append(cpa_z - v[-1])
                        # print(f"dz {dz_v}")
                        profile_dz_vertex_eta.Fill(cpa_z - v[-1], cp.cp_eta.values[0])
                        if(dz_v <= 0.1):
                            dz_wanted.append(dz_v)
                            dxy_wanted.append(dxyCP)
                            angle_w = m.acos(abs(dot_product(direction, vec)))
                            angle_wanted.append(180*angle_w/m.pi)
                            # print(f"DOT {dot_product(direction, vec)}")
        dxy_cp_mean = np.array(dxy_cp_LC_Best).mean()
        dxy_cp_std = np.array(dxy_cp_LC_Best).std()
        dxy_caloparticle_w0.append(dxy_cp_mean)
        dxy_caloparticle_w0_std.append(dxy_cp_std)
                    # if(dz_v < 1):
                    #     print(f"cpa {cpa}, vertex {v}")

    for i in range(4):
        
        fig = plt.figure()
        ax = fig.add_subplot(2, 2, 1, projection="3d")
        cm = plt.cm.get_cmap("winter")
        ax.scatter(INT_x[i], INT_y[i], INT_z[i], c="red")
        ax.scatter(LCs_x[i], LCs_y[i], LCs_z[i], c=LCs_E[i], cmap=cm)
        ax.plot(CP_x[i], CP_y[i], CP_z[i], marker = '+')
        ax.scatter(TO_CP_x[i], TO_CP_y[i], TO_CP_z[i])
        ax.scatter(FROM_CP_x[i], FROM_CP_y[i], FROM_CP_z[i])
        ax.scatter(V_x[i], V_y[i], V_z[i], marker = "+")
        ax.scatter(bs[0], bs[1], bs[2], marker = "+")
        
        # print(bad_dz_vert[i])
        # ax.set_zlim(0,10)
        s_v = str(bad_dot_vert[i])
        ax.set_title("dot cp-pca" + s_v) 
        ax = fig.add_subplot(2, 2, 2)
        ax.scatter(INT_x[i], INT_y[i], c="red")
        ax.scatter(LCs_x[i], LCs_y[i], c=LCs_E[i], cmap=cm)
        ax.scatter(CP_x[i], CP_y[i], marker = "+", label = 'CP')
        ax.plot(TO_CP_x[i], TO_CP_y[i], label = 'PCA')
        ax.plot(FROM_CP_x[i], FROM_CP_y[i], label = 'CPDIR')
        ax.scatter(V_x[i], V_y[i],  marker = "+",  label = 'vertex')
        ax.scatter(bs[0], bs[1], marker = "+",  label = 'beam spot')
        ax.set_xlabel('$X$')
        ax.set_ylabel('$Y$')
        ax.set_title("X-Y")
        plt.legend(loc = 'upper left')

        ax = fig.add_subplot(2, 2, 3)
        ax.scatter(INT_x[i], INT_z[i], c="red")
        ax.scatter(LCs_x[i], LCs_z[i], c=LCs_E[i], cmap=cm)
        ax.scatter(CP_x[i], CP_z[i], marker = "+")
        ax.plot(TO_CP_x[i], TO_CP_z[i])
        ax.plot(FROM_CP_x[i], FROM_CP_z[i])
        ax.scatter(V_x[i], V_z[i],  marker = "+")
        ax.scatter(bs[0], bs[2], marker = "+")
        ax.set_xlabel('$X$')
        ax.set_ylabel('$Z$')
        ax.set_title("X-Z")
        ax = fig.add_subplot(2, 2, 4)
        ax.scatter(INT_y[i], INT_z[i], c="red")
        ax.scatter(LCs_y[i], LCs_z[i], c=LCs_E[i], cmap=cm)
        ax.scatter(CP_y[i], CP_z[i], marker = "+")
        ax.plot(TO_CP_y[i], TO_CP_z[i])
        ax.plot(FROM_CP_y[i], FROM_CP_z[i])
        ax.scatter(V_y[i], V_z[i],  marker = "+")
        ax.scatter(bs[1], bs[2], marker = "+")
        ax.set_xlabel('$Y$')
        ax.set_ylabel('$Z$')
        ax.set_title("Y-Z")
        # plt.show()
        plt.savefig(outputPathExample + f"example_trackster{i}.png")

   
    if plotflag == True:

        mean_dxy_layer = []
        std_dxy_layer = []

        mean_dx_layer = []
        std_dx_layer = []

        mean_dy_layer = []
        std_dy_layer = []

        plt.figure(figsize = (15,15))
        plt.scatter(dxy_wanted, dz_wanted)
        plt.xlabel("dxy [cm]")
        plt.ylabel("dz [cm]")
        plt.title("Requiring dz <= 1mm")
        plt.savefig(outputPath + "wanted_dz_dxy.png")

        plt.figure(figsize = (15,15))
        plt.hist(dxy_wanted, bins = 100)
        plt.xlabel("dxy [cm]")
        plt.ylabel("Entries")
        plt.title("Requiring dz <= 1mm")
        plt.savefig(outputPath + "wanted_dxy_hist.png") 

        plt.figure(figsize = (15,15))
        plt.hist(angle_wanted)
        plt.xlabel("Cosine")
        plt.ylabel("Entries ")
        plt.title("Requiring dz <= 1mm")
        plt.savefig(outputPath + "wanted_angle_cp.png")

        plt.figure(figsize = (15,15))
        plt.scatter(dz_wanted, angle_wanted, marker = 'o')
        plt.xlabel("dz [cm]")
        plt.ylabel("Angle [$\circ$]")
        plt.title("Requiring dz <= 1mm")
        plt.savefig(outputPath + "wanted_angle_cp_dz.png")


        plt.figure(figsize = (15,15))
        plt.hist(lcs_energy, bins = 100)
        plt.xlabel("Energies [GeV]")
        plt.savefig(outputPath + "Energies_LCs.png")
        # plt.show()

        # lc_E_log = [(lcs_energy[i]/tot_energies[i]) + (lcs_energy[i]/tot_energies[i])**2 + (lcs_energy[i]/tot_energies[i])**3 for i in range(len(lcs_energy/tot_energies[i]))]
        # lc_E_log = [m.exp(100*lcs_energy[i]/tot_energies[i]) for i in range(len(lcs_energy))]
        lc_E_log = [lcs_energy[i]/tot_energies[i]  for i in range(len(lcs_energy)) ]
        # lc_E_log = [m.log(lcs_energy[i]/tot_energies[i]) + 5 for i in range(len(lcs_energy))]
        # lc_E_log = [(lcs_energy[i]/tot_energies[i])**3 for i in range(len(lcs_energy))]
        # plt.figure(figsize = (15,15))
        plt.hist(lc_E_log, bins = 100)
        # plt.hist(lc_E_log_w0, bins = 100)
        plt.xlabel("Weight")
        plt.savefig(outputPath + "LogWeighting_Energy.png")
        # plt.show()

        plt.figure(figsize = (15,15))
        x_axis = [lcs_energy[i]/tot_energies[i] for i in range(len(lcs_energy))]
        plt.scatter(x_axis, lc_E_log)
        # plt.scatter(lcs_energy, lc_E_log_w0)
        plt.xlabel("$LC_{E}/Trackster_{E}$")
        plt.ylabel("Weights")
        plt.savefig(outputPath + "LogWeightingVSEnergy.png")
        # plt.show()

        # print(f"Wos {len(w0s)} dxy {len(dxy_cp_mean)}")

        plt.figure(figsize = (15,15))
        plt.errorbar(w0s, dxy_caloparticle_w0, yerr=dxy_caloparticle_w0_std, fmt = 'o')
        plt.xlabel("w0")
        plt.ylabel("dxy CP")
        plt.title("dxy CaloParticle - PCA-Direction vs W0")
        plt.savefig(outputPath + "dxy_CP_WO.png")


        plt.figure(figsize = (15,15))
        plt.hist(dz_vertex, bins = 100)
        plt.xlabel("dz [cm]")
        plt.ylabel("Entries")
        plt.title("dz PCA-Direction and Primary Vertex")
        plt.savefig(outputPath + "dz_primary_vertex.png")
        plt.close()

        fig = plt.figure(figsize =(15,15))
        cm = plt.cm.get_cmap("winter")
        im = plt.scatter(dxy_cp_LC_Best, dz_vertex, c = E_map_dxy_dz, cmap = cm)
        plt.xlim(0, 1.5)
        plt.ylim(-30,30)
        plt.xlabel("dxy CP [cm]")
        plt.ylabel("dz vertex [cm]")
        plt.title("dz Vertex - dxy CaloParticle")
        cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax, label = 'CP Energy')
        cbar.set_label("CP Energy", fontsize = 30)
        plt.savefig(outputPath + "scatter_dxy_dz_vertex.png",  bbox_inches='tight')

        range_d = (-5,5)
        fitHisto(
            dz_vertex,
            "dz PCA-Direction and Primary Vertex",
            "dz [cm]",
            "Entries",
            "dz_Primary_Vertex_Fit",
            bins=100,
            range=range_d,
            PU = pileup,
            eta = eta_val,
            outputdir=outputPath,
        )
        range_d = (-1,1)
        


        for d in list_of_layer_list:
            mean_dxy_layer.append(np.array(d).mean())
            std_dxy_layer.append(np.array(d).std())

        for d in list_of_layer_list_dx:
            mean_dx_layer.append(np.array(d).mean())
            std_dx_layer.append(np.array(d).std())

        for d in list_of_layer_list_dy:
            mean_dy_layer.append(np.array(d).mean())
            std_dy_layer.append(np.array(d).std())



        
        # print(f"std X {std_dx_layer}")
        # print(f"std Y {std_dy_layer}")

        plt.figure(figsize=(15, 15))
        plt.hist(E_bad_LC, bins=100)
        plt.xlabel("Energy-Fraction")
        plt.ylabel("Entries")
        plt.savefig(outputPath + "E_bad_lc.png")
        # plt.show()

        plt.figure(figsize=(15, 15))
        plt.errorbar(layers_number, mean_dxy_layer, yerr=std_dxy_layer, fmt="o")
        plt.xlabel("Layer number")
        plt.ylabel("dxy [cm]")
        plt.title("dxy PCA-Direction position @ LC-z - LCs Position")
        plt.savefig(outputPath + "layers_dxy_fromPCA.png")
        # plt.savefig( outputPath + "pickle/" + "layers_dxy_fromPCA.png")
        # plt.show()
        plt.figure(figsize=(15, 15))
        plt.errorbar(layers_number, mean_dx_layer, yerr=std_dx_layer, fmt="o")
        plt.xlabel("Layer number")
        plt.ylabel("$(x_{PCA} - x_{LCs})$")
        plt.title("dx PCA-Direction position @ LC-z - LCs Position")
        plt.savefig(outputPath + "layers_dx_fromPCA.png")
        # plt.savefig( outputPath + "pickle/" + "layers_dx_fromPCA.pickle")
        # plt.show()
        plt.figure(figsize=(15, 15))
        plt.errorbar(layers_number, mean_dy_layer, yerr=std_dy_layer, fmt="o")
        plt.xlabel("Layer number")
        plt.ylabel("$(y_{PCA} - y_{LCs})$")
        plt.title("dy PCA-Direction position @ LC-z - LCs Position")
        plt.savefig(outputPath + "layers_dy_fromPCA.png")
        # plt.savefig( outputPath + "pickle/" + "layers_dy_fromPCA.picle")

        plt.figure(figsize=(15, 15))
        plt.scatter(dxLayers, dyLayers)
        plt.xlabel("dx [cm]")
        plt.ylabel("dy [cm]")
        plt.title("dy vs dx PCA-Direction position @ LC-z - LCs Position")
        plt.xlim(-max(dxLayers), max(dxLayers))
        plt.ylim(-max(dxLayers), max(dxLayers))
        # plt.show()
        plt.savefig(outputPath + "scatter_layers_dx_dy_fromPCA.png")
        # plt.savefig( outputPath + "pickle/" + "scatter_layers_dx_dy_fromPCA.pickle")

        plt.figure(figsize=(15, 15))
        plt.scatter(dxLayers, dxyLayers)
        plt.xlabel("dx [cm]")
        plt.ylabel("dxy [cm]")
        plt.title("dxy vs dx PCA-Direction position @ LC-z - LCs Position")
        # plt.show()
        plt.savefig(outputPath + "scatter_layers_dx_dxy_fromPCA.png")
        # plt.savefig( outputPath + "pickle/" + "scatter_layers_dx_dxy_fromPCA.pickle")

        plt.figure(figsize=(15, 15))
        plt.scatter(dyLayers, dxyLayers)
        plt.xlabel("dy [cm]")
        plt.ylabel("dxy [cm]")
        plt.title("dxy vs dy PCA-Direction position @ LC-z - LCs Position")
        # plt.show()
        plt.savefig(outputPath + "scatter_layers_dy_dxy_fromPCA.png")
        # plt.savefig(outputPath + "pickle/" + "scatter_layers_dy_dxy_fromPCA.pickle")



        plotHisto(
            dxy_cp_LC_Best,
            "dxy CaloParticle and PCA-Direction",
            "dxy [cm]",
            "Entries",
            "PCA-dxy_CP-PCA_Best",
            range=(0, 5),
            outputdir=outputPath,
        )

        fitHisto(
            dx_cp_LC_Best,
            "dx CaloParticle and PCA-Direction",
            "dx [cm]",
            "Entries",
            "PCA-dx_CP-PCA_Best",
            bins=150,
            range=range_d,
            PU = pileup,
            eta = eta_val,
            outputdir=outputPath,
        )

        fitHisto(
            dy_cp_LC_Best,
            "dy CaloParticle and PCA-Direction",
            "dy [cm]",
            "Entries",
            "PCA-dy_CP-PCA_Best",
            bins=150,
            range=range_d,
            PU = pileup,
            eta = eta_val,
            outputdir=outputPath,
        )
        c2 = ROOT.TCanvas("c", "c")
        leg = ROOT.TLegend(0.7, 0.7, 1, 0.9)
        profile_dz_vertex_eta.SetTitle(
            "dz PCA-Direction - Primary Vertex"
        )
        profile_dz_vertex_eta.SetStats(0000)
        profile_dz_vertex_eta.SetLineWidth(3)
        profile_dz_vertex_eta.Draw()
        c2.SaveAs(outputPath + "profile_dzPV_eta.png", "png")
        # a = input()

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


        ranges = (0.99, 1)
        c = ROOT.TCanvas("c", "c", 1000, 1000)
        profile_bestlc_pca_cp_pt.GetYaxis().SetRangeUser(0.995, 1.005)

        singlePlot(
            c,
            profile_bestlc_pca_cp_pt,
            "Angle between PCA-LC direction wrt CP direction VS CP momentum",
            "pT [GeV]",
            "Cosine",
            "Profile_BestPCA-CP_pT",
            output_dir=outputPath,
        )

        c2 = ROOT.TCanvas("c", "c")
        leg = ROOT.TLegend(0.7, 0.7, 1, 0.9)
        profile_bestlc_pca_cp_pt.SetTitle(
            "Angle between PCA-LC direction wrt CP direction vs CP momentum"
        )

        profile_bestlc_pca_cp_pt.SetLineColor(ROOT.kBlue)
        profile_bestlc_pca_cp_pt.SetStats(0000)
        profile_bestlc_pca_cp_pt.GetYaxis().SetRangeUser(0.9999, 1)

        profile_bestlc_pca_cp_pt.SetLineWidth(3)

        leg.AddEntry(profile_bestlc_pca_cp_pt, "PCA-LC-Best")

        profile_bestlc_pca_cp_pt.Draw("")

        leg.Draw("SAME")
        c2.SaveAs(outputPath + "Profile_multi_pT.png", "png")
        # a = input()

        c3 = ROOT.TCanvas("c", "c")
        leg23 = ROOT.TLegend(0.7, 0.7, 1, 0.9)
        profile_bestlc_pca_cp_eta.SetLineColor(ROOT.kBlue)
        profile_bestlc_pca_cp_eta.SetStats(0000)

        profile_bestlc_pca_cp_eta.SetLineWidth(3)
        leg23.AddEntry(profile_bestlc_pca_cp_eta, "PCA-LC-Best")
        profile_bestlc_pca_cp_eta.SetTitle(
            "Angle between PCA-LC direction wrt CP direction vs CP #eta"
        )
        profile_bestlc_pca_cp_eta.GetYaxis().SetRangeUser(0.9999, 1)
        profile_bestlc_pca_cp_eta.GetXaxis().SetTitle("#eta")
        profile_bestlc_pca_cp_eta.Draw()
        leg23.Draw("SAME")
        c3.SaveAs(outputPath + "Profile_multi_eta.png", "png")
        # a = input()
        # i = input("a")
        ranges = (0.9995, 1)
        plotHisto(
            dxy_cp_LC_Best,
            "x-y distance from CaloParticle",
            "dxy [cm]",
            "Entries",
            "dxy_CaloParticle_Best",
            range=(0, 5),
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
            cp_vs_tracksterLC_BestLCs_super,
            "Angle between PCA-LC direction and CaloParticle direction, cosine > 0.99996",
            "Cosine",
            "Entries",
            "PCA-LC_BestLCs_CaloParticle_super",
            range=ranges,
            outputdir=outputPath,
        )
        ranges = (0.99990, 1)
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
            cp_vs_tracksterLC_BestLCs_super,
            "Angle between PCA-LC direction and CaloParticle direction, cosine > 0.99996",
            "Cosine",
            "Entries",
            "PCA-LC_BestLCs_CaloParticle_super",
            range=(0.99995, 1),
            outputdir=outputZoomSubDir,
        )

        plotHisto(
            dxy_cp_LC_Best,
            "x-y distance from CaloParticle",
            "dxy [cm]",
            "Entries",
            "dxy_CaloParticle_Best",
            range=(0, 3),
            outputdir=outputZoomSubDir,
        )
