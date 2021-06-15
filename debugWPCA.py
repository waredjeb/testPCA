#!/usr/bin/env python
#coding: utf-8

#In[ ]:





#In[1]:


#from plotsUtils import *
#dataDir = "./data/"
#rootFileName = "output_SingleGamma.root"
#fileRoot = ROOT.TFile(dataDir + rootFileName)
#tree = fileRoot.Get("trackstersAnalyzer/tree")
#getBranches(tree)
#N = tree.GetEntries()

#for i in range(1):
    #tree.GetEntry(i)
    #vertices = []
    #print(v_x.size())
    #for j in range(v_x.size()):
        #vertex = [0,0,0]
        #vertex[0] = v_x.at(j)
        #vertex[1] = v_y.at(j)
        #vertex[2] = v_z.at(j
        #vertices.append(vertex)

#vtxs = np.array(vertices).T    
#fig = plt.figure(figsize = (10,10))
#ax = fig.add_subplot(2,2,1, projection = '3d')   
#ax.scatter(vtxs[2], vtxs[1], vtxs[0], marker='o', alpha=1)
#plt.show()


#In[2]:


#for evt in range(N):
    #print(f"Evt {evt}")
    #tree.GetEntry(evt)
    #print(N_trackster.at(0))
    #for j in range(N_lcs.size()):
        #print(N_lcs.at(j))


#In[3]:


#vertex_positions = []
#CaloParticleCollection = []
#TrackstersCollection = []
#for evt in range(N):
    #tree.GetEntry(evt)

    #print(f"\nEvent number: {evt}")

    #Loop over number of trackster
    #t_id = 0
    #for i in range(1):  # loop over tracksters
        #t = Trackster(evt, t_id, [])
        #tmp_hit_idx = 0
        #if N_lcs.size() > 0:
            #for j in range(N_lcs.at(0)):  # loop over layerclusters
                #x = lc_x.at(j)
                #y = lc_y.at(j)
                #z = lc_z.at(j)
                #E = lc_E.at(j)
                #lc = LayerCluster(x, y, z, E, t_id,[])
                #for k in range(N_hit_per_lc.at(j)):  # loop over hit
                    #be carefull with the index, is a single array containing the rechits for all the layer clusters
                    #new_k = tmp_hit_idx + k
                    #hit_x = rh_x.at(new_k)
                    #hit_y = rh_y.at(new_k)
                    #hit_z = rh_z.at(new_k)
                    #hit_E = rh_E.at(new_k)
                    #hit_fraction = rh_fraction.at(new_k)
                    #hit = RecHit(hit_x, hit_y, hit_z, hit_E, hit_fraction)
                    #lc.addRechit(hit)
                #tmp_hit_idx += N_hit_per_lc.at(j)
                #t.addLayerCluster(lc)
            #print(cp_px.at(i))
            #cp = CaloParticle(
                #cp_px.at(i),
                #cp_py.at(i),
                #cp_pz.at(i),
                #cp_pt.at(i),
                #cp_E.at(i),
                #cp_eta.at(i),
                #cp_phi.at(i),
                #evt,
                #i,
            #)
            #CaloParticleCollection.append(cp)
            #TrackstersCollection.append(t)
            #vertex_pos = [v_x.at(0), v_y.at(0), v_z.at(0)]
            #vertex_positions.append(vertex_pos)
        #t_id += 1

#print(len(CaloParticleCollection), len(TrackstersCollection))


#In[ ]:





#In[4]:


#def findIntersectionYplane(v1, barycenter, yplane):
    #k = (yplane-barycenter[0])/(v1[0])
    #x = v1[0] * k + barycenter[0]
    #y = v1[1] * k +  barycenter[1]
    #z = v1[2] * k +  barycenter[2]
    #return [x,y,z]


#In[5]:


#caloparticles = makeCaloParticlesDF(CaloParticleCollection, save=False)
#cp_vs_tracksterLC = []
#dz = []
#for evt in tqdm(range(len(TrackstersCollection))):
    #vertex = vertex_positions[evt]
    #t = TrackstersCollection[evt]
    #evt_id = t.event_id
    #cp = caloparticles[caloparticles.loc[:, "evt_id"] == evt_id]
    #Assign PCA on LayerClusters
    #t.assignPCAToTracksterLC(
        #3,
        #True,
        #filter_percentage=0,
        #scaling=False,
        #nonlinearweight=False,
        #bestLCS=False,
    #)  # weighted
    #vec = cp.loc[:, ["cp_px", "cp_py", "cp_pz"]].values[0]
    #if t.valid == True:
        #dot = vec.dot(t.eigenvectors[0])
        #if m.acos(dot) > m.pi / 2:
            #angle = m.pi - m.acos(dot)
        #else:
            #angle = m.acos(dot)
        #cp_vs_tracksterLC.append(angle)
    #intersection  = findIntersectionYplane(t.eigenvectors[0],t.barycenter, 0)
    #diff = abs(intersection[-1] - vertex[-1])
    #dz.append(diff)
    #print(f"intersection {findIntersectionYplane(t.eigenvectors[0],t.barycenter, 0)}")
    #print(f"Vertex {vertex}")
    #print(f"Eigen {t.eigenvectors[0]}")
    
    #print(t.barycenter)
    
#eigen = ( -7.28842209,  98.08794081, 335.69000961) + t*(0.00513576, 0.29493745, 0.95550276)


#In[6]:


#t = TrackstersCollection[0]
#t.assignPCAToTracksterLC(
    #3,
    #True,
    #filter_percentage=0,
    #scaling=False,
    #nonlinearweight=False,
    #bestLCS=False,
#)  # weighted


#In[7]:


#df = t.dataframe


#In[8]:


#df


#In[9]:

#b = t.barycenter
#%matplotlib notebook
#plots(df.iloc[:,1:5].values, [b[2],b[0],b[1]], t.eigenvectors,[400,400,400], noise = 0, save = "PROVA", arrow_scale = 6)

#plt.show()
#In[ ]:




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import mplhep as hep
from pcaUtils import *
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import warnings
warnings.filterwarnings('ignore')

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.rcParams.label.data = True
hep.rcParams.label.paper = False

def plots3DwithProjection(data, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (20,25))
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(data[:,0], data[:,1], data[:,2], marker='o', alpha=.3)
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)
    ax = fig.add_subplot(2,2,2)
    ax.scatter(data[:,0], data[:,1])#, c = data[:,3], alpha = .5, cmap=cm)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    ax.scatter(data[:,0], data[:,2])##,c = data[:,3], alpha = .5)#,cmap=cm)

    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)
    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(data[:,1], data[:,2])#, c = data[:,3], alpha = .5)#,cmap=cm)

    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)

    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label = 'Rhits energy')
    cbar.set_label("Rechits energy", fontsize = 20)
    
def plots2(data,barycenter, components, length, noise = 0, save = "", arrow_scale = 10):
    '''
    Function to plots the eigenvector obtained by the PCA
    - data: - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    - barycenter: barycenter obtained by the PCA
    - components: eigenvectors obtained by the PCA
    - length = eigenvalues obtained by the PCA
    - noise: if noise > 0 takes the last <noise> rows and plot them with an X instead of a circle (only if you add noise)
    - save = plot outputfile name
    - arrow_scale = scaling the length of the arrows in the plots
    TODO: Remove noise, implement a more general way to scale the arrows
    '''
    #quiver plots
    fig = plt.figure(figsize = (20,25))
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('hot')
    noise_idx_l = -noise -1 
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,1], data[:noise_idx_l,2])#, marker='o', alpha=.3, c = data[:noise_idx_l, 3], cmap = cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,1], data[noise_idx:,2])#, marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    col = ['r','g','b']
    colors = col
    x0 = barycenter[0]
    y0 = barycenter[1]
    z0 = barycenter[2]
    scale = 1
    # new_length = [i/length[0]*scale for i in length]
    new_length = length
    vec1 = components[0]*(new_length[0]) * scale
    vec2 = components[1]*(new_length[1]) * scale
    vec3 = components[2]*(new_length[2]) * scale
    ax.quiver(x0,y0,z0,vec1[0],vec1[1],vec1[2], color = col[0])
    ax.quiver(x0,y0,z0,vec2[0],vec2[1],vec2[2], color = col[1])
    ax.quiver(x0,y0,z0,vec3[0],vec3[1],vec3[2], color = col[2])

    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)
    
    ax = fig.add_subplot(2,2,2)
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,1])#, c = data[:noise_idx_l,3], alpha = .5, cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,1])#, marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(x0,y0,vec1[0], vec1[1], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,y0,vec2[0], vec2[1], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,y0,vec3[0], vec3[1], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)

    
    ax = fig.add_subplot(2,2,3)
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,2])#,c = data[:noise_idx_l,3], alpha = .5,cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,2])#, marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(x0,z0,vec1[0], vec1[2], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,z0,vec2[0], vec2[2], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,z0,vec3[0], vec3[2], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)
    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(data[:noise_idx_l,1], data[:noise_idx_l,2])#, c = data[:noise_idx_l,3], alpha = .5,cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,1], data[noise_idx:,2])#, marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(y0,z0,vec1[1], vec1[2], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(y0,z0,vec2[1], vec2[2], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(y0,z0,vec3[1], vec3[2], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")
    plt.xlim(-15, 15)
    plt.ylim(-15, 15)

    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs energy')
    cbar.set_label("Lcs energy", fontsize = 20)
    # fig.subplots_adjust(wspace=0.05, hspace=0.4, right=0.8)
    
angles = []
for i in range(10000):
    x = np.random.normal(0, .01, 30)
    y = np.random.normal(0, .01, 30)
    z = np.random.normal(0, 8, 30)
    d = {"x":x, "y":y, "z":z}
    df = pd.DataFrame(d)
    #if(i == 0):
        #plots3DwithProjection(df.values)
        #plt.show()
    pca, barycenter, valid = runPCA(3,df.values)
    #if(i == 0):
        #plots2(df.values, barycenter,pca.components_, pca.explained_variance_, noise = 0, save = "", arrow_scale = 10)
    dot = dot_product([0,0,1],pca.components_[0])
    angle_sign = dot_productSign([0,0,1],pca.components_[0])
    if m.acos(dot) > m.pi / 2:
        angle = m.pi - m.acos(dot)
    else:
        angle = m.acos(dot)
    angles.append(angle_sign)


plt.hist(angles, bins = 500, range = (-0.005, 0.005))

plt.savefig("./tmpPlot/angles.png")
