import pandas as pd
import matplotlib.pyplot as plt
import math as m
import numpy as np
from sklearn.preprocessing import StandardScaler

import random
from sklearn.decomposition import PCA, FactorAnalysis, IncrementalPCA
# Enable 3D plotting
from mpl_toolkits import mplot3d 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.proj3d import proj_transform
import seaborn as sns
import mplhep as hep
from sklearn.preprocessing import StandardScaler
import os
from random import gauss

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.rcParams.label.data = True
hep.rcParams.label.paper = False

# matplotlib.rcParams.update({'font.size': 22})



#Create new directory
def mkdir_p(mypath):
    '''Function to create a new directory, if it not already exist
        - mypath : directory path
    '''
    from errno import EEXIST
    from os import makedirs,path
    try:
        makedirs(mypath)
    except OSError as exc: 
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

def makeCaloParticlesDF(CaloParticlesCollection, output_name = "cp_collection_df", output_dir = './postData/', save = False):
    list_of_df = []
    for i in range(len(CaloParticlesCollection)):
        cp = CaloParticlesCollection[i]
        evt_id = cp.event_id
        trackster_id = cp.trackster_id
        d = {"evt_id": [evt_id] * 1,
            "trackster_id": [trackster_id] * 1,
            "cp_px": cp.px,
            "cp_py": cp.py,
            "cp_pz": cp.pz,
            "cp_eta": cp.eta,
            "cp_phi": cp.phi,
            "cp_x": cp.x,
            "cp_y": cp.y,
            "cp_z": cp.z,
            "cp_pt": cp.pt,
            "cp_E": cp.E
            }
        df = pd.DataFrame(d)
        list_of_df.append(df)
    df_merged = pd.concat(list_of_df)
    if(save == True):
        df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged


def makeCaloParticleDFSingle(cp, output_name = "cp_df", output_dir = './postData/', save = False):
    evt_id = cp.event_id
    trackster_id = cp.trackster_id
    d = {"evt_id": [evt_id] * 1,
        "trackster_id": [trackster_id] * 1,
        "cp_px": cp.px,
        "cp_py": cp.py,
        "cp_pz": cp.pz,
        "cp_eta": cp.eta,
        "cp_phi": cp.phi,
        "cp_pt": cp.pt,
        "cp_E": cp.E
        }
    df = pd.DataFrame(d)
    if(save == True):
        df.to_csv(output_dir + output_name + ".csv")
    return df


def makeLayerClusterDF(TracksterCollection, output_name ="lcs_df", output_dir = './postData/', save = False):
    '''Function to create a dataframe starting from a list of Tracksters. 
       The output dataframe contains all the layerclusters of each trackster:
        - evt_id: Trackster event
        - lc_x, lc_y, lc_z: Trackster's layercluster coordinates
        - lc_E: Trackster's layercluster energy
        
    - TracksterCollection: list of Tracksters
    - output_name: output file name
    - output_dir: directory in which the .csv is saved
    '''
    list_of_df = []
    for i in range(len(TracksterCollection)):
        t = TracksterCollection[i]
        evt_id = t.event_id
        lcs = t.layerclusters
        d = {"evt_id": [evt_id] * len(lcs),
            "lc_x": [lc.x for lc in lcs],
            "lc_y": [lc.y for lc in lcs],
            "lc_z": [lc.z for lc in lcs],
            "lc_E": [lc.E for lc in lcs]}
        df = pd.DataFrame(d)
        df['tot_energy'] = [df['lc_E'].sum()] * df.shape[0]
        list_of_df.append(df)
    df_merged = pd.concat(list_of_df)
    if(save == True):
        df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged

def makeLayerClusterDFSingleTrackster(t, output_name =' lcs_single_trackster_df', output_dir = './postData/', save = False):
    '''Function to create a dataframe starting from a single trackster. 
       The output dataframe contains all the layerclusters of each trackster:
        - evt_id: Trackster event
        - lc_x, lc_y, lc_z: Trackster's layercluster coordinates
        - lc_E: Trackster's layercluster energy
        
    - TracksterCollection: list of Tracksters
    - output_name: output file name
    - output_dir: directory in which the .csv is saved
    - tot_energy: Trackster's total energy
    TODO: Merge this function in makeLayerClusterDF
    '''
    list_of_df = []
    evt_id = t.event_id
    lcs = t.layerclusters
    d = {"evt_id": [evt_id] * len(lcs),
        "lc_x": [lc.x for lc in lcs],
        "lc_y": [lc.y for lc in lcs],
        "lc_z": [lc.z for lc in lcs],
        "lc_E": [lc.E for lc in lcs]}
    df = pd.DataFrame(d)
    df['tot_energy'] = [df['lc_E'].sum()] * df.shape[0]
    list_of_df.append(df)
    df_merged = pd.concat(list_of_df)
    if(save == True):
        df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged

def makeRecHitDF(TracksterCollection, output_name = "rh_df", output_dir = './postData/'):
    '''Function to create a dataframe starting from a list of Tracksters
       The output dataframe contains all the rechits of each layerclusters of each trackster:
        - evt_id: Trackster event
        - lc_id: Layercluster number
        - rh_x, rh_y, rh_z: Trackster's rechits coordinates
        - rh_E: Trackster's layercluster energy
        - tot_energy: Trackster's total energy
    - TracksterCollection: list of Tracksters
    - output_name: output file name
    - output_dir: directory in which the .csv is saved
    TODO: Merge this function in makeLayerClusterDF
    '''
    list_of_df = []
    for i in range(len(TracksterCollection)):
        t = TracksterCollection[i]
        evt_id = t.event_id
        lcs = t.layerclusters
        xs = []
        ys = []
        zs = []
        Es = []
        fractions = []
        evt_ids = []
        lc_ids = []
        for j in range(len(lcs)):
            rhits = lcs[j].rechits
            for k in range(len(rhits)):
                h = rhits[k]
                xs.append(h.x)
                ys.append(h.y)
                zs.append(h.z)
                Es.append(h.E)
                fractions.append(h.fraction)
                evt_ids.append(evt_id)
                lc_ids.append(j)


        d = {"evt_id": evt_ids,
            "lc_id" : lc_ids,
            "rh_x": xs,
            "rh_y": ys,
            "rh_z": zs,
            "rh_E": Es,
            "fraction": fractions}

        df = pd.DataFrame(d)
        df['tot_energy'] = [df['rh_E'].sum()] * df.shape[0]
        list_of_df.append(df)
    
    df_merged = pd.concat(list_of_df)
    df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged


def plots3DwithProjectionList(data, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(data[:,0], data[:,1], data[:,2], marker='o', alpha=.3, c = data[:, 3], cmap = cm)
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    ax.scatter(data[:,0], data[:,1], c = data[:,3], alpha = .5, cmap=cm)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    ax.scatter(data[:,0], data[:,2],c = data[:,3], alpha = .5,cmap=cm)

    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(data[:,1], data[:,2], c = data[:,3], alpha = .5,cmap=cm)

    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label = 'Rhits energy')
    cbar.set_label("Rechits energy", fontsize = 50)

def plots3DwithProjection(data, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(data[:,0], data[:,1], data[:,2], marker='o', alpha=.3, c = data[:, 3], cmap = cm)
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    ax.scatter(data[:,0], data[:,1], c = data[:,3], alpha = .5, cmap=cm)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    ax.scatter(data[:,0], data[:,2],c = data[:,3], alpha = .5,cmap=cm)

    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(data[:,1], data[:,2], c = data[:,3], alpha = .5,cmap=cm)

    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label = 'Rhits energy')
    cbar.set_label("Rechits energy", fontsize = 50)



def plots(data,barycenter, components, length, noise = 0, save = "", arrow_scale = 10):
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
    fig = plt.figure(figsize = (30,25))
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('hot')
    noise_idx_l = -noise -1 
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,1], data[:noise_idx_l,2], marker='o', alpha=.3, c = data[:noise_idx_l, 3], cmap = cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,1], data[noise_idx:,2], marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
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
    
    ax = fig.add_subplot(2,2,2)
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,1], c = data[:noise_idx_l,3], alpha = .5, cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,1], marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(x0,y0,vec1[0], vec1[1], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,y0,vec2[0], vec2[1], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,y0,vec3[0], vec3[1], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")

    
    ax = fig.add_subplot(2,2,3)
    ax.scatter(data[:noise_idx_l,0], data[:noise_idx_l,2],c = data[:noise_idx_l,3], alpha = .5,cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,0], data[noise_idx:,2], marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(x0,z0,vec1[0], vec1[2], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,z0,vec2[0], vec2[2], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(x0,z0,vec3[0], vec3[2], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")
    
    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(data[:noise_idx_l,1], data[:noise_idx_l,2], c = data[:noise_idx_l,3], alpha = .5,cmap=cm)
    if(noise != 0):
        noise_idx = -noise
        ax.scatter(data[noise_idx:,1], data[noise_idx:,2], marker='x', alpha=.3, c = data[noise_idx:, 3], cmap = cm)
    ax.quiver(y0,z0,vec1[1], vec1[2], color = col[0],angles='xy', scale_units='xy', scale=1)
    ax.quiver(y0,z0,vec2[1], vec2[2], color = col[1],angles='xy', scale_units='xy', scale=1)
    ax.quiver(y0,z0,vec3[1], vec3[2], color = col[2],angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs energy')
    cbar.set_label("Lcs energy", fontsize = 50)
    # fig.subplots_adjust(wspace=0.05, hspace=0.4, right=0.8)
    plt.savefig("./plots/" + save + ".png", bbox_inches='tight')

def getXYZ(df):
    return df["x"],df["y"],df["z"],df["E"]


def runPCA(n_components, data, incremental = False):
        data = data.T
        assert data.shape[0] == 3, 'Wrong shape'
        
        if(data.shape[1] < 3):            
            return False, [0,0,0], False
        else:
            # print(data.shape)
            if(incremental == False):
                pca = PCA(n_components = n_components)
                pca.fit(data.T)
            else:
                pca = IncrementalPCA(n_components = n_components)
                pca.fit(data.T)            

            barycenter = [data[0].mean(),data[1].mean(), data[2].mean()]
            return pca, barycenter, True

def runFA(n_components, data):
        data = data.T
        assert data.shape[0] == 3, 'Wrong shape'
        
        if(data.shape[1] < 3):            
            return False, [0,0,0], False
        else:
            # print(data.shape)
            # print("FA")
            fa = FactorAnalysis(n_components = n_components, rotation = 'varimax')
            fa.fit(data.T)
            

            barycenter = [data[0].mean(),data[1].mean(), data[2].mean()]
            return fa, barycenter, True


def runWPCA(n_components, xyz_T, energy, tot_energy, nonlinear = False, w0 = 0):
    xyz = xyz_T.T
    assert xyz.shape[0] == 3, "Input data shape should be 3, check your input shape"
    N = len(energy)
    barycenter = [0.,0.,0.]
    covM = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
    covMW = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
    sorted_index = [0,0,0]
    sorted_eigenvalue = [0,0,0]
    sorted_eigenvectors = [0,0,0]
    if(N > 1):
        weight_sum = 0
        weight_sum2 = 0
        # Compute the barycenter
        for i in range(N):
            if(nonlinear == True):
                constant = 1.5
                # weight = (energy[i])**1.5
                # weight = m.exp(energy[i]/1) 
                weight = energy[i]
            else:
                weight = energy[i]
            weight = energy[i] #computing the weight
            weight_sum += weight
            barycenter[0] += xyz[0][i] * weight
            barycenter[1] += xyz[1][i] * weight
            barycenter[2] += xyz[2][i] * weight
        # print(barycenter)
        barycenter /= weight_sum    
        
        #Compute the covariance matrix
        covM = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        for i in range(N):
            if(nonlinear == True):                
                # weight = m.log(energy[i]/tot_energy) + w0 #log
                # weight = energy[i]/tot_energy #linear
                weight = (energy[i]/tot_energy)
                # weight = (energy[i]/tot_energy) + (energy[i]/tot_energy)**2 + (energy[i]/tot_energy)**3 #poly3
                # weight = m.exp(100*energy[i]/tot_energy) #Exp
            else:
                weight = energy[i]/tot_energy
            if(weight > 0):
                # print(weight)
                weight_sum2 += weight * weight            
                for x in range(3):
                    for y in range(3):
                        covM[x][y] += weight * (xyz[x][i] - barycenter[x]) * (xyz[y][i] - barycenter[y])
        
        covMW = covM * 1./(1-weight_sum2) #weighted covariance matrix
        #get eigenvalues and eigenvectors and sort them
        try:
            eigen_values , eigen_vectors = np.linalg.eigh(covMW)
            norm = np.linalg.norm(eigen_values)
            vectors = eigen_values/norm
            sorted_index = np.argsort(vectors)[::-1]
            sorted_eigenvalue = vectors[sorted_index]
            sorted_eigenvectors = eigen_vectors[:,sorted_index].T
            valid = True
        except:
            print("Not Converged")
            sorted_index = [0,0,0]
            sorted_eigenvalue = [0,0,0]
            sorted_eigenvectors = [0,0,0]
            barycenter = []
            valid = False
    else:
        valid = False

    return sorted_eigenvectors, sorted_eigenvalue, covMW, covM, barycenter, valid 

def runWPCA2D(n_components, xyz_T, energy, tot_energy, nonlinear = False, w0 = 0):
    xyz = xyz_T.T
    assert xyz.shape[0] == 3, "Input data shape should be 3, check your input shape"
    N = len(energy)
    barycenter = [0.,0.,0.]
    covM = np.array([[0.,0.],[0.,0.]])
    covMW = np.array([[0.,0.],[0.,0.]])
    sorted_index = [0,0,0]
    sorted_eigenvalue = [0,0,0]
    sorted_eigenvectors = [0,0,0]
    if(N > 1):
        weight_sum = 0
        weight_sum2 = 0
        # Compute the barycenter
        for i in range(N):
            if(nonlinear == True):
                constant = 1.5
                weight = (energy[i])**1.5
            else:
                weight = energy[i]
            weight = energy[i] #computing the weight
            weight_sum += weight
            barycenter[0] += xyz[0][i] * weight
            barycenter[1] += xyz[1][i] * weight
            barycenter[2] += xyz[2][i] * weight
        # print(barycenter)
        barycenter /= weight_sum    
        
        #Compute the covariance matrix
        covM = np.array([[0.,0.],[0.,0.]])
        for i in range(N):
            if(nonlinear == True):                
                weight = m.log(energy[i]/tot_energy) + w0
            else:
                weight = energy[i]/tot_energy
            if(weight > 0):
                weight_sum2 += weight * weight            
                for x in [0,2]:
                    for y in [0,2]:
                        if(x == 2):
                            xm = 1
                        else:
                            xm = x
                        if(y == 2):
                            ym = 1
                        else:
                            ym = y
                        covM[xm][ym] += weight * (xyz[x][i] - barycenter[x]) * (xyz[y][i] - barycenter[y])
        
        covMW = covM * 1./(1-weight_sum2) #weighted covariance matrix
        #get eigenvalues and eigenvectors and sort them
        try:
            eigen_values , eigen_vectors = np.linalg.eigh(covMW)
            norm = np.linalg.norm(eigen_values)

            
            vectors = eigen_values/norm
            sorted_index = np.argsort(vectors)[::-1]
            sorted_eigenvalue = vectors[sorted_index]
            sorted_eigenvectors = eigen_vectors[:,sorted_index].T
            e_z = [sorted_eigenvectors[0][0],0,sorted_eigenvectors[0][1]]
            e_x = [sorted_eigenvectors[1][0],0,sorted_eigenvectors[1][1]]
            eig_2d = [e_z, [0,0,0], e_x]
            valid = True
        except:
            print("Not Converged")
            sorted_index = [0,0,0]
            sorted_eigenvalue = [0,0,0]
            sorted_eigenvectors = [0,0,0]
            barycenter = []
            valid = False
    else:
        valid = False

    return np.array(eig_2d), sorted_eigenvalue, covMW, covM, barycenter, valid 




def runUAPCAGauss(n_components, xyz_T, energy, tot_energy, nonlinear = False, sigmas = [1,1,1], Nsamples = 50):
    xyz =  xyz_T.T
    sigmax = sigmas[0]
    sigmay = sigmas[1]
    sigmaz = sigmas[2]

    l = []
    for i in range(len(xyz[0])):
        l_temp = [[],[],[]]
        l_temp[0] = [gauss(xyz[0][i], sigmax) for j in range(Nsamples)]
        l_temp[1] = [gauss(xyz[1][i], sigmay) for j in range(Nsamples)]
        l_temp[2] = [gauss(xyz[2][i], sigmaz) for j in range(Nsamples)]
        l.append(l_temp)
    d = {"VAR": l}
    df = pd.DataFrame(d)
    

    mean = lambda l : sum(l)/len(l)
    ut = [0,0,0]
    for i in range(df.shape[0]):
        tmp_mean = [0,0,0]
        tmp_mean[0] = mean(df.iloc[i,0][0])
        tmp_mean[1] = mean(df.iloc[i,0][1])
        tmp_mean[2] = mean(df.iloc[i,0][2])
        ut[0] += tmp_mean[0]
        ut[1] += tmp_mean[1]
        ut[2] += tmp_mean[2]
    Ut = [u/df.shape[0] for u in ut] #In this case it represents the barycenter
    barycenter = Ut #Just to avoid any confusion
    kt = [[0,0,0], [0,0,0], [0,0,0]]
    w_sum2 = 0
    for i in range(df.shape[0]):
        m = np.array([0,0,0])
        m[0] = mean(df.iloc[i,0][0])
        m[1] = mean(df.iloc[i,0][1])
        m[2] = mean(df.iloc[i,0][2])
        mmT = np.outer(m,m)
        UtUtT = np.outer(Ut, Ut)
        # weight = dataLC.iloc[i, 3]/dataLC.iloc[i, -1]
        weight = 1
        w_sum2 += weight*weight
        kt += (mmT + np.cov(np.array(df.iloc[i,0])) - UtUtT) * weight
        # KT = kt/(df.shape[0] * (1 -w_sum2))
        KT = kt/(df.shape[0])
        eigen_values , eigen_vectors = np.linalg.eigh(KT)
        sorted_index = np.argsort(eigen_values)[::-1]
        sorted_eigenvalue = eigen_values[sorted_index]
        sorted_eigenvectors = eigen_vectors[:,sorted_index].T

        return sorted_eigenvectors, sorted_eigenvalue, KT, barycenter
    return 0


def pcaSummary(pca):
    print("Covariance Matrix \n {} \n ".format(pca.get_covariance()))
    print("Eigenvectors \n {} \n".format(pca.components_))
    print("Eigenvalues \n {} \n".format(pca.explained_variance_))


def computeCovarianceMatrix(data):
    cov_mat = np.cov(data, rowvar = False)
    return cov_mat


def dot_productSign(v1,v2):
    v1_np = np.array(v1)
    v2_np = np.array(v2)
    v1_norm = v1_np.dot(v1_np)
    v2_norm = v2_np.dot(v2_np)
    dot = v1_np.dot(v2_np)
    angle = m.acos(dot)
    v3 = np.cross(v1_np,v2_np)
    dot_sign = v3.dot([0,1,0])       
    if(dot_sign < 0):
        angle = - angle
    return angle



def dot_product(v1,v2):
    v1_np = np.array(v1)
    v2_np = np.array(v2)
    v1_norm = v1_np.dot(v1_np)
    v2_norm = v2_np.dot(v2_np)
    dot = v1_np.dot(v2_np)
    # print(dot)
    return dot/(m.sqrt(v1_norm * v2_norm))

def dot_productOffset(v1,v2):
    dot = dot_product(v1,v2) - 0.5
    return dot

def findAngle(v1,v2):
    dot = dot_product(v1,v2)
    if m.acos(dot) > m.pi / 2:
        angle = m.pi - m.acos(dot)
    else:
        angle = m.acos(dot)
    return angle


def findXYAngle(v1,v2):
    v1_new = [v1[0], v1[1]]
    v2_new = [v2[0], v2[1]]
    return findAngle(v1_new, v2_new)

def findXZAngle(v1,v2):
    v1_new = [v1[0], v1[2]]
    v2_new = [v2[0], v2[2]]
    return findAngle(v1_new, v2_new)

def findYZAngle(v1,v2):
    v1_new = [v1[1], v1[2]]
    v2_new = [v2[1], v2[2]]
    return findAngle(v1_new, v2_new)    

import scipy as sci
def distance(a, b):
    """Calculate a distance between two points."""
    if(len(a)  == 2):
        return m.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
    else:
        return m.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)
    # return np.linalg.norm(a-b)

def line_to_point_distance(direction, q, r):
    """Calculate a distance between point r and a line crossing p and q."""
    def foo(t: float):
        # x is point on line, depends on t 
        x = t * direction + q
#         print(t)
        # we return a distance, which also depends on t      
        x_2d = x[:2]
        r_2d = r[:2]
        return distance(x, r)
    # which t minimizes distance?
    t0 = sci.optimize.minimize(foo, 1).x[0]
    point = t0*direction + q 
    
#     print(sci.optimize.minimize(foo, 0.1).x)
#     print(point)
#     return foo(t0)
    return point #return z

def intersectionPlane(vec, bary, z):
    k = ((z-bary[2])/vec[2])
    x = (vec[0]*k) + bary[0]
    y = (vec[1]*k) + bary[1]
    z = (vec[2]*k) + bary[2]
    return [x,y,z]

def pointDistance(v1,v2):
    v1_n = np.array(v1)
    v2_n = np.array(v2)
    return np.linalg.norm(v1_n-v2_n)


def perpendicular( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b
def intersect(x0a,y0a,dxa,dya,x0b,y0b,dxb,dyb):
    t = (dyb*(x0b-x0a)-dxb*(y0b-y0a))/(dxa*dyb-dxb*dya)
    return (x0a+dxa*t,y0a+dya*t,t)

def intersection(bary2d, vec2d, bs2d, perp2d):
    t = (perp2d[1]*(bs2d[0]-bary2d[0])-perp2d[0]*(bs2d[1]-bary2d[1]))/(vec2d[0]*perp2d[1]-perp2d[0]*vec2d[1])
    return (bary2d[0] + vec2d[0] * t, bary2d[1]*vec2d[1]*t,t)

def findCPA(vec, bary, beamspot):
    vec2d = vec[:2]
    
    perp2d = perpendicular(vec2d)
#     print(f"PERP {perp2d}")
    xint,yint,t = intersection(bary[:2], vec2d, beamspot[:2], perp2d)
    return vec * t + bary, t

def findCPAYZ(vec, bary, beamspot):
    vec2d = vec[1:]
    
    perp2d = perpendicular(vec2d)
#     print(f"PERP {perp2d}")
    xint,yint,t = intersection(bary[1:], vec2d, beamspot[1:], perp2d)
    return vec * t + bary, t

# def findCPAYZ(vec, bary, beamspot):
#     vec2d = vec[1:]