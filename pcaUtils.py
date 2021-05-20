import pandas as pd
import matplotlib.pyplot as plt
import math as m
import numpy as np

import random
from sklearn.decomposition import PCA
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

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.rcParams.label.data = True
hep.rcParams.label.paper = False



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

def makeLayerClusterDF(TracksterCollection, output_name, output_dir = './postData/'):
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
    df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged

def makeLayerClusterDFSingleTrackster(t, output_name, output_dir = './postData/'):
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
    df_merged.to_csv(output_dir + output_name + ".csv")
    return df_merged

def makeRecHitDF(TracksterCollection, output_name, output_dir = './postData/'):
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


def plots3DwithProjection(data, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (20,20))
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
    cbar.set_label("Rechits energy", fontsize = 20)



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
    fig = plt.figure(figsize = (25,20))
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
    scale = arrow_scale
    new_length = [i/length[0]*scale for i in length]
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
    cbar.set_label("Lcs energy", fontsize = 20)

    plt.savefig("./plots/" + save + ".png")

def getXYZ(df):
    return df["x"],df["y"],df["z"],df["E"]


def runPCA(n_components, data):
    data = data.T
    assert data.shape[0] == 3, 'Wrong shape'
    print(data.shape)
    pca = PCA(n_components = n_components)
    pca.fit(data.T)

    barycenter = [data[0].mean(),data[1].mean(), data[2].mean()]
    return pca, barycenter

def runWPCA(n_components, xyz_T, energy, tot_energy):
    xyz = xyz_T.T
    assert xyz.shape[0] == 3, "Input data shape should be 3, check your input shape"
    # tot_energy = energy.sum()
    N = len(energy)
    weight_sum = 0
    weight_sum2 = 0
    barycenter = [0.,0.,0.]
    covM = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])

    # Compute the barycenter
    for i in range(N):
        weight = energy[i] #computing the weight
        weight_sum += weight
        barycenter[0] += xyz[0][i] * weight
        barycenter[1] += xyz[1][i] * weight
        barycenter[2] += xyz[2][i] * weight
    barycenter /= weight_sum    
    
    #Compute the covariance matrix
    covM = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
    for i in range(N):
        weight = energy[i]/tot_energy
        weight_sum2 += weight * weight
        for x in range(3):
            for y in range(3):
                covM[x][y] += weight * (xyz[x][i] - barycenter[x]) * (xyz[y][i] - barycenter[y])
    
    covMW = covM * 1./(1-weight_sum2) #weighted covariance matrix
    #get eigenvalues and eigenvectors and sort them
    eigen_values , eigen_vectors = np.linalg.eigh(covMW)
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    return sorted_eigenvectors.T, sorted_eigenvalue, covMW, covM, barycenter

def pcaSummary(pca):
    print("Covariance Matrix \n {} \n ".format(pca.get_covariance()))
    print("Eigenvectors \n {} \n".format(pca.components_))
    print("Eigenvalues \n {} \n".format(pca.explained_variance_))


def computeCovarianceMatrix(data):
    cov_mat = np.cov(data, rowvar = False)
    return cov_mat
