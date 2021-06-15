from numpy.core.einsumfunc import _find_contraction
import ROOT
from pcaUtils import *

class CaloParticle:
    def __init__(self, px,py,pz,pt,E,eta,phi,x,y,z,event_id, trackster_id):
        self.px = px
        self.py = py
        self.pz = pz
        self.pt = pt
        self.E = E
        self.eta = eta
        self.phi = phi
        self.x = x
        self.y = y
        self.z = z
        self.event_id = event_id     
        self.trackster_id = trackster_id   

class RecHit:
    def __init__(self, x,y,z,E, fraction):
        self.x = x
        self.y = y
        self.z = z
        self.E = E
        self.fraction = fraction

class LayerCluster:
    def __init__(self, x, y, z, E, trackster_id, rechits = []):
        self.x = x 
        self.y = y
        self.z = z
        self.E = E
        self.trackster_id = trackster_id
        self.rechits = rechits
        self.Nhits = len(rechits)

    def addRechit(self, rh):
        return self.rechits.append(rh)

class Trackster:
    # barycenter, eigenvectors, eigenvalues, dataframe, covM = 0
    def __init__(self, event_id, id, layerclusters = []):
        self.event_id = event_id
        self.layerclusters = layerclusters
        self.id = id
        # self.barycenter
        # self.eigenvectors
        # self.eigenvalues
        # self.dataframe
        # self.covM
    
    def checkPCA(self):
        if(self.barycenter == []):
            self.valid = False
        else:
            self.valid = True

    def selectBestLCS(self, startingLC = 0, endLC = -1):
        energy_cut = 2
        df = self.dataframe
        tot_energy = df.loc[0,'tot_energy']
        idx = df.groupby(["lc_z"], sort=False)['lc_E'].transform(max) == df['lc_E']
        df_tmp = df[idx]
        max_energy = max(df_tmp['lc_E'].values)
        fraction_cut = 0.005 #GeV
        idx_filter = df_tmp['lc_E'] > fraction_cut
        df_tmp2 = df_tmp[idx_filter]
        # print(df_tmp2.iloc[startingLC:endLC, :])
        
        self.bestLCSdf = df_tmp2.iloc[startingLC:endLC, :]
        bestLCS = []
        for i in range(len(idx)):
            if(idx[i] == True):
                # if(len(self.layerclusters[i].rechits) > 1):
                # if(abs(self.layerclusters[i].z) > 330.0127258300781 and abs(self.layerclusters[i].z) < 340.86724853515625):
                if(self.layerclusters[i].E/tot_energy > fraction_cut):
                    bestLCS.append(self.layerclusters[i])
                # else:
                    # continue
            else:
                continue
        bestCentralLCs = bestLCS[startingLC:endLC]
        self.bestLCS = bestCentralLCs

    def selectBestLCSEnergy(self, fraction = 0.05):
        df = self.dataframe
        # idx = df.groupby(["lc_z"], sort=False)['lc_E'].transform(max) == df['lc_E']
        idx = df.loc[:, "lc_E"] / df.loc[:, "tot_energy"] > fraction
        self.bestLCSdf = df[idx]
        bestLCSFilter = []
        for i in range(len(idx)):
            if(idx[i] == True):
                # if(len(self.layerclusters[i].rechits) > 1):
                bestLCSFilter.append(self.layerclusters[i])
                # else:
                    # continue
            else:
                continue
        self.bestLCSFilter = bestLCSFilter

    def addLayerCluster(self,lc):
        return (self.layerclusters).append(lc)
    
    def Nlayerclusters(self):
        return len(self.layerclusters)
    
    def makeDataFrameLC(self, startLC = 0, endLC = -1):
        evt_id = self.event_id
        lcs = self.layerclusters
        d = {"evt_id": [evt_id] * len(lcs),
                "lc_x": [lc.x for lc in lcs],
                "lc_y": [lc.y for lc in lcs],
                "lc_z": [lc.z for lc in lcs],
                "lc_E": [lc.E for lc in lcs]}
        df = pd.DataFrame(d)
        df['tot_energy'] = [df['lc_E'].sum()] * df.shape[0]
        self.dataframe = df
        self.selectBestLCS(startingLC=startLC, endLC=endLC)
        return df

    def makeDataFrameRH(self, fraction = 0, bestLCS = False):
        evt_id = self.event_id
        if(bestLCS == True):
            lcs = self.bestLCS
        else:
            if(fraction == 0):
                lcs = self.layerclusters
            else:
                self.selectBestLCSEnergy(fraction  = fraction)
                lcs = self.bestLCSFilter
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
        self.dataframe = df
        return df

    def assignFAToTracksterLC(self, n_components, weighted = False, filter_percentage = 0, scaling = False,  nonlinearweight = False, bestLCS = False, startingLC = 0, endLC = -1, w0 = 0):
        self.makeDataFrameLC(startLC=startingLC, endLC=endLC)
        if(bestLCS == True):
            df = self.bestLCSdf
        else:
            df = self.dataframe
        filter = df.loc[:, "tot_energy"].values[0] * filter_percentage
        df_filtered = df[df.loc[:, "lc_E"] > filter]
        # print(df_filtered.loc[:, "lc_E"])
        xyz_tmp = df_filtered.loc[:, ['lc_x', 'lc_y', 'lc_z']].values
        if(xyz_tmp.shape[0] < 2):
            self.valid = False
            return False
        else:
            if(scaling == True):        
                scaler = StandardScaler()
                xyz = scaler.fit_transform(xyz_tmp)
                self.valid = False
                return False
            else:
                xyz = xyz_tmp
            if(weighted == False):            
                pca, barycenter, valid = runFA(n_components, xyz)       
                if(valid == False):
                    self.valid = valid
                    self.barycenter = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.covM = [0,0,0]
                else:
                    self.barycenter = np.array(barycenter)
                    self.eigenvectors = pca.components_
                    self.eigenvalues = [1,1,1]
                    self.covM = [0,0,0]
                    self.valid = valid
            else:
                energy = df_filtered.loc[:, 'lc_E'].values
                sorted_eigenvectors, sorted_eigenvalue, covMW, covM, barycenter, valid = runWPCA(n_components, xyz, energy, energy.sum(), nonlinear = nonlinearweight, w0 = w0)          
                self.valid = valid
                _, _, _, _, barycenter_real, _ = runWPCA(n_components, xyz_tmp, energy, energy.sum(), w0 = w0)
                self.barycenter = barycenter_real
                self.eigenvectors = sorted_eigenvectors
                self.eigenvalues = sorted_eigenvalue
                self.covM = covMW

    def assignPCAToTracksterLC(self, n_components, weighted = False, filter_percentage = 0,  scaling = False,  nonlinearweight = False, bestLCS = False, startingLC = 0, endLC = -1, incremental = False,w0 = 0):
        self.makeDataFrameLC(startLC=startingLC, endLC=endLC)
        if(bestLCS == True):
            df = self.bestLCSdf
        else:
            df = self.dataframe
        filter = df.loc[:, "tot_energy"].values[0] * filter_percentage
        df_filtered = df[df.loc[:, "lc_E"] > filter]
        
        # print(df_filtered.loc[:, "lc_E"])
        xyz_tmp = df_filtered.loc[:, ['lc_x', 'lc_y', 'lc_z']].values
        if(xyz_tmp.shape[0] < 2):
            self.valid = False
            return False
        else:
            if(scaling == True):        
                scaler = StandardScaler()
                xyz = scaler.fit_transform(xyz_tmp)
                self.valid = False
                return False
            else:
                xyz = xyz_tmp
            if(weighted == False):            
                pca, barycenter, valid = runPCA(n_components, xyz, incremental = incremental)       
                if(valid == False):
                    self.valid = valid
                    self.barycenter = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.covM = [0,0,0]
                else:
                    self.barycenter = np.array(barycenter)
                    self.eigenvectors = pca.components_
                    self.eigenvalues = pca.explained_variance_
                    self.covM = pca.get_covariance()
                    self.valid = valid
            else:
                # print(xyz)
                energy = df_filtered.loc[:, 'lc_E'].values
                sorted_eigenvectors, sorted_eigenvalue, covMW, covM, barycenter, valid = runWPCA(n_components, xyz, energy, energy.sum(), nonlinear = nonlinearweight, w0 = w0)          
                self.valid = valid
                _, _, _, _, barycenter_real, _ = runWPCA(n_components, xyz, energy, energy.sum(), w0 = w0)
                self.barycenter = barycenter_real
                self.eigenvectors = sorted_eigenvectors
                self.eigenvalues = sorted_eigenvalue
                self.covM = covMW

    def assignPCAToTracksterLC2D(self, n_components, weighted = False, filter_percentage = 0, scaling = False,  nonlinearweight = False, bestLCS = False, w0 = 0):
        self.makeDataFrameLC()
        if(bestLCS == True):
            df = self.bestLCSdf
        else:
            df = self.dataframe
        filter = df.loc[:, "tot_energy"].values[0] * filter_percentage
        df_filtered = df[df.loc[:, "lc_E"] > filter]
        # print(df_filtered.loc[:, "lc_E"])
        xyz_tmp = df_filtered.loc[:, ['lc_x', 'lc_y', 'lc_z']].values
        if(xyz_tmp.shape[0] < 2):
            self.valid = False
            return False
        else:
            if(scaling == True):        
                scaler = StandardScaler()
                xyz = scaler.fit_transform(xyz_tmp)
                self.valid = False
                return False
            else:
                xyz = xyz_tmp
            if(weighted == False):            
                pca, barycenter, valid = runPCA(n_components, xyz)       
                if(valid == False):
                    self.valid = valid
                    self.barycenter = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.eigenvalues = [0,0,0]
                    self.covM = [0,0,0]
                else:
                    self.barycenter = np.array(barycenter)
                    self.eigenvectors = pca.components_
                    self.eigenvalues = pca.explained_variance_
                    self.covM = pca.get_covariance()
                    self.valid = valid
            else:
                energy = df_filtered.loc[:, 'lc_E'].values
                sorted_eigenvectors, sorted_eigenvalue, covMW, covM, barycenter, valid = runWPCA2D(n_components, xyz, energy, energy.sum(), nonlinear = nonlinearweight, w0 = w0)          
                self.valid = valid
                _, _, _, _, barycenter_real, _ = runWPCA2D(n_components, xyz_tmp, energy, energy.sum(), w0 = w0)
                self.barycenter = barycenter_real
                self.eigenvectors = sorted_eigenvectors
                self.eigenvalues = sorted_eigenvalue
                self.covM = covMW


    def assignPCAToTracksterRH(self, n_components, weighted = False, filter_percentage = 0, scaling = False, nonlinearweight = False, bestLCS = False):
        self.makeDataFrameRH(bestLCS = bestLCS, fraction = filter_percentage)
        if(bestLCS == True):
            df = self.dataframe
        else:
            df = self.dataframe
        # filter = df.loc[:, "rh_E"].sum() * filter_percentage
        # df_filtered = df[df.loc[:, "rh_E"] > filter]
        xyz_tmp = df.loc[:, ['rh_x', 'rh_y', 'rh_z']].values
        if(scaling == True):
            scaler = StandardScaler()
            xyz = scaler.fit_transform(xyz_tmp)
        else:
            xyz = xyz_tmp

        if(weighted == False):            
            pca, barycenter, valid = runPCA(n_components, xyz)       
            if(valid == False):
                self.valid = valid
                self.barycenter = [0,0,0]
                self.eigenvalues = [0,0,0]
                self.eigenvalues = [0,0,0]
                self.covM = [0,0,0]
            else:
                self.barycenter = np.array(barycenter)
                self.eigenvectors = pca.components_
                self.eigenvalues = pca.explained_variance_
                self.covM = pca.get_covariance()
                self.valid = valid
        else:
            energy = df.loc[:, 'rh_E'].values
            sorted_eigenvectors, sorted_eigenvalue, covMW, covM, barycenter, valid  = runWPCA(n_components, xyz, energy, energy.sum(), nonlinear = nonlinearweight)
            # if(valid == False):
            self.valid = valid
            _, _, _, _, barycenter_real, _ = runWPCA(n_components, xyz_tmp, energy, energy.sum())
            self.barycenter = barycenter_real
            self.eigenvectors = sorted_eigenvectors
            self.eigenvalues = sorted_eigenvalue
            self.covM = covMW

    

