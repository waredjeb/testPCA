#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from pcaUtils import *


# In[2]:


mkdir_p("./plots")


# In[3]:


#import data rechits
df = pd.read_csv("./postData/rhits_df.csv")
df = df.drop(df.columns[0], axis = 1)
event_number = 0
df_single_event = df[df.loc[:,"evt_id"] == event_number] #taking only one event --> one trackster per event so far
dataRH =df_single_event.iloc[:, 2:]


# In[4]:


dataRH.head()


# ## Plotting trackster's Rechits

# In[5]:


plots3DwithProjection(dataRH.values)


# In[2]:


df = pd.read_csv("./postData/lc_df.csv")
# df = df.drop(df.columns[0], axis = 1)
event_number = 0
df_single_event = df[df.loc[:,"evt_id"] == event_number] #taking only one event
dataLC =df_single_event.iloc[:, 2:]


# In[7]:


dataLC.head()
print(dataLC.shape)


# ## Plotting Trackster's LayerClusters (each point is the barycenter of a layercluster)

# In[8]:


plots3DwithProjection(dataLC.values)
# plt.show()

# # Classic PCA - No reweighting - RecHits

# In[9]:


dataRH.columns


# In[10]:


pcaClassicRH, barycenterClassicRH = runPCA(3, dataRH.iloc[:, :3].values)


# In[11]:


pcaSummary(pcaClassicRH)
print(barycenterClassicRH)


# In[12]:


plots(dataRH.values, barycenterClassicRH, pcaClassicRH.components_, pcaClassicRH.explained_variance_, noise = 0, save = "Rechits_Noweight", arrow_scale = 6)


# # Classic PCA - No reweighting - LayerClusters

# In[13]:


dataLC.columns


# In[14]:


pcaClassicLC, barycenterClassicLC = runPCA(3, dataLC.iloc[:, :3].values)


# In[15]:


pcaSummary(pcaClassicLC)


# In[16]:


plots(dataLC.values, barycenterClassicLC, pcaClassicLC.components_, pcaClassicLC.explained_variance_, noise = 0, save = "LCS_Noweight", arrow_scale = 6)


# # Comparison between the Principal Directions - RecHits vs LayerClusters

# In[17]:


vec1_RH = pcaClassicRH.components_[0] #principal component PCA on RecHits
vec1_LC = pcaClassicLC.components_[0] #principal component PCA on RecHits
dot = vec1_RH.dot(vec1_LC)
print(f'Radians {m.acos(dot)}, {np.rad2deg(m.acos(dot))}')


# ## Filtering both the datasets  by energy selection

# In[18]:


dataRH_filtered = dataRH[dataRH.iloc[:, 3] > 0.2]
dataLC_filtered = dataLC[dataLC.iloc[:, 3] > 1]


# #### RecHits

# In[19]:


plots3DwithProjection(dataRH_filtered.values)


# In[20]:


pcaClassicRHFilter, barycenterClassicRHFilter = runPCA(3, dataRH_filtered.iloc[:, :3].values)
pcaSummary(pcaClassicRHFilter)
print(barycenterClassicRHFilter)


# In[21]:


plots(dataRH_filtered.values, barycenterClassicRHFilter, pcaClassicRHFilter.components_, pcaClassicRHFilter.explained_variance_, noise = 0, save = "RHS_Noweight_filtered", arrow_scale = 6)


# In[22]:


plots(dataRH.values, barycenterClassicRHFilter, pcaClassicRHFilter.components_, pcaClassicRHFilter.explained_variance_, noise = 0, save = "RHS_Noweight_filtered_total", arrow_scale = 6)


# ### Layerclusters

# In[23]:


plots3DwithProjection(dataLC_filtered.values)


# In[24]:


pcaClassicLCFilter, barycenterClassicLCFilter = runPCA(3, dataLC_filtered.iloc[:, :3].values)
pcaSummary(pcaClassicLCFilter)


# In[25]:


plots(dataLC.values, barycenterClassicLCFilter, pcaClassicLCFilter.components_, pcaClassicLCFilter.explained_variance_, noise = 0, save = "LCS_Noweight_filtered", arrow_scale = 6)


# In[26]:


plots(dataLC_filtered.values, barycenterClassicLCFilter, pcaClassicLCFilter.components_, pcaClassicLCFilter.explained_variance_, noise = 0, save = "LCS_Noweight_filtered_total", arrow_scale = 6)


# ## Comparison 

# In[27]:


vec1_RHFiltered = pcaClassicRHFilter.components_[0]
vec1_LCFiltered = pcaClassicLCFilter.components_[0]
dot = vec1_RHFiltered.dot(vec1_LCFiltered)
print(f'Radians {m.acos(dot)}, {np.rad2deg(m.acos(dot))}')


# # Weighted PCA - Rechits

# In[28]:


dataRH.columns


# In[29]:


e_vecWRH, e_valWRH, covMWRH, covMRH, barycenterWRH, valid = runWPCA(3, dataRH.iloc[:, :3].values, dataRH.iloc[:,3].values, dataRH.iloc[:,3].values.sum())


# In[30]:


plots(dataRH.values, barycenterWRH, e_vecWRH, e_valWRH, noise = 0, save = "RH_Weighted", arrow_scale = 6)


# # Weighted PCA - Layer Clusters

# In[31]:


dataLC.columns


# In[32]:


e_vecWLC, e_valWLC, covMWLC, covMLC, barycenterWLC, valid = runWPCA(3, dataLC.iloc[:, :3].values, dataLC.iloc[:,3].values, dataLC.iloc[:,3].values.sum())


# In[33]:


plots(dataLC.values, barycenterWLC, e_vecWLC, e_valWLC, noise = 0, save = "LCS_Weighted", arrow_scale = 6)
plt.show()

# ### Comparison

# In[34]:


vec1_RHW= e_vecWRH[0]
vec1_LCW = e_vecWLC[0]
dot = vec1_RHW.dot(vec1_LCW)
print(f'Radians {m.acos(dot)}, {np.rad2deg(m.acos(dot))}')


# # Filtering

# In[35]:


e_vecWRHF, e_valWRHF, covMWRHF, covMRHF, barycenterWRHF, valid = runWPCA(3, dataRH_filtered.iloc[:, :3].values, dataRH_filtered.iloc[:,3].values, dataRH.iloc[:,3].values.sum())
plots(dataRH_filtered.values, barycenterWRHF, e_vecWRHF, e_valWRHF, noise = 0, save = "RH_Weighted_Filtered", arrow_scale = 6)


# In[36]:


e_vecWLCF, e_valWLCF, covMWLCF, covMLCF, barycenterWLCF,valid = runWPCA(3, dataLC_filtered.iloc[:, :3].values, dataLC_filtered.iloc[:,3].values, dataLC.iloc[:,3].values.sum())
plots(dataLC_filtered.values, barycenterWLCF, e_vecWLCF, e_valWLCF, noise = 0, save = "LCS_Weighted_Filtered", arrow_scale = 6)


# In[37]:


vec1_RHWFiltered = e_vecWRHF[0]
vec1_LCWFiltered = e_vecWLCF[0]
dot = vec1_RHWFiltered.dot(vec1_LCWFiltered)
print(f'Radians {m.acos(dot)}, {np.rad2deg(m.acos(dot))}')


# # Scaling

# In[39]:


dataLC.head()


# In[46]:


#import data rechits
df = pd.read_csv("./postData/caloparticle.csv")
df = df.drop(df.columns[0], axis = 1)
event_number = 0
df_single_event = df[df.loc[:,"evt_id"] == event_number] #taking only one event --> one trackster per event so far
dataCP =df_single_event.iloc[:, 2:]


# In[93]:


vecCP = dataCP.iloc[:,:3].values[0]
posCP = [88.482277,-19.239126,293.500000]


# In[94]:


length = vecCP.dot(vecCP) * 30


# In[95]:


vecCP = vecCP*length


# In[96]:


# get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure(figsize = (20,20))
ax = fig.add_subplot(2,2,1, projection = '3d')    
cm = plt.cm.get_cmap('hot')
noise_idx_l = -1
ax.scatter(dataLC_filtered.iloc[:,0], dataLC_filtered.iloc[:,1], dataLC_filtered.iloc[:,2], marker='o', alpha=.3, c = dataLC_filtered.iloc[:, 3], cmap = cm)
ax.quiver(posCP[0],posCP[1],posCP[2],vecCP[0],vecCP[1],vecCP[2], color = 'r')


# In[3]:


dataLC


# In[29]:


idx = dataLC.groupby(["lc_z"], sort=False)['lc_E'].transform(max) == dataLC['lc_E']


# In[31]:


a = dataLC[idx]


# In[27]:


a


# In[ ]:




