import ROOT
import cppyy
import pandas as pd
from utils import *
import random

rootFile = "./data/output_singlephoton.root"
nameSave = "_CloseByGamma"
plotDir = "./plots_new/"
pullPlotsDir = plotDir+"pulls/"
saveFormat = ".png"
fileRoot = ROOT.TFile(rootFile)
# treeName = ROOT.TFile("tree")

tree = fileRoot.Get("trackstersAnalyzer/tree")
angles_pca_momentum = ROOT.std.vector('double')()
angles_bary_momentum = ROOT.std.vector('double')()
N_best_lcs = ROOT.std.vector('double')()
deta_lcs = ROOT.std.vector('double')()
dphi_lcs = ROOT.std.vector('double')()
deta_lcs_maxE = ROOT.std.vector('double')()
dphi_lcs_maxE = ROOT.std.vector('double')()
deta_lcs_even = ROOT.std.vector('double')()
dphi_lcs_even = ROOT.std.vector('double')()
deta_lcs_odd = ROOT.std.vector('double')()
dphi_lcs_odd = ROOT.std.vector('double')()
v_eta = ROOT.std.vector('double')()
v_phi = ROOT.std.vector('double')()
# deltaR = ROOT.std.vector('double')()
pt = ROOT.std.vector('double')()
# angles_lcs = ROOT.std.vector('double')()

tree.SetBranchAddress("angles_pca_momentum", angles_pca_momentum)
tree.SetBranchAddress("angles_bary_momentum", angles_bary_momentum)
tree.SetBranchAddress("deta_lcs", deta_lcs)
tree.SetBranchAddress("dphi_lcs", dphi_lcs)
tree.SetBranchAddress("deta_lcs_maxE", deta_lcs_maxE)
tree.SetBranchAddress("dphi_lcs_maxE", dphi_lcs_maxE)
tree.SetBranchAddress("deta_lcs_even", deta_lcs_even)
tree.SetBranchAddress("dphi_lcs_even", dphi_lcs_even)
tree.SetBranchAddress("deta_lcs_odd", deta_lcs_odd)
tree.SetBranchAddress("dphi_lcs_odd", dphi_lcs_odd)
# tree.SetBranchAddress("N_best_lcs", N_best_lcs)
tree.SetBranchAddress("vertex_eta", v_eta)
tree.SetBranchAddress("vertex_phi", v_phi)
tree.SetBranchAddress("calopt", pt)

angles_bary_momentum_list = []
angles_pca_momentum_list = []
pt_list = []
eta_list = []
phi_list = []

bin_pulls = 100
angles_bary_momentum_hist = ROOT.TH1F("Angles", "Angles", 50, 0, 0.05)
angles_pca_momentum_hist = ROOT.TH1F("Angles", "Angles", 50, 0, 0.05)
deta_hist = ROOT.TH1F("#Delta#eta", "#Delta#eta", 100, -0.05, 0.05)
dphi_hist = ROOT.TH1F("#Delta#phi", "#Delta#phi", 100, -0.05, 0.05)
deta_hist_maxE = ROOT.TH1F("#Delta#eta", "#Delta#eta", 100, -0.05, 0.05)
dphi_hist_maxE = ROOT.TH1F("#Delta#phi", "#Delta#phi", 100, -0.05, 0.05)
deta_hist_even = ROOT.TH1F("#Delta#eta", "#Delta#eta", 100, -0.05, 0.05)
dphi_hist_even = ROOT.TH1F("#Delta#phi", "#Delta#phi", 100, -0.05, 0.05)
deta_hist_odd = ROOT.TH1F("#Delta#eta", "#Delta#eta", 100, -0.05, 0.05)
dphi_hist_odd = ROOT.TH1F("#Delta#phi", "#Delta#phi", 100, -0.05, 0.05)
angles_pca_momentum_hist_vs_veta = ROOT.TH2F("Angles", "Angles", 20, 0, 0.1, 20, 1.3, 2.5)
angles_pca_momentum_hist_vs_cppt = ROOT.TH2F("Angles", "Angles", 20, 0, 0.1, 20, 0, 100)
vtx_phi = ROOT.TH1F("#Delta#eta", "#Delta#eta", 100, 1.9, 2.1)

angles_bary_momentum_pull_hist = ROOT.TH1F("Angles", "Angles", 50, -1,1)
angles_pca_momentum_pull_hist = ROOT.TH1F("Angles", "Angles", 50, -1, 1)
pt_vs_angles_bary_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 0, 100, -1,1)
pt_vs_angles_pca_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 0, 100, -1, 1)
eta_vs_angles_bary_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 1.4, 2.4, -1,1)
eta_vs_angles_pca_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 1.4, 2.4, -1, 1)
phi_vs_angles_bary_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 1.9, 2.1, -1,1)
phi_vs_angles_pca_momentum_pull_hist = ROOT.TProfile("Angles", "Angles",bin_pulls, 1.9, 2.1, -1, 1)

N = tree.GetEntries()
print(N)
for evt in range(0,N):
    tree.GetEntry(evt) 
    for i in range(angles_bary_momentum.size()):
        angles_bary_momentum_hist.Fill(angles_bary_momentum.at(i))
        angles_bary_momentum_list.append(angles_bary_momentum.at(i))
        angles_pca_momentum_list.append(angles_pca_momentum.at(i))
        angles_pca_momentum_hist.Fill(angles_pca_momentum.at(i))
        angles_pca_momentum_hist_vs_veta.Fill(angles_pca_momentum.at(i), v_eta.at(i))
        angles_pca_momentum_hist_vs_cppt.Fill(angles_pca_momentum.at(i), pt.at(i))
        pt_list.append(pt.at(i))
        eta_list.append(v_eta.at(i))
        phi_list.append(v_phi.at(i))
        
    for i in range(deta_lcs.size()):
        deta_hist.Fill(deta_lcs.at(i))  
        dphi_hist.Fill(dphi_lcs.at(i))

    for i in range(deta_lcs_even.size()):
        deta_hist_even.Fill(deta_lcs_even.at(i))
        dphi_hist_even.Fill(dphi_lcs_even.at(i))
    
    for i in range(deta_lcs_odd.size()):
        deta_hist_odd.Fill(deta_lcs_odd.at(i))
        dphi_hist_odd.Fill(dphi_lcs_odd.at(i))
    
    for i in range(deta_lcs_maxE.size()):
        deta_hist_maxE.Fill(deta_lcs_maxE.at(i))
        dphi_hist_maxE.Fill(dphi_lcs_maxE.at(i))
    


#Pull plots
pull_bary_momentum = pull(angles_bary_momentum_list)
pull_pca_momentum = pull(angles_pca_momentum_list)

for i in range(len(pull_bary_momentum)):
    angles_bary_momentum_pull_hist.Fill(pull_bary_momentum[i])
    angles_pca_momentum_pull_hist.Fill(pull_pca_momentum[i])
    pt_vs_angles_bary_momentum_pull_hist.Fill(pt_list[i], pull_bary_momentum[i])
    pt_vs_angles_pca_momentum_pull_hist.Fill(pt_list[i], pull_pca_momentum[i])
    eta_vs_angles_bary_momentum_pull_hist.Fill(eta_list[i], pull_bary_momentum[i])
    eta_vs_angles_pca_momentum_pull_hist.Fill(eta_list[i], pull_pca_momentum[i])
    phi_vs_angles_bary_momentum_pull_hist.Fill(phi_list[i], pull_bary_momentum[i])
    phi_vs_angles_pca_momentum_pull_hist.Fill(phi_list[i], pull_pca_momentum[i])   
    vtx_phi.Fill(phi_list[i]) 

c = ROOT.TCanvas("c", "c")

singlePlot(c, angles_bary_momentum_hist, "Angle between: Barycenter w.r.t Vertex Position VS CP Momentum", "Angles[rad]", "Entries", "Angles_Bary_Vertex"+nameSave )
singlePlot(c, angles_pca_momentum_hist, "Angle between: PCA direction VS CP Momentum", "Angles[rad]", "Entries", "Angles_PCA"+nameSave )
singlePlot(c, deta_hist, "#Delta#eta between layerclusters and vertex", "#Delta#eta", "Entries", "DEta_LCS_vertex"+nameSave )
singlePlot(c, dphi_hist, "#Delta#phi between layerclusters and vertex", "#Delta#phi", "Entries", "DPhi_LCS_vertex"+nameSave )
singlePlot(c, deta_hist_maxE, "#Delta#eta between Highest Energy LC and vertex", "#Delta#eta", "Entries", "DEta_LCS_HighE_vertex"+nameSave )
singlePlot(c, dphi_hist_maxE, "#Delta#phi between Highest Energy LC and vertex", "#Delta#phi", "Entries", "DPhi_LCS_HighE_vertex"+nameSave )
singlePlot(c, deta_hist_even, "#Delta#eta between Even layerclusters and vertex", "#Delta#eta", "Entries", "DEta_LCS_Even_vertex"+nameSave )
singlePlot(c, dphi_hist_even, "#Delta#phi between Even layerclusters and vertex", "#Delta#phi", "Entries", "DPhi_LCS_Even_vertex"+nameSave )
singlePlot(c, deta_hist_odd, "#Delta#eta between Odd layerclusters and vertex", "#Delta#eta", "Entries", "DEta_LCS_Odd_vertex"+nameSave )
singlePlot(c, dphi_hist_odd, "#Delta#phi between Odd layerclusters and vertex", "#Delta#phi", "Entries", "DPhi_LCS_Odd_vertex"+nameSave )
singlePlot(c, vtx_phi, "Vertex #phi", "#phi", "Entries", "Vertex_phi"+nameSave )
singlePlot(c, angles_bary_momentum_pull_hist, "Pull plot angle between: Barycenter wrt Vertex VS CP momentum", "Pull Angles", "Entries", "Pull_Angles_Bary_Vertex"+nameSave )
singlePlot(c, angles_pca_momentum_pull_hist, "Pull plot angle between: PCA direction VS CP momentum", "Pull Angles", "Entries", "Pull_Angles_PCA"+nameSave )
singlePlot(c, angles_pca_momentum_hist_vs_veta,"Angle between:  Vertex #eta VS (PCA Direction VS CP Momentum) ", "Angles[rad]", "#eta", "Angles_PCAwrtMomentum_vs_Vertex_eta"+nameSave, is2D=True )
singlePlot(c, angles_pca_momentum_hist_vs_veta,"Angle between: CP pT VS (PCA Direction VS CP Momentum)", "Angles[rad]", "pT[GeV]", "Angles_PCAwrtMomentum_vs_CP_Pt"+nameSave, is2D=True )
singlePlot(c, pt_vs_angles_bary_momentum_pull_hist,"Pull plot: Angle between Bary wrt Vertex and CP momentum VS CP pT ", "pT[Gev]", "pull angle", "Pull_Angles_BarywrtMomentum_vs_CP_pT"+nameSave, is2D=False, output_dir = pullPlotsDir )
singlePlot(c, pt_vs_angles_pca_momentum_pull_hist,"Pull plot: Angle between PCA direction and CP momentum VS CP pT ", "pT[Gev]", "pull angle", "Pull_Angles_PCAwrtMomentum_vs_CP_pT"+nameSave, is2D=False, output_dir = pullPlotsDir )
singlePlot(c, eta_vs_angles_bary_momentum_pull_hist,"Pull plot: Angle between Bary wrt Vertex and CP momentum VS Vertex #eta ", "#eta", "pull angle", "Pull_Angles_BarywrtMomentum_vs_Vtx_eta"+nameSave, is2D=False, output_dir = pullPlotsDir )
singlePlot(c, eta_vs_angles_pca_momentum_pull_hist,"Pull plot: Angle between PCA direction and CP momentum VS Vertex #eta ", "#eta", "pull angle", "Pull_Angles_PCAwrtMomentum_vs_Vtx_eta"+nameSave, is2D=False, output_dir = pullPlotsDir )
multiPlot(c,deta_hist_even, deta_hist_odd, "#Delta#eta between LayerCluster and Vertex", "#Delta#eta", "Entries", "Even", "Odd","DEta_LCS_EvenAndOdd_vertex"+nameSave )
multiPlot(c,dphi_hist_even, dphi_hist_odd, "#Delta#phi between LayerCluster and Vertex", "#Delta#phi", "Entries", "Even", "Odd","Dphi_LCS_EvenAndOdd_vertex"+nameSave )
#Useless, phi is fixed at one value
# singlePlot(c, phi_vs_angles_bary_momentum_pull_hist,"Pull plot: Angle between Bary wrt Vertex and CP momentum VS Vertex #phi ", "#phi", "pull angle", "Pull_Angles_BarywrtMomentum_vs_Vtx_phi"+nameSave, is2D=False, output_dir = pullPlotsDir )
# singlePlot(c, phi_vs_angles_pca_momentum_pull_hist,"Pull plot: Angle between PCA direction and CP momentum VS Vertex #phi ", "#phi", "pull angle", "Pull_Angles_PCAwrtMomentum_vs_Vtx_phi"+nameSave, is2D=False, output_dir = pullPlotsDir )




