import ROOT
import cppyy
import pandas as pd
from utils import *
import random
from treeVar import *


def loadTree(rootFile, subdir = "trackstersAnalyzer", tree_name = "tree"):
    fileRoot = ROOT.TFile(rootFile)
    tree = fileRoot.Get(subdir + "/" + tree_name)
    return tree

def getBranches(tree):
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
    tree.SetBranchAddress("vertex_x", v_x)
    tree.SetBranchAddress("vertex_y", v_y)
    tree.SetBranchAddress("vertex_z", v_z)
    tree.SetBranchAddress("calopt", pt)
    tree.SetBranchAddress("rh_x", rh_x)
    tree.SetBranchAddress("rh_y", rh_y)
    tree.SetBranchAddress("rh_z", rh_z)
    tree.SetBranchAddress("rh_E", rh_E)
    tree.SetBranchAddress("rh_fraction", rh_fraction)
    tree.SetBranchAddress("lc_x", lc_x)
    tree.SetBranchAddress("lc_y", lc_y)
    tree.SetBranchAddress("lc_z", lc_z)
    tree.SetBranchAddress("lc_E", lc_E)
    tree.SetBranchAddress("lc_eta", lc_eta)
    tree.SetBranchAddress("lc_phi", lc_phi)
    tree.SetBranchAddress("N_lcs", N_lcs)
    tree.SetBranchAddress("cp_x", cp_x)
    tree.SetBranchAddress("cp_y", cp_y)
    tree.SetBranchAddress("cp_z", cp_z)
    tree.SetBranchAddress("vtxs_x", vtxs_x)
    tree.SetBranchAddress("vtxs_y", vtxs_y)
    tree.SetBranchAddress("vtxs_z", vtxs_z)
    tree.SetBranchAddress("vtxs_N", vtxs_N)
    tree.SetBranchAddress("bs_x", bs_x)
    tree.SetBranchAddress("bs_y", bs_y)
    tree.SetBranchAddress("bs_z", bs_z)

    tree.SetBranchAddress("cp_px", cp_px)
    tree.SetBranchAddress("cp_py", cp_py)
    tree.SetBranchAddress("cp_pz", cp_pz)
    tree.SetBranchAddress("cp_E", cp_E)
    tree.SetBranchAddress("calopt", cp_pt)
    tree.SetBranchAddress("caloeta", cp_eta)
    tree.SetBranchAddress("calophi", cp_phi)
    tree.SetBranchAddress("N_hit_per_lc", N_hit_per_lc)
    tree.SetBranchAddress("N_trackster", N_trackster)

    return tree
