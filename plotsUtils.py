from matplotlib.pyplot import bar, title
from classes import CaloParticle, LayerCluster, RecHit, Trackster
from readTree import *
from tqdm import tqdm
import pickle
plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.rcParams.label.data = True
hep.rcParams.label.paper = False

def singlePlot(c, histo, title, xtitle, ytitle, savetitle, format = '.png', is2D = False, output_dir = "./plotsAnalysis/"):
    histo.SetTitle(title)
    histo.GetXaxis().SetTitle(xtitle)
    histo.GetYaxis().SetTitle(ytitle)
    if(is2D == True):
        histo.Draw('COLZ')
    else:
        histo.Draw()
    c.SaveAs(output_dir + savetitle + format, format)
    c.GetListOfPrimitives().Remove(histo)
    c.Modified()
    return c

    


def plotHisto(list, title, x_title, y_title, saveName, bins = 100, range = None, outputdir = './plotsAnalysis/'):    
    fig = plt.figure(figsize = (20,15))
    if(range != None):
        range_h = range
    else:
        range_h = (min(list), max(list))
    
    plt.hist(list, bins = bins, range = range_h, histtype='step', linewidth = 3)
    plt.title(title)
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)
    plt.xlabel(x_title, fontsize = 30)
    plt.ylabel(y_title, fontsize = 30)
    plt.savefig(outputdir + saveName + ".png")
    pickle.dump(fig, open(outputdir + "pickle/" + saveName +'.fig.pickle', 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
    plt.clf()


def fitHisto(list, title, x_title, y_title, saveName, bins = 100, range = None, PU = 0, eta = 0, outputdir = './plotsAnalysis/'):
    c = ROOT.TCanvas("c","c")
    if(range == None):
        h = ROOT.TH1F("h", "h", bins , min(list), max(list))
    else:
        h = ROOT.TH1F("h", "h", bins , range[0], range[1])
    gaussFit = ROOT.TF1("gaussfit","gaus" ,-5 ,5)
    for i in list:
        h.Fill(i)
    
    h.Fit(gaussFit,'E')
    h.SetStats(0)
    h.SetTitle(title)
    h.GetXaxis().SetTitle(x_title)
    h.GetYaxis().SetTitle(y_title)
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    chi2   = gaussFit.GetChisquare()
    ndof   = gaussFit.GetNDF()
    mean   = gaussFit.GetParameter(1)
    width = gaussFit.GetParameter(2)
    h.Draw()
    latex.DrawLatex(0.7 ,0.80 ,"Mean = %.4f"%(mean))
    latex.DrawLatex(0.7 ,0.75 ,"Width = %.4f"%(width ))
    latex.DrawLatex(0.7 ,0.70 ,"PU = %.0f"%(PU))
    if(eta > 0):
        latex.DrawLatex(0.7 ,0.65 ,"eta = %.1f"%(eta))
    else:
        latex.DrawLatex(0.7 ,0.65 ,"eta = [1.7, 2.7]")
    # latex.DrawLatex(0.7,0.7 ,"chi2/ndof = %.2f/%d = %.2f"%(chi2 ,ndof ,chi2/ndof))
    # latex.Draw()
    s = outputdir+saveName
    c.Show()
    c.SaveAs(s + ".png", "png")

def plotMultiHisto(list_h, title, x_title,y_title, saveName, labels = ["h1", "h2", "h3"], bins = 100, ranges  = None, density = False, outputdir = './plotsAnalysis/'):   
    fig = plt.figure(figsize = (20,15))
    if(range != None):
        range_h = ranges
    else:
        range_h = (min(list_h[0]), max(list_h[0]))
    
    for i in range(len(list_h)):
        h = list_h[i]
        l = labels[i]
        plt.hist(h, bins = bins, range = range_h, label = l, alpha = 1, histtype='step', linewidth = 3, density = density)
        plt.title(title)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.xlabel(x_title, fontsize = 30)
        plt.ylabel(y_title, fontsize = 30)
        plt.legend(loc= 'upper left')
    plt.savefig(outputdir + saveName + ".png")
    pickle.dump(fig, open(outputdir + "pickle/" + saveName +'.fig.pickle', 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
    plt.clf()


def plotMultiHistoDifferent(list_h, title, x_title,y_title, saveName, diff_h = 3, labels = ["h1", "h2", "h3"], bins = 100, ranges  = None, density = False, outputdir = './plotsAnalysis/'):   
    fig = plt.figure(figsize = (20,15))
    if(range != None):
        range_h = ranges
    else:
        range_h = (min(list_h[0]), max(list_h[0]))
    colors = ['cyan', 'fuchsia', 'black']
    j = 0
    for i in range(len(list_h)):
        h = list_h[i]
        l = labels[i]
        if(i < diff_h):            
            plt.hist(h, bins = bins, range = range_h, label = l, alpha = 1, histtype='step', linewidth = 3, density = density)
        else:
            plt.hist(h, bins = bins, range = range_h, label = l, alpha = 0.5, color=colors[j],  histtype='stepfilled', linewidth = 3, density = density)
            j+=1
        plt.title(title)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.xlabel(x_title, fontsize = 30)
        plt.ylabel(y_title, fontsize = 30)
        plt.legend(loc= 'upper left')
    plt.savefig(outputdir + saveName + ".png")
    pickle.dump(fig, open(outputdir + "pickle/" + saveName +'.fig.pickle', 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
    plt.clf()