import numpy as np
import ROOT
from pcaUtils import *

def pull(l):
    list_np = np.array(l)
    # print(l)
    mu = np.nanmean(l)
    std = np.nanstd(l)
    pull = [(i-mu)/std for i in list_np]
    return pull

def singlePlot(c, histo, title, xtitle, ytitle, savetitle, format = '.png', is2D = False, output_dir = "plots_new/"):
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

def multiPlot(c, h1, h2, title, xtitle, ytitle, legtitle1, legtitle2, savetitle, format = '.png', output_dir ="plots_new/"):
    leg = ROOT.TLegend(0.8,0.8,1,0.9)
    h1.SetTitle(title)
    h1.GetXaxis().SetTitle(xtitle)
    h1.GetYaxis().SetTitle(ytitle)
    h1.SetFillColor(ROOT.kWhite)
    h2.SetFillColor(ROOT.kWhite)
    h1.SetLineColor(ROOT.kBlue)
    h2.SetLineColor(ROOT.kRed)
    h1.SetLineWidth(2)
    h1.SetStats(0000)
    h2.SetStats(0000)
    h2.SetLineWidth(2)
    leg.AddEntry(h1, legtitle1)
    leg.AddEntry(h2, legtitle2)
    h1.Draw()
    h2.Draw("same")
    leg.Draw("same")
    c.SaveAs(output_dir + savetitle + format, format)
    c.GetListOfPrimitives().Remove(h1)
    c.GetListOfPrimitives().Remove(h2)
    c.Modified()
    return c

def PlotRatio(h1, h2, title_1 = "New Electrons", title_2 = "Legacy Electron", output_dir = "./plot/", max_y = 600, normalize = True, plot_format = "png", save = True, extra = "", title = "Electrons"):
    saveH = False
    if(h1.Integral() > 0 and h2.Integral() > 0):
        saveH = True
        c3 = ROOT.TCanvas("c3","c3",800,800)
        
        if(normalize == True):
            h1.Scale(1/h1.Integral())
            h2.Scale(1/h2.Integral())
        pad1 = ROOT.TPad("pad1","pad1", 0,0.3,1,1.0)
        pad1.SetBottomMargin(0.1)
        # pad1.SetGridx()
        pad1.Draw()
        pad1.cd()
        leg = ROOT.TLegend(0.55, 0.75, 0.7,0.9)
        leg.AddEntry(h1, title_1 + "[" + str(round(h1.Integral(),4)) + "]")
        leg.AddEntry(h2, title_2 + "[" + str(round(h2.Integral(),4)) + "]")
        leg.SetTextSize(0.03)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        h1.SetStats(0)
        h1.SetTitle(title)
        h2.SetTitle(title)
        xaxis = h1.GetXaxis()
        xmin = xaxis.GetBinCenter(xaxis.GetFirst())
        xmax = xaxis.GetBinCenter(xaxis.GetLast())
        h1.GetYaxis().SetRangeUser(0,max_y)
        h1.Draw("histo")
        h2.Draw("histo same")
        h1.GetYaxis().SetLabelSize(0.04)
        # axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
        # axis.SetLabelFont(1) # Absolute font size in pixel (precision 3)
        # axis.SetLabelSize(15)
        # axis.Draw()
        leg.Draw()
        # T = ROOT.TLatex()
        # text = "#scale[1.]{ #font[62]{CMS} #font[52]{Preliminary} }"
        # T.DrawLatexNDC(.1, .92, text)
        c3.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.1)
        pad2.SetBottomMargin(0.1)
        pad2.SetGridx()
        pad2.Draw()
        pad2.cd()
        line = ROOT.TLine(xmin,1,xmax,1)
        line.SetLineColorAlpha(ROOT.kRed, 0.6)
        line.SetLineWidth(3)
        line.SetLineStyle(7)
        

        h3 = h1.Clone("h1")
        h3.SetLineColor(ROOT.kBlack)
        h3.SetMinimum(0.1)
        h3.SetMaximum(2)
        h3.Sumw2()
        h3.SetStats(0)
        h3.Divide(h2)
        h3.SetMarkerStyle(21)
        h3.Draw("EP")
        line.Draw()

        h1.SetLineColor(ROOT.kBlue + 1)
        h1.SetFillColorAlpha(ROOT.kBlue + 1, 0.11)
        h1.SetFillStyle(1001)
        h1.SetLineWidth(2)
        h1.GetYaxis().SetTitleSize(23)
        h1.GetYaxis().SetTitleFont(43)
        h1.GetXaxis().SetTitleSize(23)
        h1.GetXaxis().SetTitleOffset(1.4)
        h1.GetXaxis().SetTitleFont(43)

        h2.SetLineColor(ROOT.kRed)
        h2.SetFillColor(ROOT.kRed)
        h2.SetFillStyle(0)
        h2.SetLineWidth(3)

        h3.SetTitle("")

        # / Y axis ratio plot settings
        h3.GetYaxis().SetTitle("#frac{"+title_1+"}{"+title_2+"} ")
        # h3.GetYaxis().SetNdivisions(505)
        h3.GetYaxis().SetTitleSize(20)
        h3.GetYaxis().SetTitleFont(43)
        h3.GetYaxis().SetTitleOffset(1.55)
        h3.GetYaxis().SetLabelFont(43) #// Absolute font size in pixel (precision 3)
        h3.GetYaxis().SetLabelSize(15)

        #X axis ratio plot settings
        h3.GetXaxis().SetTitleSize(20)
        h3.GetXaxis().SetTitleFont(43)
        h3.GetXaxis().SetTitleOffset(4.)
        h3.GetXaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
        h3.GetXaxis().SetLabelSize(15)
        c3.Draw()
        c3.Show()
    else:
        print("One of the two histos is empty")
        return 0
        
        
    output_name = ""
    output_name += output_dir
    output_name += h1.GetName() + "_Ratio"    
    output_name_root = output_name + extra + ".root"    
    output_name += extra + "."+plot_format
    
    if(save):
        c3.SaveAs(output_name, plot_format)
        c3.SaveAs(output_name_root, 'root')
        print("Histos {}_Comparison saved in {}".format(h1.GetName(), output_dir))
    else:
        print("Histos {}_Comparison drawn but not saved!".format(h1.GetName()))