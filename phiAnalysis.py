import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math
import time
import sys
from ROOT import TColor
#from ROOT import gStyle
ROOT.gStyle.SetPalette(51)
#TColor.InvertPalette()


prefixMC = 'July13_3files_6binPhi_10binPt2_'
prefix = 'July5_allfiles_6binPhi_10binPt2_'
today = 'July17_'

fInMC=ROOT.TFile(prefixMC+'RadLength2D_MC.root')
fIn=ROOT.TFile(prefix+'RadLength2D_DATA.root')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptFit()

start_time = time.time()
# loop over histograms
allKey = [x.GetName() for x in fIn.GetListOfKeys()]

ikey = 0
diff_max = 0.
fOut = ROOT.TFile(today+'diffDATAMC_RadLength.root','RECREATE')
for key in allKey:
    triplet = key.replace('radlength','')
    #if triplet != 'pixel2': continue
    th2diff = ROOT.TH2D("diffDATAMC_"+triplet, "diffDATAMC_"+triplet, 41, -4, 4, 6, -3.15, 3.15)
    
    th2 = fIn.Get(key)
    th2MC = fInMC.Get(key)

    for i in range(1, th2.GetNbinsX()+1):
        # Check for empty bins
        nFilledBins = 0
        nFilledBinsMC = 0
        for j in range(1,th2.GetNbinsY()+1):
            if abs(th2.GetBinContent(i,j)) > 1.e-6: nFilledBins +=1
            if abs(th2MC.GetBinContent(i,j)) > 1.e-6: nFilledBinsMC +=1
        if nFilledBins < 5 or nFilledBinsMC < 5: continue

        for j in range(1, th2.GetNbinsY()+1):
            data = th2.GetBinContent(i,j)
            dataerr = th2.GetBinError(i,j)
            mc = th2MC.GetBinContent(i,j)
            mcerr = th2MC.GetBinError(i,j)
            if data*mc != 0:
                diff = abs(1. - mc/data) 
                th2diff.SetBinContent(i,j,diff )
                th2diff.SetBinError(i,j, np.sqrt( (mc/data**2 * dataerr)**2 + (1./data * mcerr)**2 ) )
                if diff_max < diff:
                    if triplet == 'TID+2':
                        print 'Known negative RadLength of TID+2, false max difference'
                    else:
                        diff_max =diff
                        print diff_max,triplet,i,j
    th2diff.GetXaxis().SetTitle("local #eta")
    th2diff.GetYaxis().SetTitle("local #phi")
    th2diff.SetTitle('(1-mc/data)_'+triplet)
    th2diff.Write()
    ikey+=1
    #if ikey>0: sys.exit()
#fOut.Close()
print 'Max difference percentage: ',diff_max

fOut = ROOT.TFile(today+'chi2.root','RECREATE')
ikey = 0
for key in allKey:
    triplet = key.replace('radlength','')
    #if triplet != 'pixel2': continue

    th2 = fIn.Get(key)
    th2MC = fInMC.Get(key)
    chi2 = ROOT.TH2D("chi2_"+triplet, "chi2_"+triplet, 41, -4, 4, 20, 0, 0.4)
    # loop over eta
    for i in range(1, th2.GetNbinsX()+1):
        # Project
        th1 = th2.ProjectionY("proj_{i}".format(i=i), i,i)
        th1.SetTitle(triplet+'_{i}'.format(i=i))
        # No data
        if th1.GetEntries() < 1.e-6: continue
        # Too few data points
        nFilledBins = 0
        for j in range(1,th2.GetNbinsY()+1):
            if abs(th1.GetBinContent(j)) > 1.e-6: nFilledBins +=1
        if nFilledBins < 5: continue

        ######## FIT ############
        # fit function
        fit = ROOT.TF1("fit","[0]*(1. + [1]*cos(x+[2]) + [3]*cos(2.*x+[4]))",-3.15,3.15);
        th1.Fit(fit,"Q")
        amp = abs(fit.GetParameter(1)) + abs(fit.GetParameter(3))
        chi2.SetBinContent(i,chi2.GetYaxis().FindBin(amp),fit.GetChisquare())

        if amp<0:
            if triplet == 'TID+2':
                print 'Known negative radlength for ', triplet, ', amp = ', amp
    ikey+=1
    chi2.GetXaxis().SetTitle("local #eta")
    chi2.GetYaxis().SetTitle("Amplitude")
    chi2.Write()
    #if ikey > 2: sys.exit()
print("--- %s seconds ---" % (time.time() - start_time))

