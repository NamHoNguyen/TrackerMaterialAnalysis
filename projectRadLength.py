'''
Author: Ho Nam Nguyen
Purpose: 
- This script projects DATA and MC 2D-RadLength onto phi-direction for fitting.
- It plots DATA and MC together and separately with the fits
- It only saves plots with high chi2 or high fluctuation amplitude or zero-ndf for further analysis
'''

import matplotlib.pyplot as plt
import ROOT
import time
import sys
from ROOT import TColor
ROOT.gStyle.SetPalette(51)

prefixMC = 'July13_3files_6binPhi_10binPt2_'
#source = 'MC'
prefix = 'July5_allfiles_6binPhi_10binPt2_'
#source = 'DATA'

fInMC=ROOT.TFile(prefixMC+'RadLength2D_MC.root')
fIn=ROOT.TFile(prefix+'RadLength2D_DATA.root')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

start_time = time.time()
# loop over histograms
allKey = [x.GetName() for x in fIn.GetListOfKeys()]
chi2max = 0.

ikey = 0
for key in allKey:
    triplet = key.replace('radlength','')
    #if triplet != 'pixel2': continue

    th2 = fIn.Get(key)
    th2MC = fInMC.Get(key)
    data = np.zeros([th2.GetNbinsX(),3])
    # loop over eta
    for i in range(1, th2.GetNbinsX()+1):
        ROOT.gStyle.SetOptFit()
        th1 = th2.ProjectionY("proj_{i}".format(i=i), i,i)
        th1.SetTitle(triplet+'_{i}'.format(i=i))
        # No data
        if th1.GetEntries() < 1.e-6: continue
        # Too few data points
        nFilledBins = 0
        for j in range(1,th2.GetNbinsY()+1):
            if abs(th1.GetBinContent(j)) > 1.e-6: nFilledBins +=1
        if nFilledBins < 5: continue

        # for MC
        th1MC = th2MC.ProjectionY("projMC_{i}".format(i=i), i,i)
        th1MC.SetTitle(triplet+'MC_{i}'.format(i=i))
        # No data
        if th1MC.GetEntries() < 1.e-6: continue
        # Too few data points
        nFilledBins = 0
        for j in range(1,th2MC.GetNbinsY()+1):
            if abs(th1MC.GetBinContent(j)) > 1.e-6: nFilledBins +=1
        if nFilledBins < 5: continue
 

        ######## FIT ############
        fit = ROOT.TF1("fit","[0]*(1. + [1]*cos(x+[2]) + [3]*cos(2.*x+[4]))",-3.15,3.15);
        th1.Fit(fit,"Q0")
        data[i-1,:] = [th2.GetXaxis().GetBinCenter(i),fit.GetChisquare(),fit.GetParameter(1)]

        # Filter only chi2/ndf > 70 or amp > 10%
        if fit.GetNDF() != 0:
            if (fit.GetChisquare()/float(fit.GetNDF()) < 70.) and (abs(fit.GetParameter(1))<0.1) and (abs(fit.GetParameter(3))<0.1):  continue

        fitMC = ROOT.TF1("fitMC","[0]*(1. + [1]*cos(x+[2]) + [3]*cos(2.*x+[4]))",-3.15,3.15);
        th1MC.Fit(fitMC,"Q0")
        #print '---------------MC-----------------> chi2/ndf = ',round(fitMC.GetChisquare(),1),'/',fitMC.GetNDF()

        # Plot DATA and MC together
        ROOT.gStyle.SetOptFit(0000)
        canv0 = ROOT.TCanvas("c0","Canvas0",1000,500)
        canv0.Divide(2,1)
        canv0.cd(1)
        th1.GetYaxis().SetRangeUser(0,fit.GetParameter(0)*1.5)
        th1MC.GetYaxis().SetRangeUser(0,fitMC.GetParameter(0)*1.5)

        th1.Draw()
        th1MC.SetLineColor(3)
        th1MC.Draw('same')

        canv0.SaveAs('July13_filtered_'+triplet+'_test'+str(i)+'_combined.png')
        canv0.Close()
        
        # Plot DATA and MC separately with fit
        ROOT.gStyle.SetOptFit(111)
        th1.GetFunction('fit').ResetBit(ROOT.TF1.kNotDraw);
        th1MC.GetFunction('fitMC').ResetBit(ROOT.TF1.kNotDraw);

        canv = ROOT.TCanvas("c1","Canvas",1000,500)
        canv.Divide(2,1)

        canv.cd(1)
        #th1.GetYaxis().SetRangeUser(0,fit.GetParameter(0)*1.5)
        #th1.GetYaxis().SetRangeUser(0,0.1)
        th1.Draw()
        
        canv.cd(2)
        #th1MC.SetLineColor(3)
        #th1MC.GetYaxis().SetRangeUser(0,fitMC.GetParameter(0)*1.5)
        th1MC.Draw()
        
        canv.SaveAs('July13_filtered_'+triplet+'_test'+str(i)+'.png')
        canv.Close()
    ikey+=1
    #if ikey > 2: sys.exit()
    continue

    # Plot Chi2
    plt.plot(data[:,0],data[:,1],'*')
    plt.xlabel('$\eta$')
    plt.ylabel('$\chi^2$')
    plt.savefig('plots/July10_Chi2VsEta_'+triplet+'.png')
    plt.gcf().clear()

    plt.plot(data[:,1],data[:,2],'*')
    plt.xlabel('$\chi^2$')
    plt.ylabel('Amplitude')
    plt.savefig('plots/July10_AmpVsChi2_'+triplet+'.png')
    plt.gcf().clear()

print("--- %s seconds ---" % (time.time() - start_time))
