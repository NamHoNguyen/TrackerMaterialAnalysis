import ROOT
import math
import time
from ROOT import TColor
#from ROOT import gStyle
canv = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(1)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#TColor.InvertPalette()
import sys

prefix = 'July4_'
source = 'DATA'

fIn=ROOT.TFile(prefix+'MinBias3D_'+source+'_4binphi_20binpt.root')
fOut=ROOT.TFile(prefix+'tSumDiff_'+source+'_4binphi_20binpt.root', 'RECREATE')


triplet_old = ''

start_time = time.time()
# loop over histograms
allKey = [x.GetName() for x in fIn.GetListOfKeys()]
for key in allKey:
    # Just loop over 'minus' regions
    if not 'minus' in key: continue
    print key
    triplet = [p for p in key.split('_') if not 'minus' in p][0]
    # When switch region or at the end of 'minus' regions
    if ((triplet_old!=triplet) and (triplet_old!='')) or (key == allKey[len(allKey)/2-1]):
        '''
        if (triplet_old!='pixel2'):
            triplet_old = triplet
            continue
        '''
        th3m = [fIn.Get(x) for x in allKey if (triplet_old in x) and ('minus' in x)]
        th3pl = [fIn.Get(x) for x in allKey if (triplet_old in x) and ('plus' in x)]
        # define TH2 to fill
        tSum=ROOT.TH2D("tSum"+triplet_old, "tSum"+triplet_old, 41, -4, 4, 4, -3.15, 3.15)
        tDiff=ROOT.TH2D("tDiff"+triplet_old, "tDiff"+triplet_old, 41, -4, 4, 4, -3.15, 3.15)

        for i in range(1, th3m[0].GetNbinsX()+1):
            for j in range(1, th3m[0].GetNbinsY()+1):
                # get the first histogram
                th1m = th3m[0].ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                th1pl = th3pl[0].ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                # add other histogram
                for k in range(1,len(th3m)):
                    th1m.Add(th3m[k].ProjectionZ("binMinus_{i}_{j}_{k}".format(i=i,j=j,k=pt), i, i, j, j))
                    th1pl.Add(th3pl[k].ProjectionZ("binPlus_{i}_{j}_{k}".format(i=i,j=j,k=pt), i, i, j, j))
                
                tsum  = th1pl.GetMean() + th1m.GetMean()
                tdiff  = th1pl.GetMean() - th1m.GetMean()
                terr = math.sqrt(th1pl.GetMeanError()*th1pl.GetMeanError() + th1m.GetMeanError()*th1m.GetMeanError())

                tSum.SetBinContent(i,j,tsum)
                tSum.SetBinError(i,j,terr)
                tDiff.SetBinContent(i,j,tdiff)
                tDiff.SetBinError(i,j,terr)

                tSum.GetXaxis().SetTitle("Detectors")
                tSum.GetYaxis().SetTitle("Sum of peaks")
                tDiff.GetXaxis().SetTitle("Detectors")
                tDiff.GetYaxis().SetTitle("Difference of peaks")
        tSum.Write()
        tDiff.Write()
    triplet_old = triplet
    #if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
