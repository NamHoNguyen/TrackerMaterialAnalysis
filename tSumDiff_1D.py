import ROOT
import math
import time
from ROOT import TColor
#from ROOT import gStyle
canv = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#TColor.InvertPalette()
import sys

prefix = 'July5_allfiles_'
source = 'DATA'

fIn=ROOT.TFile(prefix+'MinBias'+source+'_sag1D.root')
fOut=ROOT.TFile(prefix+'tSumDiff_'+source+'_1D.root', 'RECREATE')

triplet_old = ''

start_time = time.time()
allKey = [x.GetName() for x in fIn.GetListOfKeys()]
tSum=ROOT.TH1D("tSum", "tSum",37, 0, 37)
tDiff=ROOT.TH1D("tDiff", "tDiff", 37, 0, 37)
k = 1
# loop over histograms
for key in allKey:
    # Just loop over 'minus' regions
    if not 'min' in key: continue
    #print key
    triplet = [p for p in key.split('_') if not 'min' in p][0]
    th1m = fIn.Get(key)
    th1pl = fIn.Get('sagpl_'+triplet)
    
    tsum  = (th1pl.GetMean() + th1m.GetMean())/(th1pl.GetMean() - th1m.GetMean())
    tdiff = th1pl.GetMean() - th1m.GetMean()
    tsumerr = 2./(tdiff**2)*math.sqrt( (th1m.GetMean()*th1pl.GetMeanError())**2 + (th1pl.GetMean()*th1m.GetMeanError())**2 )
    tdifferr  = math.sqrt(th1pl.GetMeanError()**2 + th1m.GetMeanError()**2)

    tSum.SetBinContent(k,tsum)
    tSum.SetBinError(k,tsumerr)
    tSum.GetXaxis().SetBinLabel(k,triplet)
    tDiff.SetBinContent(k,tdiff)
    tDiff.SetBinError(k,tdifferr)
    tDiff.GetXaxis().SetBinLabel(k,triplet)
    k += 1
    #print tdiff,tdifferr
    
tSum.GetXaxis().SetTitle("Detectors")
tSum.GetYaxis().SetTitle("Sum/Difference of peaks")
tDiff.GetXaxis().SetTitle("Detectors")
tDiff.GetYaxis().SetTitle("Difference of peaks")
tDiff.GetYaxis().SetRangeUser(0,0.012)

tSum.Draw()
canv.SaveAs(prefix+'tSum_1D.png')
tDiff.Draw()
canv.SaveAs(prefix+'tDiff_1D.png')

#tSum.Write()
#tDiff.Write()
#if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
