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
k = 0
# loop over histograms
for key in allKey:
    # Just loop over 'minus' regions
    if not 'min' in key: continue
    print key
    triplet = [p for p in key.split('_') if not 'min' in p][0]
    th1m = fIn.Get(key)
    th1pl = fIn.Get('sagpl_'+triplet)
    
    tsum  = th1pl.GetMean() + th1m.GetMean()
    tdiff = th1pl.GetMean() - th1m.GetMean()
    terr  = math.sqrt(th1pl.GetMeanError()*th1pl.GetMeanError() + th1m.GetMeanError()*th1m.GetMeanError())

    tSum.SetBinContent(k,tsum)
    tSum.SetBinError(k,terr)
    tDiff.SetBinContent(k,tdiff)
    tDiff.SetBinError(k,terr)
    k += 1
    
tSum.GetXaxis().SetTitle("Detectors")
tSum.GetYaxis().SetTitle("Sum of peaks")
tDiff.GetXaxis().SetTitle("Detectors")
tDiff.GetYaxis().SetTitle("Difference of peaks")

tSum.Draw()
canv.SaveAs(prefix+'tSum_1D.pdf')
tDiff.Draw()
canv.SaveAs(prefix+'tDiff_1D.pdf')

tSum.Write()
tDiff.Write()
#if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
