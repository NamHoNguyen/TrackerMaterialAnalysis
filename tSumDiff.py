import ROOT
import math
import time
ROOT.gStyle.SetPalette(51)
ROOT.gStyle.SetOptFit()
#TColor.InvertPalette()

prefix = 'July4_'
source = 'DATA'

fIn=ROOT.TFile(prefix+'MinBias3D_'+source+'_4binphi_20binpt.root')
fOut=ROOT.TFile(prefix+'tSumDiff_'+source+'_4binphi_20binpt.root', 'RECREATE')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

triplet_old = ''

start_time = time.time()
# loop over histograms
#allKey = fIn.GetListOfKeys()
allKey = [x.GetName() for x in fIn.GetListOfKeys()]
pt2_string = [fIn.Get(x).GetName().split('_')[2] for x in allKey if ('TEC+1' in x) and ('minus' in x)]

for pt2 in pt2_string: 
    th3m = [fIn.Get(x) for x in allKey if (pt2 in x) and ('minus' in x)]
    th3pl = [fIn.Get(x) for x in allKey if (pt2 in x) and ('plus' in x)]
    for i in range(1, th3m[0].GetNbinsX()+1):
        for j in range(1, th3m[0].GetNbinsY()+1):
            tSum=ROOT.TH1D("tSum_{i}_{j}_{k}".format(i=i,j=j,k=float(pt2)),"tSum_{i}_{j}_{k}".format(i=i,j=j,k=float(pt2)), 37, 0, 37)
            tDiff=ROOT.TH1D("tDiff_{i}_{j}_{k}".format(i=i,j=j,k=float(pt2)),"tDiff_{i}_{j}_{k}".format(i=i,j=j,k=float(pt2)), 37, 0, 37)
            for k in range(len(th3m)):
                th1m = th3m[k].ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                th1pl = th3pl[k].ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

                tsum = th1pl.GetMean() + th1m.GetMean()
                tdiff = th1pl.GetMean() - th1m.GetMean()
                terr = math.sqrt(th1pl.GetMeanError()*th1pl.GetMeanError() + th1m.GetMeanError()*th1m.GetMeanError())
                    
                tSum.SetBinContent(k,tsum)
                tSum.SetBinError(k,terr)
                tDiff.SetBinContent(k,tdiff)
                tDiff.SetBinError(k,terr)

                tSum.GetXaxis().SetTitle("Detectors")
                tSum.GetYaxis().SetTitle("Sum of peaks")
                tDiff.GetXaxis().SetTitle("Detectors")
                tDiff.GetYaxis().SetTitle("Difference of peaks")
            tSum.Write()
            tDiff.Write()
    #if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
