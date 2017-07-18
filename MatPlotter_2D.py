'''
Author: Ho Nam Nguyen, Elisabetta Manca
Purpose:
- This script takes in a 3D-sagitta file and calculates corresponding radiation length
- The result is binned in local eta and local phi
- It plots out fits with negative radiaton lengths
'''
import ROOT
import math
import time
import sys
from ROOT import TColor
#from ROOT import gStyle
canv = ROOT.TCanvas()
ROOT.gStyle.SetPalette(51)
ROOT.gStyle.SetOptFit()
#TColor.InvertPalette()

#prefix = 'July13_3files_6binPhi_10binPt2_'
#source = 'MC'
#prefix = 'July11_allfiles_6binPhi_10binPt2_'
prefix = 'July5_allfiles_6binPhi_10binPt2_'
source = 'DATA'

# Important variables
nphi,phimin,phimax = [6,-3.15,3.15]
npt2,pt2min,pt2max = [10,0.5,2.25]


fIn=ROOT.TFile(prefix+'MinBias3D_'+source+'.root')
fOut=ROOT.TFile(prefix+'RadLength2D_'+source+'.root', 'RECREATE')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

triplet_old = ''

start_time = time.time()
# loop over histograms
allKey = [x.GetName() for x in fIn.GetListOfKeys()]
for key in allKey:
    # Just loop over 'minus' regions
    if not 'minus' in key: continue
    print key
    keyArray = [p for p in key.split('_') if not 'minus' in p]
    triplet = '_'.join(keyArray[:-1])
    #print triplet
    # When switch region or at the end of 'minus' regions
    if ((triplet_old!=triplet) and (triplet_old!='')) or (key == allKey[len(allKey)/2-1]):
        '''
        try:
            if (triplet_old.split('_')[1]!='pixeldisk+1'):
                triplet_old = triplet
                continue
        except:
            if (triplet_old!='pixeldisk+1'):
                triplet_old = triplet
                continue
        '''
        
        th3m = [fIn.Get(x) for x in allKey if (triplet_old in x) and ('minus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        th3pl = [fIn.Get(x) for x in allKey if (triplet_old in x) and ('plus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        th3m_key = [x for x in allKey if (triplet_old in x) and ('minus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        th3pl_key = [x for x in allKey if (triplet_old in x) and ('plus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        # define TH2 to fill
        slope=ROOT.TH2D("slope"+triplet_old, "slope"+triplet_old, 41, -4, 4, nphi, phimin, phimax)
        radlength=ROOT.TH2D("radlength_"+triplet_old, "radlength_"+triplet_old, 41, -4, 4, nphi, phimin, phimax)

        for i in range(1, th3m[0].GetNbinsX()+1):
            for j in range(1, th3m[0].GetNbinsY()+1):
                rms=ROOT.TH1D("{i}_{j}".format(i=i,j=j),"{i}_{j}".format(i=i,j=j), npt2, pt2min, pt2max);
                #rms=ROOT.TH1D("{i}_{j}".format(i=i,j=j),"{i}_{j}".format(i=i,j=j), 20, 0.5, 2.25);
                # loop over pt
                for k in range(len(th3m)):
                    th1m = ROOT.TH1
                    th1pl = ROOT.TH1
                    th1m = th3m[k].ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                    # check if the corresponding plus histogram exist
                    try: kpl = th3pl_key.index(th3m_key[k].replace('sag3Dminus_','sag3Dplus_'))
                    except:
                        print 'No corresponding plus histogram'
                        continue
                    th1pl = th3pl[kpl].ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

                    pt = float([p for p in th3m[k].GetName().split('_') if not 'minus' in p][-1])

                    mean=(1./2*(th1m.GetRMS()+th1pl.GetRMS()))
                    mean_err=mean*math.sqrt(th1m.GetRMSError()*th1m.GetRMSError()+th1pl.GetRMSError()*th1pl.GetRMSError())
                    
                    if abs(mean*mean)>1.e-10:
                    #if abs(mean*mean)>1.e-6:
                        rms.SetBinContent(rms.FindBin(pt), mean*mean)
                        rms.SetBinError(rms.FindBin(pt), mean_err)
                        if abs(mean*mean)<1.e-6:
                            print triplet,i,j

                if(rms.GetEntries()<7):  continue
                fit=ROOT.TF1("fit", "pol1");
                rms.Fit(fit, "Q");

                slope.SetBinContent(i,j,fit.GetParameter(1))
                slope.SetBinError(i,j,fit.GetParError(1))
                
                if (fit.GetParameter(0) < 0): 
                    eta = radlength.GetXaxis().GetBinCenter(i)
                    phi = radlength.GetYaxis().GetBinCenter(j)
                    print 'Negative RadLength: ',eta,phi,fit.GetParameter(0),rms.GetEntries()
                    rms.Draw()
                    #rms.GetYaxis().SetRangeUser(0,12e-6)
                    canv.SaveAs(prefix+'negativeTest_'+triplet_old+'_eta'+str(round(eta,2))+'_phi'+str(round(phi,2))+'.png')
                '''
                else:
                    eta = radlength.GetXaxis().GetBinCenter(i)
                    phi = radlength.GetYaxis().GetBinCenter(j)
                    #if abs(eta+0.2)<1.e-2:
                    if True:
                        rms.Draw()
                        #rms.GetYaxis().SetRangeUser(0,20e-6)
                        canv.SaveAs(prefix+'positiveTest_'+triplet_old+'_eta'+str(round(eta,2))+'_phi'+str(round(phi,2))+'.png')
                '''     
                radlength.SetBinContent(i,j,fit.GetParameter(0)/0.013/0.013)
                radlength.SetBinError(i,j,fit.GetParError(0)/0.013/0.013)
                radlength.GetXaxis().SetTitle("local #eta")
                radlength.GetYaxis().SetTitle("local #phi")
        #radlength.Draw('colz')                                      
        #canv.SaveAs('plots/'+prefix+'RadLength2D_'+source+'_'+triplet_old+'.png')
        radlength.Write()
    triplet_old = triplet
    #if (key == allKey[10]): break
print("--- %s seconds ---" % (time.time() - start_time))
