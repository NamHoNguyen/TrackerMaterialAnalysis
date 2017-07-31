'''
Author: Ho Nam Nguyen, Elisabetta Manca
Purpose:
- This script takes in a 3D-sagitta file and calculates corresponding radiation length
- The result is binned in local eta, local phi, and track eta
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
prefix = 'July5_allfiles_6binPhi_10binPt2_41binTrackEta_'
source = 'DATA'

# Important variables
nphi,phimin,phimax = [6,-3.15,3.15]
npt2,pt2min,pt2max = [10,0.5,2.25]
ntrackEta,trackEtamin,trackEtamax = [41,-4.,4.]


fIn=ROOT.TFile(prefix+'MinBias3D_'+source+'.root')
fOut=ROOT.TFile(prefix+'RadLength3D_'+source+'.root', 'RECREATE')

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
    triplet = '_'.join(keyArray[:-2])
    #print triplet
    # When switch region or at the end of 'minus' regions
    if ((triplet_old!=triplet) and (triplet_old!='')) or (key == allKey[len(allKey)/2-1]):
        '''
        try:
            if (triplet_old.split('_')[1]!='pixel2'):
                triplet_old = triplet
                continue
        except:
            if (triplet_old!='pixel2'):
                triplet_old = triplet
                continue
        '''
        
        th3m_key = [x for x in allKey if (triplet_old in x) and ('minus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        th3pl_key = [x for x in allKey if (triplet_old in x) and ('plus' in x) and (x.split('_')[-1].replace('.','').isdigit())]
        # define TH3 to fill
        slope=ROOT.TH3D("slope"+triplet_old, "slope"+triplet_old, 41, -4, 4, nphi, phimin, phimax, ntrackEta, trackEtamin, trackEtamax)
        radlength=ROOT.TH3D("radlength_"+triplet_old,"radlength_"+triplet_old, 41,-4,4, nphi,phimin,phimax, ntrackEta,trackEtamin,trackEtamax)

        TrackEta = []
        for y in th3m_key:
            tEta = y.split('_')[-2]
            if tEta not in TrackEta:
                TrackEta.append(tEta)
        th3m_key_new = {}
        th3pl_key_new = {}
        th3m = {}
        th3pl = {}
        for tEta in TrackEta:
            th3m[tEta] = [fIn.Get(x) for x in th3m_key if (tEta in x)]
            th3pl[tEta] = [fIn.Get(x) for x in th3pl_key if (tEta in x)]
            th3m_key_new[tEta] = [x for x in th3m_key if (tEta in x)]
            th3pl_key_new[tEta] = [x for x in th3pl_key if (tEta in x)]
        for i in range(1, th3m[list(th3m)[0]][0].GetNbinsX()+1):
            for j in range(1, th3m[list(th3m)[0]][0].GetNbinsY()+1):
                for tEta in th3m:
                    rms=ROOT.TH1D("{i}_{j}_{k}".format(i=i,j=j,k=float(tEta)),"{i}_{j}_{k}".format(i=i,j=j,k=float(tEta)),npt2,pt2min, pt2max)
                    # loop over pt
                    for k in range(len(th3m[tEta])):
                        th1m = ROOT.TH1
                        th1pl = ROOT.TH1
                        th1m = th3m[tEta][k].ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                        # check if the corresponding plus histogram exist
                        try:
                            kpl = th3pl_key_new[tEta].index(th3m_key_new[tEta][k].replace('sag3Dminus_','sag3Dplus_'))
                            #if k!=kpl: print 'Diff k--------------------------------------'
                        except:
                            #print '----------------Plus missing---------',th3m_key_new[tEta][k]
                            continue
                        th1pl = th3pl[tEta][kpl].ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

                        pt = float([p for p in th3m[tEta][k].GetName().split('_') if not 'minus' in p][-1])

                        mean=(1./2*(th1m.GetRMS()+th1pl.GetRMS()))
                        mean_err=mean*math.sqrt(th1m.GetRMSError()*th1m.GetRMSError()+th1pl.GetRMSError()*th1pl.GetRMSError())
                        
                        if abs(mean*mean)>1.e-10:
                        #if abs(mean*mean)>1.e-6:
                            rms.SetBinContent(rms.FindBin(pt), mean*mean)
                            rms.SetBinError(rms.FindBin(pt), mean_err)
                            #if abs(mean*mean)<1.e-6:
                                #print triplet,i,j
                                
                    if(rms.GetEntries()<7):  continue
                    fit=ROOT.TF1("fit", "pol1");
                    rms.Fit(fit, "Q");

                    #slope.SetBinContent(i,j,fit.GetParameter(1))
                    #slope.SetBinError(i,j,fit.GetParError(1))
                    
                    if (fit.GetParameter(0) < 0): 
                        eta = radlength.GetXaxis().GetBinCenter(i)
                        phi = radlength.GetYaxis().GetBinCenter(j)
                        print 'Negative RadLength: ',eta,phi,fit.GetParameter(0),rms.GetEntries()
                        #rms.Draw()
                        #rms.GetYaxis().SetRangeUser(0,12e-6)
                        #canv.SaveAs(prefix+'negativeTest_'+triplet_old+'_eta'+str(round(eta,2))+'_phi'+str(round(phi,2))+'.png')
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
                    radlength.SetBinContent(i,j,radlength.GetZaxis().FindBin(float(tEta)),fit.GetParameter(0)/0.013/0.013)
                    radlength.SetBinError(i,j,radlength.GetZaxis().FindBin(float(tEta)),fit.GetParError(0)/0.013/0.013)
                    radlength.GetXaxis().SetTitle("local #eta")
                    radlength.GetYaxis().SetTitle("local #phi")
                    radlength.GetZaxis().SetTitle("track #eta")
        #radlength.Draw('colz')                                      
        #canv.SaveAs('plots/'+prefix+'RadLength2D_'+source+'_'+triplet_old+'.png')
        radlength.Write()
    triplet_old = triplet
    #if (key == allKey[10]): break
print("--- %s seconds ---" % (time.time() - start_time))
