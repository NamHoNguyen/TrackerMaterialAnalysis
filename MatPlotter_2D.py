import ROOT
import math
import time
from ROOT import TColor
#from ROOT import gStyle
canv = ROOT.TCanvas()
ROOT.gStyle.SetPalette(51)
ROOT.gStyle.SetOptFit()
#TColor.InvertPalette()

prefix = 'July4_'
source = 'DATA'

fIn=ROOT.TFile(prefix+'MinBias3D_'+source+'_4binphi_20binpt.root')
fOut=ROOT.TFile(prefix+'RadLength_'+source+'_2D_4binphi_20binpt.root', 'RECREATE')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

triplet_old = ''

start_time = time.time()
# loop over histograms
allKey = fIn.GetListOfKeys()
for key in allKey:
    # Just loop over 'minus' regions
    if not 'minus' in key.GetName(): continue
    print key.GetName()
    keyArray = [p for p in key.GetName().split('_') if not 'minus' in p]
    triplet = keyArray[0]
    # When switch region or at the end of 'minus' regions
    if ((triplet_old!=triplet) and (triplet_old!='')) or (key.GetName() == allKey[len(allKey)/2-1].GetName()):
        '''
        if (triplet_old!='pixel2'):
            triplet_old = triplet
            continue
        '''
        th3m = [fIn.Get(x.GetName()) for x in fIn.GetListOfKeys() if (triplet_old in x.GetName()) and ('minus' in x.GetName())]
        th3pl = [fIn.Get(x.GetName()) for x in fIn.GetListOfKeys() if (triplet_old in x.GetName()) and ('plus' in x.GetName())]
        #print th3m,th3pl
        # define TH2 to fill
        slope=ROOT.TH2D("slope"+triplet_old, "slope"+triplet_old, 41, -4, 4, 4, -3.15, 3.15)
        radlength=ROOT.TH2D("radlength"+triplet_old, "radlength"+triplet_old, 41, -4, 4, 4, -3.15, 3.15)

        for i in range(1, th3m[0].GetNbinsX()+1):
            for j in range(1, th3m[0].GetNbinsY()+1):
                #rms=ROOT.TH1D("{i}_{j}".format(i=i,j=j),"{i}_{j}".format(i=i,j=j), 10, 0.5, 2.25);
                rms=ROOT.TH1D("{i}_{j}".format(i=i,j=j),"{i}_{j}".format(i=i,j=j), 20, 0.5, 2.25);
                # loop over pt
                for k in range(len(th3m)):
                    pt = float([p for p in th3m[k].GetName().split('_') if not 'minus' in p][1])

                    th1m = ROOT.TH1
                    th1pl = ROOT.TH1
                    th1m = th3m[k].ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                    th1pl = th3pl[k].ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

                    mean=(1./2*(th1m.GetRMS()+th1pl.GetRMS()))
                    mean_err=mean*math.sqrt(th1m.GetRMSError()*th1m.GetRMSError()+th1pl.GetRMSError()*th1pl.GetRMSError())
                    
                    #if mean!=0.: print mean,mean_err
                    if abs(mean*mean)>1.e-10:
                    #if abs(mean*mean)>1.e-6:
                        rms.SetBinContent(rms.FindBin(pt), mean*mean)
                        rms.SetBinError(rms.FindBin(pt), mean_err)

                if(rms.GetEntries()<7): continue
                fit=ROOT.TF1("fit", "pol1");
                rms.Fit(fit, "Q");

                slope.SetBinContent(i,j,fit.GetParameter(1))
                slope.SetBinError(i,j,fit.GetParError(1))
                
                radlength.SetBinContent(i,j,fit.GetParameter(0)/0.013/0.013)
                if (fit.GetParameter(0) < 0): 
                    eta = radlength.GetXaxis().GetBinCenter(i)
                    phi = radlength.GetYaxis().GetBinCenter(j)
                    print 'Negative RadLength: ',eta,phi,fit.GetParameter(0),rms.GetEntries()
                    rms.Draw()
                    #rms.GetYaxis().SetRangeUser(0,12e-6)
                    canv.SaveAs('plots/'+prefix+'negativeTest_'+triplet_old+'_eta'+str(round(eta,2))+'_phi'+str(round(phi,2))+'.png')
                '''
                else:
                    eta = radlength.GetXaxis().GetBinCenter(i)
                    phi = radlength.GetYaxis().GetBinCenter(j)
                    if abs(eta+0.2)<1.e-2:
                    #if True:
                        rms.Draw()
                        #rms.GetYaxis().SetRangeUser(0,20e-6)
                        canv.SaveAs('plots/'+prefix+'positiveTest_'+triplet_old+'_eta'+str(round(eta,2))+'_phi'+str(round(phi,2))+'.png')
                '''
                radlength.SetBinError(i,j,fit.GetParError(0)/0.013/0.013)
                radlength.GetXaxis().SetTitle("local #eta")
                radlength.GetYaxis().SetTitle("local #phi")
        radlength.Draw('colz')                                      
        canv.SaveAs('plots/'+prefix+'RadLength2D_'+source+'_'+triplet_old+'.png')
        radlength.Write()
    triplet_old = triplet
    #if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
'''
start_time = time.time()
# loop over histograms
allKey = fIn.GetListOfKeys()
print allKey[len(allKey)/2-1].GetName()
for key in allKey:

    if not 'minus' in key.GetName(): continue
    print key.GetName()

    # take only the minus histograms fro the list
    th3=ROOT.TH3
    th3=fIn.Get(key.GetName())
    keyArray = [p for p in th3.GetName().split('_') if not 'minus' in p]
    triplet = keyArray[0]
    
    if ((triplet_old!=triplet) and (triplet_old!='')) or (key.GetName() == allKey[len(allKey)/2-1].GetName()):
        print [x.GetName() for x in fIn.GetListOfKeys() if (triplet_old in x.GetName()) and ('minus' in x.GetName())]
        # define TH2 to fill
        slope=ROOT.TH2D("slope"+triplet_old, "slope"+triplet_old, 41, -4, 4, 65, -3.2, 3.2)
        radlength=ROOT.TH2D("radlength"+triplet_old, "radlength"+triplet_old, 41, -4, 4, 65, -3.2, 3.2)

        for i in range(1, th3.GetNbinsX()+1):
            for j in range(1, th3.GetNbinsY()+1):
                rms=ROOT.TH1D("{i}_{j}".format(i=i,j=j),"{i}_{j}".format(i=i,j=j), 10, 0.5, 2.25);
                #print [x.GetName() for x in fIn.GetListOfKeys() if (triplet_old in x.GetName()) and ('minus' in x.GetName())]
                #continue
                for name in [x.GetName() for x in fIn.GetListOfKeys() if (triplet_old in x.GetName()) and ('minus' in x.GetName())]:
                    # take only the minus histograms fro the list
                    th3m = ROOT.TH3
                    th3m = fIn.Get(name)
                    keyArray = [p for p in th3.GetName().split('_') if not 'minus' in p]
                    # find the paired plus
                    th3pl=ROOT.TH3
                    th3pl=fIn.Get('sag3Dplus_'+'_'.join(keyArray))
                    
                    #print th3.GetName(), th3pl.GetName()
                    #continue
                    pt = float(keyArray[1])

                    th1m = ROOT.TH1
                    th1pl = ROOT.TH1
                    th1m = th3m.ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
                    th1pl = th3pl.ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

                    mean=(1./2*(th1m.GetRMS()+th1pl.GetRMS()))
                    mean_err=mean*math.sqrt(th1m.GetRMSError()*th1m.GetRMSError()+th1pl.GetRMSError()*th1pl.GetRMSError())

                    rms.SetBinContent(rms.FindBin(pt), mean*mean)
                    rms.SetBinError(rms.FindBin(pt), mean_err)

                if(rms.GetEntries()<7): continue
                fit=ROOT.TF1("fit", "pol1");
                rms.Fit(fit, "Q");
                
                slope.SetBinContent(i,j,fit.GetParameter(1))
                slope.SetBinError(i,j,fit.GetParError(1))
                
                radlength.SetBinContent(i,j,fit.GetParameter(0)/0.013/0.013)
                #radlength.SetBinError(i,j,fit.GetParError(0)/0.013/0.013)
            
                #radlength.GetXaxis().SetTitle("local #eta")
                #radlength.GetYaxis().SetTitle("x/X_{0}")
        radlength.Draw('colz')    
        canv.SaveAs('out.png')                       
        radlength.Write()
    triplet_old = triplet
    if (key.GetName() == allKey[30].GetName()): break
print("--- %s seconds ---" % (time.time() - start_time))
'''
