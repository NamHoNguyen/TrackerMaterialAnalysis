#ifndef __CINT__
#endif
#include <string>
#include <sstream>
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TLine.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <vector>
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRint.h"
#include "TString.h"
#include "TSystem.h"
#include "Math/IFunction.h"
#include <cmath>
#include "TStopwatch.h"
#include "THnSparse.h"

int countIndex;
// Global definitions
std::map<std::string, TH1F*> h_1d;
std::map<std::string, TH1F*> h_1dm;
std::map<std::string, TH1F*> h_1dpl;
std::map<std::string, TH2F*> h_2d;
std::map<std::string, TH2F*> h_2dpl;
std::map<std::string, TH2F*> h_2dm;
std::map<std::string, TH2F*> h_2dptpl;
std::map<std::string, TH2F*> h_2dptm;
std::map<std::string, TH2F*> h_2dxpl;
std::map<std::string, TH2F*> h_2dxm;
std::map<std::string, TH2F*> h_2dxptpl;
std::map<std::string, TH2F*> h_2dxptm;
std::map<std::string, TH2F*> h_2dp;
std::map<std::string, TH2F*> h_2dp2;
//std::map<std::string, TH3F*> h_3dpl;
//std::map<std::string, TH3F*> h_3dm;
//std::map<std::string, TH3F*> h_3dpl_phi;
//std::map<std::string, TH3F*> h_3dm_phi;
//std::map<std::string, TH3F*> h_3dbm;
//std::map<std::string, TH3F*> h_3dbpl;
std::map<std::string, THnSparse*> h_ndm;
std::map<std::string, THnSparse*> h_ndpl;

// Compute the sagitta given the radii and the phi coordinate of three point in the tracker
float calculSag(float r1, float r2, float r3, float phi1, float phi2, float phi3){

	float a[2]={0};
	float b[2]={0};

	a[0]=r2*TMath::Cos(phi2)-r1*TMath::Cos(phi1);
	a[1]=r2*TMath::Sin(phi2)-r1*TMath::Sin(phi1);

	b[0]=r3*TMath::Cos(phi3)-r1*TMath::Cos(phi1);
	b[1]=r3*TMath::Sin(phi3)-r1*TMath::Sin(phi1);

	float sin_theta=(a[0]*b[1]-a[1]*b[0])/TMath::Sqrt(b[0]*b[0]+b[1]*b[1])/TMath::Sqrt(a[0]*a[0]+a[1]*a[1]);
	float a_par=TMath::Sqrt(a[0]*a[0]+a[1]*a[1])*TMath::Cos(TMath::ASin(sin_theta));

	if(a[0]*a[1]!=0&&b[0]*b[1]!=0) return sin_theta/(TMath::Sqrt(b[0]*b[0]+b[1]*b[1])-a_par);
	else return -99;
}

// Compute the squared distance from point 1 to point 3
float calcb2(float r1, float r2, float r3, float phi1, float phi2, float phi3){

	float b[2]={0};

	b[0]=r3*TMath::Cos(phi3)-r1*TMath::Cos(phi1);
	b[1]=r3*TMath::Sin(phi3)-r1*TMath::Sin(phi1);

	float cst=(b[0]*TMath::Cos(phi2)+b[1]*TMath::Sin(phi2))/TMath::Sqrt(b[0]*b[0]+b[1]*b[1]);
	return (b[0]*b[0]+b[1]*b[1])*cst;
}

// Compute sin given eta
float calcSin(float eta){

	float exp=TMath::Exp(-eta);

	return 2*exp/(1+exp*exp);
}

// Compute cos given eta
float calcCos(float eta){

	float exp=TMath::Exp(-eta);

	return (1-exp*exp)/(1+exp*exp);
}

// Plotting functions
void plot1D(std::string title, float xval, double weight, std::map<std::string, TH1F*> &allhistos,
		int numbinsx, float xmin, float xmax)
{
	std::map<std::string, TH1F*>::iterator iter= allhistos.find(title);
	if(iter == allhistos.end()) //no histo for this yet, so make a new one

	{
		TH1F* currentHisto= new TH1F(title.c_str(), title.c_str(), numbinsx, xmin, xmax);
		currentHisto->Fill(xval, weight);
		allhistos.insert(std::pair<std::string, TH1F*> (title,currentHisto) );
	}
	else // exists already, so just fill it

	{
		(*iter).second->Fill(xval, weight);

	}
}

void plot2D(std::string title, float xval, float yval, double weight, std::map<string, TH2F*> &allhistos,
		int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax)
{

	std::map<std::string, TH2F*>::iterator iter= allhistos.find(title);
	if(iter == allhistos.end()) //no histo for this yet, so make a new one

	{
		TH2F* currentHisto= new TH2F(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax);
		currentHisto->Fill(xval, yval, weight);
		allhistos.insert(std::pair<std::string, TH2F*> (title,currentHisto) );
	}
	else // exists already, so just fill it

	{
		(*iter).second->Fill(xval, yval, weight);
	}

	return;

}

void plot3D(std::string title, float xval, float yval, float zval, double weight, std::map<std::string, TH3F*> &allhistos,
		int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax, int numbinsz, float zmin, float zmax)
{

	std::map<std::string, TH3F*>::iterator iter= allhistos.find(title);
	if(iter == allhistos.end()) //no histo for this yet, so make a new one

	{
	  //cout << "Write new histo: " << title << endl;
		TH3F* currentHisto= new TH3F(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax, numbinsz, zmin, zmax);
		currentHisto->Fill(xval, yval, zval, weight);
		allhistos.insert(std::pair<std::string, TH3F*> (title,currentHisto) );
	}
	else // exists already, so just fill it

	{
		(*iter).second->Fill(xval, yval, zval, weight);
	}

	return;

}
void plotND(std::string title, Double_t val[], double weight, std::map<std::string, THnSparse*> &allhistos,
	    Int_t bin[], Double_t valmin[], Double_t valmax[], int dim)
{
  //cout << "In plotND " << endl;
  //for (std::map<std::string,THnSparse*>::iterator it=allhistos.begin(); it!=allhistos.end(); ++it) std::cout << it->first << ",";
  std::map<std::string, THnSparse*>::iterator iter= allhistos.find(title);
  if(iter == allhistos.end()) //no histo for this yet, so make a new one
    
    {
      //for (std::map<std::string,THnSparse*>::iterator it=allhistos.begin(); it!=allhistos.end(); ++it) std::cout << it->first << ",";
      //cout << endl;
      cout << "Write new histo: " << title << endl;
      //THnSparse* currentHisto;
      //currentHisto= new THnSparseF(title.c_str(), title.c_str(), dim, bin, valmin, valmax);
      THnSparse* currentHisto= new THnSparseF(title.c_str(), title.c_str(), dim, bin, valmin, valmax);
      currentHisto->Sumw2();
      currentHisto->Fill(val, weight);
      allhistos.insert(std::pair<std::string, THnSparse*> (title,currentHisto) );
    }
  else // exists already, so just fill it
    
    {
      //cout << "filling in " << title << endl;
      (*iter).second->Fill(val, weight);
    }
  
  int countIndex = 0;
  cout << "------------------------------------------------In plot ND " << endl;
  for (std::map<std::string,THnSparse*>::iterator it=allhistos.begin(); it!=allhistos.end(); ++it) countIndex+=1;
  cout << countIndex << " elements, " << allhistos.begin()->first << endl;
  countIndex = 0;
  for (std::map<std::string,THnSparse*>::iterator it=h_ndpl.begin(); it!=h_ndpl.end(); ++it) countIndex+=1;
  if (countIndex!=0) cout << countIndex << " elements, " << h_ndpl.begin()->first << endl;
  return;
  
}
