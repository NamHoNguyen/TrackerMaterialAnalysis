#include "MinBias.h"

// Define the observables that we need in the analysis
int n=0;
int stereo[35]={0};
float zUnc[35]={0};
float rUnc[35]={0};
float phiUnc[35]={0};
float z[35]={0};
float eta[35]={0};
float phi[35]={0};
float r[35]={0};
float pt[35]={0};
float trackEta=0;
float trackPhi=0;
float trackZ0=0;
float trackY0=0;
float trackX0=0;
float trackPt=-99;
float trackPtErr=-99;
int layer[35]={-99};
int detector[35]={0};
float sag=-99;
int n_pixel=0;
string prefix="June22_";

void MinBiasAnalysis(TTree *Tree, bool isFirstRun, bool isMC){

	int nEntries = Tree->GetEntries();

	string position[100]={""};
	string region[100]={""};
	string local_eta[100]={""};
	string triplet="";
	string tripl_name="";
	int i_next=0;
	int i_prev=0;
	int number[100]={0};
	double eta_max[100]={-99};
	double eta_min[100]={-99};
	float sin_track=-99;
	int sel_tracks=0;
	float mean_min=0;
	float mean_pl=0;
	int in_plots=0;

	TFile *fileIN;
	if(!isFirstRun){
		if(isMC) fileIN = TFile::Open((prefix+"MinBiasMC_sag1D.root").c_str());
		else fileIN = TFile::Open((prefix+"MinBiasDATA_sag1D.root").c_str());
	}

	// Loop over tracks
	for(int i = 1; i < nEntries; i++) {

	  if(i%100000==0) cout <<"Analyzed entry "<< i <<"/"<< nEntries <<" selected tracks "<< sel_tracks<<" entering in plots "<<in_plots<<endl;

		Tree->GetEntry(i);
		n_pixel=0;

		// Quality cuts
		if(trackPt>1.5) continue;
		if(trackPt<0.75) continue;
		if(n>100) continue;
		if(n<14) continue;

		trackPtErr=trackPtErr*trackPt*trackPt; // Since Mike defined it wrong in the ntuples :)
		if(trackPtErr>0.01) continue;

		sel_tracks++;

		// Loop over hits
		for(int j = 0; j < n; j++) {

			if(stereo[j]==1) continue; // avoid stereo modules

			// identify layers with strings instead of numbers
			if(detector[j]==0){
				position[j]="_pixel"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j];
			}
			if(detector[j]==1){
				position[j]="_pixeldisk"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_pixeldisk+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j]);
			}
			if(detector[j]==2){
				position[j]="_TIB"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j]+3;
			}
			if(detector[j]==3){
				position[j]="_TOB"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j]+7;
			}
			if(detector[j]==4){
				position[j]="_TID"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_TID+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j])+2;
			}
			if(detector[j]==5){
				position[j]="_TEC"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_TEC+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j])+5;
			}

			if(position[j]=="_pixel1"||position[j]=="_pixel2"||position[j]=="_pixel3"||position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk+2"||position[j]=="_pixeldisk-1"||position[j]=="_pixeldisk-2")
			n_pixel++;
		}

		if(n_pixel<2||n_pixel>3) continue;

		// Loop again to form the triplet strings excluding first and last hits
		for (int j=1; j<n-1; j++){

			if(stereo[j]==1) continue;
			i_next=0;
			i_prev=0;

			for(int k=0; k<n; k++){

				if(stereo[k]==1) continue;
				if(i_next*i_prev!=0) break;

				if(region[k]==region[j]){
					if(number[j]==number[k]+1) i_prev=k;
					else if(number[j]==number[k]-1) i_next=k;
				}

				if(position[j]=="_TIB1"||position[j]=="_TIB2"||position[j]=="_TIB3"||position[j]=="_TIB4"){
					if(position[k]=="_TID+1"||position[k]=="_TID-1") i_next=k;
				}
				if(position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk-1"){
					if(position[k]=="_pixel2"||position[k]=="_pixel1"||position[k]=="_pixel3") i_prev=k;
				}

				if(position[j]=="_TID+1"||position[j]=="_TID-1"){
					if(position[k]=="_TIB1"||position[k]=="_TIB2"||position[k]=="_TIB3"||position[k]=="_TIB4") i_prev=k;
					//if(position[k]=="_TIB1"||position[k]=="_TIB2"||position[k]=="_TIB3") i_prev=k;
				}
				if(position[j]=="_TEC+1"||position[j]=="_TEC-1"){
					if(position[k]=="_TOB1"||position[k]=="_TOB2"||position[k]=="_TOB3"||position[k]=="_TOB4"||position[k]=="_TOB5") i_prev=k;
				}

				if(position[j]=="_pixeldisk+2"||position[j]=="_pixeldisk-2"){
					if(position[k]=="_TIB1") i_next=k;
				}
				if(position[j]=="_TID+2"||position[j]=="_TID-2"){
					if(position[k]=="_TOB1") i_next=k;
				}
				if(position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk-1"){
					if(position[k]=="_TIB1") i_next=k;
				}
				if(position[j]=="_TID+1"||position[j]=="_TID-1"){
					if(position[k]=="_TOB1") i_next=k;
				}
			}

			if(i_next*i_prev==0) continue;

			tripl_name=position[i_prev]+position[j]+position[i_next];
			triplet=position[j];
			
			if(position[j]=="_pixel2") in_plots++;

			// A bunch of plots
			plot1D("pt", trackPt, 1, h_1d, 1000, 0.5, 2);
			plot1D("local pt"+triplet, pt[j], 1, h_1d, 1000, 0.5, 2);
			plot1D("eta", trackEta, 1, h_1d, 1000, -5, 5);
			plot1D("local eta"+triplet, eta[j], 1, h_1d, 160, -2.0, 2.0);

			// Correct for incidence angle
			if(region[j]=="_barrel") sin_track=calcSin(trackEta);
			else sin_track=TMath::Abs(calcCos(trackEta));

			float b2=calcb2(r[i_prev],r[j],r[i_next], phi[i_prev],phi[j],phi[i_next])*sin_track;
			sag=calculSag(r[i_prev],r[j],r[i_next], phi[i_prev],phi[j],phi[i_next]);

			plot2D("b2vspt", trackPt*trackPt, sqrt(b2), 1, h_2d, 100, 0.5, 2.25, 500, 0., 50);

			TString sagmin_name="sagmin"+triplet;
			TString sagpl_name="sagpl"+triplet;

			if(!isFirstRun){
				TH1 *sagmean_min=(TH1*)fileIN->Get(sagmin_name);
				TH1 *sagmean_pl=(TH1*)fileIN->Get(sagpl_name);

				if(!sagmean_min) continue;
				if(!sagmean_pl) continue;

				mean_min=sagmean_min->GetMean();
				mean_pl=sagmean_pl->GetMean();
			}
			if(triplet!=""){
				if(isFirstRun){
					if(sag<0) plot1D("sagmin"+triplet, pt[j]*sag, 1, h_1dm, 600, -0.01, 0);
					else plot1D("sagpl"+triplet, pt[j]*sag, 1, h_1dpl, 600, 0, 0.01);
				}
				else{
					if(sag<0) plot1D("sagminus_new"+triplet, ((pt[j]*sag)-mean_min)*sqrt(b2), 1, h_1dm, 600, -0.1, 0.1);
					else plot1D("sagplus_new"+triplet, ((pt[j]*sag)-mean_pl)*sqrt(b2), 1, h_1dpl, 600, -0.1, 0.1);

					if(sag<0) plot3D("sag3Dminus"+triplet, eta[j], pt[j]*pt[j], ((pt[j]*sag)-mean_min)*sqrt(b2), 1, h_3dm, 41, -4, 4, 10, 0.5, 2.25, 600, -0.01, 0.01);
					else plot3D("sag3Dplus"+triplet, eta[j], pt[j]*pt[j], ((pt[j]*sag)-mean_pl)*sqrt(b2), 1, h_3dpl, 41, -4, 4, 10, 0.5, 2.25, 600, -0.01, 0.01);
				}

				triplet="";
				tripl_name="";
				sagmin_name="sagminus";
				sagpl_name="sagplus";
			}
		}
		for (int j=0; j<n; j++){
			position[j]="";
			r[j]=0;
			phi[j]=0;
			sag=0;
		}

	}
}

void MinBias(bool isFirstRun = true, bool isMC = true, TString dirname="root://eoscms//eos/cms//store/cmst3/user/bachtis/ZeroBias/crab_PointsDATA/160606_201244/0000/", const char *ext=".root")
{
	if(isMC) dirname="root://eoscms//eos/cms//store/cmst3/user/bachtis/MinBias_TuneMBReps08_13TeV-pythia8/crab_PointsMC/160607_144658/0000/";
	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		int i=0;
		TChain *Tree = new TChain("tree");
		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				if(!isMC&&i>2) break; //set the number of files to run
				i++;
				//cout << fname.Data() << endl;
				TString filename;
				if(!isMC) filename ="root://eoscms//eos/cms//store/cmst3/user/bachtis/ZeroBias/crab_PointsDATA/160606_201244/0000/";
				else filename="root://eoscms//eos/cms//store/cmst3/user/bachtis/MinBias_TuneMBReps08_13TeV-pythia8/crab_PointsMC/160607_144658/0000/";
				filename=filename+fname.Data();
				cout << filename<< endl;
				Tree->AddFile(filename);

			}
		}

		Tree->SetBranchAddress("n" ,&n);
		Tree->SetBranchAddress("zUnc" ,zUnc);
		Tree->SetBranchAddress("rUnc" ,rUnc);
		Tree->SetBranchAddress("phiUnc" ,phiUnc);
		Tree->SetBranchAddress("z" ,z);
		Tree->SetBranchAddress("eta" ,eta);
		Tree->SetBranchAddress("phi" ,phi);
		Tree->SetBranchAddress("r" ,r);
		Tree->SetBranchAddress("pt" ,pt);
		Tree->SetBranchAddress("stereo" , stereo);
		Tree->SetBranchAddress("layer" , layer);
		Tree->SetBranchAddress("detector" , detector);
		Tree->SetBranchAddress("trackEta" ,&trackEta);
		Tree->SetBranchAddress("trackPhi" ,&trackPhi);
		Tree->SetBranchAddress("trackPt" ,&trackPt);
		Tree->SetBranchAddress("trackPtErr" ,&trackPtErr);
		Tree->SetBranchAddress("trackZ0" ,&trackZ0);
		Tree->SetBranchAddress("trackX0" ,&trackX0);
		Tree->SetBranchAddress("trackY0" ,&trackY0);


		MinBiasAnalysis(Tree, isFirstRun, isMC);


	}


	std::map<std::string, TH1F*>::iterator it1d;
	std::map<std::string, TH1F*>::iterator it1dm;
	std::map<std::string, TH1F*>::iterator it1dpl;
	std::map<std::string, TH3F*>::iterator it3d;
	std::map<std::string, TH3F*>::iterator it3dpl;
	std::map<std::string, TH3F*>::iterator it3d2;

	// Writes 1D histos
	if(isFirstRun){
		TFile *fileOUT;
		if(isMC) fileOUT = new TFile((prefix+"MinBiasMC_sag1D.root").c_str(),"RECREATE");
		else fileOUT = new TFile((prefix+"MinBiasDATA_sag1D.root").c_str(),"RECREATE");
	}
	else{
	TFile *fileOUT;
	if(isMC) fileOUT = new TFile((prefix+"MinBiasMC_sag1D_new.root").c_str(),"RECREATE");
	else fileOUT = new TFile((prefix+"MinBiasDATA_sag1D_new.root").c_str(),"RECREATE");
}

for(it1dm=h_1dm.begin(); it1dm!=h_1dm.end(); it1dm++) {

	it1dm->second->Write();

	delete it1dm->second;
}

for(it1dpl=h_1dpl.begin(); it1dpl!=h_1dpl.end(); it1dpl++) {

	it1dpl->second->Write();
	delete it1dpl->second;
}

// Writes 3D histos

if(!isFirstRun){
	TFile *file3D;
	if(isMC) file3D= new TFile((prefix+"MinBias3D_MC.root").c_str(),"RECREATE");
	else file3D=new TFile((prefix+"MinBias3D_DATA.root").c_str(),"RECREATE");
}
for(it3d=h_3dm.begin(); it3d!=h_3dm.end(); it3d++) {

	it3d->second->Write();
	delete it3d->second;
}

for(it3d2=h_3dpl.begin(); it3d2!=h_3dpl.end(); it3d2++) {

	it3d2->second->Write();
	delete it3d2->second;
}

}
