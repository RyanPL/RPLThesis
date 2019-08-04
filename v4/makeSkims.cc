#include<stdlib.h>
#include<stdio.h>

#define makeSkims_cxx

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFormula.h"

#include "makeSkims.h"

#include "/home/ryanpl/Thesis/mycuts/allEg4KineCorsPass2noMemLeak.h"
#include "/home/ryanpl/Thesis/mycuts/ecCutsDevPass2.h"
#include "/home/ryanpl/Thesis/mycuts/vzCutsPass2.h"
//#include "/home/ryanpl/Thesis/mycuts/nwTranslnOsiCutPass2.h"
#include "/home/ryanpl/Thesis/mycuts/fiducialCutsPass2.h"
#include "/home/ryanpl/Thesis/mycuts/moreFidCuts.h"

using namespace std;

int main(int argc, char *argv[]){
	if(argc < 3){
		cout << "Usage: " << argv[0] << " <ENERGY LEVEL>" << " <METHOD NUMBER>";
		cout << "<RUN NUMBER>: singleRunLoop()" << endl;
		return 1;
	}
	makeSkims m;
	cout << "Initializing singleRunLoop(" << argv[1] << "," << argv[2] << ")..." << endl;
	m.singleRunLoop(argv[1], argv[2]);

}

void makeSkims::singleRunLoop(Char_t energylvl[10], Char_t runNum[10])
{
	Int_t Ebi; //0,1,2,3,4 for 1.0539, 1.339, 1.9889,2.2256,2.999 respectively
	if(strcmp(energylvl,"1p1")==0)
		Ebi = 0;
	else if(strcmp(energylvl,"1p3")==0)
		Ebi = 1;
	else if(strcmp(energylvl,"1p5")==0)
		Ebi = 1;
	else if(strcmp(energylvl,"2")==0)
		Ebi = 2;
	else if(strcmp(energylvl,"2p2")==0)
		Ebi = 3;	
	else if(strcmp(energylvl,"3")==0)
		Ebi = 4;
	else
	{
		cout << "ERROR: Invalid energy level" << endl;
		return;
	}

	ifstream binfile;
	binfile.open("/home/ryanpl/Thesis/common/Q2bins.dat");
	Char_t data[100];
	Int_t numlines = 0;
	string reader;
	while(getline(binfile, reader)) //reading number of lines
		numlines++;
	const Int_t size = numlines;
	const Int_t sizeMinusOne = size - 1;
	binfile.close();
	binfile.open("/home/ryanpl/Thesis/common/Q2bins.dat");
	while(true){ //skipping the title line
		binfile >> data;
		if(isdigit(data[0]))
			break;
	}
	Float_t edges[size];
	binfile >> data;
	edges[0] = atof (data);
	for(Int_t i = 1; i <= size; i++) {
		binfile >> data;
		edges[i] = atof (data);
		for(Int_t j = 0; j < 4; j++) { //skip averages, next line number, and min
			binfile >> data;
		}
	}
	binfile.close();

	Float_t qSq, wSq, w, v, phi, theta, Eb;
	Float_t EbgNominal[5] = {1.0539, 1.339, 1.9889, 2.256, 2.999};
	Char_t name[100], filename[100];
	const Float_t M = 0.938;
	const Float_t PI = 3.14159;
	const Int_t NUM_BINS = 75; //NOTE THAT THIS HAS BEEN CHANGED
	//cut stuff//
	Float_t thDclPosRad, phDclPosRad, phDclPosDeg, thVtxDeg;
	/////////////
	Float_t a[sizeMinusOne][NUM_BINS], sig[sizeMinusOne][NUM_BINS];//[q2 bin][w bin] 
	Float_t nplus[sizeMinusOne][NUM_BINS], nminus[sizeMinusOne][NUM_BINS], fcplus[sizeMinusOne], fcminus[sizeMinusOne], fcpcount[sizeMinusOne], fcmcount[sizeMinusOne];
	for(Int_t i = 0; i<NUM_BINS; i++){
		for(Int_t j = 0; j<sizeMinusOne; j++){
			a[j][i] = 0.0;
			sig[j][i] = 0.0;
			nplus[j][i] = 0.0;
			nminus[j][i] = 0.0;
		}
	}

	ifstream infile;
	sprintf(name, "/home/ryanpl/Thesis/common/fileList_nh3_%s.dat", energylvl);
	infile.open(name);
	
	while(getline(infile,reader)){
		if(reader.find(runNum) == string::npos)//skip if it's not the correct run
			continue;
		strcpy(filename,reader.c_str());//convert string to char array cuz sprintf is fickle
		sprintf(name,"/cache/clas/eg4a/production/pass2/%sgev/nh3/root/%s",energylvl,filename);
		cout << name << endl; //For testing + progress tracking
		TFile *f = new TFile(name);
		if(!f){
			cout << "Error: Unable to load " << name << endl;
			continue;
		}
		TTree* tree = (TTree*)gDirectory->Get("h10");
		Init(tree);
		Int_t nentries = (Int_t)tree->GetEntriesFast();
		for(Int_t i = 0; i<nentries; i++) {
			tree->GetEntry(i);
			//Int_t Ebi = 4; //0,1,2,3,4 for 1.0539,1.9889,2.2256,2.999 respectively
			Eb = EbgNominal[Ebi] - 0.002;
			//kinematic corrections//
			Float_t ppC=p[0],thetaRadC = TMath::ACos(cz[0])*180.0/PI,ZrcorC=vz[0], KineCorrectedThetaVertex = TMath::ACos(tl1_cz[0])*180.0/PI;//TODO: verify
			Int_t Dc = int(dc[0]);
			Double_t corCxCyCzVzP[5];
			for(Int_t i = 0; i<5; i++){corCxCyCzVzP[i] = 0.0;}
			Int_t PID = 0; //0 or 1 for e- and hadrons, respectively
			applyAllKineCorPass2(PID, Ebi, (int)q[0], raster_x, raster_y, p[0], cx[0], cy[0], cz[0], vz[0], tl1_cx[Dc-1], tl1_cy[Dc-1], tl1_x[Dc-1], tl1_y[Dc-1], tl1_z[Dc-1], corCxCyCzVzP);
			/////////////////////////
			phi = TMath::ATan(cy[0]/cx[0]);
			if(cx[0] < 0)
				phi += PI;
			if(phi < 0)
				phi += 2*PI;
			theta = TMath::ACos(cz[0]);
			phDclPosRad = TMath::ATan(tl1_cy[0]/tl1_cx[0]);
			if(tl1_cx[0] < 0)
				phDclPosRad += PI;
			if(phDclPosRad < 0)
				phDclPosRad += 2*PI;
			phDclPosDeg = phDclPosRad * 180.0 / PI;
			thDclPosRad = TMath::ACos(tl1_cz[0]);
			qSq = 2*Eb*p[0]*(1-cz[0]);
			v = Eb - p[0];
			wSq = M*M+2*M*v - qSq;
			if(wSq > 0){
				w = sqrt(wSq);
			}
			else{
				w = 0;
			}
			////////////////cuts go here////////////////
			Bool_t stCtNoCc = 0; //Detector status cuts
			Int_t Cc = int(cc[0]), Sc = int(sc[0]), Ec = int(ec[0]), Qq = int(q[0]);
			if(Dc>0 && Sc>0 && Ec>0 && dc_part>0 && sc_part>0 && ec_part>0 && Dc<=dc_part 
     && Sc<=sc_part && Ec<=ec_part && int(dc_stat[Dc-1])>0 
     && gpart>0 && Qq<0 && p[0]>0.3) stCtNoCc=true;
			Bool_t ecInQ2Ct = ecCuts4ElPass2(Ebi,p[0],cz[0],ec_ei[0],etot[0]);//EC cut
			//Bool_t osiCutTruth = OsipenkoCutsPass2Test(Ebi); //Osipenko cut
			Bool_t vzInQ2Ct = vzCutsInQ2Bins(Ebi,ppC,thetaRadC,ZrcorC);//vz cuts
			Bool_t pMinCt = true; //pMin cut
			if(ppC < 0.4 || ppC < 0.2*Eb)
				pMinCt = false;
			Bool_t sectCt = false; //sector cut
			if(Dc>0 && dc_sect[Dc-1]==6)
				sectCt = true;
			Bool_t npheMinCt = true; //nphe minimum cut
			if(nphe[0]<25.0)
				npheMinCt=false;
			Bool_t fidCtReg2EC = Pass2FidCutLatestFromRegECcomparison(Ebi,p[0],thDclPosRad); //fiducial cut 1
			Bool_t fidCtExp2Sim = Pass2FidCutLatest(phDclPosRad,theta);
			Bool_t fidCtV1 = false;
			if(fidCtExp2Sim==true && fidCtReg2EC==true)
				fidCtV1 = true;
			Bool_t fidCtExpBySimAndEConlyInInvPvsThVtx = Pass2FidCutOnInvPvsThVtx(Ebi, ppC, thetaRadC); //fiducial cut 2
			Double_t phDc1DegPlusMinus30 = phDclPosDeg, thVtxDeg= KineCorrectedThetaVertex; //fiducial cut 3
  			Bool_t fidMoreMIPB = false;
			fidMoreMIPB = moreFidCutsWithMoreInvPBins(Ebi, ppC, phDc1DegPlusMinus30, thVtxDeg);
			if(!ecInQ2Ct || !stCtNoCc || w!=w || /*!osiCutTruth ||*/ !vzInQ2Ct || !sectCt || !npheMinCt || !fidCtV1 || !fidCtExpBySimAndEConlyInInvPvsThVtx || !fidMoreMIPB) continue;
			////////////////////////////////////////////
			Int_t qBin = 0;
			for(Int_t j = 1; j<size; j++) { //find q^2 bin
				if(qSq < edges[j]) {
					qBin = j;
					break;
				}
			}
			Int_t wBin = (int)(w*(NUM_BINS/3));//quickly determine bin
			if((hel_ofl == 0)||(hel_ofl == 32512)){ //determine helicity
				nminus[qBin][wBin]++; //and increment n for that bin
				fcminus[qBin] += (float)fcup_g; //and add to fcup average
				fcmcount[qBin]++;
			}
			else if((hel_ofl == 1)||(hel_ofl == 32513)){
				nplus[qBin][wBin]++;
				fcplus[qBin] += (float)fcup_g;
				fcpcount[qBin]++;
			}
		}//end tree loop
	}//end overall loop
	//asym section//
	for(Int_t j = 0; j<sizeMinusOne; j++){//j = q2 bin (j comes first cuz I say so)
		for(Int_t i = 0; i<NUM_BINS; i++) {//i = w bin
			if(((nplus[j][i]==0)&&(nminus[j][i]==0))||(fcminus[j]==0)||(fcplus[j]==0)){
				a[j][i] = 0;
				sig[j][i] = 0;
			}
			else{
				sig[j][i] = 1.0/sqrt(nplus[j][i] + nminus[j][i]);//finding sig before fcup, is that fine?
				nminus[j][i] = nminus[j][i] / fcminus[j];
				nplus[j][i] = nplus[j][i] / fcplus[j];//fcup normalization
				a[j][i] = (nplus[j][i] - nminus[j][i])/(nplus[j][i] + nminus[j][i]);
			}
		}
	}
	////////////////
	/**
	 *Summary consists of sets of 2 lines
	 *Asymmetry values for all W bins in that Q^2 bin
	 *Sigma values for W all bins in that Q^2
	 **/
	Char_t outName[100];
	sprintf(outName,"/home/ryanpl/Thesis/common/skimsV4/%s_%s_skim.dat",energylvl,runNum);
	ofstream outfile;
	outfile.open(outName);
	outfile << size << endl;
	for(Int_t j=0;j<sizeMinusOne;j++){
		for(Int_t i=0;i<NUM_BINS;i++){
			outfile << a[j][i] << " ";
		}
		outfile << endl;
		for(Int_t i=0;i<NUM_BINS;i++){
			outfile << sig[j][i] << " ";
		}
		outfile << endl;
	}
	outfile.close();
}

