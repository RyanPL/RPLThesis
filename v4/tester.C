#include<stdlib.h>
#include<stdio.h>

#define tester_cxx

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TROOT.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFormula.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "tester.h"

using namespace std;

//#define E 2.0;
//#define M 0.938;
//#define PI 3.14159

void tester::lastgraphs()
{
	Float_t edges[56] = {0, 0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};
	const Int_t nbins = 55;
	Float_t w, qSq, v, phi, az;
	const Double_t E = 2.0;
	const Double_t PI = 3.14159;
	const Double_t M = 0.938;
	TH1F *hqsq = new TH1F("hqsq", "Q^{2} - Virtual Photon 4-Momentum", nbins, edges);
	TH1F *hw = new TH1F("hw","W - Lorentz-Invariant Missing Mass", 150, 0, 3);
	TH2F *hqvw = new TH2F("hqvw","Q^{2} (strength of interaction) vs W (mass of scattered products; Q^2; W)", 400, -1, 3, 300, 0, 3);
	Int_t nentries = (Int_t)h10->GetEntries();
	for(Int_t i = 0; i<nentries; i++) {
		h10->GetEntry(i);
		phi = TMath::ATan(cy[0]/cx[0]);
		if(cx[0] < 0)
			phi += PI;
		if(phi < 0)
			phi += 2*PI;
		az = TMath::ACos(cz[0]);
		qSq = 2*E*p[0]*(1-cz[0]);
		v = E - p[0];
		w = TMath::Sqrt((M*M)+(2*M*v) - qSq);
		hqsq->Fill(qSq);
		hw->Fill(w);
		hqvw->Fill(qSq,w);
	}

	TCanvas *c1 = new TCanvas("c1", "Q^{2} - Virtual Photon 4-Momentum", 40, 40, 1200, 800);
	hqsq->GetXaxis()->SetTitle("Q^{2} (GeV)");
	hqsq->GetYaxis()->SetTitle("Event Count");
	hqsq->GetXaxis()->CenterTitle();
	hqsq->GetYaxis()->CenterTitle();
	gStyle->SetOptStat(0);
	hqsq->Draw();
	TCanvas *c2 = new TCanvas("c2", "W - Lorentz-Invariant Missing Mass", 80, 80, 1200, 800);
	hw->GetXaxis()->SetTitle("W (GeV)");
	hw->GetYaxis()->SetTitle("Event Count");
	hw->GetXaxis()->CenterTitle();
	hw->GetYaxis()->CenterTitle();
	gStyle->SetOptStat(0);
	hw->Draw();
}




/*
 * Dismantled the old depricated(sp?) comparison method. Now using it to look at individual runs to nail own a flipcheck zone
 */
void tester::onegraph()
{
	char energy[10], run[10], filename[100], data[100];
	bool combine;
	cout << "Enter energy" << endl;
	scanf("%s",&energy);
	cout << "Enter run number" << endl;
	scanf("%s",&run);
	cout << "Combine?" << endl;
	scanf("%s",&data);
	if(strcmp(data,"Yes")==0)
		combine = true;
	else
		combine = false;
	sprintf(filename,"../common/skimsV4/%s_%s_skim.dat",energy,run);
	
	ifstream file;
	file.open(filename);
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	Float_t xaxis[numWBins];
	for(int i=0; i<numWBins; i++){xaxis[i]=(float)i*0.04;}//create basic x-axis
	Float_t asym[numQBins][numWBins], error[numQBins][numWBins];

	for(Int_t qBin = 0; qBin < numQBins; qBin++){
		for(Int_t wBin = 0; wBin < numWBins; wBin++) {
			file >> data;
			asym[qBin][wBin] = atof(data);
		}
		for(Int_t wBin = 0; wBin < numWBins; wBin++) {
			file >> data;
			if(atof(data) == 0)
				error[qBin][wBin] = 0.0;
			else
				error[qBin][wBin] = atof(data);
		}
	}
	file.close();

	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};
	
	Float_t zero[numWBins];
	for(int i = 0; i < numWBins; i++){zero[i] = 0.0;}

	if(combine)
	{
		Float_t y[numWBins];
		zero1DArray(y,numWBins);
		for(int qBin = 0; qBin < numQBins; qBin++)
		{
			for(int wBin = 0; wBin < numWBins; wBin++)
			{
				y[wBin] += asym[qBin][wBin];
			}
		}
		TCanvas *c = new TCanvas("c","Combined run asym",40,40,800,600);
		TH2D *blank = new TH2D("Combined graph", "Combined graph", 300, 0, 3, 500, -5.0, 5.0);
		TGraph *g = new TGraph(numWBins,xaxis,y);
		blank->Draw();
		g->Draw("same");
		gStyle->SetOptStat(0); //remove annoying stat box
		return;
	}

	
	//graphing loop
	for(Int_t i = 0; i < 9; i++)
	{
		Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "A run with a name", 20*i, 20*i, 1200, 600);
		c->Divide(3,2);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			char blankName[10];
			sprintf(blankName,"%f < Q^2 < %f GeV",edges[graphNum],edges[graphNum+1]);
			TH2D *blank = new TH2D(blankName, blankName, 300, 0, 3, 500, -1.0, 1.0);
			Float_t y[numWBins], yerror[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				y[bin] = asym[graphNum][bin];
				yerror[bin] = error[graphNum][bin];
			}
			TGraphErrors *gasym = new TGraphErrors(numWBins, xaxis, y, zero, yerror);
			c->cd(j+1);
			blank->Draw();
			gasym->Draw("same");
			gStyle->SetOptStat(0); //remove annoying stat box
		}
	}
}

void tester::countcheck()
{
	char energy[10], run[10], target[10], filename1[100], filename2[100];
	cout << "Enter target" << endl;
	scanf("%s",&target);
	cout << "Enter energy" << endl;
	scanf("%s",&energy);
	cout << "Enter run number" << endl;
	scanf("%s",&run);
	sprintf(filename1,"../common/n_skimsV2/%s/%s_%s_%s_np.dat",target,target,energy,run);
	sprintf(filename2,"../common/n_skimsV2/%s/%s_%s_%s_nm.dat",target,target,energy,run);
	
	ifstream file1, file2;
	file1.open(filename1);
	file2.open(filename2);
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	Char_t data[100];
	Float_t xaxis[numWBins], zero[numWBins];
	for(int i=0; i<numWBins; i++){
		xaxis[i]=(float)i*0.04;//create basic x-axis
		zero[i] = 0.0;}
	Float_t pos[numQBins][numWBins], neg[numQBins][numWBins];

	for(Int_t qBin = 0; qBin < numQBins; qBin++){
		for(Int_t wBin = 0; wBin < numWBins; wBin++) {
			file1 >> data;
			pos[qBin][wBin] = atof(data);
			file2 >> data;
			neg[qBin][wBin] = atof(data);
		}
		file1 >> data; //skipping fcup
		file2 >> data;
	}
	file1.close();
	file2.close();

	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};

	
	//graphing loop
	for(Int_t i = 0; i < 9; i++)
	{
		Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "A run with a name", 20*i, 20*i, 1600, 800);
		c->Divide(3,2);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			char blankName[10];
			sprintf(blankName,"%f < Q^2 < %f GeV",edges[graphNum],edges[graphNum+1]);
			TH2D *blank = new TH2D(blankName, blankName, 300, 0, 3, 500, -1.0, 1.0);
			Float_t y1[numWBins], y2[numWBins], y1e[numWBins], y2e[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				y1[bin] = pos[graphNum][bin];
				y2[bin] = neg[graphNum][bin];
				y1e[bin] = pow(pos[graphNum][bin], -0.5);
				y2e[bin] = pow(neg[graphNum][bin], -0.5);
			}
			TGraphErrors *g1 = new TGraphErrors(numWBins, xaxis, y1, zero, y1e);
			TGraphErrors *g2 = new TGraphErrors(numWBins, xaxis, y2, zero, y2e);
			g1->SetLineColor(kRed);
			g2->SetLineColor(kBlue);
			c->cd(j+1);
			g1->Draw();
			g2->Draw("same");
			leg = new TLegend(0.6, 0.7, 0.89, 0.89);
			leg->AddEntry(g1,"Positive","l");
			leg->AddEntry(g2,"Negative","l");
			leg->Draw();
		}
	}
}

void tester::zero1DArray(Float_t a[], Int_t size)
{
	for(int i = 0; i < size; i++)
	{
		a[i] = 0.0;
	}
}

void tester::zero2DArray(Float_t a[][75], Int_t s1, Int_t s2)
{
	for(int i = 0; i < s1; i++)
	{
		for(int j = 0; j < s2; j++)
		{
			a[i][j] = 0.0;
		}
	}
}

void tester::parsegraph()
{
	Char_t data[100], filename[100], cutname[100];
	string fileToRead;	
	ifstream fileList, infile;
	const Int_t numWBins = 75;
	const Int_t numQBins = 53;

	//carb first
	system("touch tempfile.dat");
	system("ls ../common/n_skimsV2/carb >> carbtempfile.dat");
	fileList.open("carbtempfile.dat");
	Float_t cpraw[numQBins][numWBins], cmraw[numQBins][numWBins];
	Float_t fccp = 0.0;//fc+ for carb
	Float_t fccm = 0.0;//fc- for carb
	zero2DArray(cpraw,numQBins,numWBins);
	zero2DArray(cmraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/carb/%s",cutname);
		infile.open(filename);
		bool pos;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					cpraw[qBin][wBin] += atof(data);
				else
					cmraw[qBin][wBin] += atof(data);
			}
			infile >> data;
			if(pos)
				fccp += atof(data);
			else
				fccm += atof(data);
		}
		
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm carbtempfile.dat");
	
	//now emp
	system("touch tempfile.dat");
	system("ls ../common/n_skimsV2/emp >> emptempfile.dat");
	fileList.open("emptempfile.dat");
	Float_t epraw[numQBins][numWBins], emraw[numQBins][numWBins];
	Float_t fcep = 0.0;//fc+ for emp
	Float_t fcem = 0.0;//fc- for emp
	zero2DArray(epraw,numQBins,numWBins);
	zero2DArray(emraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/emp/%s",cutname);
		infile.open(filename);
		bool pos;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					epraw[qBin][wBin] += atof(data);
				else
					emraw[qBin][wBin] += atof(data);
			}
			infile >> data;
			if(pos)
				fcep += atof(data);
			else
				fcem += atof(data);
		}
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm emptempfile.dat");

	//now nh3
	system("touch tempfile.dat");
	system("ls ../common/n_skimsV2/nh3 >> nh3tempfile.dat");
	fileList.open("nh3tempfile.dat");
	Float_t npraw[numQBins][numWBins], nmraw[numQBins][numWBins];
	Float_t fcnp = 0.0;//fc+ for nh3
	Float_t fcnm = 0.0;//fc- for nh3
	zero2DArray(npraw,numQBins,numWBins);
	zero2DArray(nmraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/nh3/%s",cutname);
		infile.open(filename);
		bool pos;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					npraw[qBin][wBin] += atof(data);
				else
					nmraw[qBin][wBin] += atof(data);
			}
			infile >> data;
			if(pos)
				fcnp += atof(data);
			else
				fcnm += atof(data);
			
		}
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm nh3tempfile.dat");

	//some numbers for this particular method
	Float_t fcntot = fcnp + fcnm;
	Float_t fcetot = fcep + fcem;
	Float_t fcctot = fccp + fccm;

	//prepare final arrays
	Float_t carb[numQBins][numWBins], emp[numQBins][numWBins], nh3[numQBins][numWBins], xaxis[numWBins];
	for(int i=0;i<75;i++)xaxis[i]=(float)i*0.01;
	Float_t maxval = 0.0;
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			carb[i][j] = (fcntot / fcctot) * (cpraw[i][j] + cmraw[i][j]);
			emp[i][j] = (fcntot / fcetot) * (epraw[i][j] + emraw[i][j]);
			nh3[i][j] = (npraw[i][j] + nmraw[i][j]);
			//carb[i][j] = (cpraw[i][j] / fccp) + (cmraw[i][j] / fccm);
			//emp[i][j] = (epraw[i][j] / fcep) + (emraw[i][j] / fcem);
			//nh3[i][j] = (npraw[i][j] / fcnp) + (nmraw[i][j] / fcnm);
			if(carb[i][j] > maxval) maxval = carb[i][j];
			if(emp[i][j] > maxval) maxval = emp[i][j];
			if(nh3[i][j] > maxval) maxval = nh3[i][j];
		}
	}
	maxval *= 1.1;
	//cout << carb[15][150] << " " << emp[15][150] << " " << nh3[15][150] << endl;
	//cout << fcntot << " " << fcetot << " " << fcctot << endl;

	//graphing loop
	for(Int_t i = 0; i < 9; i++)
	{
		Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "n spectra", 20*i, 20*i, 1200, 600);
		c->Divide(3,2);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			char blankName[10];
			sprintf(blankName,"Bin #%d",graphNum);
			TH2D *blank = new TH2D(blankName, blankName, 75, 0, 3, 500, 0, maxval);
			Float_t ye[numWBins], yc[numWBins], yn[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				ye[bin] = emp[graphNum][bin];
				yc[bin] = carb[graphNum][bin];
				yn[bin] = nh3[graphNum][bin];
			}
			TGraph *gemp = new TGraph(numWBins, xaxis, ye);
			TGraph *gcarb = new TGraph(numWBins, xaxis, yc);
			gcarb->SetLineColor(kBlue);
			TGraph *gnh3 = new TGraph(numWBins, xaxis, yn);
			gnh3->SetLineColor(kRed);
			c->cd(j+1);
			blank->Draw();
			gemp->Draw("same");
			gcarb->Draw("same");
			gnh3->Draw("same");
		}
	}
	
}

void tester::cutcheck()
{
	ifstream file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12;
	file1.open("../common/cutChecks/all.dat");
	file2.open("../common/cutChecks/allbutfid.dat");
	file3.open("../common/cutChecks/ecUbQ2Ct.dat");
	file4.open("../common/cutChecks/fidCtExpBySimAndEConlyInInvPvsThVtx.dat");
	file5.open("../common/cutChecks/fidCtV1.dat");
	file6.open("../common/cutChecks/fidMoreMIPB.dat");
	file7.open("../common/cutChecks/none.dat");
	file8.open("../common/cutChecks/npheMinCt.dat");
	file9.open("../common/cutChecks/sectCt.dat");
	file10.open("../common/cutChecks/ctCtNoCc.dat");
	file11.open("../common/cutChecks/vzInQ2Ct.dat");
	file12.open("../common/cutChecks/ycut.dat");
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	Char_t data[100];
	Float_t xaxis[numWBins];
	for(int i=0; i<numWBins; i++){xaxis[i]=(float)i*0.04;}//create basic x-axis

	Float_t c1[numQBins][numWBins], c2[numQBins][numWBins], c3[numQBins][numWBins], c4[numQBins][numWBins], c5[numQBins][numWBins], c6[numQBins][numWBins], c7[numQBins][numWBins], c8[numQBins][numWBins], c9[numQBins][numWBins], c10[numQBins][numWBins], c11[numQBins][numWBins], c12[numQBins][numWBins];
	Float_t fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8, fc9, fc10, fc11, fc12, maxval[numQBins];
	zero1DArray(maxval,numQBins);
	for(int qBin = 0; qBin < numQBins; qBin++)
	{
		for(int wBin = 0; wBin < numWBins; wBin++)
		{
			file1 >> data;
			c1[qBin][wBin] = atof(data);
			file2 >> data;
			c2[qBin][wBin] = atof(data);
			file3 >> data;
			c3[qBin][wBin] = atof(data);
			file4 >> data;
			c4[qBin][wBin] = atof(data);
			file5 >> data;
			c5[qBin][wBin] = atof(data);
			file6 >> data;
			c6[qBin][wBin] = atof(data);
			file7 >> data;
			c7[qBin][wBin] = atof(data);
			if(c7[qBin][wBin] > maxval[qBin])
				maxval[qBin] = c7[qBin][wBin];
			file8 >> data;
			c8[qBin][wBin] = atof(data);
			file9 >> data;
			c9[qBin][wBin] = atof(data);
			file10 >> data;
			c10[qBin][wBin] = atof(data);
			file11 >> data;
			c11[qBin][wBin] = atof(data);
			file12 >> data;
			c12[qBin][wBin] = atof(data);
		}
		maxval[qBin] = maxval[qBin] * 1.1;

		file1 >> data;
		fc1 += atof(data);
		file2 >> data;
		fc2 += atof(data);
		file3 >> data;
		fc3 += atof(data);
		file4 >> data;
		fc4 += atof(data);
		file5 >> data;
		fc5 += atof(data);
		file6 >> data;
		fc6 += atof(data);
		file7 >> data;
		fc7 += atof(data);
		file8 >> data;
		fc8 += atof(data);
		file9 >> data;
		fc9 += atof(data);
		file10 >> data;
		fc10 += atof(data);
		file11 >> data;
		fc11 += atof(data);
		file12 >> data;
		fc12 += atof(data);
	}


	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};

	
	Float_t zero[numWBins];
	for(int i = 0; i < numWBins; i++)
	{
		zero[i] = 0.0;
	}

	//graphing loop
	for(Int_t i = 0; i < 9; i++)
	{
		Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "Cut check", 20*i, 20*i, 1600, 800);
		c->Divide(3,2);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			char blankName[10];
			sprintf(blankName,"%f < Q^2 < %f GeV",edges[graphNum],edges[graphNum+1]);
			TH2D *blank = new TH2D(blankName, blankName, 75, 0, 3, 500, 0, maxval[graphNum]);
			Float_t y1[numWBins], y2[numWBins], y3[numWBins], y4[numWBins], y5[numWBins], y6[numWBins], y7[numWBins], y8[numWBins], y9[numWBins], y10[numWBins], y11[numWBins], y12[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				y1[bin] = c1[graphNum][bin];
				y2[bin] = c2[graphNum][bin];
				y3[bin] = c3[graphNum][bin];
				y4[bin] = c4[graphNum][bin];
				y5[bin] = c5[graphNum][bin];
				y6[bin] = c6[graphNum][bin];
				y7[bin] = c7[graphNum][bin];
				y8[bin] = c8[graphNum][bin];
				y9[bin] = c9[graphNum][bin];
				y10[bin] = c10[graphNum][bin];
				y11[bin] = c11[graphNum][bin];
				y12[bin] = c12[graphNum][bin];
			}
			TGraph *g1 = new TGraph(numWBins, xaxis, y1);
			TGraph *g2 = new TGraph(numWBins, xaxis, y2);
			g2->SetLineColor(2);
			TGraph *g3 = new TGraph(numWBins, xaxis, y3);
			g3->SetLineColor(3);
			TGraph *g4 = new TGraph(numWBins, xaxis, y4);
			g4->SetLineColor(4);
			TGraph *g5 = new TGraph(numWBins, xaxis, y5);
			g5->SetLineColor(5);
			TGraph *g6 = new TGraph(numWBins, xaxis, y6);
			g6->SetLineColor(6);
			TGraph *g7 = new TGraph(numWBins, xaxis, y7);
			g7->SetLineColor(7);
			TGraph *g8 = new TGraph(numWBins, xaxis, y8);
			g8->SetLineColor(8);
			TGraph *g9 = new TGraph(numWBins, xaxis, y9);
			g9->SetLineColor(9);
			TGraph *g10 = new TGraph(numWBins, xaxis, y10);
			g10->SetLineColor(kOrange);
			TGraph *g11 = new TGraph(numWBins, xaxis, y11);
			g11->SetLineColor(21);
			TGraph *g12 = new TGraph(numWBins, xaxis, y12);
			g12->SetLineColor(28);
			c->cd(j+1);
			blank->Draw();
			g1->Draw("same");
			g2->Draw("same");
			g3->Draw("same");
			g4->Draw("same");
			g5->Draw("same");
			g6->Draw("same");
			g7->Draw("same");
			g8->Draw("same");
			g9->Draw("same");
			g10->Draw("same");
			g11->Draw("same");
			g12->Draw("same");
			if(graphNum == (i*6)){
				leg = new TLegend(0.4,0.6,0.89,0.89);
				leg->AddEntry(g1,"all","l");
				leg->AddEntry(g2,"allbutfid","l");
				leg->AddEntry(g3,"ecUbQ2Ct","l");
				leg->AddEntry(g4,"fidCtExpBySimAndEConlyInInvPvsThVtx","l");
				leg->AddEntry(g5,"fidCtV1","l");
				leg->AddEntry(g6,"fidMoreMIPB","l");
				leg->AddEntry(g7,"none","l");
				leg->AddEntry(g8,"npheMinCt","l");
				leg->AddEntry(g9,"sectCt","l");
				leg->AddEntry(g10,"ctCtNoCc","l");
				leg->AddEntry(g11,"vzInQ2Ct","l");
				leg->AddEntry(g12,"ycut","l");
				leg->Draw();
			}
		}
	}
}

/*
 * Dismantled the old depricated(sp?) comparison method. Now using it to look at individual runs to nail own a flipcheck zone
 */
void tester::systerr()
{
	cout << "お前はもう死んでいる" << endl;
	char filename2[100], data[100];
	cout << "What error are we looking at? (pbpt,fdf)" << endl;
	scanf("%s",&data);
	if((strcmp(data,"pbpt")!=0)&&(strcmp(data,"fdf")!=0))
	{
		if((strcmp(data,"nani?")==0))
		{
			cout << "SSKSKKKKKRRRRRRRRT" << endl;
			return;
		}
		cout << "Enter something correct ya dumbass" << endl;
		return;
	}
	sprintf(filename2,"../common/systerr/asym_%s.dat",data);
	
	ifstream file1, file2;
	file1.open("../common/asym_final.dat");
	file2.open(filename2);
	const Int_t numELevels = 3;
	const Int_t numQBins = 54;
	const Int_t numWBins = 75;
	Float_t asym[numELevels][numQBins][numWBins], errorasym[numELevels][numQBins][numWBins],
			ave[numELevels][numQBins][numWBins], diff[numELevels][numQBins][numWBins];

	for(Int_t eBin = 0; eBin < numELevels; eBin++)
	{
		for(Int_t qBin = 0; qBin < numQBins; qBin++){
			for(Int_t wBin = 0; wBin < numWBins; wBin++) {
				file1 >> data;
				float num1 = atof(data); //I swear there's a reason I'm being so roundabout
				asym[eBin][qBin][wBin] = num1;
				file2 >> data;
				float num2 = atof(data); //it's coming
				errorasym[eBin][qBin][wBin] = num2;
				ave[eBin][qBin][wBin] = (num1 + num2) / 2; //told you
				diff[eBin][qBin][wBin] = abs((num1 - num2) / 2); //error bar is half the difference
				//if(abs(num1 - num2) < 0.02) //if there is no/low error bar
				//	diff[eBin][qBin][wBin] = 0.01; //set a basic, visible error bar
				//IMPORTANT: DELETE the above segment once you find a better option!
			}
		}
	}
	file1.close();
	file2.close();

	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};
	
	Float_t zero[numWBins];
	for(int i = 0; i < numWBins; i++){zero[i] = 0.0;}
	Float_t xaxis[numWBins];
	for(int i=0; i<numWBins; i++){xaxis[i]=(float)i*0.04;}//create basic x-axis

	ofstream lerr1, lerr2, lerr3;
	lerr1.open("../common/systerr/last_fdf1.dat"); //i just don't fucking care anymore
	lerr2.open("../common/systerr/last_fdf2.dat");
	lerr3.open("../common/systerr/last_fdf3.dat");


	//graphing loop
	for(Int_t i = 0; i < 9; i++)
	{
		Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "A run with a name", 20*i, 20*i, 1600, 800);
		c->Divide(3,2);
		c->SetFillColor(42);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			char blankName[10];
			sprintf(blankName,"%f < Q^2 < %f GeV",edges[graphNum],edges[graphNum+1]);
			TH2D *blank = new TH2D(blankName, blankName, 300, 0, 3, 500, -1.0, 1.0);
			Float_t y1a[numWBins], y1b[numWBins],
					y2a[numWBins], y2b[numWBins],
					y3a[numWBins], y3b[numWBins], //a is original, b is error
					b1[numWBins], b1e[numWBins],
					b2[numWBins], b2e[numWBins],
					b3[numWBins], b3e[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				y1a[bin] = asym[0][graphNum][bin];
				y1b[bin] = errorasym[0][graphNum][bin];
				y2a[bin] = asym[1][graphNum][bin];
				y2b[bin] = errorasym[1][graphNum][bin];
				y3a[bin] = asym[2][graphNum][bin];
				y3b[bin] = errorasym[2][graphNum][bin];
				b1[bin] = ave[0][graphNum][bin];
				b1e[bin] = diff[0][graphNum][bin];
				b2[bin] = ave[1][graphNum][bin];
				b2e[bin] = diff[1][graphNum][bin];
				b3[bin] = ave[2][graphNum][bin];
				b3e[bin] = diff[2][graphNum][bin];
			}
			TGraph *g1a = new TGraph(numWBins, xaxis, y1a);
			TGraph *g1b = new TGraph(numWBins, xaxis, y1b);
			TGraph *g2a = new TGraph(numWBins, xaxis, y2a);
			TGraph *g2b = new TGraph(numWBins, xaxis, y2b);
			TGraph *g3a = new TGraph(numWBins, xaxis, y3a);
			TGraph *g3b = new TGraph(numWBins, xaxis, y3b);
			TGraphErrors *bar1 = new TGraphErrors(numWBins, xaxis, b1, zero, b1e);
			TGraphErrors *bar2 = new TGraphErrors(numWBins, xaxis, b2, zero, b2e);
			TGraphErrors *bar3 = new TGraphErrors(numWBins, xaxis, b3, zero, b3e);
			/*TH1F *err1 = new TH1F(blankName, blankName, numWBins, b1e);
			TH1F *err2 = new TH1F(blankName, blankName, numWBins, b2e);
			TH1F *err3 = new TH1F(blankName, blankName, numWBins, b3e);*/
			g1a->SetLineColor(3);
			g2a->SetLineColor(8);
			g3a->SetLineColor(30);
			g1b->SetLineColor(4);
			g2b->SetLineColor(9);
			g3b->SetLineColor(38);
			bar1->SetLineColor(kRed);
			bar2->SetLineColor(kGreen);
			bar3->SetLineColor(kBlue);
			/*err1->SetFillColor(kRed);
			err2->SetFillColor(kGreen);
			err3->SetFillColor(kBlue);*/
			c->cd(j+1);
			gPad->SetFillColor(42);
			blank->Draw();
			/*g1a->Draw("same");
			g1b->Draw("same");
			g2a->Draw("same");
			g2b->Draw("same");
			g3a->Draw("same");
			g3b->Draw("same");*/
			bar1->Draw("P");
			bar2->Draw("P");
			bar3->Draw("P");
			/*err1->Draw("same");
			err2->Draw("same");
			err3->Draw("same");*/
			gStyle->SetOptStat(0); //remove annoying stat box
			gStyle->SetEndErrorSize(3);
			leg = new TLegend(0.4, 0.6, 0.89, 0.89);
			leg->AddEntry(bar1,"1.1 and 1.3 GeV","l");
			leg->AddEntry(bar2,"2 and 2.2 GeV","l");
			leg->AddEntry(bar3,"3 GeV","l");
			if(j==1)
				leg->Draw(); //keep the legend to one window
			//quick errors
			/*if((graphNum < 10)||(graphNum > 40))
				continue;
			cout << blankName << "    ";*/
			float num[6];
			zero1DArray(num,6);
			//for(int g = 20; g < 35; g++)
			for(int g = 0; g < numWBins; g++)
			{
				if((asym[0][graphNum][g] != 0)&&(b1e[g]!=0))
				{
					//num[0] += (b1e[g] * 2)/b1[g];
					float val = abs((b1e[g] * 2)/asym[0][graphNum][g]);
					if(val <= 0.25){
						num[0] += val;
						num[1] += 1.0;
					}
				}
				if((asym[1][graphNum][g] != 0)&&(b2e[g]!=0))
				{
					//num[2] += (b2e[g] * 2)/b2[g];
					float val = abs((b2e[g] * 2)/asym[1][graphNum][g]);
					if(val <= 0.25){
						num[2] += val;
						num[3] += 1.0;
					}
				}
				if((asym[2][graphNum][g] != 0)&&(b3e[g]!=0))
				{
					//num[4] += (b3e[g] * 2)/b3[g];
					float val = abs((b3e[g] * 2)/asym[2][graphNum][g]);
					if(val <= 0.25){
						num[4] += val;
						num[5] += 1.0;
					}
				}
				lerr1 << (b1e[g] * 2) << " ";
				lerr2 << (b2e[g] * 2) << " ";
				lerr3 << (b3e[g] * 2) << " ";
			}
			lerr1 << endl;
			lerr2 << endl;
			lerr3 << endl;
			float perc1 = /*100 */ (num[0]/num[1]);
			float perc2 = /*100 */ (num[2]/num[3]);
			float perc3 = /*100 */ (num[4]/num[5]);
			//cout << perc1 << "% : " << perc2 << "% : " << perc3 << "%" << endl;
		}
	}
}

void tester::errorcheck()
{
	Char_t data[100], filename[100], cutname[100];
	string fileToRead;	
	ifstream fileList, infile;
	const Int_t numELevels = 3;
	const Int_t numWBins = 75;
	const Int_t numQBins = 54;

	system("touch tempfile.dat");
	system("ls ../common/n_skimsV2/nh3 >> nh3tempfile.dat");
	fileList.open("nh3tempfile.dat");
	Float_t npraw[numQBins][numWBins], nmraw[numQBins][numWBins];
	Float_t fcnp = 0.0;//fc+ for nh3
	Float_t fcnm = 0.0;//fc- for nh3
	zero2DArray(npraw,numQBins,numWBins);
	zero2DArray(nmraw,numQBins,numWBins);

	int tick = 0;
	while(getline(fileList,fileToRead)){
		if(tick == 4)
			break;
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/nh3/%s",cutname);
		Int_t energyIndex = getEnergy(fileToRead.substr(4,3));
		if(energyIndex != 2)
			continue;
		tick++;
		infile.open(filename);
		bool pos;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					npraw[qBin][wBin] += atof(data);
				else
					nmraw[qBin][wBin] += atof(data);
			}
			infile >> data;
			if(pos)
				fcnp += atof(data);
			else
				fcnm += atof(data);
			
		}
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm nh3tempfile.dat");

	ifstream dfstream;
	Double_t df[numELevels][numQBins][numWBins];
	dfstream.open("../common/dfapprox.dat");
	for(int e = 0; e < numELevels; e++)
	{
		for(int i = 0; i < numQBins; i++)
		{
			for(int j = 0; j < numWBins; j++)
			{
				dfstream >> data;
				df[e][i][j] = atof(data);
			}
		}
	}

	//some numbers for this particular method
	Float_t fcntot = fcnp + fcnm;

	//prepare final arrays
	Float_t n[numQBins][numWBins], nnorm[numQBins][numWBins], subback[numQBins][numWBins],
		    a[numQBins][numWBins], apar[numQBins][numWBins], xaxis[numWBins], zero[numWBins];
	for(int i=0;i<75;i++){
		xaxis[i]=(float)i*0.04;
		zero[i]=0.0;}
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			n[i][j] = npraw[i][j] + nmraw[i][j];
			nnorm[i][j] = n[i][j] / fcntot;
			subback[i][j] = (npraw[i][j] - nmraw[i][j])/fcntot;
			a[i][j] = subback[i][j] / nnorm[i][j];
			apar[i][j] = a[i][j] / (df[2][i][j] * 0.614);
		}
	}

	Float_t y1[numWBins], y1e[numWBins],
			y2[numWBins], y2e[numWBins],
			y3a[numWBins], y3ae[numWBins],
			y3b[numWBins], y3be[numWBins],
			y4[numWBins], y4e[numWBins],
			y5[numWBins], y5e[numWBins],
			y6[numWBins], y6e[numWBins];
	for(int bin = 0; bin<numWBins; bin++)
	{
		y1[bin] = n[30][bin];
		if(n[30][bin]!=0)
			y1e[bin] = pow(n[30][bin],-0.5);
		if(y1e[bin] > 100)
			y1e[bin] = 0.0;
		y2[bin] = nnorm[30][bin];
		if(nnorm[30][bin]!=0)
			y2e[bin] = pow(n[30][bin],-0.5)/fcntot;
		if(y2e[bin] > 100)
			y2e[bin] = 0.0;
		y3a[bin] = npraw[30][bin]/fcnp;
		if(y3a[bin]!=0)
			y3ae[bin] = pow(npraw[30][bin],-0.5)/fcnp;
		if(y3ae[bin] > 100)
			y3ae[bin] = 0.0;
		y3b[bin] = nmraw[30][bin]/fcnp;
		if(y3b[bin]!=0)
			y3be[bin] = pow(nmraw[30][bin],-0.5)/fcnm;
		if(y3ae[bin] > 100)
			y3ae[bin] = 0.0;
		y4[bin] = subback[30][bin];
		if(subback[30][bin]!=0)
			y4e[bin] = sqrt((y3ae[bin]*y3ae[bin])+(y3be[bin]*y3be[bin]));
		if(y4e[bin] > 100)
			y4e[bin] = 0.0;
		y5[bin] = a[30][bin];
		if(n[30][bin]!=0)
			y5e[bin] = pow(n[30][bin],-0.5);
		if(y5e[bin] > 100)
			y5e[bin] = 0.0;
		y6[bin] = apar[30][bin];
		if(df[2][30][bin]!=0)
			y6e[bin] = y5e[bin]/(df[2][30][bin]*0.614);
		if(y6e[bin] > 100)
			y6e[bin] = 0.0;
	}

	//graphing segment
	TCanvas *c = new TCanvas("Errors", "Errors and Stuff", 20*i, 20*i, 1600, 800);
	c->Divide(3,2);
	TGraphErrors *g1 = new TGraphErrors(numWBins, xaxis, y1, zero, y1e);
	TGraphErrors *g2 = new TGraphErrors(numWBins, xaxis, y2, zero, y2e);
	TGraphErrors *g3a = new TGraphErrors(numWBins, xaxis, y3a, zero, y3ae);
	TGraphErrors *g3b = new TGraphErrors(numWBins, xaxis, y3b, zero, y3be);
	TGraphErrors *g4 = new TGraphErrors(numWBins, xaxis, y4, zero, y4e);
	TGraphErrors *g5 = new TGraphErrors(numWBins, xaxis, y5, zero, y5e);
	TGraphErrors *g6 = new TGraphErrors(numWBins, xaxis, y6, zero, y6e);
	g3a->SetLineColor(kRed);
	g3b->SetLineColor(kBlue);
	c->cd(1);
	TH2D *blank1 = new TH2D("","Raw counts",300,0,3,200,0,35000);
	blank1->Draw();
	g1->Draw("P");
	c->cd(2);
	TH2D *blank2 = new TH2D("","Normalized counts",300,0,3,200,0,0.000009);
	blank2->Draw();
	g2->Draw("P");
	c->cd(3);
	TH2D *blank3 = new TH2D("","n+ vs n-",300,0,3,200,0,0.000009);
	blank3->Draw();
	g3a->Draw("P");
	g3b->Draw("P");
	c->cd(4);
	TH2D *blank4 = new TH2D("","n+ - n-",300,0,3,200,-0.0000001,0.0000001);
	blank4->Draw();
	g4->Draw("P");
	c->cd(5);
	TH2D *blank5 = new TH2D("","Raw Asymmetry",300,1,2,200,-0.05,0.05);
	blank5->Draw();
	g5->Draw("P");
	c->cd(6);
	TH2D *blank6 = new TH2D("","A_||",300,1,2,200,-1,1);
	blank6->Draw();
	g6->Draw("P");
	
}

int getEnergy(string inp) //simple method for getting which E bin to put things in
{
	char e[3];
	strcpy(e,inp.c_str()); //what is even going on anymore
	if(strcmp(e,"1p1")==0)
		return 0;
	else if(strcmp(e,"1p3")==0)
		return 0;
	else if(strcmp(e,"1p5")==0) //1p5 appears to be bad. Removing for now
		return 6;
	else if(strcmp(e,"2_5")==0) //this relies on the run number starting with 5 cuz I'm a bad programmer
		return 1;
	else if(strcmp(e,"2p2")==0) //cutting this one too cuz its taking results from bad to worse
		return 6;
	else if(strcmp(e,"3_5")==0) 
		return 2;
	else
	{
		cout << "何?" << endl;
		return 6;
	}
}