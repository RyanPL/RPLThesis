#include<stdlib.h>
#include<stdio.h>

#define dfapprox_cxx

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

#include "dfapprox.h"

using namespace std;

int main(int argc, char *argv[]){
	
}

/*
 * The ultimate all in one DF approximator
 */
void dfapprox::parsegraph(Bool_t outputinstead = false, Bool_t fulldivide = true)
{
	if(!outputinstead)
		cout << "Parsing to graph (like the method says but whatever I'm not making a new one)" << endl;
	else if(outputinstead)
		cout << "Parsing to output" << endl;
	else{
		cout << "I don't know what the hell you	did but you really broke it" << endl;
		return;}

	if(fulldivide)
		cout << "Dividing EVERY energy level" << endl;
	else
		cout << "Merging similar energy levels" << endl;

	/*
	 * Deceptively small yet important bit where the master lists are filled
	 */
	int masterList[1080][3]; //first index is run, second index is run num, target type, torus current
	getMasterList(masterList);
	int masterListCarb[72][3];
	getMasterListCarb(masterListCarb);

	Char_t data[100], filename[100], cutname[100];
	string fileToRead;	
	ifstream fileList, infile;
	const Int_t numWBins = 75; //it DOES have 75 W bins, right?
	const Int_t numQBins = 53;
	const Int_t numELevels = 9; //rather than futz around with if-else everywhere,
								//I'll just pretend the top six levels don't exist if !fulldivide

	//carb first
	system("ls ../common/n_skimsV2/carb >> carbtempfile.dat");
	fileList.open("carbtempfile.dat");
	Float_t cpraw[numELevels][numQBins][numWBins], cmraw[numELevels][numQBins][numWBins];
	Float_t fccp[numELevels], fccm[numELevels];//fc+, fc- for carb
	zero3DArray(cpraw,numELevels,numQBins,numWBins);
	zero3DArray(cmraw,numELevels,numQBins,numWBins);
	zero1DArray(fccp,numELevels);
	zero1DArray(fccm,numELevels);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/carb/%s",cutname);
		infile.open(filename);
		bool pos;
		int index = findIndexCarb(fileToRead.substr(5,3),fileToRead.substr(7,7), masterListCarb);
		if(index == -1)
			continue;
		Int_t dataSet = getSetCarb(fileToRead.substr(5,3), masterListCarb[index][1], masterListCarb[index][2], fulldivide);
		if(dataSet == 9)
			continue;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					cpraw[dataSet][qBin][wBin] += atof(data);
				else
					cmraw[dataSet][qBin][wBin] += atof(data);				
			}
			infile >> data;
			if(pos)
				fccp[dataSet] += atof(data);
			else
				fccm[dataSet] += atof(data);
		}
		
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm carbtempfile.dat");

	for(int i = 0; i < numQBins; i++) //copying to other relevant sets
	{
		for(int j = 0; j < numWBins; j++)
		{
			/*cpraw[2][i][j] = cpraw[1][i][j];
			cpraw[3][i][j] = cpraw[1][i][j];
			cpraw[5][i][j] = cpraw[4][i][j];
			cpraw[7][i][j] = cpraw[6][i][j];

			cmraw[2][i][j] = cmraw[1][i][j];
			cmraw[3][i][j] = cmraw[1][i][j];
			cmraw[5][i][j] = cmraw[4][i][j];
			cmraw[7][i][j] = cmraw[6][i][j];

			fccp[2] = fccp[1];
			fccp[3] = fccp[1];
			fccp[5] = fccp[4];
			fccp[7] = fccp[6];

			fccm[2] = fccm[1];
			fccm[3] = fccm[1];
			fccm[5] = fccm[4];
			fccm[7] = fccm[6];*/
			//above part is for the old non-set split carbon. Keeping it just in case

			cpraw[3][i][j] = cpraw[1][i][j];
			cpraw[5][i][j] = cpraw[4][i][j];
			cpraw[7][i][j] = cpraw[4][i][j]; //new bit to account for lack of 2p2 data

			cmraw[3][i][j] = cmraw[1][i][j];
			cmraw[5][i][j] = cmraw[4][i][j];
			cmraw[7][i][j] = cmraw[4][i][j];

			fccp[3] = fccp[1];
			fccp[5] = fccp[4];
			fccp[7] = fccp[4];

			fccm[3] = fccm[1];
			fccm[5] = fccm[4];
			fccm[7] = fccm[4];
		}
	}

	//now nh3
	system("ls ../common/n_skimsV2/nh3 >> nh3tempfile.dat");
	fileList.open("nh3tempfile.dat");
	Float_t npraw[numELevels][numQBins][numWBins], nmraw[numELevels][numQBins][numWBins];
	Float_t fcnp[numELevels], fcnm[numELevels] = 0.0;//fc+ and fc- for nh3
	zero3DArray(npraw,numELevels,numQBins,numWBins);
	zero3DArray(nmraw,numELevels,numQBins,numWBins);
	zero1DArray(fcnp,numELevels);
	zero1DArray(fcnm,numELevels);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skimsV2/nh3/%s",cutname);
		infile.open(filename);
		bool pos;
		int index = findIndex(fileToRead.substr(4,3), fileToRead.substr(6,7));
		if(index == -1)
			continue;
		int dataSet = getSet(fileToRead.substr(4,3), masterList[index][1], masterList[index][2], fulldivide);
		if(dataSet == 9)
			continue;
		if(fileToRead.find("np") != string::npos)
			pos = true;
		else
			pos = false;
		
		for(Int_t qBin = 0; qBin < numQBins; qBin++)
		{
			for(Int_t wBin = 0; wBin < numWBins; wBin++){
				infile >> data;
				if(pos)
					npraw[dataSet][qBin][wBin] += atof(data);
				else
					nmraw[dataSet][qBin][wBin] += atof(data);
			}
			infile >> data;
			if(pos)
				fcnp[dataSet] += atof(data);
			else
				fcnm[dataSet] += atof(data);
			
		}
		infile.close();
		infile.clear();
	}
	fileList.close();
	fileList.clear();
	system("rm nh3tempfile.dat");

	//some numbers for this particular method
	Float_t fcntot[numELevels];
	Float_t fcctot[numELevels];
	Float_t ntot[numELevels][numQBins][numWBins], ctot[numELevels][numQBins][numWBins];

	//normalized counts
	Float_t cnorm[numELevels][numQBins][numWBins], nnorm[numELevels][numQBins][numWBins];

	for(Int_t e = 0; e < numELevels; e++)
	{
		fcntot[e] = fcnp[e] + fcnm[e];
		fcctot[e] = fccp[e] + fccm[e];
		for(Int_t i = 0; i < numQBins; i++)
		{
			for(Int_t j = 0; j < numWBins; j++)
			{
				ntot[e][i][j] = npraw[e][i][j] + nmraw[e][i][j];
				ctot[e][i][j] = cpraw[e][i][j] + cmraw[e][i][j];
				//cnorm[e][i][j] = (fcntot[e] / fcctot[e]) * (cpraw[e][i][j] + cmraw[e][i][j]);
				//nnorm[e][i][j] = npraw[e][i][j] + nmraw[e][i][j];
				//trying out a new, more normal method that should hopefully fix error bars
				cnorm[e][i][j] = ctot[e][i][j] / fcctot[e];
				nnorm[e][i][j] = ntot[e][i][j] / fcntot[e];
			}
		}
	}

	//now we find scaling factor
	
	cout << "Calculating multiplier..." << endl;
	Float_t multiplier[numELevels];
	for(Int_t e = 0; e < numELevels; e++){
		if(!fulldivide && (e > 2))
			break;
		Float_t bestdiff = 9999999999.99;
		for(Float_t mult = 0.5; mult < 1.5; mult += 0.001)
		{
			Float_t diff = fitfunction(cnorm, nnorm, mult, e);
			if(diff < bestdiff)
			{
				bestdiff = diff;
				multiplier[e] = mult;
			}
		}
		cout << "Using a multiplier of " << multiplier[e] << " for level " << e << endl;
	}

	//load up actual DF
	ifstream jlabdf;
	jlabdf.open("../common/F_DFfinal_set6.txt.0");
	Float_t realdf[numQBins][300];
	Float_t maxval = 0.0; //we'll need this later to scale the graph
	for(int i=0;i<numQBins;i++){for(int j=0;j<300;j++){realdf[i][j]=0.0;}} //quick custom zeroing

	for(int i = 0; i<301; i++)
	{
		for(int j = 0; j<42; j++)
		{
			jlabdf >> data;
			if((i == 0)  //first line
			|| (j == 0)) //first column
			{
				continue;
			}
			realdf[j+11][i-1] = atof(data);
			if(maxval < realdf[j+11][i-1])
				maxval = realdf[j+11][i-1];
		}
	}

	//new section: find "desired" multiplier
	/*Float_t desiredmult[3];
	for(int e = 0; e < numELevels; e++){
	Float_t average = 0.0;
	Float_t avecounter = 0.0;
	for(int i = 10; i < 40; i++)
	{
		for(int j = 18; j < 60; j++)
		{
			if((cnorm[e][i][j] != 0.0)&&(nnorm[e][i][j] != 0.0)){
			average += ((1-realdf[i][j*4])*nnorm[e][i][j])/cnorm[e][i][j];
			avecounter++;}
		}
	}
	desiredmult[e] = average / avecounter;
	cout << "Desired level " << e << " multiplier: " << desiredmult[e] << endl;}*/

	//calculate DF approx
	Float_t fdf[numELevels][numQBins][numWBins];
	for(int e = 0; e < numELevels; e++){
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			if(nnorm[e][i][j] != 0.0)
				fdf[e][i][j] = (nnorm[e][i][j] - (cnorm[e][i][j] * multiplier[e]))/nnorm[e][i][j];
			else
				fdf[e][i][j] = 0.0;
			if(fdf[e][i][j] > maxval)
				maxval = fdf[e][i][j];
		}
	}}

	maxval *= 1.10;
	cout << "Max = " << maxval << endl;

	//this annoying thing
	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};

	//utility arrays
	Float_t zero[numWBins], xaxis[numWBins], realxaxis[300];
	for(int i = 0; i < numWBins; i++)
	{
		zero[i] = 0.0;
		xaxis[i] = i * 0.04;
	}
	for(int i = 0; i < 300; i++)
	{
		realxaxis[i] = i * 0.01;
	}

	//I'M DOING THESE DAMN ERROR BARS HERE AND NOW
	Float_t fdferr[numELevels][numQBins][numWBins];
	for(int e = 0; e < numELevels; e++){
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			if(nnorm[e][i][j]==0)
			{
				fdferr[e][i][j] = 0.0;
				continue;
			}
			Float_t nh3errsq = ntot[e][i][j] / (fcntot[e] * fcntot[e]);
			Float_t carberrsq = (ctot[e][i][j] * multiplier[e] * multiplier[e]) / (fcctot[e] * fcctot[e]);
			Float_t toperrsq = nh3errsq + carberrsq;
			Float_t term1 = toperrsq / pow((nnorm[e][i][j] - multiplier[e]*cnorm[e][i][j]),2.0);
			Float_t term2 = nh3errsq / (nnorm[e][i][j] * nnorm[e][i][j]);
			fdferr[e][i][j] = (term1 + term2) * (fdf[e][i][j] * fdf[e][i][j]);
			if(fdferr[e][i][j] > 2.0)
				fdferr[e][i][j] = 0.0;
		}
	}}

	if(outputinstead)
	{
		output(fdf, fulldivide);
		outputerror(fdferr, fulldivide);
		cout << "Results saved!" << endl;
		return;
	}

	
	//let's graph this mf
	for(Int_t i = 0; i < 9; i++)
	{
		/*Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "Dilution factor approximation", 20*i, 20*i, 1600, 800);
		c->Divide(3,2);*/
		//c->SetFillColor(42);
		for(Int_t j = 0; j < 6; j++)
		{
			Int_t graphNum = (i*6) + j;
			if(graphNum >= numQBins)
				break;
			if((graphNum <= 16)||(graphNum >= 38))
				break;
			Char_t cName[5];
			sprintf(cName,"c%d",graphNum);
			TCanvas *cc = new TCanvas(cName, "Dilution Factor Approximation", 10*graphNum, 10*graphNum, 1200, 800);
			char blankName[10];
			sprintf(blankName,"%.3f < Q^2 < %.3f GeV",edges[graphNum],edges[graphNum+1]);
			TH2D *blank = new TH2D(blankName, blankName, 75, 1.0, 2.2, 500, 0, 0.4);
			Float_t y0[numWBins], y1[numWBins], y2[numWBins], y3[numWBins], y4[numWBins], y5[numWBins],
					y6[numWBins], y7[numWBins], y8[numWBins],
			        y0error[numWBins], y1error[numWBins], y2error[numWBins], y3error[numWBins],
			        y4error[numWBins], y5error[numWBins], y6error[numWBins], y7error[numWBins],
			        y8error[numWBins], yreal[300];
			for(int bin = 0; bin<numWBins; bin++)
			{
				y0[bin] = fdf[0][graphNum][bin];
				y0error[bin] = pow(fdferr[0][graphNum][bin],0.5);
				y1[bin] = fdf[1][graphNum][bin];
				y1error[bin] = pow(fdferr[1][graphNum][bin],0.5);
				if(y1error[bin] >= 0.1)
					y1error[bin] = 0.0;
				y2[bin] = fdf[2][graphNum][bin];
				y2error[bin] = pow(fdferr[2][graphNum][bin],0.5);
				y3[bin] = fdf[3][graphNum][bin];
				y3error[bin] = pow(fdferr[3][graphNum][bin],0.5);
				y4[bin] = fdf[4][graphNum][bin];
				y4error[bin] = pow(fdferr[4][graphNum][bin],0.5);
				if(y4error[bin] > 0.1)
					y4error[bin] = 0.0;
				y5[bin] = fdf[5][graphNum][bin];
				y5error[bin] = pow(fdferr[5][graphNum][bin],0.5);
				y6[bin] = fdf[6][graphNum][bin];
				y6error[bin] = pow(fdferr[6][graphNum][bin],0.5);
				y7[bin] = fdf[7][graphNum][bin];
				y7error[bin] = pow(fdferr[7][graphNum][bin],0.5);
				y8[bin] = fdf[8][graphNum][bin];
				y8error[bin] = pow(fdferr[8][graphNum][bin],0.5);
				if(y8error[bin] > 0.1)
					y8error[bin] = 0.0;
			}
			for(int bin = 0; bin < 300; bin++)
			{
				yreal[bin] = realdf[graphNum][bin];
			}
			TGraphErrors *gfdf0 = new TGraphErrors(numWBins, xaxis, y0, zero, y0error);
			TGraphErrors *gfdf1 = new TGraphErrors(numWBins, xaxis, y1, zero, y1error);
			TGraphErrors *gfdf2 = new TGraphErrors(numWBins, xaxis, y2, zero, y2error);
			TGraphErrors *gfdf3 = new TGraphErrors(numWBins, xaxis, y3, zero, y3error);
			TGraphErrors *gfdf4 = new TGraphErrors(numWBins, xaxis, y4, zero, y4error);
			TGraphErrors *gfdf5 = new TGraphErrors(numWBins, xaxis, y5, zero, y5error);
			TGraphErrors *gfdf6 = new TGraphErrors(numWBins, xaxis, y6, zero, y6error);
			TGraphErrors *gfdf7 = new TGraphErrors(numWBins, xaxis, y7, zero, y7error);
			TGraphErrors *gfdf8 = new TGraphErrors(numWBins, xaxis, y8, zero, y8error);
			TGraph *grealdf = new TGraph(300, realxaxis, yreal);
			//gfdf1->setTitle(blankName);
			/*gfdf1->SetLineColor(3);
			gfdf2->SetLineColor(4);
			gfdf3->SetLineColor(6);
			gfdf4->SetLineColor(7);
			gfdf5->SetLineColor(8);
			gfdf6->SetLineColor(9);
			gfdf7->SetLineColor(28);
			gfdf8->SetLineColor(40);*/
			gfdf4->SetLineColor(kBlue);
			gfdf8->SetLineColor(kGreen);
			grealdf->SetLineColor(kRed);
			//c->cd(j+1);
			//gPad->SetFillColor(42);
			blank->Draw();
			//gfdf0->Draw("P");
			gfdf1->Draw("P");
			//gfdf2->Draw("P");
			if(fulldivide)
			{
				//gfdf3->Draw("P");
				gfdf4->Draw("P");
				//gfdf5->Draw("P");
				//gfdf6->Draw("P");
				//gfdf7->Draw("P");
				gfdf8->Draw("P");
			}
			grealdf->Draw("L"); //this is the red one
			leg = new TLegend(0.7, 0.77, 0.89, 0.89);
			if(fulldivide)
			{
				//leg->AddEntry(gfdf0,"1p1 GeV, Target 11","l");
				leg->AddEntry(gfdf1,"1p3 GeV, Target 1","l");
				//leg->AddEntry(gfdf2,"1p3 GeV, Target 5","l");
				//leg->AddEntry(gfdf3,"1p3 GeV, Target 11","l");
				leg->AddEntry(gfdf4,"2 GeV, Target 1","l");
				//leg->AddEntry(gfdf5,"2 GeV, Target 11","l");
				//leg->AddEntry(gfdf6,"2p2 GeV, Target 1","l");
				//leg->AddEntry(gfdf7,"2p2 GeV, Target 5","l");
				leg->AddEntry(gfdf8,"3 GeV, Target 1","l");
			}
			else
			{
				leg->AddEntry(gfdf0,"1p1 and 1p3 GeV","l");
				leg->AddEntry(gfdf1,"2 GeV","l");
				leg->AddEntry(gfdf2,"3 GeV","l");
			}
			leg->AddEntry(grealdf,"JLab DF","l");
			/*if(j == 1)
				leg->Draw();*/
			gStyle->SetOptStat(0); //remove annoying stat box
			leg->Draw();
			gStyle->SetEndErrorSize(5);
			//cout << edges[graphNum] << " -> " << edges[graphNum+1] << endl;
		}
	}
	//outputerror(fdferr, fulldivide);
}

void zero1DArray(Float_t a[], Int_t size)
{
	for(int i = 0; i < size; i++)
	{
		a[i] = 0.0;
	}
}

void zero2DArray(Float_t a[][75], Int_t s1, Int_t s2)
{
	for(int i = 0; i < s1; i++)
	{
		for(int j = 0; j < s2; j++)
		{
			a[i][j] = 0.0;
		}
	}
}

void zero3DArray(Float_t x[][53][75], Int_t d1, Int_t d2, Int_t d3) //this is utter madness
{
	for(Int_t i=0; i < d1; i++) {
		for(Int_t j=0; j < d2; j++) {
			for(Int_t k=0; k < d3; k++) {
				x[i][j][k] = 0.0;
			}
		}
	}
}

int getEnergyCarb(string inp, bool fulldivide) //used ONLY for carbon
{
	char e[3];
	strcpy(e,inp.c_str()); //what is even going on anymore
	if(!fulldivide)
	{
		if(strcmp(e,"1p1")==0)
			return 0;
		else if(strcmp(e,"1p3")==0)
			return 0;
		else if(strcmp(e,"1p5")==0) //1p5 appears to be bad. Removing for now
			return 9;
		else if(strcmp(e,"2_5")==0) //this relies on the run number starting with 5 cuz I'm a bad programmer
			return 1;
		else if(strcmp(e,"2p2")==0) //cutting this one too cuz its taking results from bad to worse
			return 9;
		else if(strcmp(e,"3_5")==0) 
			return 2;
		else
		{
			cout << "ERROR: Invalid energy level " << e << endl;
			cout << "(Does the run number start with something besides 5?)" << endl;
			return 6;
		}
	}
	else if (fulldivide)
	{
		if(strcmp(e,"1p1")==0)
			return 0;
		else if(strcmp(e,"1p3")==0)
			return 1;
		else if(strcmp(e,"1p5")==0) //1p5 appears to be bad. Removing for now
			return 9;
		else if(strcmp(e,"2_5")==0) //this relies on the run number starting with 5 cuz I'm a bad programmer
			return 4;
		else if(strcmp(e,"2p2")==0) //cutting this one too cuz its taking results from bad to worse
			return 6;
		else if(strcmp(e,"3_5")==0) 
			return 8;
		else
		{
			cout << "ERROR: Invalid energy level " << e << endl;
			cout << "(Does the run number start with something besides 5?)" << endl;
			return 9;
		}
	}
	else
	{
		cout << "何?" << endl;
	}
}

int getEnergy(string inp) //nh3 only
{
	char e[3];
	strcpy(e,inp.c_str()); //what is even going on anymore
	if(strcmp(e,"1p1")==0)
		return 0;
	else if(strcmp(e,"1p3")==0)
		return 0;
	else if(strcmp(e,"1p5")==0)
		return 9;
	else if(strcmp(e,"2_5")==0) //this relies on the run number starting with 5
		return 1;
	else if(strcmp(e,"2p2")==0)
		return 1;
	else if(strcmp(e,"3_5")==0) //because I am a very bad programmer
		return 2;
	else
	{
		cout << "ERROR: Invalid energy level " << e << endl;
		cout << "(Does the run number start with something besides 5?)" << endl;
		return 9;
	}
}

int getSet(string energy, int targ, int tor, bool fulldivide) //used ONLY for nh3
{
	int e = getEnergy(energy);
	if(!fulldivide)
		return e;
	char derp[3];
	strcpy(derp,energy.c_str());
	if(tor != -1) //there are no meaningful groups of positive tor, plus zero is the error one
		return 9; //-1 = "skip this one"
	if(e == 0)
	{
		if(targ==1)
			return 1;
		else if(targ==5)
			return 2;
		else if(targ==11){
			if(strcmp(derp,"1p1")==0)
				return 0;
			else //1p3
				return 3;
		}
	}
	else if(e == 1)
	{
		if(targ==1){
			if(strcmp(derp,"2_5")==0)
				return 4;
			else //2p2
				return 6;
		}
		else if(targ==11)
			return 5;
		else if(targ==5)
			return 7;
		else
			return 9;
	}
	else if(e == 2)
		return 8;
	return 9; //catch-all "didn't fit another box"
}

int getSetCarb(string energy, int targ, int tor, bool fulldivide)
{ 
	int e = getEnergyCarb(energy, fulldivide);
	if(!fulldivide)
		return e;

	char derp[3];
	strcpy(derp,energy.c_str());
	if(tor != -1) //there are no meaningful groups of positive tor, plus zero is the error one
		return 9; //-1 = "skip this one"
	if(e == 0)//1p1
	{
		if(targ==4)
			return 0;
	}
	if(e == 1)//1p3
	{
		if(targ==4)
			return 1; //duplicate to 3
		else if(targ==6){
			return 2;
		}
	}
	if(e == 4)//2
	{
		if(targ==4)
			return 4;//duplicate to 5
	}
	if(e == 6)//2p2
	{
		if(targ==4)
			return 6;
		else if(targ==6)
			return 7;
	}
	if(e == 8)//3
	{
		return 8;
	}
	return 9; //catch-all "didn't fit another box"
}

//finds the index of a run on the master list. Relies on the master list not changing
int findIndex(string energy, string runChunk)
{
	char e[3], bit[1];
	int run;
	strcpy(e,energy.c_str()); //what is even going on anymore
	strcpy(bit,energy.substr(2,1).c_str());//just to identify how to align the run
	if(strcmp(bit,"5")==0)//single character energy
		istringstream(runChunk.substr(0,5)) >> run;
	else //three character energy
		istringstream(runChunk.substr(2,5)) >> run;

	if(strcmp(e,"1p1")==0)
		return (run - 52051);
	else if(strcmp(e,"1p3")==0)
		return 255 + (run - 51158);
	else if(strcmp(e,"1p5")==0)
		return -1;
	else if(strcmp(e,"2_5")==0) //this relies on the run number starting with 5
		return 574 + (run - 51487);
	else if(strcmp(e,"2p2")==0)
		return 669 + (run - 50991);
	else if(strcmp(e,"3_5")==0) //because I am a very bad programmer
		return 825 + (run - 50744);
	else
	{
		cout << "ERROR: Invalid energy level " << e << endl;
		cout << "(Does the run number start with something besides 5?)" << endl;
		return -1;
	}
}

int findIndexCarb(string energy, string runChunk, int masterListCarb[72][3])
{
	char e[3], bit[1];
	int run;
	strcpy(e,energy.c_str()); //what is even going on anymore
	strcpy(bit,energy.substr(2,1).c_str());//just to identify how to align the run
	if(strcmp(bit,"5")==0)//single character energy
		istringstream(runChunk.substr(0,5)) >> run;
	else //three character energy
		istringstream(runChunk.substr(2,5)) >> run;
	for(int i = 0; i < 72; i++)
	{
		if(masterListCarb[i][0] == run)
			return i;
	}
	cout << "ERROR: Invalid energy level " << e << endl;
	cout << "(Does the run number start with something besides 5?)" << endl;
	return -1; //can't find it
}

int energyFromSet(int s, bool sepBySet)
{
	if(!sepBySet)
	{
		if((s >= 0)&&(s <= 3))
			return 0;
		else if((s >= 4)&&(s <= 7))
			return 1;
		else if(s == 8)
			return 2;
		else
			return -1;
	}
	else
	{
		if(s == 0)
			return 0;
		else if(s <= 3)
			return 1;
		else if(s <= 5)
			return 3; //no 1p5 forever I guess
		else if(s <= 7)
			return 4;
		else if (s == 8)
			return 5;
		else
			return -1;
	}
}

//test how close carb and nh3 are when a multiplier is applied to carb
Float_t dfapprox::fitfunction(Float_t c[][53][75], Float_t n[][53][75], Float_t mult, Int_t e_index)
{
	Float_t sum = 0.0;
	int max;
	if(e_index == 2)
		max = 22;
	else
		max = 21;
	for(int i = 0; i < 53; i++)
	{
		for(int j = 12; j < max; j++)           //originally 12 and 21
		{
			sum += sqrt(pow(mult*c[e_index][i][j] - n[e_index][i][j],2));
		}
	}
	return sum;
}

/*Float_t dfapprox::errorfunction(Float_t calculated[75], Float_t actual[300])
{
	Float_t average = 0.0;
	for(int i = 0; i < 75; i++)
	{
		Float_t got = calculated[i];
		int j = i*4; //just making the next line easier
		Float_t expected = (actual[j]+actual[j+1]+actual[j+2]+actual[j+3])/4.0;
		avarage += ((got - expected)/expected * 100.0
	}
	average /= 75.0
	return average;
}*/

void getMasterList(int masterInfo[1080][3])
{
	ifstream nh3list;
	nh3list.open("../common/masterLists/masterList_nh3.dat");
	char data[100];
	for(int i=0;i<9;i++){nh3list >> data;} //skip first line
	for(int i = 0; i < 1080; i++)
	{
		nh3list >> data; //skip """run""" column
		nh3list >> data; //skip calib
		nh3list >> data; //this is the run
		masterInfo[i][0] = atoi(data);
		nh3list >> data; //this is the target
		masterInfo[i][1] = atoi(data); //some testing shows the "NA" and "4+7" cases shouldn't be a problem
		nh3list >> data; //skip ebeam
		nh3list >> data; //this is torus;
		float num = atof(data);
		if(num > 0.0)
			masterInfo[i][2] = 1;
		else if(num < 0.0)
			masterInfo[i][2] = -1;
		else
			masterInfo[i][2] = 0;
		nh3list >> data; //who knows what this is but we're skipping it
		nh3list >> data; //skip tgt_spin
		nh3list >> data; //skip #_total_triggers
	}
}

void getMasterListCarb(int masterInfoCarb[72][3])
{
	ifstream carblist;
	carblist.open("../common/masterLists/masterList_carb.dat");
	char data[100];
	for(int i = 0; i < 72; i++)
	{
		carblist >> data; //this is the run
		masterInfoCarb[i][0] = atoi(data);
		carblist >> data; //this is the target
		masterInfoCarb[i][1] = atoi(data); //again, the plus ones shouldn't be a problem
		carblist >> data;
		float num = atof(data);
		if(num > 0.0)
			masterInfoCarb[i][2] = 1;
		else if(num < 0.0)
			masterInfoCarb[i][2] = -1;
		else
			masterInfoCarb[i][2] = 0;
	}
}

/*
 * To be used instead of graphing, to output the resultant approximated DF to a text file to be used in asymmetry
 */
void dfapprox::output(Float_t fdf[9][53][75], Bool_t fulldivide)
{
	ofstream outfile;
	int max;
	if(!fulldivide)
	{
		outfile.open("../common/dfapprox.dat");
		max = 3;
	}
	else
	{
		outfile.open("../common/dfapprox_fulldivide.dat");
		max = 9;
	}
	for(int e = 0; e < max; e++)
	{
		for(int i = 0; i < 53; i++)
		{
			for(int j = 0; j < 75; j++)
			{
				outfile << fdf[e][i][j] << " ";
			}
			outfile << endl;
		}
	}
}

/*
 * For use in systerr
 */
void outputerror(Float_t fdferr[9][54][75], Bool_t fulldivide)
{
	ofstream outfile;
	int max;
	if(!fulldivide){
		outfile.open("../common/fdferr.dat");
		max = 3;
	}
	else{
		outfile.open("../common/fdferr_fulldivide.dat");
		max = 9;
	}
	for(int e = 0; e < 9; e++)
	{
		for(int i = 0; i < 53; i++)
		{
			for(int j = 0; j < 75; j++)
			{
				Float_t num = pow(fdferr[e][i][j],0.5);
				outfile << num << " ";
			}
			outfile << endl;
		}
	}
}
