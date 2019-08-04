#include<stdlib.h>
#include<stdio.h>

#define ngraph_cxx

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

#include "ngraph.h"

using namespace std;

//#define E 2.0;
//#define M 0.938;
//#define PI 3.14159

int main(int argc, char *argv[]){
	
}

void ngraph::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout << jentry << endl;
      // if (Cut(ientry) < 0) continue;
   }
}

/*
 * Doesn't actually graph, but crunches all the numbers and creates a single file for the grapher to read
 */

void ngraph::parse()
{
	Char_t data[100], filename[100], cutname[100];
	string fileToRead;	
	ifstream fileList, infile;
	const Int_t numWBins = 300;
	const Int_t numQBins = 53;

	//carb first
	system("touch tempfile.dat");
	system("ls ../common/n_skims/carb >> carbtempfile.dat");
	fileList.open("carbtempfile.dat");
	Float_t cpraw[numQBins][numWBins], cmraw[numQBins][numWBins];
	Float_t fccp = 0.0;//fc+ for carb
	Float_t fccm = 0.0;//fc- for carb
	zero2DArray(cpraw,numQBins,numWBins);
	zero2DArray(cmraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skims/carb/%s",cutname);
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
	system("ls ../common/n_skims/emp >> emptempfile.dat");
	fileList.open("emptempfile.dat");
	Float_t epraw[numQBins][numWBins], emraw[numQBins][numWBins];
	Float_t fcep = 0.0;//fc+ for emp
	Float_t fcem = 0.0;//fc- for emp
	zero2DArray(epraw,numQBins,numWBins);
	zero2DArray(emraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skims/emp/%s",cutname);
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
	system("ls ../common/n_skims/nh3 >> nh3tempfile.dat");
	fileList.open("nh3tempfile.dat");
	Float_t npraw[numQBins][numWBins], nmraw[numQBins][numWBins];
	Float_t fcnp = 0.0;//fc+ for nh3
	Float_t fcnm = 0.0;//fc- for nh3
	zero2DArray(npraw,numQBins,numWBins);
	zero2DArray(nmraw,numQBins,numWBins);

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/n_skims/nh3/%s",cutname);
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

	//////////////////////////////EMERGENCY EXTRA SECTION/////////////////////////////////
	const Int_t reducedWBins = 75;
	const int six = 6;
	Float_t smallAverage[six];
	for(int i = 0; i < six; i++){smallAverage[i] = 0.0;}
	Int_t counter = 0;
	Float_t cprawR[numQBins][reducedWBins], cmrawR[numQBins][reducedWBins],
		eprawR[numQBins][reducedWBins], emrawR[numQBins][reducedWBins],
		nprawR[numQBins][reducedWBins], nmrawR[numQBins][reducedWBins];
	for(Int_t i = 0; i < numQBins; i++)
	{
		for(Int_t j = 0; j < numWBins; j++)
		{
			counter++;
			smallAverage[0] += cpraw[i][j];
			smallAverage[1] += cmraw[i][j];
			smallAverage[2] += epraw[i][j];
			smallAverage[3] += emraw[i][j];
			smallAverage[4] += npraw[i][j];
			smallAverage[5] += nmraw[i][j];
			if(counter == 4)
			{
				counter = 0;
				int index;
				if(j <= 1)
					index = 0;
				else
					index = j/4;
				cprawR[i][index] = smallAverage[0] / 4.0;
				cmrawR[i][index] = smallAverage[1] / 4.0;
				eprawR[i][index] = smallAverage[2] / 4.0;
				emrawR[i][index] = smallAverage[3] / 4.0;
				nprawR[i][index] = smallAverage[4] / 4.0;
				nmrawR[i][index] = smallAverage[5] / 4.0;
				
				for(int q=0;q<six;q++){smallAverage[q]=0.0;}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////

	//some numbers for this particular method
	Float_t fcntot = fcnp + fcnm;
	Float_t fcetot = fcep + fcem;
	Float_t fcctot = fccp + fccm;

	Char_t outName[100] = "../common/n_data_final.dat";
	ofstream outfile;
	outfile.open(outName);
	//first carb
	for(Int_t i = 0; i < numQBins; i++)
	{
		for(Int_t j = 0; j < reducedWBins; j++)
		{
			Float_t num = (fcntot / fcctot) * (cprawR[i][j] + cmrawR[i][j]);
			outfile << num << " ";
		}
		outfile << endl;
	}

	//then emp
	for(Int_t i = 0; i < numQBins; i++)
	{
		for(Int_t j = 0; j < reducedWBins; j++)
		{
			Float_t num = (fcntot / fcetot) * (eprawR[i][j] + emrawR[i][j]);
			outfile << num << " ";
		}
		outfile << endl;
	}

	//then nh3
	for(Int_t i = 0; i < numQBins; i++)
	{
		for(Int_t j = 0; j < reducedWBins; j++)
		{
			Float_t num = (nprawR[i][j] + nmrawR[i][j]);
			outfile << num << " ";
		}
		outfile << endl;
	}
	outfile.close();
}

void ngraph::zero1DArray(Float_t a[], Int_t size)
{
	for(int i = 0; i < size; i++)
	{
		a[i] = 0.0;
	}
}

void ngraph::zero2DArray(Float_t a[][300], Int_t s1, Int_t s2)
{
	for(int i = 0; i < s1; i++)
	{
		for(int j = 0; j < s2; j++)
		{
			a[i][j] = 0.0;
		}
	}
}


/*
 * Where the graphing happens
 */
void ngraph::graph()
{
	ifstream file;
	file.open("../common/n_data_final.dat");
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	const Float_t multiplier = 1.059;//Use findmultiplier() to verify
	Char_t data[100];
	Float_t xaxis[numWBins];
	Float_t carb[numQBins][numWBins], emp[numQBins][numWBins], nh3[numQBins][numWBins], mcarb[numQBins][numWBins];

	for(int i=0;i<75;i++)xaxis[i]=(float)i*0.01;//create x-axis

	Float_t maxval[numQBins];

	//first carb
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			carb[i][j] = atof(data);
			mcarb[i][j] = multiplier * atof(data);
		}
	}
	//then emp
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			emp[i][j] = atof(data);
		}
	}
	//then nh3
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < 300; j++)
		{
			file >> data;
			nh3[i][j] = atof(data);
			if(atof(data) > maxval[i]) maxval[i] = atof(data);
		}
	}
	for(int i=0;i<numQBins;i++)maxval[i] *= 1.1;

	

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
			Float_t thismax = maxval[graphNum];
			TH2D *blank = new TH2D(blankName, blankName, 75, 0, 3, 500, 0, thismax);
			Float_t ye[numWBins], yc[numWBins], yn[numWBins], ym[numWBins];
			for(int bin = 0; bin<numWBins; bin++)
			{
				ye[bin] = emp[graphNum][bin];
				yc[bin] = carb[graphNum][bin];
				yn[bin] = nh3[graphNum][bin];
				ym[bin] = mcarb[graphNum][bin];
			}
			TGraph *gemp = new TGraph(numWBins, xaxis, ye);
			TGraph *gcarb = new TGraph(numWBins, xaxis, yc);
			gcarb->SetLineColor(kBlue);
			TGraph *gnh3 = new TGraph(numWBins, xaxis, yn);
			gnh3->SetLineColor(kRed);
			TGraph *gmcarb = new TGraph(numWBins, xaxis, ym);
			gmcarb->SetLineColor(kGreen);
			c->cd(j+1);
			blank->Draw();
			gemp->Draw("same");
			gcarb->Draw("same");
			gnh3->Draw("same");
			gmcarb->Draw("same");
		}
	}
}

//test how close carb and nh3 are when a multiplier is applied to carb
Float_t ngraph::fitfunction(Float_t c[][75], Float_t n[][75], Float_t mult)
{
	Float_t sum = 0.0;
	for(int i = 0; i < 53; i++)
	{
		for(int j = 13; j < 22; j++)
		{
			sum += sqrt(pow(mult*c[i][j] - n[i][j],2));
		}
	}
	cout << mult << ":" << sum << endl;
	return sum;
}

void ngraph::findMultiplier()
{
	ifstream file;
	file.open("../common/n_data_final.dat");
	Char_t data[100];
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	Float_t carb[numQBins][numWBins], emp[numQBins][numWBins], nh3[numQBins][numWBins];
	//first carb
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			carb[i][j] = atof(data);
		}
	}
	//then emp
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			emp[i][j] = atof(data);
		}
	}
	//then nh3
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			nh3[i][j] = atof(data);
		}
	}

	Float_t bestdiff = 9999999999.99;
	Float_t bestmult = 0.0;
	for(Float_t mult = 0.5; mult < 1.5; mult += 0.001)
	{
		Float_t diff = fitfunction(carb, nh3, mult);
		if(diff < bestdiff)
		{
			bestdiff = diff;
			bestmult = mult;
		}
	}
	cout << bestmult << endl;
}

void ngraph::prelimdiff(Float_t multiplier)
{
	ifstream file;
	file.open("../common/n_data_final.dat");
	Char_t data[100];
	const Int_t numQBins = 53;
	const Int_t numWBins = 75;
	Float_t carb[numQBins][numWBins], emp[numQBins][numWBins], nh3[numQBins][numWBins], mcarb[numQBins][numWBins];
	//first carb
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			carb[i][j] = atof(data);
			mcarb[i][j] = multiplier * atof(data);
		}
	}
	//then emp
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			emp[i][j] = atof(data);
		}
	}
	//then nh3
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			file >> data;
			nh3[i][j] = atof(data);
		}
	}

	ofstream outfile;
	outfile.open("../common/prelim_difference.dat");
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			Float_t num;
			if(nh3[i][j] != 0.0)
				num = (nh3[i][j] - mcarb[i][j])/nh3[i][j];
			else
				num = 0.0;
			outfile << num << " ";
		}
		outfile << endl;
	}
}
