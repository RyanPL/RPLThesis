#define makeGraphs_cxx
#include "makeGraphs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>


void zero2DArray(Double_t x[][75], Int_t d1, Int_t d2) //sets all positions on a 2d array to zero
{
	for(Int_t i=0; i < d1; i++) {
		for(Int_t j=0; j < d2; j++) {
			x[i][j] = 0.0;
		}
	}
}

void zero1DArray(Double_t x[75], Int_t d2) //sets all positions on a 2d array to zero
{
	for(Int_t j=0; j < d2; j++) {
		x[j] = 0.0;
	}
}

void zero3DArray(Double_t x[][54][75], Int_t d1, Int_t d2, Int_t d3) //this is utter madness
{
	for(Int_t i=0; i < d1; i++) {
		for(Int_t j=0; j < d2; j++) {
			for(Int_t k=0; k < d3; k++) {
				x[i][j][k] = 0.0;
			}
		}
	}
}

int getEnergy(string inp) //simple method for getting which E bin to put things in
{
	char e[3];
	strcpy(e,inp.c_str()); //what is even going on anymore
	if(strcmp(e,"1p1")==0)
		return 0;
	else if(strcmp(e,"1p3")==0)
		return 0;
	else if(strcmp(e,"1p5")==0)
		return 3;
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
		return 3;
	}
}

int getSet(string energy, int targ, int tor)
{
	int e = getEnergy(energy);
	char derp[3];
	strcpy(derp,energy.c_str());
	if(tor != -1) //there are no meaningful groups of positive tor, plus zero is the error one
		return -1; //-1 = "skip this one"
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
			return -1;
	}
	else if(e == 2)
		return 8;
	return -1; //catch-all "didn't fit another box"
}

void makeGraphs::asymDS(Bool_t useapprox = true)
{
	if(useapprox)
		cout << "Using approximated DF..." << endl;
	else
		cout << "Using JLab DF..." << endl;

	//cout << "Do you want to enable switches? Only answer 'Yes' if you have run switchScanner with the current batch of results: ";
	Char_t data[100], filename[100], cutname[100];
	Bool_t flipThem = true;

	/*scanf("%s",&data);
	if(strcmp(data,"Yes")==0)
	{
		cout << "Roger that. Switching appropriate skims..." << endl;
		flipThem = true;
	}
	else
		cout << "Roger that. NOT switching appropriate skims..." << endl;*/
	ifstream flips;
	flips.open("../common/switchesV4.dat");

	Bool_t raw = false;
	/*cout << "Display raw asymmetry? Enter 'Yes' if so, anything else if no: ";
	scanf("%s",&data);
	if(strcmp(data,"Yes")==0)
	{
		cout << "Roger that. Displaying raw asymmetry..." << endl;
		raw = true;
	}
	else
		cout << "Roger that. Displaying double-spin asymmetry..." << endl;*/

	Bool_t sepBySet = true;
	/*cout << "Do you want to apply separated DF? Enter 'Yes' if so, anything else for energy level separation: ";
	scanf("%s",&data);
	if(strcmp(data,"Yes")==0)
	{
		cout << "Separating by set and recombining..." << endl;
		sepBySet = true;
	}
	else
		cout << "Separating by energy..." << endl;*/




	/*
	 * Deceptively small yet important bit where the master list is filled
	 */
	int masterList[1080][3]; //first index is run, second index is run num, target type, torus current
	getMasterList(masterList);

	string fileToRead;
	system("ls ../common/skimsV4 >> tempfile.dat");
	ifstream fileList;
	fileList.open("tempfile.dat");
	ifstream infile;

	//doing first file manually to get array sizes
	//this is a super clunky way of doing this but I can't think of anything better
	getline(fileList,fileToRead);
	strcpy(cutname,fileToRead.c_str()); //converting to char array cuz problems
	sprintf(filename,"../common/skimsV4/%s",cutname);
	//cout << filename << endl;
	infile.open(filename);
	infile >> data;
	Int_t numLines = atoi(data);
	const Int_t numSets = 9; //THIS USED TO BE NUMELEVELS see below for definitions
	const Int_t numQBins = numLines;
	const Int_t numWBins = 75;
	Double_t numerator[numSets][numQBins][numWBins], denominator[numSets][numQBins][numWBins], truedf[numQBins][300], real[numQBins][300], realCount[numQBins][300];
	Double_t pbpt[numSets] = {0.568, 0.571, 0.552, 0.535, 0.605, 0.636, 0.560, 0.597, 0.614};
	Int_t energyIndex;
	zero3DArray(numerator,numSets,numQBins,numWBins);
	zero3DArray(denominator,numSets,numQBins,numWBins);
	zero2DArray(truedf,numQBins,numWBins);
	zero2DArray(real,numQBins,numWBins);
	zero2DArray(realCount,numQBins,numWBins);

	/* SET DEFINITIONS - the new E levels
	0: 1p1, target 11, outbending
	1: 1p3, target 1, outbending
	2: 1p3, target 5, outbending
	3: 1p3, target 11, outbending
	4: 2, target 1, outbending
	5: 2, target 11,  outbending
	6: 2p2, target 1, outbending
	7: 2p2, target 5, outbending
	8: 3, target 1, outbending
	*/

	cout << "Getting secondary quantities (dilution factor, etc.)" << endl;

	//getting dilution factor
	ifstream dilution;
	dilution.open("../common/F_DFfinal_set6.txt.0");
	for(int i = 0; i<301; i++)
	{
		for(int j = 0; j<42; j++)
		{
			dilution >> data;
			if((i == 0)  //first line
			|| (j == 0)) //first column
			{
				continue;
			}
			truedf[j+11][i-1] = atof(data);
		}
	}
	dilution.close();
	//done with dilution factor

	//so it ends up using a 300 W bin DF doesn't work well with 75 W bin results
	//so I'm gonna average it into 75 bins until the new DF is working
	Double_t fdf[numSets][numQBins][numWBins];
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			Double_t ave = 0.0;
			for(int k = 0; k < 4; k++)
			{
				ave += truedf[i][(j*4)+k];
			}
			for(int s = 0; s < numSets; s++){
				fdf[s][i][j] = ave / 4.0;}
		}
	}
	//////////////////////////////////////////////////////////
	//    IN CASE YOU DIDN'T NOTICE THAT MEANS DELETE       //
	//     THE ABOVE SECTION AND REPLACE FDFR WITH FDF      //
	//      ONCE THE APPROXIMATED DF IS HOOKED UP           //
	//////////////////////////////////////////////////////////

	//getting approx DF
	ifstream dfstream;
	Double_t dfapprox[numSets][numQBins][numWBins];
	int max = 3;
	if(!sepBySet)
		dfstream.open("../common/dfapprox.dat");
	else{
		dfstream.open("../common/dfapprox_fulldivide.dat");
		max = 6;
	}
	for(int e = 0; e < max; e++)
	{
		for(int i = 0; i < 54; i++)
		{
			for(int j = 0; j < 75; j++)
			{
				dfstream >> data;
				dfapprox[e][i][j] = atof(data);
			}
		}
	}

	//swap em out if "useapprox" option is chose
	if(useapprox)
	{
		for(int e = 0; e < max; e++)
		{
			for(int i = 0; i < 54; i++)
			{
				for(int j = 0; j < 75; j++)
				{
					fdf[e][i][j] = dfapprox[e][i][j];
				}
			}
		}
	}

	cout << "Getting JLab results..." << endl;;

	//getting real results
	ifstream realResults;
	realResults.open("../common/Aparmod.set6.txt");
	for(int i = 0; i<301; i++)
	{
		for(int j = 0; j<42; j++)
		{
			realResults >> data;
			if((i == 0)  //first line
			|| (j == 0)) //first column
				continue;
			real[j+11][i-1] = atof(data);
			realCount[j+11][i-1]++;
		}
	}
	realResults.close();
	//done with real results

	/*for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			if(realCount[i][j] != 0)
				real[i][j] /= realCount[i][j]; //WHAT DOES THIS PART EVEN DO
		}
	}*/

	///Q^2 bins
	Float_t edges[55] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97, 13.1};
	//done getting bins

	//now back to that first file we're doing manually
	cout << "Collecting data...";
	energyIndex = getEnergy(fileToRead.substr(0,3));
	//cout << "Energy class: " << energyIndex << endl;
	if(energyIndex == 3)cout << "uh oh spaghetti o's" << endl;
	else{
		bool thisOneIsFlipped = false;
		flips >> data;
		if((strcmp(data,"F")==0)&&(flipThem))
			thisOneIsFlipped = true;

		int index = findIndex(fileToRead.substr(0,3),fileToRead.substr(2,7));
		if(index == -1)
			continue;

		int set = getSet(fileToRead.substr(0,3), masterList[index][1], masterList[index][2]);
		//cout << "Set = " << set << endl;
		if(set == -1)
			continue;

		for(Int_t qBin = 0; qBin < numQBins; qBin++){
			Double_t points_temp[numWBins];
			Double_t sig_temp[numWBins];
			zero1DArray(points_temp,numWBins);
			zero1DArray(sig_temp,numWBins);
			for(Int_t wBin = 0; wBin < numWBins; wBin++) {
				infile >> data;
				Float_t readNum = atof(data);
				if(thisOneIsFlipped)
				{
					readNum *= -1.0;
				}
				points_temp[wBin] = readNum;
			}
			for(Int_t wBin = 0; wBin < numWBins; wBin++) {
				infile >> data;
				if(atof(data) == 0)
					sig_temp[wBin] = 0.0;
				else
					sig_temp[wBin] = pow(atof(data),-2);//actually n
				numerator[set][qBin][wBin] += (points_temp[wBin]*sig_temp[wBin]);
				denominator[set][qBin][wBin] += sig_temp[wBin];
			}
		}
	}
	infile.close();

	while(getline(fileList,fileToRead)){
		strcpy(cutname,fileToRead.c_str());
		sprintf(filename,"../common/skimsV4/%s",cutname);
//		cout << filename << endl; //just so's we knows
		infile.open(filename);
		infile >> data; //skip numlines
		energyIndex = getEnergy(fileToRead.substr(0,3));

		bool thisOneIsFlipped = false;
		flips >> data;
		if((strcmp(data,"F")==0)&&(flipThem))
			thisOneIsFlipped = true;

		int index = findIndex(fileToRead.substr(0,3),fileToRead.substr(2,7));
		if(index == -1)
			continue;

		int set = getSet(fileToRead.substr(0,3), masterList[index][1], masterList[index][2]);
		//is energyindex still relevant now due to the changes in the above line?
		if(set == -1)
			continue;

		for(Int_t qBin = 0; qBin < numQBins; qBin++){
			Double_t points_temp[numWBins];
			Double_t sig_temp[numWBins];
			zero1DArray(points_temp,numWBins);
			zero1DArray(sig_temp,numWBins);
			for(Int_t wBin = 0; wBin < numWBins; wBin++) {
				infile >> data;
				Float_t readNum = atof(data);
				if(thisOneIsFlipped)
				{
					readNum *= -1.0;
				}
				points_temp[wBin] = readNum;
			}
			for(Int_t wBin = 0; wBin < numWBins; wBin++) {
				infile >> data;
				if(atof(data) == 0)
					sig_temp[wBin] = 0.0;
				else
					sig_temp[wBin] = pow(atof(data),-2);//actually n
				numerator[set][qBin][wBin] += (points_temp[wBin]*sig_temp[wBin]);
				denominator[set][qBin][wBin] += sig_temp[wBin];
			}
		}
		infile.close();
	}
	cout << " Done!" << endl;
	fileList.close();
	system("rm tempfile.dat");

	Float_t x[numWBins], zero[numWBins], realx[300]; //utility arrays
	for(Int_t i = 0; i<numWBins; i++) {
		x[i] = (((Float_t)i)*3)/numWBins;
		zero[i] = 0;
	}
	for(Int_t i = 0; i<300; i++) {
		realx[i] = i * 0.01;
	}

	//big decisions: are we separating by set or by energy?
	//so basically, since I don't want to rebuild everything for this option
	//if it's by energy, I'm just gonna load it all into the first three bins of the first index
	//then the grapher will ignore the rest
	//im a good programmer dammit
	if(!sepBySet)
	{
		for(int i = 0; i < numQBins; i++)
		{
			for(int j = 0; j < numWBins; j++)
			{
				numerator[0][i][j] += numerator[1][i][j] + numerator[2][i][j] + numerator[3][i][j];
				denominator[0][i][j] += denominator[1][i][j] + denominator[2][i][j] + denominator[3][i][j];
				numerator[1][i][j] += numerator[4][i][j] + numerator[5][i][j] + numerator[6][i][j] + numerator[7][i][j];
				denominator[1][i][j] += denominator[4][i][j] + denominator[5][i][j] + denominator[6][i][j] + denominator[7][i][j];
				numerator[2][i][j] = numerator[8][i][j];
				denominator[2][i][j] = denominator[8][i][j];
			}
		}
		pbpt[0] = 0.5565; //these are all averages
		pbpt[1] = 0.5995;
		pbpt[2] = 0.614;
	}

	//surprise bitch, systerr is back
	Float_t ferr1[numQBins][numWBins], ferr2[numQBins][numWBins], ferr3[numQBins][numWBins],
			perr1[numQBins][numWBins], perr2[numQBins][numWBins], perr3[numQBins][numWBins];
	ifstream fe1, fe2, fe3, pe1, pe2, pe3;
	fe1.open("../common/systerr/last_fdf1.dat");
	fe2.open("../common/systerr/last_fdf2.dat");
	fe3.open("../common/systerr/last_fdf3.dat");
	pe1.open("../common/systerr/last_pbpt1.dat");
	pe2.open("../common/systerr/last_pbpt2.dat");
	pe3.open("../common/systerr/last_pbpt3.dat");
	for(int i = 0; i < numQBins; i++)
	{
		for(int j = 0; j < numWBins; j++)
		{
			fe1 >> data;
			ferr1[i][j] = atof(data);
			fe2 >> data;
			ferr2[i][j] = atof(data);
			fe3 >> data;
			ferr3[i][j] = atof(data);
			pe1 >> data;
			perr1[i][j] = atof(data);
			pe2 >> data;
			perr2[i][j] = atof(data);
			pe3 >> data;
			perr3[i][j] = atof(data);
		}
	}


	//this part is for the systerr analysis. Should be commented out normally
	/*ifstream fdferr;
	fdferr.open("../common/fdferr_fulldivide.dat");
	if(!sepBySet)
		return;
	for(int e = 0; e < 9; e++)//don't use this without fulldivide!
	{
		for(int i = 0; i < numQBins; i++)
		{
			for(int j = 0; j < numWBins; j++)
			{
				fdferr >> data;
				float errbit = atof(data);
				fdf[e][i][j] += errbit;
			}
		}
	}
	fdferr.close();/**/

	Float_t finalasym[3][numQBins][numWBins]; //for outputting
	Float_t finalerror[3][numQBins][numWBins];

	

	for(Int_t i = 0; i<9; i++) { //graphing loop
		/*Char_t cName[5];
		sprintf(cName,"c%d",i);
		TCanvas *c = new TCanvas(cName, "Double-spin Asymmetry", 20*i, 20*i, 1600, 800);
		c->Divide(3,2);*/ //no more 6-groupings!
		//c->SetFillColor(42);
		for(Int_t j = 0; j<6; j++) {
			Int_t histNum = (i*6) + j;
			//histNum += 4;
			if(histNum >= numWBins)
				break;
			if((histNum <= 16)||(histNum >= 38))
				break;
			Char_t cName[5];
			sprintf(cName,"c%d",histNum);
			TCanvas *cc = new TCanvas(cName, "Double-Spin Asymmetry", 10*histNum, 10*histNum, 1200, 800);
			char blankName[10];
			sprintf(blankName,"%.3f < Q^{2} < %.3f GeV",edges[histNum],edges[histNum+1]);
			TH2D *blank = new TH2D(blankName,blankName,300,0.9,1.6,200,-0.2,0.2); //is this still necessary?
			Float_t y0[numWBins], y0Error[numWBins],
					y1[numWBins], y1Error[numWBins],
					y2[numWBins], y2Error[numWBins],
					y3[numWBins], y3Error[numWBins],
					y4[numWBins], y4Error[numWBins],
					y5[numWBins], y5Error[numWBins],
					y6[numWBins], y6Error[numWBins],
					y7[numWBins], y7Error[numWBins],
					y8[numWBins], y8Error[numWBins], yReal[300];

			getPoints(y0,numerator,denominator,fdf,pbpt[0],0,histNum,raw,sepBySet);
			getErrors(y0Error,denominator,fdf,pbpt[0],0,histNum,raw,sepBySet);
			getPoints(y1,numerator,denominator,fdf,pbpt[1],1,histNum,raw,sepBySet);
			getErrors(y1Error,denominator,fdf,pbpt[1],1,histNum,raw,sepBySet);
			getPoints(y2,numerator,denominator,fdf,pbpt[2],2,histNum,raw,sepBySet);
			getErrors(y2Error,denominator,fdf,pbpt[2],2,histNum,raw,sepBySet);
			getPoints(y3,numerator,denominator,fdf,pbpt[3],3,histNum,raw,sepBySet);
			getErrors(y3Error,denominator,fdf,pbpt[3],3,histNum,raw,sepBySet);
			getPoints(y4,numerator,denominator,fdf,pbpt[4],4,histNum,raw,sepBySet);
			getErrors(y4Error,denominator,fdf,pbpt[4],4,histNum,raw,sepBySet);
			getPoints(y5,numerator,denominator,fdf,pbpt[5],5,histNum,raw,sepBySet);
			getErrors(y5Error,denominator,fdf,pbpt[5],5,histNum,raw,sepBySet);
			getPoints(y6,numerator,denominator,fdf,pbpt[6],6,histNum,raw,sepBySet);
			getErrors(y6Error,denominator,fdf,pbpt[6],6,histNum,raw,sepBySet);
			getPoints(y7,numerator,denominator,fdf,pbpt[7],7,histNum,raw,sepBySet);
			getErrors(y7Error,denominator,fdf,pbpt[7],7,histNum,raw,sepBySet);
			getPoints(y8,numerator,denominator,fdf,pbpt[8],8,histNum,raw,sepBySet);
			getErrors(y8Error,denominator,fdf,pbpt[8],8,histNum,raw,sepBySet);

			Float_t yy0[numWBins], yy0e[numWBins],
					yy1[numWBins], yy1e[numWBins];
			if(sepBySet) //statistically combining sets
			{
				for(int z = 0; z < numWBins; z++)
				{
					float n0 = pow(y0Error[z],-2); //for error bar purposes
					float n1 = pow(y1Error[z],-2);
					float n2 = pow(y2Error[z],-2);
					n2 = 0.0;
					float n3 = pow(y3Error[z],-2);
					float n4 = pow(y4Error[z],-2);
					float n5 = pow(y5Error[z],-2);
					float n6 = pow(y6Error[z],-2);
					float n7 = pow(y7Error[z],-2);
					n7 = 0.0;
					//now do the ol' numerator/denominator method

					//but first some raw numbers
					float n0r = pow(y0Error[z]*pbpt[0]*fdf[0][histNum][z],-2); //raw numbers
					float n1r = pow(y1Error[z]*pbpt[1]*fdf[1][histNum][z],-2);
					float n2r = pow(y2Error[z]*pbpt[2]*fdf[2][histNum][z],-2);
					n2r = 0.0;
					float n3r = pow(y3Error[z]*pbpt[3]*fdf[3][histNum][z],-2);
					float n4r = pow(y4Error[z]*pbpt[4]*fdf[4][histNum][z],-2);
					float n5r = pow(y5Error[z]*pbpt[5]*fdf[5][histNum][z],-2);
					float n6r = pow(y6Error[z]*pbpt[6]*fdf[6][histNum][z],-2);
					float n7r = pow(y7Error[z]*pbpt[7]*fdf[7][histNum][z],-2);
					n7r = 0.0;

					yy0[z] = ((y0[z]*n0r)+(y1[z]*n1r)+(y2[z]*n2r)+(y3[z]*n3r))/(n0r+n1r+n2r+n3r);
					yy0e[z] = pow((n0+n1+n2+n3),-0.5);
					if(yy0[z] != yy0[z])
					{
						yy0[z] = 0.0;
						yy0e[z] = 0.0;
					}
					yy1[z] = ((y4[z]*n4r)+(y5[z]*n5r)+(y6[z]*n6r)+(y7[z]*n7r))/(n4r+n5r+n6r+n7r);
					yy1e[z] = pow((n4+n5+n6+n7),-0.5);
					if(yy1[z] != yy1[z])
					{
						yy1[z] = 0.0;
						yy1e[z] = 0.0;
					}
					//note that the short targets have been cut
				}
				//now that the new ones have been produced, plug them into the old
				//so I can still use the old graph objects
				//top level is unchanged so I just move it over
				for(int z = 0; z < numWBins; z++)
				{
					y0[z] = yy0[z];
					y1[z] = yy1[z];
					y2[z] = y8[z];
					y0Error[z] = yy0e[z];
					y1Error[z] = yy1e[z];
					y2Error[z] = y8Error[z];
				}
			}

			//I apologize to anyone that has to read this
			//It's the night before my final corrections are due and I have to redo these graphs
			//And I'm about to lose my damn mind
			Float_t e1f[numWBins], e2f[numWBins], e3f[numWBins],
					e1p[numWBins], e2p[numWBins], e3p[numWBins], nuy[numWBins];
			for(int t = 0; t < numWBins; t++)
			{
				if(ferr1[histNum][t] < 0.2)
					e1f[t] = ferr1[histNum][t];
				else
					e1f[t] = 0.0;
				if(ferr2[histNum][t] < 0.2)
					e2f[t] = ferr2[histNum][t];
				else
					e2f[t] = 0.0;
				if(ferr3[histNum][t] < 0.2)
					e3f[t] = ferr3[histNum][t];
				else
					e3f[t] = 0.0;
				if(perr1[histNum][t] < 0.2)
					e1p[t] = perr1[histNum][t];
				else
					e1p[t] = 0.0;
				if(perr2[histNum][t] < 0.2)
					e2p[t] = perr2[histNum][t];
				else
					e2p[t] = 0.0;
				if(perr3[histNum][t] < 0.2)
					e3p[t] = perr3[histNum][t];
				else
					e3p[t] = 0.0;
				nuy[t] = -0.2;
			}
			
			for(int bin=0; bin<300; bin++){
				yReal[bin] = real[histNum][bin];
			}

			TGraphErrors *gasym0 = new TGraphErrors(numWBins,x,y0,zero,y0Error);
			TGraphErrors *gasym1 = new TGraphErrors(numWBins,x,y1,zero,y1Error);
			TGraphErrors *gasym2 = new TGraphErrors(numWBins,x,y2,zero,y2Error);

			TGraph *gef1 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e1f);
			TGraph *gef2 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e2f);
			TGraph *gef3 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e3f);
			TGraph *gep1 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e1p);
			TGraph *gep2 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e2p);
			TGraph *gep3 = new TGraphAsymmErrors(numWBins, x, nuy, zero, zero, zero, e3p);

			gef1->SetFillColorAlpha(kRed, 0.5);
			gep1->SetFillColorAlpha(kBlue, 0.5);
			gef2->SetFillColorAlpha(kRed, 0.5);
			gep2->SetFillColorAlpha(kBlue, 0.5);
			gef3->SetFillColorAlpha(kRed, 0.5);
			gep3->SetFillColorAlpha(kBlue, 0.5);

			leg = new TLegend(0.7, 0.77, 0.89, 0.89);

			/*TGraph *gasym0 = new TGraph(numWBins,x,y0);
			TGraph *gasym1 = new TGraph(numWBins,x,y1);
			TGraph *gasym2 = new TGraph(numWBins,x,y2);*/
			
			blank->GetXaxis()->SetTitle("W - Lorentz-invariant missing mass (GeV)");
			blank->GetYaxis()->SetTitle("Asymmetry");
			blank->GetXaxis()->CenterTitle();
			blank->GetYaxis()->CenterTitle();
			//gasym1->SetLineColor(kBlue);
			//gasym2->SetLineColor(kRed);

			gasym0->SetMarkerStyle(7);
			gasym1->SetMarkerStyle(7);
			gasym2->SetMarkerStyle(7);
			
			//c->cd(j+1); //no longer valid!
			//gPad->SetFillColor(42);
			blank->Draw();
			gasym0->Draw("P");
			//gasym1->Draw("P");
			//gasym2->Draw("P");
			gef1->Draw("3");
			gep1->Draw("3");

			TGraph *greal = new TGraph(300,realx,yReal);
			greal->SetLineColor(kGreen);
			//greal->Draw("L");
			//leg->AddEntry(greal,"JLab Asymmetry fit","l");
			leg->AddEntry(gasym0,"1.1 and 1.3 GeV","l");
			//leg->AddEntry(gasym1,"2 and 2.2 GeV","l");
			//leg->AddEntry(gasym2,"3 GeV","l");	
			leg->AddEntry(gef1,"FDF Error", "f");
			leg->AddEntry(gep1,"PbPt Error", "f");
			//if(j==1)
			//	leg->Draw(); //keep the legend to one window
			leg->Draw();
			gStyle->SetOptStat(0); //remove annoying stat box
			gStyle->SetEndErrorSize(5);
			
			for(int f = 0; f < numWBins; f++)
			{
				finalasym[0][histNum][f] = y0[f];
				finalasym[1][histNum][f] = y1[f];
				finalasym[2][histNum][f] = y2[f];
				finalerror[0][histNum][f] = y0Error[f];
				finalerror[1][histNum][f] = y1Error[f];
				finalerror[2][histNum][f] = y2Error[f];
			}
		}
	}

	

	/*Float_t f1[numWBins], f2[numWBins], f3[numWBins],
			f1e[numWBins], f2e[numWBins], f3e[numWBins], count[6];
	Float_t realf[300], rcount;
	for(int t = 0; t < numWBins; t++)
	{
		f1[t] = 0;
		f2[t] = 0;
		f3[t] = 0;
		f1e[t] = 0;
		f2e[t] = 0;
		f3e[t] = 0;
	}
	for(int u = 0; u < 6; u++)
		count[u] = 0;
	for(int y = 0; y < 300; y++)
		realf[y] = 0;
	for(int qBin = 0; qBin < numQBins; qBin++)
	{
		for(int wBin = 0; wBin < numWBins; wBin++)
		{
			if(finalasym[0][qBin][wBin] != 0.0)
			{
				f1[wBin] += finalasym[0][qBin][wBin];
				count[0] += 1;
			}
			if(finalasym[1][qBin][wBin] != 0.0)
			{
				f2[wBin] += finalasym[1][qBin][wBin];
				count[1] += 1;
			}
			if(finalasym[2][qBin][wBin] != 0.0)
			{
				f3[wBin] += finalasym[2][qBin][wBin];
				count[2] += 1;
			}
			if(finalerror[0][qBin][wBin] != 0.0)
			{
				f1e[wBin] += finalerror[0][qBin][wBin];
				count[3] += 1;
			}
			if(finalerror[1][qBin][wBin] != 0.0)
			{
				f2e[wBin] += finalerror[1][qBin][wBin];
				count[4] += 1;
			}
			if(finalerror[2][qBin][wBin] != 0.0)
			{
				f3e[wBin] += finalerror[2][qBin][wBin];
				count[5] += 1;
			}
		}
		for(int w; w < 300; w++)
		{
			if(real[qBin][w] != 0.0)
			{
				realf[w] += real[qBin][w];
				rcount += 1;
			}
		}
	}
	for(int i = 0; i < numWBins; i++)
	{
		f1[wBin] /= count[0];
		f2[wBin] /= count[1];
		f3[wBin] /= count[2];
		f1e[wBin] /= count[3];
		f2e[wBin] /= count[4];
		f3e[wBin] /= count[5];
	}

	for(int i = 0; i < 300; i++)
	{
		realf[i] /= rcount;
	}

	TCanvas *cccc = new TCanvas("cccc","Combined Asymmetry",40,40,800,600);
	TH2D *blank = new TH2D("Combined graph", "Combined graph", 300, 0, 3, 500, -5.0, 5.0);
	TGraphErrors *g1 = new TGraphErrors(numWBins,x,f1,zero,f1e);
	TGraphErrors *g2 = new TGraphErrors(numWBins,x,f2,zero,f2e);
	TGraphErrors *g3 = new TGraphErrors(numWBins,x,f3,zero,f3e);
	TGraph *gr = new TGraph(numWBins, realx, realf);
	blank->Draw();
	g1->Draw("P");
	g2->Draw("P");
	g3->Draw("P");
	gr->Draw("L");
	gStyle->SetOptStat(0); //remove annoying stat box
	leg2 = new TLegend(0.55, 0.65, 0.89, 0.89);
	leg2->AddEntry(gr,"1.1 and 1.3 GeV","l");
	leg2->AddEntry(gr,"2.0 and 2.2 GeV","l");
	leg2->AddEntry(gr,"3.0 GeV","l");
	leg2->AddEntry(gr,"JLab Asymmetry fit","l");
	leg2->Draw();*/


	//gonna finish off with a cheeky little save here. Cause screw asking to do so each time
	//output(finalasym);
	cout << "Done!" << endl;
}

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

void getPoints(Float_t points[75], Double_t numerator[9][54][75], Double_t denominator[9][54][75], Double_t fdf[9][54][75], Double_t pbpt, int set, int histNum, bool raw, bool sepBySet)
{
	int e = energyFromSet(set, sepBySet);
	Double_t pbptErr[9] = {0.080, 0.034, 0.068, 0.030, 0.034, 0.040, 0.072, 0.036, 0.048};
	//only use above when looking at pbpt systerr
	if(!sepBySet)
	{
		pbptErr[0] = 0.06875;
		pbptErr[1] = 0.0685;
		pbptErr[2] = 0.066;
	}
	for(int bin = 0; bin < 75; bin++)
	{
		if(denominator[set][histNum][bin] == 0)
		{
			points[bin] = 0;
		}
		else{//if everything is fine, statistically combine
			points[bin] = numerator[set][histNum][bin]/denominator[set][histNum][bin];
			points[bin] /= (pbpt/* + pbptErr[set]/**/);//only use when looking at pbpt error
			if(fdf[e][histNum][bin] != 0.0)
			{
				if(!raw)
					points[bin] /= fdf[e][histNum][bin];
			}
			else
			{
				points[bin] = 0.0; //no results if fdf is zero
			}
		}
	}
}

void getErrors(Float_t errors[75], Double_t denominator[9][54][75], Double_t fdf[9][54][75], Double_t pbpt, int set, int histNum, bool raw, bool sepBySet)
{
	int e = energyFromSet(set, sepBySet);
	for(int bin = 0; bin < 75; bin++)
	{
		if(denominator[set][histNum][bin] == 0)
		{
			errors[bin] = 0;
		}
		else{//if everything is fine, statistically combine
			if(fdf[e][histNum][bin] == 0.0)
			{
				errors[bin] = 0.0;
				continue;
			}
			errors[bin] = pow(denominator[set][histNum][bin],-0.5);
			errors[bin] /= pbpt;
			errors[bin] /= fdf[e][histNum][bin];
		}
	}
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

/*
 * For use in comparing to adjusted asyms for systematic error analysis
 */
void output(Float_t final[3][53][75])
{
	ofstream outfile;
	int max;
	outfile.open("../common/asym_final.dat");
	max = 3;
	for(int e = 0; e < max; e++)
	{
		for(int i = 0; i < 53; i++)
		{
			for(int j = 0; j < 75; j++)
			{
				outfile << final[e][i][j] << " ";
			}
			outfile << endl;
		}
	}
}

/*/I'm getting sick of this not working so I'm doing it manually. F the police.
	Float_t edges[54] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223,
			     0.00266, 0.00317, 0.00379, 0.00452, 0.00540, 0.00645,
			     0.00770, 0.00919, 0.0110, 0.0131, 0.0156, 0.0187,
			     0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540,
			     0.0645, 0.0770, 0.0919, 0.110, 0.131, 0.156,
			     0.187, 0.223, 0.226, 0.317, 0.379, 0.452,
			     0.540, 0.645, 0.770, 0.919, 1.10, 1.31,
			     1.56, 1.87, 2.23, 2.66, 3.17, 3.79,
			     4.52, 5.40, 6.45, 7.70, 9.19, 10.97};//*/
