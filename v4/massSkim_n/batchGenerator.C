#include<stdlib.h>
#include<stdio.h>

#define makeSkims_cxx

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

int batchGenerator()
{
	ifstream infile;
	char filename[100], target[10], energy[10], run[10];
	cout << "Just enter the filename, I'm so sick of this" << endl;
	scanf("%s",&filename);
	infile.open(filename);
	string line;
	while(infile.peek() != EOF)
	{
		infile >> target;
		infile >> energy;
		infile >> run;
		system("touch tempfile");
		ofstream outfile;
		outfile.open("tempfile");
		outfile << "PROJECT: eg1b" << endl;
		outfile << "TRACK: analysis" << endl;
		//outfile << "OS: centos7" << endl; //no longer seems necessary
		outfile << "JOBNAME: n_fullcutcount_" << target << "_" << energy << "_" << run << endl;
		outfile << "INPUT_FILES: /home/ryanpl/Thesis/v4/AnalysisScript" << endl;
		outfile << "COMMAND: AnalysisScript " << target << " " << run << " " << energy << endl;
		outfile << "MEMORY: 9500 MB" << endl;
		outfile.close();
		system("jsub tempfile");
		system("rm tempfile");
	}
	return 0;
}
