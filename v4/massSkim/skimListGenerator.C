int skimListGenerator()
{
	ifstream f;
	char filename[100], target[10], energy[10], newfilename[100], command[100];
	cout << "Target?" << endl;
	scanf("%s", &target);
	cout << "Energy?" << endl;
	scanf("%s", &energy);
	sprintf(filename,"../../common/fileList_%s_%s.dat",target,energy);
	cout << filename << endl;
	f.open(filename);
	string line, currentRun;
	sprintf(newfilename,"batchList_%s_%s.dat",target,energy);
	sprintf(command,"touch %s",newfilename);
	system(command);
	ofstream outfile;
	outfile.open(newfilename);
	while(getline(f,line))
	{
		string thisRun = line.substr(5,5);
		if(thisRun.compare(currentRun)==0)
			continue;
		else
			currentRun = thisRun;
		cout << "<TARGET> = " << target << " | <ENERGY> = " << energy << " | <RUN> = " << currentRun  << endl;
		outfile << target << " " << energy << " " << currentRun << endl;
	}
	f.close();
}
