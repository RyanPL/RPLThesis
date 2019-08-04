#ifndef EC_CUTS_DEV_PASS2_H
#define EC_CUTS_DEV_PASS2_H
/*******
 *  Since there were bins which had not enough or no data/statistics to do
 *     good Gaussian fittings, so I have picked up the last good value for
 *     all the higher bins (or lowest good value for all the lower bins, 
 *     if any)
 *  For example, the 0th Q2-bin is ignored but it still has some data in it
 *    and I copied the lowest value of both Mean and Sigma and assigned that
 *    to the 0th bin to make the code consistent and have no issue during
 *    the data analysis run. 
 *
 * Update: //0.06 changed to 0.05 on 1/14/16
***/

//SF = etot/p = sampling fraction
//Numbers collected from files such as https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/EC/Plots/ecSF_InQ2binsWithCts06Pass2Ebi0GausFits.txt corresponding to plots such as https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/EC/Plots/ecSF_InQ2binsWithCts06Pass2Ebi0.gif
//Unfortunately, due to a stupid mistake today (1/2/16), I deleted the corresponding source file. Only a 4 day older (not-up-to-date) file exists in /home/adhikari/AnaEG4Pass2/Cuts/EC/ (Some idea about what updates that was lost can be obtained from emails that I wrote to Dr. Kuhn with subject 'EC-cuts (Gaussian fits on etot/p)' and dated 12/16/15 and 12/17/15.)
double MeanSF[5][41] = {  //[5] for Eb indices, 40 for Q2-bin indices (1 to 40), 0th bin added by hand
  {0.24276,   0.24276,0.248848,0.249858,0.246993,0.242765,0.241928,0.247511,0.249595,0.249591,0.249571,0.248884,0.249028,0.249379,0.246809,0.245065,0.246434,0.249373,0.250194,0.251312,0.252599,   0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,   0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599,0.252599}, //1.0 GeV
  {0.232695,  0.232695,0.236073,0.240385,0.24383,0.248251,0.251552,0.251159,0.250325,0.249258,0.24689,0.24316,0.245959,0.24904,0.249968,0.249104,0.249759,0.248349,0.24546,0.246441,0.24949,0.250711,0.252504,0.252466,   0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,   0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,0.252466,0.252466}, //1.3 GeV
  {0.22243,   0.22243,0.224964,0.227952,0.230887,0.234114,0.236967,0.240539,0.244537,0.247268,0.250544,0.249655,0.248938,0.247588,0.246531,0.244812,0.246701,0.249161,0.249893,0.249758,0.249477,0.2491,0.246875,0.247365,0.250075,0.251202,0.252496,    0.252496,0.252496,0.252496,0.252496,    0.252496,0.252496,0.252496,0.252496,0.252496,0.252496,0.252496,0.252496,0.252496,0.252496}, //2.0 GeV
  {0.2301,   0.2301,0.231506,0.23335,0.235287,0.238292,0.241129,0.243258,0.246061,0.248772,0.251087,0.253089,0.256559,0.256424,0.255576,0.254968,0.253833,0.252023,0.252113,0.254608,0.25578,0.255908,0.25579,0.254253,0.252898,0.25591,0.257061,0.259403,   0.259403,0.259403,0.259403,   0.259403,0.259403,0.259403,0.259403,0.259403,0.259403,0.259403,0.259403,0.259403,0.259403}, //2.3 GeV
  {0.231399,   0.231399,0.232021,0.232852,0.234764,0.237536,0.239851,0.242787,0.245055,0.247086,0.249647,0.251708,0.253433,0.255389,0.257532,0.259772,0.262021,0.264115,0.264177,0.263614,0.26279,0.260057,0.258497,0.260809,0.262144,0.262119,0.261286,0.25973,0.261159,0.262895,0.264891    ,0.264891,0.264891,0.264891,0.264891,0.264891,0.264891,0.264891,0.264891,0.264891,0.264891} //3.0 GeV
};

double SigmaSF[5][41] = {  //[5] for Eb indices, 40 for Q2-bin indices (1 to 40), 0th bin added by hand
  {0.0294313,   0.0294313,0.0245648,0.0240581,0.0238775,0.0234791,0.0236234,0.024215,0.0242233,0.0248571,0.0252079,0.0256451,0.0263732,0.0267299,0.0263702,0.0260845,0.0262313,0.0259479,0.025613,0.0256279,0.025947,   0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,   0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947,0.025947},  //1.0 GeV
  {0.0327678,   0.0327678,0.0307928,0.0286355,0.0265059,0.0246985,0.0219915,0.0218816,0.0219649,0.0219435,0.022242,0.0220096,0.0225479,0.0227154,0.0229373,0.0235889,0.0242705,0.0242589,0.0234786,0.0235435,0.0234042,0.022955,0.0226406,0.023475,      0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,      0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,0.023475,0.023475}, //1.3 GeV
  {0.0326168,  0.0326168,0.031814,0.0302171,0.029416,0.0286307,0.0270238,0.0249259,0.0225275,0.0211732,0.0192195,0.0187138,0.0188568,0.0187289,0.0188228,0.0188836,0.019296,0.0194162,0.0196555,0.0199324,0.0207492,0.0206905,0.0202058,0.0202842,0.0202476,0.0198083,0.0194148,    0.0194148,0.0194148,0.0194148,0.0194148,   0.0194148,0.0194148,0.0194148,0.0194148,0.0194148,0.0194148,0.0194148,0.0194148,0.0194148,0.0194148}, //2.0 GeV
  {0.0323514,  0.0323514,0.0324079,0.0312327,0.0299641,0.0294516,0.0285603,0.0266752,0.0249825,0.0231347,0.0218469,0.0206324,0.019149,0.0183763,0.01837,0.0185419,0.0184001,0.018643,0.0187495,0.0191608,0.0195656,0.0200402,0.020652,0.0205402,0.0200238,0.0203742,0.01979,0.0198384,   0.0198384,0.0198384,0.0198384,  0.0198384,0.0198384,0.0198384,0.0198384,0.0198384,0.0198384,0.0198384,0.0198384,0.0198384,0.0198384}, //2.3 GeV
  {0.0356647,  0.0356647,0.0349906,0.0346109,0.0332861,0.0315748,0.030388,0.0297525,0.0286662,0.0273775,0.026132,0.0251068,0.024093,0.0230616,0.0216324,0.0205016,0.0192276,0.0180918,0.0173286,0.0178533,0.0182804,0.0177496,0.017818,0.0179902,0.0187981,0.0195226,0.0196121,0.018765,0.0191855,0.0183465,0.0202183   ,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183,0.0202183} //3.0 GeV
};

double Bins[53] = {0.000919, 0.00110, 0.00131, 0.00156, 0.00187, 0.00223, 0.00266, 0.00317,
		   0.00379, 0.00452, 0.00540, 0.00645, 0.00770, 0.00919, 0.0110, 0.0131,
		   0.0156, 0.187, 0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 0.0540, 0.0645,
		   0.0770, 0.0919, 0.110, 0.131, 0.156, 0.187, 0.223, 0.266, 0.317, 0.379,
		   0.452, 0.540, 0.645, 0.770, 0.919, 1.10, 1.31, 1.56, 1.87, 2.23, 2.66,
		   3.17, 3.79, 4.52, 5.40, 6.45, 7.70, 9.19};
bool ecSFcut(int Ebi, int Q2bin, double SF);
int GetQ2Index(double Q2, double Q2min, double Q2max);

bool ecSFcut(int Ebi, int Q2bin, double SF)
{
  double SFcut = MeanSF[Ebi][Q2bin] - 3*SigmaSF[Ebi][Q2bin];
  if(SF>SFcut) return true; else return false;
}

/**/int GetQ2Index(double Q2, double Q2min, double Q2max){
  for(int i = 0; i<53; i++){
    if(Q2 > Bins[i]){
      if(i == 52){
        Q2max = 10.97;
      }
      else
        Q2max = Bins[i+1];
      Q2min = Bins[i];
      return i;
    }
  }
  return 0;
}

//p0 and cz0 are p[0] and cz[0] i.e the uncorrected p and z-direction cosine
//Chose p0 and cz0 as inputs to ensure, there will be no confusion about whether to use
//kinematic corrected p and theta or uncorrected ones.
bool ecCuts4ElPass2(int Ebi,double p0,double cz0, double ecIn, double ecTot)
{
  /**/double Ebg[5] = {1.0539, 1.339, 1.9889, 2.256, 2.999};
  double Eb=Ebg[Ebi]-0.002; //incoming energy-loss corrected beam energy
  double theta=acos(cz0);
  double Q2min=0.0,Q2max=0.0, Q2=4*Eb*p0*pow(sin(theta/2),2.0);  
  int Q2bin=-1; bool sfCt=false, eciCt=false;
  Q2bin = GetQ2Index(Q2, Q2min, Q2max); //Input is Q2 and outputs are Q2index, Q2min and Q2max 
  if(Q2bin>-1 && Q2bin<41)
    {
      //if(ecIn>0.05) eciCt=true;//0.06 changed to 0.05 on 1/14/16
      if(ecIn>0.06) eciCt=true;
      sfCt = ecSFcut(Ebi,Q2bin,ecTot/p0);
    }
  //if(ecIn>0.05) eciCt=true; else eciCt=false;//0.06 changed to 0.05 on 3/20/16
  if(eciCt==true && sfCt==true) return true;
  else return false;
}

#endif