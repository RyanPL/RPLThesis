#ifndef MORE_FID_CUTS  //6/30/16
#define MORE_FID_CUTS

//================ Cuts on PhDc1Deg vs ThVtxDeg dists using the new inv-P bins with less # of bins =======
//==== Such as https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cors/CCinef/New/eventsUsed4AvgNpheReg_Eb2.gif
double thVtxMaxCuts[8] = {12.0,30.0, 28.5, 27.0,30.0, 25.0,25.0,25.0};

bool moreFidCutsWithLessInvPBins(int Ebi, double pp, double phDc1DegPlusMinus30, double thVtxDeg)//phDc1Deg must be between -30.0, 30.0
{
  double phDc1Deg = phDc1DegPlusMinus30;
  double invPbinEdge[]={0.3, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5, 2.9, 3.4};
  //Calculate invP & invPBin based on Ebi & p. 
  int ipBin=-1; double invP = ItorP[Ebi]/(pp*2250.0); 
  for(int i=0;i<8;i++) {if(invP > invPbinEdge[i] && invP <= invPbinEdge[i+1]) ipBin=i; if(ipBin>-1) break;}


  bool thMaxCt = true; if(thVtxDeg > thVtxMaxCuts[ipBin]) thMaxCt = false;
  bool patch1Ct = true; if(ipBin == 1 && phDc1Deg < 0.0 && thVtxDeg > 23.5) patch1Ct = false;
  //Upper dark patch on ip=4th invP-bin: area above line joining points (23,5) and (30,-2) which 
  //    gives the equation (using 2-point form):    y = ((-2.0 - 5.0)/(30.0 - 23.0)) * (x - 23.0) - 5.0
  bool patch2Ct = true; 
  if(ipBin == 4 && thVtxDeg > 23.0 
     &&  (phDc1Deg > ( ((-2.0 - 5.0)/(30.0 - 23.0)) * (thVtxDeg - 23.0) + 5.0) )) patch2Ct = false;

  //Upper dark patch on ip=5th invP-bin: area above line joining points (18,5) and (25,-1) which 
  //    gives the equation (using 2-point form):    y = ((-1.0 - 5.0)/(25.0 - 18.0)) * (x - 18.0) - 1.0
  bool patch3Ct = true; 
  if(ipBin == 5 && thVtxDeg > 18.0 
     &&  (phDc1Deg >  ( ((-1.0 - 5.0)/(25.0 - 18.0)) * (thVtxDeg - 18.0) + 5.0) )) patch3Ct = false;
  if(thMaxCt == false || patch1Ct == false || patch2Ct == false || patch3Ct == false) return false;
  else return true;
}

//
//  7/2/16:
//  Since moreFidCutsWithLessInvPBins(..) above seemed to make the ratios better, I am remaking a second version of it
//    refining the cuts or adding new ones (particularly at lower angles as well)
//
//double thVtxMaxCuts2[8] = {12.0, 30.0, 30.0, 31.5, 30.0, 25.0, 25.0, 25.0};
double thVtxMaxCuts2[8] = {12.0, 30.0, 30.0, 31.5, 29.5, 25.0, 25.0, 25.0};//7/5/16
//See: https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cors/CCinef/New/eventsUsed4AvgNpheReg_Eb2.gif
//   & https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/MoreCts/moreFiducialCuts.gif 
bool moreFidCutsWithLessInvPBins2(int Ebi, double pp, double phDc1DegPlusMinus30, double thVtxDeg)//phDc1Deg must be between -30.0, 30.0
{
  double phDc1Deg = phDc1DegPlusMinus30;
  double invPbinEdge[]={0.3, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5, 2.9, 3.4};
  //Calculate invP & invPBin based on Ebi & p. 
  int ipBin=-1; double invP = ItorP[Ebi]/(pp*2250.0); 
  for(int i=0;i<8;i++) {if(invP > invPbinEdge[i] && invP <= invPbinEdge[i+1]) ipBin=i; if(ipBin>-1) break;}


  bool thMaxCt = true; if(thVtxDeg > thVtxMaxCuts2[ipBin]) thMaxCt = false;
  //bool patch1Ct = true; if(ipBin == 1 && phDc1Deg < -1.5 && thVtxDeg > 30.0) patch1Ct = false;
  bool patch1Ct  = true; if(ipBin == 1 && phDc1Deg > 12.5 && thVtxDeg > 12.0) patch1Ct = false;
  bool patch21Ct = true; if(ipBin == 2 && phDc1Deg > 3.0 && thVtxDeg > 28.5) patch21Ct = false;
  bool patch22Ct = true; if(ipBin == 2 && phDc1Deg < -6.0 && thVtxDeg > 28.5) patch22Ct = false;
  bool patch23Ct = true; if(ipBin == 2 && thVtxDeg < 7.0) patch23Ct = false;
  bool patch31Ct = true; if(ipBin == 3 && phDc1Deg > 6.0 && thVtxDeg > 27.0) patch31Ct = false;
  bool patch32Ct = true; if(ipBin == 3 && phDc1Deg > 1.5 && thVtxDeg > 29.0) patch32Ct = false;
  bool patch41Ct = true; if(ipBin == 4 && phDc1Deg > 4.0 && thVtxDeg > 23.5) patch41Ct = false;
  bool patch42Ct = true; if(ipBin == 4 && phDc1Deg > -1.5 && thVtxDeg > 27.3) patch41Ct = false;
 
  bool patch51Ct = true, patch52Ct = true, patch53Ct = true; 
  //Upper left dark patch on ip=5th invP-bin: area above line joining points (8,13) and (10,4.0) which 
  //    gives the equation (using 2-point form):    y = ((4.0 - 13.0)/(10.0 - 8.0)) * (x - 8.0) - 4.0
  if(ipBin == 5 && thVtxDeg < 10.0 
     &&  (phDc1Deg >  ( ((4.0 - 13.0)/(10.0 - 8.0)) * (thVtxDeg - 8.0) + 13.0) )) patch51Ct = false;
  //Upper dark patch on ip=5th invP-bin: area above line joining points (18.5,4) and (25,-0.5) which 
  //    gives the equation (using 2-point form):    y = ((-0.5 - 4.0)/(25.0 - 18.5)) * (x - 8.0) - 0.5
  if(ipBin == 5 && thVtxDeg > 18.5 
     &&  (phDc1Deg >  ( ((-0.5 - 4.0)/(25.0 - 18.5)) * (thVtxDeg - 18.5) + 4.0) )) patch52Ct = false;

  //Lower dark patch on ip=5th invP-bin: area below line joining points (16.0,-10.0),(21.5,-5.0)
  //    and to the left of line joining (21.5,-5.0) & (24.0,-13.0)
  if(ipBin == 5 && thVtxDeg > 16.0 //&& thVtxDeg < 22.5 
     &&  (phDc1Deg <  ( ((-5.0 + 10.0)/(21.5 - 16.0)) * (thVtxDeg - 16.0) - 10.0) )
     &&  (phDc1Deg <  ( ((-13.0 + 5.0)/(24.0 - 21.5)) * (thVtxDeg - 21.5) - 5.0) )) patch53Ct = false;


  //if(thMaxCt == false || patch1Ct == false || patch2Ct == false || patch31Ct == false) return false;
  if(thMaxCt == true && patch1Ct == true && patch21Ct == true &&  patch22Ct == true &&  patch23Ct == true 
     &&  patch31Ct == true &&  patch32Ct == true &&  patch41Ct == true &&  patch42Ct == true &&  patch51Ct == true 
     &&  patch52Ct == true &&  patch53Ct == true) return true;
  else return false;
}

//7/4/16: For the triangular patch at low theta end on ip=2 bin as seen in the following plot:
//       https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cors/CCinef/New/eventsUsed4AvgNpheReg_Eb2.gif  
bool moreFidCutsWithLessInvPBins2extra(int Ebi, double pp, double phDc1DegPlusMinus30, double thVtxDeg)//phDc1Deg must be between -30.0, 30.0
{
  double phDc1Deg = phDc1DegPlusMinus30;
  double invPbinEdge[]={0.3, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5, 2.9, 3.4};
  //Calculate invP & invPBin based on Ebi & p. 
  int ipBin=-1; double invP = ItorP[Ebi]/(pp*2250.0); 
  for(int i=0;i<8;i++) {if(invP > invPbinEdge[i] && invP <= invPbinEdge[i+1]) ipBin=i; if(ipBin>-1) break;}

  //defining the patch by two st. lines joining the following pairs of points: 
  // Line1:   joining (7.0,8.5) & (13.0,1.0)
  // Line2:   joining (7.0,-2.5) & (13.0,1.0)  //(8.0,1.0) (7.0,8.0)
  bool patch2Ct = true;
  if(ipBin == 2   &&  (phDc1Deg <  ( ((1.0 - 8.5)/(13.0 - 7.0)) * (thVtxDeg - 7.0) + 8.5) )
     &&  (phDc1Deg >  ( ((1.0 + 2.5)/(13.0 - 7.0)) * (thVtxDeg - 7.0) - 2.5) )   ) patch2Ct = false;
  return patch2Ct;
}



//
// 7/2-10/16: I wanted to make similar set of cuts based on dark-spots (low cc-efficiency regions) in the following plots
//      https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cors/CCinef/Newr/eventsUsed4AvgNpheReg_Eb2.gif
//      But, for the moment, I decided to not use it. I may come back to it later.
//

double thVtxMaxCutsOldBins[13] = {11.0, 12.0, 30.0, 35.0, 34.0, 34.0, 33.0, 32.0, 30.0, 29.0, 25.0, 24.0, 25.0};
//End points of the Inclinded straight lines defining some cuts:
const int stLnNum=34; //# of straight lines defining the cuts (will be used to make and draw corresponding TLines as well)
double ctPt[stLnNum][4] = { //[4] is for cut points (x1,y1), (x2, y2) (points defining the cut-line)
  {12.0,12.5, 12.0,14.0},  {12.0,12.5,  35.0,12.5},//invPBin=2: Two lines for the patch cut on the top region 
  {29.0,0.0,  29.0,-14.0}, {29.0,0.0,   35.0,0.0}, //invPBin=2: Two lines for the patch cut on bottom-right region 
  {12.0,12.5, 12.0,14.0},  {12.0,12.5,  35.0,12.5},//invPBin=3: Two lines for the patch cut on the top region 
  {30.5,1.0,  30.5,14.0},  {30.5,1.0,   35.0,5.0}, //invPBin=3: Two lines for the patch cut on top right region
  {30.0,-2.0, 30.0,-14.0}, {30.0,-2.0,  35.0,-7.0},//invPBin=3: Two lines for the patch cut on bottom-right region 
  {30.0,3.0,  30.0,14.0},  {30.0,3.0,   36.0,3.0}, //invPBin=4: Two lines for the patch cut on top right region
  {30.0,-5.0, 30.0,-14.0}, {30.0,-5.0,  35.0,-8.0},//invPBin=4: Two lines for the patch cut on bottom-right region
  {29.5,3.5,  29.5,14.0},  {29.5,3.5,   35.0,3.5}, //invPBin=5: Two lines for the patch cut on top right region
  {29.5,-5.5, 29.5,-14.0}, {29.5,-5.5,  35.0,-7.5},//invPBin=5: Two lines for the patch cut on bottom-right region
  {28.5,1.0,  28.5,14.0},  {28.5,4.0,   33.0,2.5}, //invPBin=6: Two lines for the patch cut on top right region
  {28.5,-7.0, 28.5,-14.0}, {28.5,-7.0,  33.0,-6.0},//invPBin=6: Two lines for the patch cut on bottom right region
  {27.5,7.0,  27.5,14.0},  {27.5,7.0,   32.0,2.0}, //invPBin=7: Line for the patch cut on top right region
  {25.0,13.0, 31.0,-1.0},                          //invPBin=8: Line for the patch cut on top right region
  {23.0,5.0,  23.0,14.0},  {23.0,5.0,   29.0,-2.0},//invPBin=9: Two lines for the patch cut on top right region
  {21.0,3.0,  21.0,14.0},  {21.0,3.0,   25.0,4.0}, //invPBin=10: Two lines for the patch cut on top right region
  {10.0,0.0,  10.0,14.0},  {10.0,-1.0,  15.0,8.5}, //invPBin=11: Two lines for the patch cut on top left region
  {15.0,8.5,  25.0,-1.0},                          //invPBin=11: Two lines for the patch cut on top right region
  {12.0,-14.0,21.0,-6.0},  {21.5,-6.0, 21.5,-14.0},//invPBin=11: Two lines for the patch cut on bottom right region
};


//https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/MoreCts/moreFiducialCutsMoreInversePBinsEbi2.gif //11/13/16
//
bool moreFidCutsWithMoreInvPBins(int Ebi, double pp, double phDc1DegPlusMinus30, double thVtxDeg)//phDc1Deg must be between -30.0, 30.0
{
  double ph = phDc1DegPlusMinus30, th = thVtxDeg;
  double invPbinEdge[]={0.3,0.4,0.5,0.6,0.7,0.85,1.0,1.25,1.5,1.75,2.0,2.25,2.53,3.4};
  //Calculate invP & invPBin based on Ebi & p. 
  int ipBin=-1; double invP = ItorP[Ebi]/(pp*2250.0); 
  for(int i=0;i<13;i++) {if(invP > invPbinEdge[i] && invP <= invPbinEdge[i+1]) ipBin=i; if(ipBin>-1) break;}

  bool thMaxCt = true; if(th > thVtxMaxCutsOldBins[ipBin]) thMaxCt = false;
  bool patchCt21=true, patchCt22=true, patchCt31=true, patchCt32=true, patchCt33=true, patchCt41=true, patchCt42=true;
  bool patchCt51=true, patchCt52=true, patchCt61=true, patchCt62=true, patchCt71=true, patchCt81=true, patchCt91=true;
  bool patchCtA1=true, patchCtB1=true, patchCtB2=true, patchCtB3=true; //Using hex-notation for 10 and 11.
  //Straight line equation (joing points (x1,y1) & (x2,y2)):  y - y1 = ((y2-y1)/(x2-x1)) * (x - x1)
  //       or           y  = (x - x1)*(y2-y1)/(x2-x1)    + y1
  if(ipBin == 2  &&  th > 12.0  &&  ph > 12.0) patchCt21 = false;
  if(ipBin == 2  &&  th > 29.0  &&  ph < 0.0)  patchCt22 = false;
  if(ipBin == 3  &&  th > 12.0  &&  ph > 12.0) patchCt31 = false;
  if(ipBin == 3  &&  th > 30.5  &&  ph > (th - ctPt[7][0]) * (ctPt[7][3] - ctPt[7][1])/(ctPt[7][2] - ctPt[7][0]) + ctPt[7][1]) patchCt32 = false;
  if(ipBin == 3  &&  th > 30.0  &&  ph < (th - ctPt[9][0]) * (ctPt[9][3] - ctPt[9][1])/(ctPt[9][2] - ctPt[9][0]) + ctPt[9][1]) patchCt33 = false;
  if(ipBin == 4  &&  th > 30.0  &&  ph > 3.0) patchCt41 = false;
  if(ipBin == 4  &&  th > 30.0  &&  ph < (th - ctPt[13][0]) * (ctPt[13][3] - ctPt[13][1])/(ctPt[13][2] - ctPt[13][0]) + ctPt[13][1]) patchCt42 = false;
  if(ipBin == 5  &&  th > 29.5  &&  ph > 3.5) patchCt51 = false;
  if(ipBin == 5  &&  th > 30.0  &&  ph < (th - ctPt[17][0]) * (ctPt[17][3] - ctPt[17][1])/(ctPt[17][2] - ctPt[17][0]) + ctPt[17][1]) patchCt52 = false;
  if(ipBin == 6  &&  th > 28.5  &&  ph > (th - ctPt[19][0]) * (ctPt[19][3] - ctPt[19][1])/(ctPt[19][2] - ctPt[19][0]) + ctPt[19][1]) patchCt61 = false;
  if(ipBin == 6  &&  th > 28.5  &&  ph < (th - ctPt[21][0]) * (ctPt[21][3] - ctPt[21][1])/(ctPt[21][2] - ctPt[21][0]) + ctPt[21][1]) patchCt62 = false;
  if(ipBin == 7  &&  th > 27.5  &&  ph > (th - ctPt[23][0]) * (ctPt[23][3] - ctPt[23][1])/(ctPt[23][2] - ctPt[23][0]) + ctPt[23][1]) patchCt71 = false;
  if(ipBin == 8                 &&  ph > (th - ctPt[24][0]) * (ctPt[24][3] - ctPt[24][1])/(ctPt[24][2] - ctPt[24][0]) + ctPt[24][1]) patchCt81 = false;
  if(ipBin == 9  &&  th > 23.0  &&  ph > (th - ctPt[26][0]) * (ctPt[26][3] - ctPt[26][1])/(ctPt[26][2] - ctPt[26][0]) + ctPt[26][1]) patchCt91 = false;
  if(ipBin == 10 &&  th > 21.0  &&  ph > (th - ctPt[28][0]) * (ctPt[28][3] - ctPt[28][1])/(ctPt[28][2] - ctPt[28][0]) + ctPt[28][1]) patchCtA1 = false;
  if(ipBin == 11 &&  th > 10.0  &&  ph > (th - ctPt[30][0]) * (ctPt[30][3] - ctPt[30][1])/(ctPt[30][2] - ctPt[30][0]) + ctPt[30][1]) patchCtB1 = false;
  if(ipBin == 11                &&  ph > (th - ctPt[31][0]) * (ctPt[31][3] - ctPt[31][1])/(ctPt[31][2] - ctPt[31][0]) + ctPt[31][1]) patchCtB2 = false;
  if(ipBin == 11 &&  th < 21.5  &&  ph < (th - ctPt[32][0]) * (ctPt[32][3] - ctPt[32][1])/(ctPt[32][2] - ctPt[32][0]) + ctPt[32][1]) patchCtB3 = false;
 
  if(thMaxCt == false || patchCt21==false || patchCt22==false || patchCt31==false || patchCt32==false || patchCt33==false || patchCt41==false || 
     patchCt42==false || patchCt51==false || patchCt52==false || patchCt61==false || patchCt62==false || patchCt71==false || patchCt81==false || 
     patchCt91==false || patchCtA1==false || patchCtB1==false || patchCtB2==false || patchCtB3==false) return false;
  else return true;
}



double thVtxMaxCutsOldBinsEbi1[13] = {19.0, 19.0, 35.0, 35.0, 35.0, 34.5, 33.5, 33.0, 31.5, 30.0, 25.0, 24.0, 25.0};
//End points of the Inclinded straight lines defining some cuts:
const int stLnNumEbi1=30; //# of straight lines defining the cuts (will be used to make and draw corresponding TLines as well)
double ctPtEbi1[stLnNum][4] = { //[4] is for cut points (x1,y1), (x2, y2) (points defining the cut-line)
  {11.0,13.0,  19.0,11.0},                          //invPBin=1: Two lines for the upper patch cut
  {14.0,-13.0, 19.0,-10.0},                         //invPBin=1: Two lines for the lower patch cut

  {15.0,14.0, 15.0,12.5},  {15.0,12.5,  31.0,12.5}, //invPBin=2: Two lines for the patch cut on the top region 
  {31.0,12.5, 31.0,1.0},   {31.0,1.0,   36.0,4.0},  //invPBin=2: Two lines for the patch cut on top-right region 
  {30.0,-3.0, 30.0,-14.0}, {30.0,-3.0,  36.0,-8.0}, //invPBin=2: Two lines for the patch cut on bottom-right region 

  {15.0,14.0, 15.0,12.5},  {15.0,12.5,  31.0,12.5}, //invPBin=3: Two lines for the patch cut on the top region 
  {31.0,12.5, 31.0,1.0},   {31.0,1.0,   36.0,4.0},  //invPBin=3: Two lines for the patch cut on top-right region 
  {30.0,-3.0, 30.0,-14.0}, {30.0,-3.0,  36.0,-8.0}, //invPBin=3: Two lines for the patch cut on bottom-right region 

  {30.0,12.5, 30.0,3.0},   {30.0,3.0,   36.0,3.0},  //invPBin=4: Two lines for the patch cut on top-right region 
  {30.0,-4.0, 30.0,-14.0}, {30.0,-4.0,  36.0,-8.0}, //invPBin=4: Two lines for the patch cut on bottom-right region 
 
  {30.0,13.5, 30.0,4.0},   {30.0,4.0,   35.0,3.0},  //invPBin=5: Two lines for the patch cut on top-right region 
  {30.0,-4.0, 30.0,-14.0}, {30.0,-4.0,  35.0,-8.0}, //invPBin=5: Two lines for the patch cut on bottom-right region 

  {29.0,13.5, 29.0,5.0},   {29.0,5.0,   35.0,2.0},  //invPBin=6: Two lines for the patch cut on top-right region 
  {29.0,-8.0, 29.0,-14.0}, {29.0,-8.0,  35.0,-5.0}, //invPBin=6: Two lines for the patch cut on bottom-right region 

  {27.0,13.0, 33,2.0},                              //invPBin=7: Two lines for the patch cut on the top region 

  {24.0,14.0, 32,-4.5},                             //invPBin=8: Two lines for the patch cut on the top region 
  {8.5,13.5, 10.0,8.0},   {10.0,8.0,    17.0,13.5}, //invPBin=8: Two lines for the patch cut on bottom-right region 
};




//https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/MoreCts/moreFiducialCutsMoreInversePBinsEbi1.gif //11/13/16
bool moreFidCutsWithMoreInvPBinsEbi1(int Ebi, double pp, double phDc1DegPlusMinus30, double thVtxDeg)//phDc1Deg must be between -30.0, 30.0
{
  double ph = phDc1DegPlusMinus30, th = thVtxDeg;
  double invPbinEdge[]={0.3,0.4,0.5,0.6,0.7,0.85,1.0,1.25,1.5,1.75,2.0,2.25,2.53,3.4};
  //Calculate invP & invPBin based on Ebi & p. 
  int ipBin=-1; double invP = ItorP[Ebi]/(pp*2250.0); 
  for(int i=0;i<13;i++) {if(invP > invPbinEdge[i] && invP <= invPbinEdge[i+1]) ipBin=i; if(ipBin>-1) break;}

  bool thMaxCt = true; if(th > thVtxMaxCutsOldBinsEbi1[ipBin]) thMaxCt = false;
  bool patchCt11=true, patchCt12=true, patchCt21=true, patchCt22=true, patchCt23=true, 
    patchCt31=true, patchCt32=true, patchCt33=true, patchCt41=true, patchCt42=true;
  bool patchCt51=true, patchCt52=true, patchCt61=true, patchCt62=true, patchCt71=true, 
    patchCt81=true, patchCt82=true; 
  //Straight line equation (joing points (x1,y1) & (x2,y2)):  y - y1 = ((y2-y1)/(x2-x1)) * (x - x1)
  //       or           y  = (x - x1)*(y2-y1)/(x2-x1)    + y1
  if(ipBin == 1  &&  ph > ((th - ctPtEbi1[0][0]) * (ctPtEbi1[0][3] - ctPtEbi1[0][1])/(ctPtEbi1[0][2] - ctPtEbi1[0][0]) + ctPtEbi1[0][1])) patchCt11 = false;
  if(ipBin == 1  &&  ph < ((th - ctPtEbi1[1][0]) * (ctPtEbi1[1][3] - ctPtEbi1[1][1])/(ctPtEbi1[1][2] - ctPtEbi1[1][0]) + ctPtEbi1[1][1])) patchCt12 = false;
  if(ipBin == 2  &&  th > 15.0  &&  ph > 12.5) patchCt21 = false;
  if(ipBin == 2  &&  th > 31.0  &&  ph > ((th - ctPtEbi1[5][0]) * (ctPtEbi1[5][3] - ctPtEbi1[5][1])/(ctPtEbi1[5][2] - ctPtEbi1[5][0]) + ctPtEbi1[5][1])) patchCt22 = false;
  if(ipBin == 2  &&  th > 30.0  &&  ph < ((th - ctPtEbi1[7][0]) * (ctPtEbi1[7][3] - ctPtEbi1[7][1])/(ctPtEbi1[7][2] - ctPtEbi1[7][0]) + ctPtEbi1[7][1])) patchCt23 = false;

  if(ipBin == 3  &&  th > 15.0  &&  ph > 12.5) patchCt21 = false;
  if(ipBin == 3  &&  th > 31.0  &&  ph > ((th - ctPtEbi1[11][0]) * (ctPtEbi1[11][3] - ctPtEbi1[11][1])/(ctPtEbi1[11][2] - ctPtEbi1[11][0]) + ctPtEbi1[11][1])) patchCt31 = false;
  if(ipBin == 3  &&  th > 30.0  &&  ph < ((th - ctPtEbi1[13][0]) * (ctPtEbi1[13][3] - ctPtEbi1[13][1])/(ctPtEbi1[13][2] - ctPtEbi1[13][0]) + ctPtEbi1[13][1])) patchCt32 = false;

  if(ipBin == 4  &&  th > 30.0  &&  ph > 3.0)  patchCt41 = false;
  if(ipBin == 4  &&  th > 30.0  &&  ph < ((th - ctPtEbi1[17][0]) * (ctPtEbi1[17][3] - ctPtEbi1[17][1])/(ctPtEbi1[17][2] - ctPtEbi1[17][0]) + ctPtEbi1[17][1])) patchCt42 = false;

  if(ipBin == 5  &&  th > 30.0  &&  ph > ((th - ctPtEbi1[19][0]) * (ctPtEbi1[19][3] - ctPtEbi1[19][1])/(ctPtEbi1[19][2] - ctPtEbi1[19][0]) + ctPtEbi1[19][1])) patchCt51 = false;
  if(ipBin == 5  &&  th > 30.0  &&  ph < ((th - ctPtEbi1[21][0]) * (ctPtEbi1[21][3] - ctPtEbi1[21][1])/(ctPtEbi1[21][2] - ctPtEbi1[21][0]) + ctPtEbi1[21][1])) patchCt52 = false;
 
  if(ipBin == 6  &&  th > 29.0  &&  ph > ((th - ctPtEbi1[23][0]) * (ctPtEbi1[23][3] - ctPtEbi1[23][1])/(ctPtEbi1[23][2] - ctPtEbi1[23][0]) + ctPtEbi1[23][1])) patchCt61 = false;
  if(ipBin == 6  &&  th > 29.0  &&  ph < ((th - ctPtEbi1[25][0]) * (ctPtEbi1[25][3] - ctPtEbi1[25][1])/(ctPtEbi1[25][2] - ctPtEbi1[25][0]) + ctPtEbi1[25][1])) patchCt62 = false;

  if(ipBin == 7  &&  ph > ((th - ctPtEbi1[26][0]) * (ctPtEbi1[26][3] - ctPtEbi1[26][1])/(ctPtEbi1[26][2] - ctPtEbi1[26][0]) + ctPtEbi1[26][1])) patchCt81 = false;
  if(ipBin == 8  &&  ph > ((th - ctPtEbi1[27][0]) * (ctPtEbi1[27][3] - ctPtEbi1[27][1])/(ctPtEbi1[27][2] - ctPtEbi1[27][0]) + ctPtEbi1[27][1])) patchCt71 = false;
  if(ipBin == 8  &&  ph > ((th - ctPtEbi1[28][0]) * (ctPtEbi1[28][3] - ctPtEbi1[28][1])/(ctPtEbi1[28][2] - ctPtEbi1[28][0]) + ctPtEbi1[28][1])
     &&  ph > ((th - ctPtEbi1[29][0]) * (ctPtEbi1[29][3] - ctPtEbi1[29][1])/(ctPtEbi1[29][2] - ctPtEbi1[29][0]) + ctPtEbi1[29][1])) patchCt82 = false;

  if(thMaxCt == false || patchCt11==false || patchCt21==false || patchCt22==false || patchCt23==false || 
     patchCt31==false || patchCt32==false || patchCt33==false || patchCt41==false || patchCt42==false || 
     patchCt51==false || patchCt52==false || patchCt61==false || patchCt62==false || patchCt71==false || 
     patchCt81==false || patchCt82==false) return false;
  else return true;
}
#endif
