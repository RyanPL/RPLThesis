#ifndef FIDUCIAL_CUTS_PASS2_H
#define FIDUCIAL_CUTS_PASS2_H

/**/double rad2deg = 180.0 / 3.14159;
/**/double ItorP[5] = {2250.0, 2250.0, 2250.0, 2250.0, 2250.0}; //TODO: Probably not the actual answer

bool aboveThisLine(double X, double Y, double x1, double y1, double x2, double y2)
{
  bool thisLineCut = true; 
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  if(Y < ( (X - x1) * (y2-y1) / (x2-x1) + y1) ) thisLineCut = false;
  return thisLineCut;
}

bool belowThisLine(double X, double Y, double x1, double y1, double x2, double y2)
{
  bool thisLineCut = true; 
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  if(Y > ( (X - x1) * (y2-y1) / (x2-x1) + y1) ) thisLineCut = false;
  return thisLineCut;
}

//A copy of /u/home/adhikari/AnaEG4Pass2/Cuts/FidCuts/fiducialCutsPass2.h  (3/14/16)
//Deleted  Pass2FidCutLatestFromPhpVsThpRegSimComp(..) which is not needed anymore

//========== thetaVtxRad - kinematic corrected theta in radians, phiDc1 is also in Rad

bool Pass2FidCutOnInvPvsThVtx(int Ebi, double pp, double thetaRadVtx) //3/8/16 (related to drawMyFidCutsOnIpVsThVtx(..) in retrieveProcessHistos.C) //11/11/16: ~/secure/Analysis/Pass2/Cuts/Fid/invMomVsThVtxPass2Ebi1RatioRegByEConlyFidCut09.gif
{
  /*
    drawMyTLine(6.0,0.3,6.0,2.6,2,lCol,lSty);//(x1,y1,x2,y2,width,col,style); th_vtx>6.0 (vertical cut)
    drawMyTLine(6.0,0.55,9.0,0.25,2,lCol,lSty);//Lower angular cut
    drawMyTLine(6.0,2.35,8.0,2.17,2,lCol,lSty); //Upper left angular cut
    drawMyTLine(8.0,2.17,9.0,2.17,2,lCol,lSty); //Upper horizontal cut
    drawMyTLine(9.0,2.17,15.0,2.6,2,lCol,lSty); //Upper right angular cut
    drawMyTLine(23.5,0.1,25.5,2.6,2,lCol,lSty); //Left cut around hole at theta~24 degree
    drawMyTLine(25.0,0.1,28.0,2.6,2,lCol,lSty); //Right cut around hole at theta~24 degree
  */
  double thDeg = rad2deg*thetaRadVtx;    if(pp==0.0) return false;
  double invP =0.0; if(pp>0.0) invP = ItorP[Ebi]/(pp*2250.0);   double X=thDeg, Y = invP;
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1

  //the near vertical cut (at about going between 5.5 and 6.0 deg in thDC1)
  bool vertCt =true; if( X < 6.0 ) vertCt=false;
  bool angLower=true; if( Y < ( (X-6.0)*(0.25-0.55)/(9.0-6.0) + 0.55 ) ) angLower=false;
  bool angUpperLeft=true; if(Y > ( (X-6.0)*(2.17-2.35)/(8.0-6.0) + 2.35 ) ) angUpperLeft=false;
  bool horizCt =true; if(Y > 2.17) horizCt=false;
  bool angUpperRight=true; if(Y > ( (X-9.0)*(2.6-2.17)/(15.0-9.0) + 2.17 ) ) angUpperRight=false;


  bool thHoleCutLeft=true; if(Y < ( (X-23.5)*(2.6-0.1)/(25.5-23.5) + 0.1 ) ) thHoleCutLeft=false;
  bool thHoleCutRight=true; if(Y > ( (X-25.0)*(2.6-0.1)/(28.0-25.0) + 0.1 ) ) thHoleCutRight=false;


  //If we simply AND all of above cuts, it may leave out the good corners such as those on the left of upper hole
  //   therefore, I am doing it in careful steps as follows:
  bool upperHoleCut=true; if(angUpperLeft==false && horizCt==false && angUpperRight==false) upperHoleCut=false;

  bool set1Cut=true; if(upperHoleCut==false || vertCt==false || angLower==false)     set1Cut=false;
  bool set2Cut=true; if(thHoleCutLeft==false && thHoleCutRight==false)               set2Cut=false;

  if(set1Cut==true && set2Cut==true) return true; else return false;
}


//========== thetaVtxRad - kinematic corrected theta in radians, phiDc1 is also in Rad
bool Pass2FidCutVersion0(double phiRadDC1, double thetaRadVtx) //2/16/16 (related to drawMyFidCutsV0(..) in retrieveProcessHistos.C)
{
  //https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/fidCutPlotsSet2_Eb2_Ratio.gif
  //https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/fidCutPlotsSet2_Eb1_Ratio.gif
  /*
    void drawMyFidCuts(int lCol)
    {
    drawMyTLine(6.0,-13.0,6.0,13.0,2,lCol,7);  //(x1,y1,x2,y2,width,col,style); th_vtx>6.0

    drawMyTLine(13,13.0,40.0,13.0,2,lCol,7);   //phDc1<13.0 for theta>12
    drawMyTLine(13,-13.0,40.0,-13.0,2,lCol,7); //phDc1<13.0 for theta>12

    drawMyTLine(6.0,13.0,8.0,18.0,2,lCol,7);   drawMyTLine(8.0,18.0,13.0,13.0,2,lCol,7);   //Upper angular cuts
    drawMyTLine(6.0,-13.0,8.0,-18.0,2,lCol,7); drawMyTLine(8.0,-18.0,13.0,-13.0,2,lCol,7); //Upper angular cuts
    }
  */
  double thDeg = rad2deg*thetaRadVtx, phDeg = rad2deg*phiRadDC1; 
  if(phDeg<0.0) phDeg = phDeg + 360.0;  phDeg = phDeg - 300.0; double yy = phDeg, xx = thDeg;

  bool lowThetaCut = true; if(thDeg<6.0) lowThetaCut = false;
  bool phiCut = true; if( thDeg>13.0 && (phDeg<-13.0 || phDeg>13.0)) phiCut = false;

  bool lowerAngularCut = true, upperAngularCut = true; 
 
  double x1=6.0, y1=13.0, x2=8.0, y2=18.0, x3=13.0, y3=13.0;//Points defining the upper angular cuts
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  if( (xx>6.0 && xx<8.0 && (yy > ((xx-x1)*(y2-y1)/(x2-x1) + y1))) //Upper angular cut with +ve slope, the following one is with -ve slope
      ||   (xx>8.0 && xx<13.0 && yy > ((xx-x2)*(y3-y2)/(x3-x2) + y2))    
      ) upperAngularCut = false; 
 
  if( (xx>6.0 && xx<8.0 &&   (yy < ((xx-x1)*((-y2)-(-y1))/(x2-x1) + (-y1)) ) ) //All the y-coordinates have -ve values for lower cuts
      ||  
      (xx>8.0 && xx<13.0 &&  (yy < ((xx-x2)*((-y3)-(-y2))/(x3-x2) + (-y2))  ))
      ) lowerAngularCut = false;
  
  if(lowThetaCut==true && phiCut==true && upperAngularCut==true && lowerAngularCut==true) return true;
  else return false;
}

bool Pass2FidCutLatest(double phiRadDC1, double thetaRadVtx) //2/17/16 (related to drawMyFidCuts(..) in retrieveProcessHistos.C)
{
  //https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/fidCutPlotsSet2_Eb2_RatioBigger.gif &.../fidCutPlotsSet2_Eb1_Ratio.gif

  double thDeg = rad2deg*thetaRadVtx, phDeg = rad2deg*phiRadDC1; 
  if(phDeg<0.0) phDeg = phDeg + 360.0;  phDeg = phDeg - 300.0; double yy = phDeg, xx = thDeg;

  bool lowThetaCut = true; if(thDeg<6.0) lowThetaCut = false;
  bool phiCut = true; if( thDeg>8.0 && (phDeg<-13.0 || phDeg>13.0)) phiCut = false;

  bool lowerAngularCut = true, upperAngularCut = true; 
  double x1=6.0, y1=11.0, x2=8.0, y2=13.0;//Points defining the upper angular cuts, (x1,-y1) & (x2,-y2) will define the lower one.
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  if( xx>6.0 && xx<8.0 && (yy > ((xx-x1)*(y2-y1)/(x2-x1) + y1))) //Upper angular cut with +ve slope, 
    upperAngularCut = false;  
  if( xx>6.0 && xx<8.0 &&   (yy < ((xx-x1)*((-y2)-(-y1))/(x2-x1) + (-y1)) ) ) //All the y-coordinates have -ve values for lower cut
    lowerAngularCut = false;
  
  if(lowThetaCut==true && phiCut==true && upperAngularCut==true && lowerAngularCut==true) return true;
  else return false;
}

bool Pass2FidCutLatestFromRegECcomparison(int Ebi, double pp, double thetaRadDc1) 
{
  //related to https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/BkUp/invMomVsThDc1Pass2Ebi4Ratio.gif
  double ItorScaledInvP = ItorP[Ebi]/(pp*2250.0);//fabs(Itor)/(pp*2250.0);
  double thDeg = rad2deg*thetaRadDc1;
  double xx=thDeg, yy=ItorScaledInvP;
  bool verticalCut = true;  if(xx<4.0) verticalCut = false;


  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  double x1=0.0, y1=1.0, x2=11.0, y2=0.0;//Points defining the upper angular cuts, (x1,-y1) & (x2,-y2) will define the lower one.
  bool lowerAngularCut = aboveThisLine(xx, yy, x1, y1, x2, y2);
  x1=0.0; y1=2.5; x2=9.0; y2=2.0; 
  bool upperAngularCut = belowThisLine(xx, yy, x1, y1, x2, y2);
  bool upperHorizontalCut=true, upperVerticalCut=true;
  if(xx>6.5 && xx<9.0 && yy>2.2)  upperHorizontalCut=false;
  if(yy>2.2 && xx<0.9) upperVerticalCut=false;

  if(verticalCut==false || lowerAngularCut==false ||
     (upperAngularCut == false && upperHorizontalCut==false && upperVerticalCut==false)) return false; 
  else return true;
}

bool Pass2FidCutLatestFromRegECcomparisonVtx(int Ebi, double pp, double thetaRadVtx) 
{
  //related to https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/F2017/invMomVsThVtxPass2Ebi4RatioN.png
  double ItorScaledInvP = ItorP[Ebi]/(pp*2250.0);//fabs(Itor)/(pp*2250.0);
  double thDeg = rad2deg*thetaRadVtx;
  double xx=thDeg, yy=ItorScaledInvP;
  bool verticalCut = true;  if(xx<6.0) verticalCut = false;

  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1
  double x1=0.0, y1=1.1, x2=11.0, y2=0.0;//Points defining the upper angular cuts, (x1,-y1) & (x2,-y2) will define the lower one.
  bool lowerAngularCut = aboveThisLine(xx, yy, x1, y1, x2, y2);
  x1=6.0; y1=2.35; x2=8.5; y2=2.2; bool upperAngularCut1 = belowThisLine(xx, yy, x1, y1, x2, y2);
  x1=8.5; y1=2.2; x2=15.0; y2=2.5; bool upperAngularCut2 = belowThisLine(xx, yy, x1, y1, x2, y2);
  bool upperHorizontalCut=true, upperVerticalCut=true;
  if(xx>6.5 && xx<9.0 && yy>2.2)  upperHorizontalCut=false;
  if(yy>2.2 && xx<0.9) upperVerticalCut=false;

  if(verticalCut==false || lowerAngularCut==false ||
     (upperAngularCut1 == false && upperAngularCut2 == false)) return false; 
  else return true;
}


/**///COMMENTED OUT BECAUSE IT WAS CAUSING PROBLEMS - Ryan L
/*bool Pass2FidCutLatestFrom_ExpBySimInInvPvsThDc1(int Ebi, double pp, double thetaDC1Rad) //2/21/16 (related to DrawAllCutLines_OnExpBySimInInvPvsThDc1(..) in retrieveProcessHistos.C)
{
  /* 
     https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Fid/invMomVsThDc1Pass2Ebi2RatioExpBySim.gif and 
     invMomVsThDc1Pass2Ebi2RatioExpBySim.gif 

     //iP definition must be the same as used in filling h2IpVthDc1 histo in  makeSimHistos2Compare.C & makeClasHistos2Compare.C
     //double ItorScaledInvP = ItorP[Ebi]/(pp*2250.0);//fabs(Itor)/(pp*2250.0); with pp=p[0] for sim and ppC for exp.
  *//*
  if(pp==0.0) return false;
  double invP =0.0; if(pp>0.0) invP = ItorP[Ebi]/(pp*2250.0);   double X=thetaDC1Rad*rad2deg, Y = invP;
  // Equation joining points (x1,y1) & (x2,y2):    y = (x-x1)*(y2-y1)/(x2-x1) + y1

  //the near vertical cut (at about going between 5.5 and 6.0 deg in thDC1)
  bool vertCt1 =true; if( Y < ( (X-6.0)*(2.9-0.1)/(5.5-6.0) + 0.1 ) ) vertCt1=false;
  bool angLower=true; if( Y < ( (X-12.0)*(1.0-0.1)/(0.0-12.0) + 0.1 ) ) angLower=false;
  bool angUpper=true; if(Y > ( (X-6.8)*(2.8-2.22)/(0.0-6.8) + 2.22 ) ) angUpper=false;
  bool horizCt =true; if(Y>2.22) horizCt=false;
  bool vertCt2 =true; if(X<9.0) vertCt2=false;
  //If we simply AND all of above cuts, it may leave out the good corners such as those on the left of upper hole
  //   therefore, I am doing it in careful steps as follows:
  bool upperHoleCut=true; if(angUpper==false && horizCt==false && vertCt2==false) upperHoleCut=false;
  if(upperHoleCut==true && vertCt1==true && angLower==true) return true; else return false;
}*/

#endif
