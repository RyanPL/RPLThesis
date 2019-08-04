#ifndef ALL_KINE_COR_PASS2
#define ALL_KINE_COR_PASS2

#include "/home/ryanpl/Thesis/mycuts/trackingCorResults4Pass18parsNew.h"
#include "/home/ryanpl/Thesis/mycuts/trackingCorResults4Pass2wd18parsNewNoMemLeak.h"
#include "/home/ryanpl/Thesis/mycuts/momCorPass2.h"




//This file created and following blocks/routines added on 1/24/16 to make a single function call (to make life easier)
//double *returnTrackingCorOutput(int PID, int Eb_index, int qq, UShort_t Xadc, UShort_t Yadc, double pp, 
void returnTrackingCorOutput(int PID, int Eb_index, int qq, UShort_t Xadc, UShort_t Yadc, double pp, 
			     double Csx, double Csy, double Csz, double vZ, //Vars upto this except PID needed by RstCor
			     //double cxdc,double cydc,double xdc,double ydc, double zdc)
			     double cxdc,double cydc,double xdc,double ydc, double zdc, double target[])
{
  double phicor=0.0, vZcor=0.0, xr=0.0, yr=0.0;//, *target = new double[4];
  double theta = acos(Csz), phiPhys = atan2(Csy,Csx);
  nRastCor21parReturnXrYrToo(Eb_index, qq, Xadc, Yadc, pp, Csx, Csy, theta, phiPhys, vZ, phicor, vZcor, xr, yr);
  //target = target_info18par(PID, xr, yr, cxdc, cydc, xdc, ydc, zdc, pp, qq);
  target_info18par(PID, xr, yr, cxdc, cydc, xdc, ydc, zdc, pp, qq/*, target*/);//TODO: commented out target cuz problems
  //return target;
} 


//double *applyAllKineCorPass2(int PID, int Eb_index, int qq, UShort_t Xadc, UShort_t Yadc, double pp, 
void applyAllKineCorPass2(int PID, int Eb_index, int qq, UShort_t Xadc, UShort_t Yadc, double pp, 
			  double Csx, double Csy, double Csz, double vZ, 
			  //double cxdc,double cydc,double xdc,double ydc, double zdc)
			  double cxdc,double cydc,double xdc,double ydc, double zdc, double targetKC[])
{
  double phicor=0.0, vZcor=0.0, xr=0.0, yr=0.0,// *target = new double[4], *targetKC = new double[5];
    target[4];
  //target = returnTrackingCorOutput(PID, Eb_index, qq, Xadc, Yadc, pp, Csx, Csy, Csz, vZ, cxdc, cydc, xdc, ydc, zdc);
  returnTrackingCorOutput(PID, Eb_index, qq, Xadc, Yadc, pp, Csx, Csy, Csz, vZ, cxdc, cydc, xdc, ydc, zdc, target);

  double thR = acos(target[2]), phR = atan2(target[1], target[0]);
  double pCor =  momCorPass2(Eb_index, thR, phR, pp);
  for(int i=0;i<4;i++) targetKC[i] = target[i];   targetKC[4] = pCor;
  //return targetKC;
} 





#endif
