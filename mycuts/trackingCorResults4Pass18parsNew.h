#ifndef TRACKING_COR_RESULTS_PASS2_18PAR
#define TRACKING_COR_RESULTS_PASS2_18PAR
/*******************
*
*  The input and outputs to the following tracking correction routine for EG4 are as follows:
*       
*  Inputs: 
*      PID, xr, yr, cxdc, cydc, xdc, ydc, zdc, p, q
*          where,
*            PID           :  PID is equal to 0, or 1 for electron or proton/hadron (respectively)
*            xr, yr        :  Beam position (in cm) returned by my raster-correction routine 
*                               (but, remember, we wont use raster corrected vz)
*            cxdc,cydc,czdc:  direction cosines at DC1 (tl1 variables)
*            xdc,ydc,zdc   :  track hit positions at DC1 (tl1 variables)
*            p, q          :  momentum and charge of the track
*       
*  Outupts: 
*      target[4]
*         where, the four elements of 'target' array carry the corrected values of cx, cy, cz and vz
*                i.e. the corrected direction cosines and z-coordinate at the vertex     
*
*   
*  Note: P. Bosted chose in his tech. note & code the following notation
*              : subscripts 0 and f used to denote angles/positions at vertex and DC1/tl1_ respectively
********************/
double par18TrkC[18] = 
  //{2.91769e-05,-0.000185976,0.000656243,-1.71564e-05,-0.00138385,0.00180966,-4.7297e-06,9.8246e-05,0.000692423,0.00216883,0.00663362,-0.000829092,-1.05155e-05,-1.28438e-06,5.5983e-08,2.02356e-08,0.133172,-0.00247964}; 
  {0.000181606,-0.000280946,-0.000722115,0.00043964,-0.000272506,0.000193008,-6.97254e-06,-2.66085e-06,0.000700284,0.00220063,0.00373464,0.000131815,2.35626e-06,9.88537e-08,1.11373e-06,2.50411e-07,-0.00301914,-0.000505913};
double *target_info18par(int PID,double xr,double yr,double cxdc,double cydc,double xdc,double ydc,
		    double zdc,double p,int q)
{
  double cxc, cyc, czc, thf, dth, phif, xp, yp, th0, phi0;//th,phi,vznom
  double rdc, fc, vz, czdc, temp, temp2, temp3, *target = new double[4];
  double targsign=1.0;  // sign of target field compared to simulations  

  bool electronID=false, protonID=false;   if(PID==0) electronID=true;   else if(PID==1) protonID=true; 
 
  // protect against garbage in, garbage out. This is needed because some tracks seem to have p=0
  cxc=0; cyc=0; czc=1; vz=0; target[0] = 0; target[1]=0; target[2]=0; target[3]=0;
  if(p>=0.001  &&  (q==1 || q==-1) && (pow(cxdc,2) + pow(cydc,2))<1 )
    {
      fc = 0.995;//fc=B/50 with B in kG overall integral B*dl correction factor. Our B=4.97T, so "field corr." is 0.995
      
      // Get nominal th and phi at DC1 in radians
      temp = 1. - pow(cxdc,2) - pow(cydc,2);
      if(temp>0) czdc = sqrt(temp); else czdc = 0.;   thf = acos(czdc);  phif = atan2(cydc, cxdc);


      if(electronID==true)   thf = thf + (par18TrkC[0] + par18TrkC[1]*phif)*cos(thf)/cos(phif) + (par18TrkC[2]+par18TrkC[3]*phif)*sin(thf);
      else if(protonID==true)thf = thf + (par18TrkC[4] + par18TrkC[5]*phif);


      // Get phi0 without raster correction yet
      phi0 = phif + targsign * fc * (0.186+par18TrkC[10] + (0.045+par18TrkC[11])*pow(thf,2) + (0.008+par18TrkC[12])*pow(thf,3) +(0.0032+par18TrkC[13]) * pow(thf,3) / pow(p,2)) * q / p; 
   
      // correction to polar angle from focusing effect. First, get focusing term for beam (x,y)=0.
      dth = fc * (0.90*thf + 1.2*pow(thf,3)) / (100.0 * pow(p,2));
       
      // displacement of beam along trajectory (xp) and perpendicular to it (yp)
      xp =  xr * cos(phi0) + yr * sin(phi0);
      yp = -(xr + par18TrkC[6]) * sin(phi0) + (yr + par18TrkC[7]) * cos(phi0);
      
      // correction to dth from radial target field, which only depends on raster x and y but not vertex z
      // Also, no effect on peak at zero!
      //dth = dth * (1. + targsign * q * p * (0.5 / thf) * (yp/0.75));
      //dth = dth * (1. + targsign * q * p * (1.0/ thf) * (yp/1.5));//0.5/0.75=1.0/1.5=0.667
      dth = dth*(1.+targsign*q*p*(1.0/thf)*yp*(par18TrkC[17]+0.667)); //8/30/15:SEK suggestion: use sigma=0.2 for p17

      
      // Now can get cz
      th0 = thf + dth;      czc = cos(th0);
      
      //Now phi0 again, this time including raster correction
      phi0 = phif + targsign * fc * (0.186+par18TrkC[10] + (0.045+par18TrkC[11])*pow(thf,2) + (0.008+par18TrkC[12])*pow(thf,3) +(0.0032+par18TrkC[13]) * pow(thf,3) / pow(p,2)) * q / p* (1. - (0.09+par18TrkC[14]) * ((0.35+par18TrkC[15])/thf) * xp ); 

      // Get cx and cy using this cz
      cxc = sin(th0) * cos(phi0);
      cyc = sin(th0) * sin(phi0);

      // renomalize czc
      //temp2=1.0-cxc*cxc-cyc*cyc;  if(temp2>0)czc=sqrt(temp2); else czc=0.0; //Disabled on 8/30/15 (c/o SEK)

      // Apply target field rotation correction     
      cxc = cxc - targsign * double(q) * czc * par18TrkC[8] / p;
      cyc = cyc + targsign * double(q) * czc * par18TrkC[9] / p;

      // renomalize czc
      temp3=1.0-cxc*cxc-cyc*cyc;  if(temp3>0)czc=sqrt(temp3); else czc=0.0;
      th0 = acos(czc);
      
      // Get vertex z in cm
      rdc = sqrt(pow(xdc - xr,2) + pow(ydc - yr,2));

      vz = zdc - (rdc -  (22.+par18TrkC[16]) * cos(th0) * (tan(th0) - tan(thf))) /tan(thf);

      target[0]= cxc;    target[1] = cyc;    target[2] = czc;    target[3] = vz;
    }
  
  return target;
}


//
//kp:  Following is the full copy of the comments from P. Bosted's orgininal correction code
//          http://www.jlab.org/Hall-B/secure/eg1-dvcs/software/targetinfo/target.C
//
//
//
// EG1-DVCS experiment
// Function: find direction cosines at target using DC1 c           anlges, taking into account raster position
//           also finds vertex z position using DC1 info
//           and beam position from raster information
// Usage: First call beam_info routine to get xr and yr
//        Can be used to over-write cx,cy,cz,vz in ntuples
// Inputs: irun: run number (integer)
//         xr: beam position in x from raster in cm (real)
//         yr: beam position in y from raster in cm (real) 
//         cxdc: cx at DC1 (trl_cx in nt10 or nt22) (real)
//         cydc: cy at DC1 (trl_cy in nt10 or nt22) (real)
//         xdc: x at DC1 (cm) (trl_x in nt10/nt22) (real)
//         ydc: y at DC1 (cm) (trl_y in nt10/nt22) (real)
//         zdc: z at DC1 (cm) (trl_z in nt10/nt22) (real)
//         p: particle momentum in GeV (real)
//         q: particle charge (-1 or 1) (integer)
// Outputs: cxc: cx at target (real)
//          cyc: cy at target (real)
//          czc: cz at target (real)
//           vz: vz at target in cm (real)
// Documentation: EG1-DVCS TN-004 (April, 2010)
// Author: Peter Bosted
// Date: April 20, 2010
////////////////////////////////////////////////////////////
// C++ translation by Erin Seder
// Define in main program:
//    #include "target.C"
//    float *target;
// during event loop get xy[0], xy[1] from raster correction subroutine and define the run_number
// then call during particle loop from nt22/root22 variables as:
//    target = target_info(run_number,xy[0],xy[1],tl1_cx[jj],tl1_cy[jj],tl1_x[jj],tl1_y[jj] ,tl1_z[jj] ,p[jj],q[jj]);
// this will return an array of 4 variables:
//    target[0] = new cx,
//    target[1] = new cy,
//    target[2] = new cz, 
//    target[3] = new vertex
// Please note this is only for CHARGED PARTICLES


#endif


/*5/30/15: Following function is copied from /u/site/www/html/Hall-B/secure/eg4/adhikari/Corrections/Mcor/Pass2/TrackingCor/getInclusiveEvents.C so that I could be using it from this common place in several programs that I expect to write in future EG4 data analysis.
Uses double nwPar21[] = {3.42931e-05, 0.000422835, -4.07979e-05, -0.000420968, 4330.72, 4440.27, -101.022, 4129.42, 4340.07, -101.064,4104.95, 4251.29, -101.041, 4014.65, 3946.27, -101.052, 4192.43, 3343.59, -101.002, 0.00096982, 0.00201097}; 
from LinkedFiles/NewEG4_MomCorNw.h (which in turn is included through LinkedFiles/myOwnHeaderFiles.h
*/
#ifndef NewEG4_MomCorNw_H
double nwPar21[] = {3.42931e-05, 0.000422835, -4.07979e-05, -0.000420968, 4330.72, 4440.27, -101.022, 4129.42, 4340.07, -101.064,4104.95, 4251.29, -101.041, 4014.65, 3946.27, -101.052, 4192.43, 3343.59, -101.002, 0.00096982, 0.00201097}; 
#endif


#ifndef __nRastCor21parReturnXrYrToo__
#define __nRastCor21parReturnXrYrToo__
//8/2/14: nRastCor21par copied from .../LinkedFiles/NewEG4_MomCorNw.h:254 and modified a little to make return xr & yr as well. 
//void nRastCor21par(int Eb_index, int q, UShort_t Xadc, UShort_t Yadc, double p, double Csx, double Csy, double theta,
void nRastCor21parReturnXrYrToo(int Eb_index, int q, UShort_t Xadc, UShort_t Yadc, double p, double Csx, double Csy, double theta,
				double phiPhys, double vZ, double & phicor, double & vZcor, double & xr, double & yr)
{
  /**/double Ebg[5] = {1.0539, 1.339, 1.9889, 2.256, 2.999};
  double rad2deg = 180.0/TMath::Pi(); double Ebeam = Ebg[Eb_index] - 0.002;  //theta = acos(cz[0]); if (theta==0 || theta<0) return;
  double phiDeg = rad2deg*phiPhys;                                         // phiDeg = rad2deg*atan2(cy[0],cx[0]); 
  if(phiDeg<-30.0) phiDeg = phiDeg + 360.0; int Sector = int((phiDeg - 30.0)/60.0 + 2);
  //cout<<rad2deg*phiPhys<<" "<<phiDeg<<" "<<Sector<<endl;//8/3/15
  double phiRadPhys = phiDeg/rad2deg;
  double phi_sec = (Sector - 1)*60.0/rad2deg;//Sector angle given in radians
  
  //Directly applying the new 21 par raster correction
  double Xgain = nwPar21[0] + nwPar21[1]/Ebeam, Ygain = nwPar21[2] + nwPar21[3]/Ebeam;
  double Xoffset = nwPar21[3*Eb_index + 4], Yoffset = nwPar21[3*Eb_index + 5]; //Zave=nwPar21[3*Eb_index + 6] (see in fcn7par21())
  double x = (Xadc*1.0 - Xoffset)*Xgain, y = (Yadc*1.0 - Yoffset)*Ygain;//Factor 1.0 is just to convert UShort_t to double
  double Xp = (x*cos(phi_sec) + y*sin(phi_sec))/cos(phiRadPhys - phi_sec);//phi here is actually (phi - phi_sec) as needed
  double factor1 = int(q)*10.0*(nwPar21[19]*Csx+ nwPar21[20]*Csy)/(p*sin(theta)*sin(theta));
  double factor2 = (1 - 0.25*tan(theta)*tan(theta))/ (sin(theta)*cos(theta));
  vZcor = vZ + Xp/tan(theta) + factor1*factor2;
 
  phicor = 0.0;

  double Pt = p*sin(theta);
  phicor = phicor - q*50.0*(Xp/100.0)/(33.356*Pt);//Phi_cor = phi - q*50.0*(Xp/100.0)/(33.356*Pt);//Total correction (the -ve part from other pars)
  xr = x; yr = y;
}
#endif

