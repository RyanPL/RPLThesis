#ifndef MOMENTUM_COR_RESULTS_PASS2
#define MOMENTUM_COR_RESULTS_PASS2

/**********
 *  Inputs:
 *      Ebi         = beam-energy index (0,1,2,3,4 for E_beam=1.0,1.3,2.0,2.3 & 3.0 GeV respectively)
 *      thR         = polar angle (theta) in radians (measured at vertex)  - after tracking correction
 *      phR         = azimuthal angle (phi) in radians (measured at vertex)- after tracking correction
 *      pOriginal   = the original/reconstructed momentum before correction
 *
 *  Output/Returns:
 *       pCorrected = the momentum (p) after the momentum correction is applied
 *
 *  The routine also uses another routine ECorrEG1bNoZave(..) for energy loss correction. The routine is defined
 *      in https://www.jlab.org/Hall-B//secure/eg4/adhikari/myHomeLinked/MyHm/LinkedFiles/outElossCorrNewRoutine.h
 *      or locally in /u/home/adhikari/LinkedFiles/outElossCorrNewRoutine.h
 *      and that file must be included for this routine to work.
 *
 *  If above statement is not clear enough about thR and phR, I want to reiterate that the input angles 
 *     'thR' and 'phR' are the tracking corrected angles expressed in radians. In other words, two of the 
 *     four outputs of the trackign correction routine go in as input in the momentum correction routine.
 *
 *  The tracking correction routine can be found at 
 *     https://www.jlab.org/Hall-B//secure/eg4/adhikari/myHomeLinked/MyHm/LinkedFiles/trackingCorResults4Pass2wd18parsNew.h
 *
 ***********/

//#include "/u/home/adhikari/LinkedFiles/outElossCorrNewRoutine.h"
double ECorrEG1bNoZave(double pp, double theta, double mass)
{
  double rad2deg = 180./TMath::Pi();  double deg2rad = TMath::Pi()/180.;
  double dE=0., DEDX=0.;  double gm=0.;
  double E = sqrt(pp*pp + mass*mass);
  double fbeta = pp/E, Beta2 = fbeta * fbeta, Gamma2 = 1.0 / ( 1.0 - Beta2 );  

  /*
    double zcenter    = -101.0; //-55.1;//#ignore
    double cR = zcenter - zave + 0.5;//#ignore
    //kp:  cR is the fraction of the distance the particle travelled i.e. cR = (zave - zcenter + 0.5)/Lt  (Lt = target length = 1.0 for EG1B)
    //     (denoted by delta_z in Nevzat's thesis), also in thesis, the order of zave -zcenter is wrong.
    
    //if( ( cR <= 1.0 ) && ( cR >= 0.0 ) )  gm = cR * 0.6 + 0.4;
    */

  /*
    kp: I think gm is total effective mass thickness traversed by the particle
    The factor 0.6 is effective mass thickness of NH3 (density of NH3 (about 1 g/cm^3) 
    multiplied by the packing fraction which is 
    roughly 0.6; See R. Fersch's thesis at page 215) and 0.4 is the sum of 0.3 
    and 0.1 where 0.4 is for mass thickness of He and 0.1 for that of window foils (see 
    Nevzat's thesis, page 158)
  */
  
  gm = 0.75;
  //Ignore all above lines with '#ignore' and use 'gm = 0.75' for all cases for a reasonable approximation (Dr. Kuhn)
    
  if( theta * rad2deg <= 45. ) gm = gm / cos( theta );
  else  gm = gm / cos( 45. * deg2rad );
   
  if( mass < 0.01) DEDX = 2.8; //For electrons
  else DEDX = 0.307 * (0.5/Beta2) * 
         ( log( 2.*511000.0*Beta2 * (Gamma2/90.) ) - Beta2 ); //For handrons
    
    
  // this is the energy that is lost in the target.
  // this much energy should be added back to measured energy to get original scattering energy
  dE =( DEDX / 1000.0 ) * gm;
  return dE;
}



double Pars[11] = {0.00361523,-0.00163235,-0.00847164,-0.0361214,-0.00132437,0.00472667,-0.0184257,0.0403751,-0.0743197,-0.0565383,0.0319114};
  //{0.00278808,-0.000449,-0.00490333,-0.0402083,-0.000561403,-0.00181541,-0.0186599,0.0424059,-0.0488417,-0.0633917,0.0233577};
double momCorPass2(int Ebi, double thR, double phR, double pOriginal)//, double *Par)//11/3/15
{
  /**/double Ebg[5] = {1.0539, 1.339, 1.9889, 2.256, 2.999};
  /**/double I_torus[5] = {1500, -1497, -2250, -2249, -2250};//TODO: may not be right
  /**/double mmPr = 0.93827;
  /**/double mmEl = 0.00051;
  double Eb=Ebg[Ebi]-0.002, I_tor = I_torus[Ebi];   //Beam energy in GeV and torus current in Amps
  double E_p = Eb/(1.0 + Eb*(1.0 - cos(thR))/mmPr); //Energy E_prime of scattered electron
  double B_torus = (0.76*I_tor*pow(sin(4*thR),2))/(3375*thR);//sin^2(4*x) = (1-cos(8*x))/2
  double r2d = 57.2957795, phD=phR*r2d, thD=thR*r2d, dp2p=0.0;
  //Both phD and phRd should be within the sector limits ie. -30 to 30 or -pi/6 to pi/6
  if(phD<-30.0) { phD=phD+360.0; phD = phD - 300.0; phR=phD/r2d; thR=thD/r2d;}


  int qq=-1; //electron charge (in units of 'e')
  double dE = ECorrEG1bNoZave(E_p,thR,mmEl);//Energy loss corr. 
  double E = sqrt(E_p*E_p+mmEl*mmEl); E=E-dE; 
  double pOEC = sqrt(E*E-mmEl*mmEl), dp2pOEL = (E_p-pOEC)/E_p;
  double pp = E_p; //Previously
  pp = pOEC;

  //Only the DC-dependent part is parameterized via Chi-sq minimization
  dp2p=((Pars[0]+Pars[1]*phR)*cos(thR)/cos(phR) 
	+(Pars[2]+Pars[3]*phR)*sin(thR))*pp/(qq*B_torus)
    + (Pars[4]*cos(thR) + Pars[5]*sin(thR)) + (Pars[6]*cos(thR) + Pars[7]*sin(thR) )*phR
    + 0.02*(Pars[8] + (Pars[9] + Pars[10]*phD/30.0 ) * pow(10.0/thD,3.0) );
  //return dp2pOEL + dp2p;
  double dpByp = dp2pOEL + dp2p;
  double pCorrected =  (1.0+dpByp)*pOriginal;
  return pCorrected; 
}

#endif
