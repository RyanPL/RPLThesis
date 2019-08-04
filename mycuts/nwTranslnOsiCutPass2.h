/******
 *  Two things changed on 9/10/15:
 *    1) initialized the ccmatch variable to true at the beginning
 *    2) Enabled the |phi_proj|<10.0 (in phi_cut3) back again (was disabled before)
 *
 *  11/14/15:
 *  Started from a copy of https://www.jlab.org/Hall-B//secure/eg4/adhikari/myHomeLinked/MyHm/LinkedFiles/ImpoRoutines/nwOsi/nwTranslnOsiCut.C
 */

using namespace std;
//#define DEBUG_HISTOS_PASS2 100 //Disable this when you are not debugging.

#ifndef MY_OSIP_CUT_4EL_PASS2_H
#define MY_OSIP_CUT_4EL_PASS2_H

//#include "/home/hkkang/Analysis/macros/osi_include.C"
#include "/u/home/adhikari/public_html/EG4/ImpoRoutines/dt_sc_cc_cut.h"
//#include "/u/home/adhikari/public_html/EG4/ImpoRoutines/cc_segm_geom.h"
#include "/u/home/adhikari/public_html/EG4/ImpoRoutines/cc_segm_geomNw.h"

/*
void CCVCRPL(float R0[3], float DIR[3], float P[3], int &IRES, float &S, float  R[3]);
void sphere_inter(int isec, float xo, float yo, float zo, float dx, float dy, float dz, float &ar, float &theta, float &phi, int &ires);
//sphere_inter(dc_sect[dc[k] -1], dc_xsc[dc[k] -1] , dc_ysc[dc[k] -1], dc_zsc[dc[k] -1], dc_cxsc[dc[k] -1], dc_cysc[dc[k] -1], dc_czsc[dc[k] -1], ar, theta, phi, ires);
*/

const int nSec=6, nSeg=18; //We've to define these array histos (in a function, called from main().)
//const double dtEcCcCt[18]={};

//#define PRINT_SOME_INFO
#ifdef DEBUG_HISTOS_PASS2
/*
//cout<<"kp debug"<<endl; //Showed some errors
TH1F *h1all_dbg = new TH1F("all_dbg", "all_nphe w/ 1st two grp cuts",310,-10.0,300.0);
TH1F *h1allcchit1_dbg = new TH1F("allcchit1_dbg", "nphe with cchit=1",310,-10.0,300.0);
TH1F *h1allcchitM_dbg = new TH1F("allcchitM_dbg", "nphe with cchit>1",310,-10.0,300.0);
TH1F *h1allosi_dbg = new TH1F("allosi_dbg", "all_nphe w/ osi cut",310,-10.0,300.0);
TH1F *h1allthP_dbg = new TH1F("allthP_dbg", "nphe with thP cut added",310,-10.0,300.0);
TH1F *h1allph1_dbg = new TH1F("allph1_dbg", "nphe with ph1 cut added",310,-10.0,300.0);
TH1F *h1allph2_dbg = new TH1F("allph2_dbg", "nphe with ph2 cut added",310,-10.0,300.0);
TH1F *h1allph3_dbg = new TH1F("allph3_dbg", "nphe with ph3 cut added",310,-10.0,300.0);
TH1F *h1allph4_dbg = new TH1F("allph4_dbg", "nphe with ph4 cut added",310,-10.0,300.0);
TH1F *h1alldtScCc_dbg = new TH1F("alldtScCc_dbg", "nphe with dtScCc cut added",310,-10.0,300.0);
TH1F *h1alldtEcCc_dbg = new TH1F("alldtEcCc_dbg", "nphe with dtEcCc cut added",310,-10.0,300.0);
TH1F *h1alldtScCcUprCt = new TH1F("alldtScCcUprCt", "nphe with dtScCc upper cut added",310,-10.0,300.0);
TH1F *h1alldtEcCcUprCt = new TH1F("alldtEcCcUprCt", "nphe with dtEcCc upper cut added",310,-10.0,300.0);
*/
//12/2/15: Want to replace above individuals with the array to minimize coding with loops
const int NphCts=13; 
TH1F *h1nphAr[4][NphCts]; /*12/2/15: [3][] for all, el-only, pi-only, [12] for above individuals 
			   Or without ec-cuts, with ecCt==1, & ecCt==0., 
			  Fourth for ecCt==0 candidates but using 'false' or -ve osi-component cuts.*/
char sNphCuts[NphCts][100]={"Grp1+2","cchit=1","cchit>1","osi","+thP i.e. #Delta#theta(seg)<3*#sigma",
			    "+ph1 i.e. |#Delta ipmt|<2","+ph2 i.e. |#phi(sec)|<30.0",
			    "+ph3 i.e. |#phi(proj)|<13.0","+ph4 i.e. |#phi(proj)|<16.0","+dtScCc(Lower)",
		     "+dtEcCc(Lower)","+dtScCc(Upper)","+dtEcCc(Upper)"};

//========== 11/14/15
TH1F *h1tl1Norm = new TH1F("tl1Norm","(norm-1.0) of DC1 dir. consines",200,-0.00007,0.00007);//>0.03
TH1F *h1iseg = new TH1F("h1iseg","iseg",80,0.0,20.0);
TH1F *h1ipmt = new TH1F("h1ipmt","ipmt",20,-2.0,2.0);
TH1F *h1dtScCcAll = new TH1F("h1dtScCcAll","sc_t - cc_t (time-walk & TOF corrected)",120,-7.0,5.0);//change binnings
TH1F *h1dtEcCcAll = new TH1F("h1dtEcCcAll","ec_t - cc_t (time-walk & TOF corrected)",120,-7.0,5.0);
TH1F *h1DphiSec = new TH1F("h1DphiSec","#phi w.r.t sector mid-plane",160,-40.0,40.0);
TH2F *h2PhsVsTh = new TH2F("h2PhsVsTh","#phi_{sec}(DC1) vs #theta(DC1)",160,0.0,40.0,100,-40.0,40.0);
TH1F *h1PhiProj = new TH1F("h1PhiProj","CC-projected #phi",200,-25.0,25.0);
TH2F *h2PhpVsThp = new TH2F("h2PhpVsThp","CC-projected #phi vs #theta",160,10.0,50.0,100,-25.0,25.0);

const int ctN=7; /* 0 - #Delta t(Ec,Cc)>2.5,    1 - with #Delta t(Cc,Cc)>2.0,  2-With  ecIn<0.05,
		  3- With ecIn<0.05 & etot/p<0.2, 4- with ecIn>0.05 & etot/p>0.2 & nphe>2.5,
		  3- With ecIn<0.05 & etot/p<0.2 & nphe<2.5, with dtScCc & dtEcCc > upper cuts*/
TH1F *h1ddPhi[nSec][nSeg][4],//*h1ddPhi[nSec][nSeg];//[4] for ipmt=all,-1,0,1
  *h1dTheta[nSec][nSeg]; //Confirm segment#, print th_mid[seg] & th_offset[sec][seg] & cuts on plots
TH1F *h1Dipmt[nSeg], *h1dtScCc[nSeg], *h1dtEcCc[nSeg], 
  *h1dtEcCcCt[ctN][nSeg], *h1dtScCcCt[ctN][nSeg], *h1npheCt[ctN][nSeg], *h1ecInCt[ctN][nSeg], 
  *h1sfCt[ctN][nSeg], *h1dtScCcNoTWc[nSeg], *h1dtEcCcNoTWc[nSeg], *h1TWc[nSeg];
/* Upper dt cuts decided based on plots such as https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Osi/Plots/partsWith_ecInGt0.05_sfGt0.2Pass2Ei0.gif which is for Ei=Eb_index=0 i.e. 1.0 GeV NH3 data */



void DefineHistogramsOsiTest(int Ebi)
{
  char hNm[100],hTtl[255]; int nBins=150; double hMn=-0.1,hMx=0.2;
  char sPMT[4][20]={"All ipmt(ph)","ipmt(ph)=-1","ipmt(ph)=0","ipmt(ph)=1"};
  for(int j=0;j<nSeg;j++)
    {
      sprintf(hNm,"DipmtSeg%02d",j+1); sprintf(hTtl,"#Delta ipmt = (ipmt(cc-segm) - ipmt(#phi)) for Seg=%d",j+1);
      h1Dipmt[j] = new TH1F(hNm,hTtl,40,-4.0,4.0);
      sprintf(hNm,"dtScCcSeg%02d",j+1); sprintf(hTtl,"sc_t-cc_t (time-walk & TOF corrected) for Seg=%d",j+1);
      h1dtScCc[j] = new TH1F(hNm,hTtl,120,-7.0,5.0);
      sprintf(hNm,"dtEcCcSeg%02d",j+1); sprintf(hTtl,"ec_t-cc_t (time-walk & TOF corrected) for Seg=%d",j+1);
      h1dtEcCc[j] = new TH1F(hNm,hTtl,120,-7.0,5.0);
      for(int i=0;i<nSec;i++)
	{
	  sprintf(hNm,"dThetaSec%dSeg%02d",i+1,j+1); sprintf(hTtl,"#Delta #theta for Sec=%d & Seg=%d",i+1,j+1);
	  h1dTheta[i][j] = new TH1F(hNm,hTtl,200,-5.0,5.0);
	  //sprintf(hNm,"ddPhiSec%dSeg%02d",i+1,j+1); sprintf(hTtl,"#Delta#phi for Sec=%d & Seg=%d",i+1,j+1);
	  //h1ddPhi[i][j] = new TH1F(hNm,hTtl,200,-10.0,10.0);
	  for(int ipmt=0;ipmt<4;ipmt++) //11/28/15
	    {
	      sprintf(hNm,"ddPhiSec%dSeg%02dPmt%d",i+1,j+1,ipmt);
	      sprintf(hTtl,"#Delta#phi for Sec=%d & Seg=%d, %s",i+1,j+1,sPMT[ipmt]);      
	      h1ddPhi[i][j][ipmt] = new TH1F(hNm,hTtl,150,-15.0,15.0); 
	      h1ddPhi[i][j][ipmt]->SetLineColor(ipmt+1);
	    }
	}
      //Histograms to see the effects of upper cuts on dtEcCc (g.t. about 2.0)
      for(int i=0;i<ctN;i++)//i for thre different cuts (one each for dtScCc,dtEcCc, & ec_i cuts)
	{
	  sprintf(hNm,"dtEcCcCt%dSeg%02d",i,j+1); sprintf(hTtl,"ec_t-cc_t (time-walk & TOF corrected) for Seg=%d",j+1);
	  h1dtEcCcCt[i][j] = new TH1F(hNm,hTtl,120,-3.0,7.0);
	  sprintf(hNm,"dtScCcCt%dSeg%02d",i,j+1); sprintf(hTtl,"sc_t-cc_t (time-walk & TOF corrected) for Seg=%d",j+1);
	  h1dtScCcCt[i][j] = new TH1F(hNm,hTtl,120,-3.0,7.0);
	  //Wanted to use h1dtCt[2] & h1dtCt[3] for dtScCc, & dtEcCc cutss, so the [3] histos below will remain unused
	  sprintf(hNm,"nphWdCt%dSeg%02d",i,j+1); sprintf(hTtl,"nphe for Seg=%d",j+1);
	  h1npheCt[i][j] = new TH1F(hNm,hTtl,310,-10.0,300.0);
	  sprintf(hNm,"ecInWdCt%dSeg%02d",i,j+1); sprintf(hTtl,"ec_in for Seg=%d",j+1);
	  h1ecInCt[i][j] = new TH1F(hNm,hTtl,200,0.0,0.3);
	  sprintf(hNm,"sfWdCt%dSeg%02d",i,j+1); sprintf(hTtl,"etot/p for Seg=%d",j+1);
	  h1sfCt[i][j] = new TH1F(hNm,hTtl,200,0.0,0.4);
	}
      //These three will be made with ec_in<0.05 and etot/p=sf<0.05
      sprintf(hNm,"dtScCcNoTWcSeg%02d",j+1); sprintf(hTtl,"sc_t-cc_t (no time-walk-cor) for Seg=%d",j+1);
      h1dtScCcNoTWc[j] = new TH1F(hNm,hTtl,160,-7.0,9.0);
      sprintf(hNm,"dtEcCcNoTWcSeg%02d",j+1); sprintf(hTtl,"ec_t-cc_t (no time-walk-cor) for Seg=%d",j+1);
      h1dtEcCcNoTWc[j] = new TH1F(hNm,hTtl,160,-7.0,9.0);
      sprintf(hNm,"TWcSeg%02d",j+1); sprintf(hTtl,"time-walk-cor = (Nphe-170)*0.005 for Seg=%d",j+1);
      h1TWc[j] = new TH1F(hNm,hTtl,120,-2.0,2.0);
    }

  for(int j=0;j<4;j++)//both,e-,pi-,pi- with -ve/false osi-cuts
    {
      for(int i=0;i<NphCts;i++)//different cuts
	{
	  sprintf(hNm,"nphP%dCts%02d",j,i); sprintf(hTtl,"Nphe with %s cuts",sNphCuts[i]);
	  h1nphAr[j][i] = new TH1F(hNm,hTtl,310,-10.0,300.0);
	}
    }
}
//========== 11/14/15
#endif



#ifdef PH_CUT2_STUDY
TH2F *h2phC_phDc1_all = new TH2F("phC_phDc1_all","phC vs phDc1 for all El, CCsect=6",320,200.0,360.0,320,200.0,360.0);
TH2F *h2phC_phDc1_rej = new TH2F("phC_phDc1_rej","phC vs phDc1 for El cut by phCut2",320,200.0,360.0,320,200.0,360.0);
TH1F *h1secDC_rej = new TH1F("secDC_rej","DC sector for El cut by phCut2", 32, 0.0, 8.0);
TH1F *h1secCC_rej = new TH1F("secCC_rej","CC sector for El cut by phCut2", 16, 0.0, 8.0);
TH1F *h1p_rej = new TH1F("p_rej","p of El cut by phCut2", 150, 0.0, 1.5);
TH1F *h1phP_rej = new TH1F("phP_rej","proj_ph of El cut by phCut2", 200, -20.0, 20.0);
TH1F *h1thP_rej = new TH1F("thP_rej","proj_th of El cut by phCut2", 200, 10.0, 50.0);
TH1F *h1nphe_rej = new TH1F("nphe_rej","nphe of El cut by phCut2", 360, -10.0, 350.0);
TH1F *h1W_rej = new TH1F("W_rej","Inv. mass of El cut by phCut2", 210, -0.1, 2.0);
TH2F *h2betaVp_rej = new TH2F("betaVp_rej","beta vs p of El cut by phCut2", 300, 0.0, 1.5, 200, 0.8, 1.2);
#endif

/*
double dtUprCtScCc[5][nSeg]= //[5] for 5 different Ebi, [nSeg] for different segments
  {
    //First 7 segments don't have data, so we're using 2.4 (max) as cuts, although they will be useless
    {2.4,2.4,2.4,2.4,2.4,2.4,2.4,   2.0,2.0,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4}, //Ebi=0
    {2.4,2.4,2.4,2.4,2.4,2.4,2.4,   2.0,2.0,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4}, //Ebi=1
    {2.4,2.4,2.4,2.4,2.4,2.4,2.4,   2.0,2.0,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4}, //Ebi=2
    {2.4,2.4,2.4,2.4,2.4,2.4,2.4,   1.6,1.6,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4}, //Ebi=3
    {2.4,2.4,2.4,2.4,2.4,2.4,2.4,   2.0,2.0,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4}  //Ebi=4
  };
double dtUprCtEcCc[5][nSeg]= //[5] for 5 different Ebi, [nSeg] for different segments
  {
    //First 7 segments don't have data, so we're using 2. (max) as cuts, although they will be useless
    {2.,2.,2.,2.,2.,2.,2.,   1.1,1.1,1.2,1.2,1.4,1.2,1.4,1.8,1.6,1.8,2.0}, //Ebi=0
    {2.,2.,2.,2.,2.,2.,2.,   1.1,1.0,1.2,1.2,1.6,1.2,1.4,1.8,1.6,1.8,2.0}, //Ebi=1
    {2.,2.,2.,2.,2.,2.,2.,   1.1,1.1,1.2,1.2,1.4,1.2,1.4,1.8,1.6,1.8,2.0}, //Ebi=2
    {2.,2.,2.,2.,2.,2.,2.,   1.1,1.1,1.2,1.2,1.4,1.2,1.6,2.0,1.4,1.8,2.0}, //Ebi=3
    {2.,2.,2.,2.,2.,2.,2.,   1.1,1.1,1.2,1.2,1.4,1.2,1.4,1.8,1.6,1.8,2.0}  //Ebi=4
  };
*/
//Now chaning all dt-upper cuts to 2.0 (upon Dr. Kuhn's suggestion)
double dtUprCtScCc[5][nSeg]= //[5] for 5 different Ebi, [nSeg] for different segments
  {
    //First 7 segments don't have data, so we're using 2.4 (max) as cuts, although they will be useless
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=0
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=1
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=2
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=3
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}  //Ebi=4
  };
double dtUprCtEcCc[5][nSeg]= //[5] for 5 different Ebi, [nSeg] for different segments
  {
    //First 7 segments don't have data, so we're using 2. (max) as cuts, although they will be useless
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=0
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=1
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=2
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}, //Ebi=3
    {2.,2.,2.,2.,2.,2.,2.,   2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.}  //Ebi=4
  };


bool cc_matchingPass2Test(int gpart, UChar_t cc, UChar_t dc, UChar_t sc, UChar_t ec, Char_t q, Char_t stat, Char_t dc_stat, 
#ifdef PH_CUT2_STUDY
			  float Eb, float pp, float ccx, float ccy, float ccz, float Tr_time,
#endif
#ifdef DEBUG_HISTOS_PASS2
			  int Ebi, double pp, double ec_in, double ec_out, double ec_tot,
#endif
			  float tlx, float tly, float tlz, UChar_t cc_sect, UShort_t cc_segm, float dc_ccth, float dc_ccph,
			  float sc_t, float sc_r, float cc_t, float cc_r, float ec_t, float ec_r, float beta, UChar_t dc_sect,
			  float dc_xsc , float dc_ysc, float dc_zsc, float dc_cxsc, float dc_cysc, float dc_czsc,
			  UShort_t Nphe, UChar_t CC_hit);
///////////////////////
////THIS IS THE ONE////actually I'm commenting it out to try to get it to compile
///////////////////////
/*bool OsipenkoCutsPass2Test(int Ebi)//Ebi added on 11/29/15
{
  return cc_matchingPass2Test(gpart, cc[0], dc[0], sc[0], ec[0], q[0], stat[0], dc_stat[dc[0]-1], 
#ifdef PH_CUT2_STUDY
			      Eb, p[0],cx[0],cy[0],cz[0], Tr_time,//This line absent before 11/22/15 (still not used)
#endif
#ifdef DEBUG_HISTOS_PASS2 
			      Ebi, p[0], ec_ei[ec[0]-1], ec_eo[ec[0]-1], etot[ec[0]-1],
#endif
			      tl1_cx[dc[0]-1], tl1_cy[dc[0]-1], tl1_cz[dc[0]-1],
			      cc_sect[cc[0]-1], cc_segm[cc[0]-1], dc_ccth[dc[0]-1], dc_ccph[dc[0]-1],
			      sc_t[sc[0]-1], sc_r[sc[0]-1], cc_t[cc[0]-1],
			      cc_r[cc[0]-1], ec_t[ec[0]-1], ec_r[ec[0]-1], b[0], dc_sect[dc[0]-1],
			      dc_xsc[dc[0]-1],dc_ysc[dc[0]-1],dc_zsc[dc[0]-1],dc_cxsc[dc[0]-1],dc_cysc[dc[0]-1],dc_czsc[dc[0]-1], 
			      nphe[cc[0]-1], cc_hit[cc[0]-1]);
}*/
//////////////////////////////////
////THIS IS THE END OF THE ONE////
//////////////////////////////////
bool
cc_matchingPass2Test(int gpart, UChar_t cc, UChar_t dc, UChar_t sc, UChar_t ec, Char_t q, Char_t stat, Char_t dc_stat,
#ifdef PH_CUT2_STUDY
		     float Eb, float pp, float ccx, float ccy, float ccz, float Tr_time,
#endif
#ifdef DEBUG_HISTOS_PASS2
		     int Ebi, double pp, double ec_in, double ec_out, double ec_tot, //added on 11/22/15
#endif
		 float tlx, float tly, float tlz, UChar_t cc_sect, UShort_t cc_segm, float dc_ccth, float dc_ccph,
		 float sc_t, float sc_r, float cc_t, float cc_r, float ec_t, float ec_r, float beta, UChar_t dc_sect,
		 float dc_xsc , float dc_ysc, float dc_zsc, float dc_cxsc, float dc_cysc, float dc_czsc, 
		 UShort_t Nphe, UChar_t CC_hit)
{
#ifdef DEBUG_HISTOS_PASS2
  double sf=ec_tot/pp; bool ecCt=false; if(sf>0.2 && ec_in>0.05) ecCt=true;
  h1nphAr[0][0]->Fill(Nphe); 
  if(ecCt==true)  h1nphAr[1][0]->Fill(Nphe); 
  else if(ecCt==false)  {h1nphAr[2][0]->Fill(Nphe); h1nphAr[3][0]->Fill(Nphe); }
#endif


  //kp: Defining constants
  //float pi=3.1415926536; float cc_thresh=2.5; 
  float cc_n=1.00153; //xc: cc_n is the index of refraction of the material(gas) inside the CC.
  //data pi/3.1415926536e0/,cc_thresh/2.5e0/,cc_n/1.00153e0/  
  float cl=29.9792458e0;//data cl/29.9792458e0/   //speed of light
  // common/cc_effn/nphe_array;   //http://home.ifi.uio.no/steikr/doc/f77/tutorial/common.html (Making variables global)
 

      
  //bool ccmatch=false;  //Initialize the ccmatch variable with the false value. (till 9/9/15)
  bool ccmatch=true;  //Initialize the ccmatch variable with the true value.//9/9/15


  //Cuts1
  //kp: All detector fired, and track welll reconstructed for the first particle //c Check the electron
  if(!(gpart>0 && cc>0 && dc>0 && sc>0 && ec>0 && int(q)==-1) ) return ccmatch = false;

  //Cuts2
  //if ( !( stat>0 && dc_stat>0 ) ) return ccmatch = false;
  if ( !( int(stat)>0 && int(dc_stat)>0 ) ) return ccmatch = false; //kp: 5/15/12


  //if(cc_hit(cc(1)).ne.1)  then return;
  /* //active till 11/29/15
  if(!(CC_hit==1))   //11/29/15
    {
#ifdef DEBUG_HISTOS_PASS2
      if(CC_hit>1)h1allcchitM_dbg->Fill(Nphe);
#endif
      return ccmatch==false;
    }
#ifdef DEBUG_HISTOS_PASS2
  else if(CC_hit==1) h1allcchit1_dbg->Fill(Nphe);
#endif
  */
  if(CC_hit<1) return ccmatch==false;
  else 
    {
#ifdef DEBUG_HISTOS_PASS2
      if(CC_hit==1) h1nphAr[0][1]->Fill(Nphe);//h1allcchit1_dbg->Fill(Nphe);
      else if(CC_hit>1)h1nphAr[0][2]->Fill(Nphe);//h1allcchitM_dbg->Fill(Nphe); //M for More or Multi?
      if(CC_hit==1) {if(ecCt==true) h1nphAr[1][1]->Fill(Nphe); 
	else if(ecCt==false) {h1nphAr[2][1]->Fill(Nphe); h1nphAr[3][1]->Fill(Nphe);}}
      else if(CC_hit>1) {if(ecCt==true) h1nphAr[1][2]->Fill(Nphe); 
	else if(ecCt==false) {h1nphAr[2][2]->Fill(Nphe);h1nphAr[3][2]->Fill(Nphe);}}
#endif
    }



  //kp: To make sure angles (at DC1) of the first particle are reasonably well constructed/normalized   //c Electron angles
  Float_t norm = (tlx*tlx + tly*tly + tlz*tlz); 
#ifdef PRINT_SOME_INFO
  cout<<"tlx^2, tly^2, tlz^2 norm "<<tlx*tlx<<" "<<tly*tly<<" "<<tlz*tlz<<" "<<norm<<endl;
#endif
#ifdef DEBUG_HISTOS_PASS2
  h1tl1Norm->Fill(norm-1.0);
#endif
  if(fabs(norm-1.0)>0.03) return ccmatch = false;

  //kp: Calculate angles at DC1
  float rad2deg = 57.2957795131e0;  //float th_El = acos(tlz/sqrt(norm))*rad2deg;  //float ph_El = atan2(tly,tlx);       //Before I had this
  float ph_El = atan2(tly,tlx)*rad2deg; //Had forgotten to convert this ph_El to deg
  if(ph_El<0.0) ph_El = ph_El + 360.0;




  //c Find out CC hit sector, hit segment and hit PMT
  int isec=cc_sect;


  /*
  //http://www.jlab.org/Hall-B/secure/eg4/xiaochao/pion/osicut/update_04262010/
  the only change to the source code is
  from
  iseg=nint(real(cc_segm(cc(1))-int(cc_segm(cc(1))/1000.)*1000)/10.)
  to
  iseg=(cc_segm(cc(1))-int(cc_segm(cc(1))/1000.)*1000)/10

  And the cuts have been updated in osipenko/cc_segm_geom.inc
  */


  /*
  //Before April 26

  int iseg = 0;  
  //iseg=nint(real(cc_segm(cc(1))-int(cc_segm(cc(1))/1000.)*1000)/10.)
  // float iseg_tmp = (cc_segm*1.0 - int(cc_segm*1.0/1000.0)*1000.0)/10.0;//Not putting .0 at the end of 1000 or 10 gave different results (different data type)
  float iseg_tmp = 1.0*(cc_segm - int(cc_segm/1000.0)*1000)/10.0;
  // http://icl.cs.utk.edu/projectsfiles/f2j/javadoc/org/netlib/util/StrictUtil.html
  //if(iseg_tmp>0) iseg = int(iseg_tmp + 0.5); else iseg = int(iseg_tmp - 0.5);//Initially it was used
  // if(iseg_tmp>=0) iseg = int(iseg_tmp + 0.5); else iseg = int(iseg_tmp - 0.5);
  if(iseg_tmp<0) iseg = int(iseg_tmp - 0.5); else iseg = int(iseg_tmp + 0.5);  //iseg = int(iseg_tmp);
  // Eqvt of nint(iseg_tmp) in fortran //kp: 02/17/10
  */



  //Now After April 26, 2010
  //Fortran: iseg=(cc_segm(cc(1))-int(cc_segm(cc(1))/1000.)*1000)/10
  int iseg = (cc_segm-int(cc_segm/1000.0)*1000)/10;
  int ipmt=int(cc_segm/1000.0)-1;      

#ifdef PRINT_SOME_INFO
  cout<<"iseg= "<<iseg<<endl;
#endif

  //kp: Convert phi_dc1  to the phi w.r.t. the sector mid plane  //c scattering phi within each sector
  float dphi_sec=0.0;
  //if(ph_El<330.0) dphi_sec=ph_El-(isec-1)*60.0;  else  dphi_sec=ph_El-360.0; 
  //Here angle is measured at DC1 and isec is for CC. This gives the loss in the 1st row of NH3n1054gof01_03.gif ( phi_cut2 )
  if(ph_El<330.0) dphi_sec=ph_El-(int(dc_sect)-1)*60.0;  else  dphi_sec=ph_El-360.0; //Replaced above line with this on June 05, 2010
  //if(int(dc_sect)<6) cout<<int(dc_sect)<<endl;


    
  //kp: Get/Calculate CC-projected angles  //c Calculate angles and the distance to CC imaginary plane 
  //call sphere_inter(dc_sect(dc(1)),dc_xsc(dc(1)),dc_ysc(dc(1)), dc_zsc(dc(1)),dc_cxsc(dc(1)),dc_cysc(dc(1)),dc_czsc(dc(1)),dist,theta,phi,idir);
  float theta =  dc_ccth;  float phi   =  dc_ccph;  //These values will be over-ridden by the values returned by the function call of 'sphere_inter()'
  // float ar =0.0; int ires = 0; sphere_inter(dc_sect, dc_xsc, dc_ysc, dc_zsc, dc_cxsc, dc_cysc, dc_czsc, ar, theta, phi, ires);//Gives slightly more events than using dc_ccth and dc_ccph from the pass1 ntuple
  float dphi_p=phi;


  double timeWalkCor=1.0*(Nphe-170)*0.005;//kp: vars defined in 11/23/15
  double ecTofCor=(ec_r - cc_r)/(cl*beta), scTofCor=(sc_r - cc_r)/(cl*beta);
  //Get three time differences  
  //Before I forgot to put cl*beta inside parentheses ()   //if(beta>1.0) beta=1.0;//Didn't make any difference.
  //float dt_ScCc = sc_t - cc_t - (sc_r - cc_r)/(cl*beta); //c Calculate dt SC-CC  //b is the beta of the particle.
  //float dt_ScCc = sc_t - cc_t - (sc_r - cc_r)/(cl*beta) - (Nphe-170)*0.005; //!time-walk correction, added 100407 (increases events slightly)
   float dt_ScCc = sc_t - cc_t - scTofCor - timeWalkCor; //!time-walk correction, added 100407 (increases events slightly) 

  //float dt_EcCc = ec_t - cc_t - (ec_r - cc_r)/(cl*beta); //c Calculate dt EC-CC
   float dt_EcCc = ec_t - cc_t - ecTofCor - timeWalkCor; //  ! time-walk correction, added 100407
#ifdef DEBUG_HISTOS_PASS2
  if(isec==6)
    {
      h1iseg->Fill(iseg);  h1ipmt->Fill(ipmt);
      h1PhiProj->Fill(dphi_p);   h2PhpVsThp->Fill(theta,dphi_p); //CC-projected phi, and theta
      h1dtScCcAll->Fill(dt_ScCc);  h1dtEcCcAll->Fill(dt_EcCc);
      h1dtScCc[iseg-1]->Fill(dt_ScCc); h1dtEcCc[iseg-1]->Fill(dt_EcCc);
      if(dt_EcCc>2.5)
	{
	  h1dtScCcCt[0][iseg-1]->Fill(dt_ScCc); h1npheCt[0][iseg-1]->Fill(Nphe);
	  h1ecInCt[0][iseg-1]->Fill(ec_in); h1sfCt[0][iseg-1]->Fill(ec_tot/pp);
	}
      if(dt_ScCc>2.0)
	{
	  h1dtEcCcCt[1][iseg-1]->Fill(dt_EcCc);  h1npheCt[1][iseg-1]->Fill(Nphe); 
	  h1ecInCt[1][iseg-1]->Fill(ec_in); h1sfCt[1][iseg-1]->Fill(ec_tot/pp);
	}
      if(ec_in<0.05)
	{
	  h1dtScCcCt[2][iseg-1]->Fill(dt_ScCc);  h1dtEcCcCt[2][iseg-1]->Fill(dt_EcCc);  
	  h1npheCt[2][iseg-1]->Fill(Nphe); h1ecInCt[2][iseg-1]->Fill(ec_in);   h1sfCt[2][iseg-1]->Fill(ec_tot/pp);
	}
      if(ec_in<0.05 && (ec_tot/pp)<0.2)
	{
	  h1dtScCcCt[3][iseg-1]->Fill(dt_ScCc);  h1dtEcCcCt[3][iseg-1]->Fill(dt_EcCc);  
	  h1npheCt[3][iseg-1]->Fill(Nphe);  h1ecInCt[3][iseg-1]->Fill(ec_in);  h1sfCt[3][iseg-1]->Fill(ec_tot/pp);
	  h1TWc[iseg-1]->Fill(timeWalkCor);
	  h1dtScCcNoTWc[iseg-1]->Fill(dt_ScCc+timeWalkCor);  h1dtEcCcNoTWc[iseg-1]->Fill(dt_EcCc+timeWalkCor);
	}
      if(ec_in>0.05 && (ec_tot/pp)>0.2 && Nphe>25.0)
	{
	  h1dtScCcCt[4][iseg-1]->Fill(dt_ScCc);  h1dtEcCcCt[4][iseg-1]->Fill(dt_EcCc);  
	  h1npheCt[4][iseg-1]->Fill(Nphe);  h1ecInCt[4][iseg-1]->Fill(ec_in);  h1sfCt[4][iseg-1]->Fill(ec_tot/pp);
	}
      if(ec_in<0.05 && (ec_tot/pp)<0.2 && Nphe<25.0)
	{
	  h1dtScCcCt[5][iseg-1]->Fill(dt_ScCc);  h1dtEcCcCt[5][iseg-1]->Fill(dt_EcCc);  
	  h1npheCt[5][iseg-1]->Fill(Nphe);  h1ecInCt[5][iseg-1]->Fill(ec_in);  h1sfCt[5][iseg-1]->Fill(ec_tot/pp);
	}
      if(dt_ScCc>dtUprCtScCc[Ebi][iseg-1] && dt_EcCc>dtUprCtEcCc[Ebi][iseg-1])
	{
	  h1dtScCcCt[6][iseg-1]->Fill(dt_ScCc);  h1dtEcCcCt[6][iseg-1]->Fill(dt_EcCc);  
	  h1npheCt[6][iseg-1]->Fill(Nphe);  h1ecInCt[6][iseg-1]->Fill(ec_in);  h1sfCt[6][iseg-1]->Fill(ec_tot/pp);
	}
    }
#endif

  //c Find geometrical PMT ID   -- xc note: pmt_hit_geom is not used
  int pmt_hit_geom = 1;
  float beta_el = beta; if(beta_el>1.0) beta_el=1.0;
  if(fabs(dphi_sec-phi_offset[isec-1][iseg-1])>(acos(1.0/(cc_n*beta_el))*rad2deg)) //What is cc_n? //kp
    {
      if((dphi_sec-phi_offset[isec-1][iseg-1])>0.0)  pmt_hit_geom=1;
      else    pmt_hit_geom=-1; 
    }
  else   pmt_hit_geom=0;




  //c     Hit PMT: left-right choice
  int pmt_hit = 0;  float constMid = 0.5;
  //if (!(isec==6)) //! sector 1-5
  if (isec<6) //! sector 1-5
    {
      if((dphi_p-phi_offset[isec-1][iseg-1])>0.0)             pmt_hit=1;
      //else if((dphi_sec-phi_offset[isec-1][iseg-1])<0.0)      pmt_hit=-1;//This was the line before: http://www.jlab.org/Hall-B/secure/eg4/xiaochao/pion/osicut/
      else if((dphi_p-phi_offset[isec-1][iseg-1])<0.0)      pmt_hit=-1;
    }
  else    //! sector 6
    {
      if((dphi_p-phi_offset[isec-1][iseg-1]) > (3.0*max(phi_wid[iseg-1],constMid)))         pmt_hit=1;
      //else if((dphi_sec-phi_offset[isec-1][iseg-1])<-3.0*max(phi_wid[iseg-1],constMid))  pmt_hit=-1;
      else if((dphi_p-phi_offset[isec-1][iseg-1]) < (-3.0*max(phi_wid[iseg-1],constMid)))  pmt_hit=-1;
#ifdef DEBUG_HISTOS_PASS2
      h1DphiSec->Fill(dphi_sec); h2PhsVsTh->Fill(acos(tlz)*rad2deg,dphi_sec);
      //h1ddPhi[isec-1][iseg-1]->Fill(dphi_p-phi_offset[isec-1][iseg-1]);
      h1ddPhi[isec-1][iseg-1][0]->Fill(dphi_p-phi_offset[isec-1][iseg-1]); //for all ipmt values
      if(pmt_hit==-1) h1ddPhi[isec-1][iseg-1][1]->Fill(dphi_p-phi_offset[isec-1][iseg-1]);
      else if(pmt_hit==0) h1ddPhi[isec-1][iseg-1][2]->Fill(dphi_p-phi_offset[isec-1][iseg-1]);
      else if(pmt_hit==1) h1ddPhi[isec-1][iseg-1][3]->Fill(dphi_p-phi_offset[isec-1][iseg-1]);
#endif 
    }

#ifdef PRINT_SOME_INFO
  cout<<"isec= "<<isec<<endl;
#endif

  //kp: Assumming that all cuts are off/false before   //c CC matching cuts
  bool theta_p_cut=false;   bool phi_cut1=false; bool phi_cut2=false;  bool phi_cut3=false;
  bool phi_cut4=false; bool phi_cut=false; bool dt_sc_cc_cut=false; bool dt_ec_cc_cut=false;
  bool dtScCcUprCt=false, dtEcCcUprCt=false;//11/29/15
  double DtScCcUprCt=2.0,DtEcCcUprCt=2.0; //12/13/15

  //#ifdef THP_CUT_NO_ADJUST
  if (isec<6)
    {
      if(fabs(theta-th_mid[iseg-1]-th_offset[isec-1][iseg-1]) < (3.0*th_wid[isec-1][iseg-1]))
	theta_p_cut= true;
      else theta_p_cut= false;
    }
  else  // !sector 6
    {
      if(fabs(theta-th_offset[isec-1][iseg-1]) < (3.0*th_wid[isec-1][iseg-1]))
	theta_p_cut= true;
      else theta_p_cut= false;
#ifdef DEBUG_HISTOS_PASS2
      h1dTheta[isec-1][iseg-1]->Fill(theta-th_offset[isec-1][iseg-1]);
      DtScCcUprCt=dtUprCtScCc[Ebi][iseg-1]; DtEcCcUprCt=dtUprCtEcCc[Ebi][iseg-1]; //12/13/15
#endif  
      if(dt_ScCc<DtScCcUprCt) dtScCcUprCt=true;
      if(dt_EcCc<DtEcCcUprCt) dtEcCcUprCt=true;
   }
  //#endif


  /*
#ifdef THP_CUT_ADJUST_M3P6
  if (isec<6)
    {
      if(fabs(theta-th_mid[iseg-1]-th_offset[isec-1][iseg-1]) < (3.0*th_wid[isec-1][iseg-1]))
	theta_p_cut= true;
      else theta_p_cut= false;
    }
  else  // !sector 6
    {
      if((theta-th_offset[isec-1][iseg-1]) > (-9.0*th_wid[isec-1][iseg-1])  && (theta-th_offset[isec-1][iseg-1]) < (9.0*th_wid[isec-1][iseg-1]) )
	theta_p_cut= true;
      else theta_p_cut= false;
    }
#endif
  */






#ifdef PRINT_SOME_INFO
  cout<<"debug2= "<<isec<<endl;
#endif

  // if (abs(pmt_hit-ipmt)<2) phi_cut1=true;  else phi_cut1=false;

  //if (iseg.eq.8) then phi_cut1=.true.   else   phi_cut1=(abs(pmt_hit-ipmt).lt.2)
  // if(iseg==8) phi_cut1=true;  else if(abs(pmt_hit-ipmt)<2) phi_cut1=true;  else phi_cut1=false;
  //if((iseg==8) && abs(pmt_hit-ipmt)<3) phi_cut1=true;  else if(abs(pmt_hit-ipmt)<2) phi_cut1=true;  else phi_cut1=false;
  if(abs(pmt_hit-ipmt)<2) phi_cut1=true;  else phi_cut1=false;//11/22/15 (<3 replaced by <2 for iseg==8)
  /*after suggestions during EG4 meeting a week ago based on 
    https://www.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/Pass2/Cuts/Osi/Plots/dpmtIdInDiffCCsegms2.gif */

  // phi_cut1=true;
#ifdef DEBUG_HISTOS_PASS2
  if(isec==6) h1Dipmt[iseg-1]->Fill(pmt_hit-ipmt);
#endif  
#ifdef THP_CUT_ADJUST_M3P6
#endif
#ifdef PRINT_SOME_INFO
  cout<<"debug3= "<<isec<<endl;
#endif


  //c Azimuthal angle cut 2 (note: this is the only one using phi_tl1)
  if (fabs(dphi_sec)>30.) 
    phi_cut2=false;       //! outside the sector
  else  phi_cut2=true;

	

  //c Azimuthal angle cut 3
  if (isec==6)
    {//! phi cuts 3 and 4 only for sector 6    /////***********
      //if (fabs(dphi_p)>10.0) phi_cut3=false;   //Reviving back in 9/10/15 (to make agree with Lamiaa's #)
      //if (fabs(dphi_p)>10.0) phi_cut3=true;    //12/21/12 (removing this cut from OsiCut, Fiducial cuts will take care of it)
      if (fabs(dphi_p)>13.0) phi_cut3=false;  //SEK recommendation (11/16/15 EG4 meeting)
      else phi_cut3=true;
      
      
      if (fabs(dphi_p)>16.0)  // ! cut 4 (a looser cut than #3, not used)
	phi_cut4=false;
      else
	phi_cut4=true;
    }
  else  //! sectors 1-5                                         /////****************
    {
      phi_cut3=true;
      phi_cut4=true;
    }                                                     /////****************

  if(phi_cut1==true && phi_cut2==true && phi_cut3 == true)
    //if( phi_cut2==true && phi_cut3 == true)  //if( phi_cut1==true && phi_cut3 == true)
    //if( phi_cut1==true && phi_cut2 == true)  //if( phi_cut1==true)
    phi_cut=true;
  



				       
  //c Timing SC-CC
   if (isec<6)
    {
      if(dt_ScCc>dt_cut[isec-1][iseg-1]) dt_sc_cc_cut=true;
      else dt_sc_cc_cut=false;
    }
  else if (isec==6)
    {
      if(dt_ScCc>-6.0) dt_sc_cc_cut=true;
      else dt_sc_cc_cut=false;
    }
  
  
  //c Timing EC-CC  ! only for sector 6 //Time difference more than -6.0 nano sec.
   if (isec<6)
    dt_ec_cc_cut=true;
  else 
    if (isec==6)
    {
      if(dt_EcCc>-6.0)  dt_ec_cc_cut=true;
      else dt_ec_cc_cut=false;
    }
  
  
   //Debug
   //theta_p_cut = true; phi_cut1 = true; phi_cut2 = true; phi_cut3 = true; dt_sc_cc_cut =true; dt_ec_cc_cut = true;
   //theta_p_cut = true; phi_cut1 = true; 
   //phi_cut2 = true; //The loss in the 1st row of NH3n1054gof01_03.gif seems to be from phi_cut2
   //phi_cut3 = true; dt_sc_cc_cut =true; dt_ec_cc_cut = true;
   //Debug

   //if(phi_cut2==false) cout<<int(dc_sect)<<" "<<dphi_sec<<" "<<ph_El<<endl;

   //if(dt_sc_cc_cut==true)  return ccmatch = true;  //Debug line 5/15/2012


  //c Combine all cuts				   
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true  && phi_cut3==true && dt_sc_cc_cut==true && dt_ec_cc_cut==true
     && dtScCcUprCt==true && dtEcCcUprCt==true) //last 2 added on 11/29/15
    //if(theta_p_cut==true && phi_cut==true && dt_sc_cc_cut==true && dt_ec_cc_cut==true) 
    /*
   return ccmatch=true;
  else return ccmatch=false;
    */
    ccmatch = true;
  else ccmatch = false;




 #ifdef PH_CUT2_STUDY
  float phC_El = atan2(ccy,ccx)*rad2deg; if(phC_El<0.0) phC_El = phC_El + 360.0;
  if(ccmatch==true) h2phC_phDc1_all->Fill(ph_El, phC_El);
  else if(theta_p_cut==true && phi_cut1==true && phi_cut2==false  && phi_cut3==true && dt_sc_cc_cut==true && dt_ec_cc_cut==true)
    {
      h2phC_phDc1_rej->Fill(ph_El, phC_El);
      h1secDC_rej->Fill(dc_sect);
      h1secCC_rej->Fill(cc_sect);
      h1p_rej->Fill(pp);
      h1phP_rej->Fill(dc_ccph); h1thP_rej->Fill(dc_ccth);
      h1nphe_rej->Fill(Nphe*1.0);
      Float_t W = sqrt(Mt_1H*Mt_1H + 2*Mt_1H*(Eb - pp) - 4*Eb*pp*pow(sin((acos(ccz))/2),2));
      h1W_rej->Fill(W);

      float c_l = 30; //speed of linght in cm/ns
      float dt = sc_t - Tr_time; float beta=sc_r/(c_l*dt);
      h2betaVp_rej->Fill(pp,beta);
    }
#endif
  


#ifdef DEBUG_HISTOS_PASS2
  if(theta_p_cut==true) 
    { h1nphAr[0][4]->Fill(Nphe);//h1allthP_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][4]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][4]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true) 
    { h1nphAr[0][5]->Fill(Nphe);//h1allph1_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][5]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][5]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true) 
    { h1nphAr[0][6]->Fill(Nphe);//h1allph2_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][6]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][6]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true) 
    { h1nphAr[0][7]->Fill(Nphe);//h1allph3_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][7]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][7]->Fill(Nphe);}
  //if(phi_cut4==true) //till 12/5/15
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true && phi_cut4==true) //From 12/5/15
    { h1nphAr[0][8]->Fill(Nphe);//  h1allph4_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][8]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][8]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true && 
     dt_sc_cc_cut==true) 
    { h1nphAr[0][9]->Fill(Nphe);//h1alldtScCc_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][9]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][9]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true && dt_sc_cc_cut==true 
     && dt_ec_cc_cut==true) 
    { h1nphAr[0][10]->Fill(Nphe);//h1alldtEcCc_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][10]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][10]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true && dt_sc_cc_cut==true 
     && dt_ec_cc_cut==true && dtScCcUprCt==true) 
    { h1nphAr[0][11]->Fill(Nphe);//h1alldtScCcUprCt->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][11]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][11]->Fill(Nphe);}
  if(theta_p_cut==true && phi_cut1==true && phi_cut2==true && phi_cut3==true && dt_sc_cc_cut==true 
     && dt_ec_cc_cut==true && dtScCcUprCt==true && dtEcCcUprCt==true) 
    { h1nphAr[0][12]->Fill(Nphe);//h1alldtEcCcUprCt->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][12]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][12]->Fill(Nphe);}
  if(ccmatch==true) 
    { h1nphAr[0][3]->Fill(Nphe); //h1allosi_dbg->Fill(Nphe); 
      if(ecCt==true) h1nphAr[1][3]->Fill(Nphe); else if(ecCt==false) h1nphAr[2][3]->Fill(Nphe);}   


  //Using -ve osipenko-component cuts on pi- candidates (selected using ecCt==false condition)
  else if(ecCt==false) //12/7/15: Although defined & filled some already, I stopped plotting them, may do so later
    {
      if(ccmatch==false) h1nphAr[3][4]->Fill(Nphe);
    }
#endif

  return ccmatch;
}



/*
//==========================================
//these subroutines are to change the physics angles into the CC projected angles.

void sphere_inter(int isec, float xo, float yo, float zo, float dx, float dy, float dz, float &ar, float &theta, float &phi, int &ires)
  // void fcn(float p, int q, float &Sum, float &Prod, float &vSin)

//The ampersand here is used to pass the value by reference instead of by value. And, in this particular case, both ways seem to work.

{
  float rad2deg = 57.2957795;
  float point[3] = {xo, yo, zo}; //1D array variable of size three defined as well as assigned values in terms of other variables.
  //cc_pln6 is for the new cherenkov counter of EG4 which was installed only in the sixth sector.
  float cc_pln6[3] = {-0.0007840784063, 0.0, -0.0017805641}; //Distance is 509.0002 cm
  float cc_pln[3] = {-0.0007754463505, 0.0, -0.001662950067};//Distance is 539.0002cm
  
  float dir[3] = {dx, dy, dz};
  float cross[3] = {0,0,0};

  if (isec<=5)
    CCVCRPL(point, dir, cc_pln, ires, ar, cross);
  else
    CCVCRPL(point, dir, cc_pln6, ires, ar, cross);
  // kp: //For goto use; for flow control; tried to 'stop' the function call if the condition tested fails.
  ar = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]); 
  if (ar!=0) 
    {
      theta = rad2deg*acos(cross[2]/ar);
      // phi = rad2deg*ata2((cross[0]/ar, (cross[1]/ar)) - 90;
      //phi = rad2deg*ata2((cross[1]/ar, (cross[0]/ar));//Alexandor Vlassov: no need to divide by ar
      phi = rad2deg*atan2(cross[0], cross[1]); //This line will be superfluous
      phi = rad2deg*atan2(cross[1], cross[0]); 
      //  if (phi<0)    phi = phi + 360.0;	
    }
  else
    {
      theta = 0;			
      phi = 0;
    }
}
*/

 /*
void CCVCRPL(float R0[3], float DIR[3], float P[3], int &IRES, float &S, float  R[3])
{
 */
  /*
 Documentation for subroutine CCVCRPL(R0,DIR,P,IRES,S,R)
c
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : crossing of the stright line(R0,d)
C-                         with a plane
C-
C-   Inputs  :   R0(3) - initial point of line
C-               D (3) - vector direction: r = R0 + s*D
C-               P (3) - plane parameters:
C-               P(1)*x + P(2)*y + P(3)*z + 1 = 0
C-
C-   Outputs :   IRES =  0 - no cross with the plane.
C-                       1 - cross in positive directionsq
C-                      -1 - cross in negative direction
C-               S    =  Distance to the cross point
C-               R(3) =  Cross point coordinates.
C-
C-   Controls:
C-
C-   Created    14-APR-1994   Alexander V. Vlassov
C-   Modified   21-AUG-1995   Alexander V. Vlassov
C-
C----------------------------------------------------------------------
c
c_end_doc
   */

/*
  float un = 1.0, a = 0.0, b = 0.0, t, c = 0.0, D[3]; //Local user defined varieables.
  S = 0;
  float vsmall = 0.000001;

  //// if((DIR[0]==0)||(DIR[1]==0)||(DIR[2]==0))
  ////    cout<<"CCVCRPL warning: DIR = "<<DIR[0]<<" "<<DIR[1]<<" "<<DIR[2]<<endl;

     for (int i=0; i<3;i++)
       c = c + un * DIR[i] * DIR[i];
     
     if(c==0)
  {
    ////// cout<<"CCVCRPL error: c = 0, DIR = "<<DIR[0]<<" "<<DIR[1]<<" "<<DIR[2]<<endl;
     //   goto kp;
     // return;
     // return;
     R[0] = 0; R[1] = 0; R[2] = 0; return;
  }
     c = sqrt(c);
     for(int i =0; i<3; i++)
     D[i] = un * DIR[i]/c;
     
     for(int i =0; i<3; i++)
  {
    a = a + un * P[i]*D[i];
    b = b + un * P[i]*R0[i];
  }
     b = b + un;
     
     if(fabs(b)<= vsmall) // point on the plane
  {
    for(int i =0; i<3; i++)
      R[i] = R0[i];
    
    IRES = 1;
  }
     else
  {
    if(fabs(a)<=vsmall) 
      {
	for(int i =0; i<3; i++)
	  R[i] = 0;
	IRES = 0;
      }
    else
      {
	t = -b/a;
	for(int i =0; i<3; i++)
	  R[i] = t*D[i] + R0[i];
	IRES = 1;
	S = t;
	if(t<0)
	  IRES = -1;
      }
  }
}
*/
#endif
