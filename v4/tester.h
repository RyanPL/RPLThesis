//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 24 11:38:22 2017 by ROOT version 5.34/30
// from TTree h10/All_out
// found on file: ../root_51577_pass1.a16.root
//////////////////////////////////////////////////////////

#ifndef tester_h
#define tester_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class tester {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UShort_t        run_num;
   UInt_t          evntid;
   Char_t          evstat;
   UShort_t        evntype;
   Short_t         evntclas;
   UShort_t        l1bit;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time;
   UInt_t          hlsc;
   Int_t           intt;
   UChar_t         helicity;
   UChar_t         sync;
   UChar_t         mps;
   UChar_t         qrt;
   Float_t         hel_flag;
   Int_t          hel_ofl;
   Int_t          hel_onl;
   Int_t          hel_r26;
   Int_t           clock_u;
   Int_t           fcup_u;
   Int_t           slm_u;
   Int_t           clock_g;
   Int_t           fcup_g;
   Int_t           slm_g;
   UShort_t        raster_x;
   UShort_t        raster_y;
   UChar_t         npart;
   Int_t           gpart;
   Short_t         id[20];   //[gpart]
   Char_t          stat[20];   //[gpart]
   UChar_t         dc[20];   //[gpart]
   UChar_t         cc[20];   //[gpart]
   UChar_t         sc[20];   //[gpart]
   UChar_t         ec[20];   //[gpart]
   Float_t         p[20];   //[gpart]
   Float_t         m[20];   //[gpart]
   Char_t          q[20];   //[gpart]
   Float_t         b[20];   //[gpart]
   Float_t         cx[20];   //[gpart]
   Float_t         cy[20];   //[gpart]
   Float_t         cz[20];   //[gpart]
   Float_t         vx[20];   //[gpart]
   Float_t         vy[20];   //[gpart]
   Float_t         vz[20];   //[gpart]
   Int_t           dc_part;
   UChar_t         dc_sect[20];   //[dc_part]
   UChar_t         dc_trk[20];   //[dc_part]
   Char_t          dc_stat[20];   //[dc_part]
   UInt_t          tb_st[20];   //[dc_part]
   Float_t         dc_xsc[20];   //[dc_part]
   Float_t         dc_ysc[20];   //[dc_part]
   Float_t         dc_zsc[20];   //[dc_part]
   Float_t         dc_cxsc[20];   //[dc_part]
   Float_t         dc_cysc[20];   //[dc_part]
   Float_t         dc_czsc[20];   //[dc_part]
   Float_t         dc_vx[20];   //[dc_part]
   Float_t         dc_vy[20];   //[dc_part]
   Float_t         dc_vz[20];   //[dc_part]
   Float_t         dc_vr[20];   //[dc_part]
   Float_t         tl1_cx[20];   //[dc_part]
   Float_t         tl1_cy[20];   //[dc_part]
   Float_t         tl1_cz[20];   //[dc_part]
   Float_t         tl1_x[20];   //[dc_part]
   Float_t         tl1_y[20];   //[dc_part]
   Float_t         tl1_z[20];   //[dc_part]
   Float_t         tl1_r[20];   //[dc_part]
   Float_t         dc_c2[20];   //[dc_part]
   Float_t         dc_ccth[20];   //[dc_part]
   Float_t         dc_ccph[20];   //[dc_part]
   Int_t           ec_part;
   UShort_t        ec_stat[20];   //[ec_part]
   UChar_t         ec_sect[20];   //[ec_part]
   UChar_t         ec_whol[20];   //[ec_part]
   UInt_t          ec_inst[20];   //[ec_part]
   UInt_t          ec_oust[20];   //[ec_part]
   Float_t         etot[20];   //[ec_part]
   Float_t         ec_ei[20];   //[ec_part]
   Float_t         ec_eo[20];   //[ec_part]
   Float_t         ec_t[20];   //[ec_part]
   Float_t         ec_r[20];   //[ec_part]
   Float_t         ech_x[20];   //[ec_part]
   Float_t         ech_y[20];   //[ec_part]
   Float_t         ech_z[20];   //[ec_part]
   Float_t         ec_m2[20];   //[ec_part]
   Float_t         ec_m3[20];   //[ec_part]
   Float_t         ec_m4[20];   //[ec_part]
   Float_t         ec_c2[20];   //[ec_part]
   Int_t           sc_part;
   UChar_t         sc_sect[20];   //[sc_part]
   UChar_t         sc_hit[20];   //[sc_part]
   UChar_t         sc_pd[20];   //[sc_part]
   UChar_t         sc_stat[20];   //[sc_part]
   Float_t         edep[20];   //[sc_part]
   Float_t         sc_t[20];   //[sc_part]
   Float_t         sc_r[20];   //[sc_part]
   Float_t         sc_c2[20];   //[sc_part]
   Int_t           cc_part;
   UChar_t         cc_sect[20];   //[cc_part]
   UChar_t         cc_hit[20];   //[cc_part]
   UShort_t        cc_segm[20];   //[cc_part]
   UShort_t        nphe[20];   //[cc_part]
   Float_t         cc_t[20];   //[cc_part]
   Float_t         cc_r[20];   //[cc_part]
   Float_t         cc_c2[20];   //[cc_part]

   // List of branches
   TBranch        *b_run_num;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_l1bit;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time;   //!
   TBranch        *b_hlsc;   //!
   TBranch        *b_intt;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_sync;   //!
   TBranch        *b_mps;   //!
   TBranch        *b_qrt;   //!
   TBranch        *b_hel_flag;   //!
   TBranch        *b_hel_ofl;   //!
   TBranch        *b_hel_onl;   //!
   TBranch        *b_hel_r26;   //!
   TBranch        *b_clock_u;   //!
   TBranch        *b_fcup_u;   //!
   TBranch        *b_slm_u;   //!
   TBranch        *b_clock_g;   //!
   TBranch        *b_fcup_g;   //!
   TBranch        *b_slm_g;   //!
   TBranch        *b_raster_x;   //!
   TBranch        *b_raster_y;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_tb_st;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_vx;   //!
   TBranch        *b_dc_vy;   //!
   TBranch        *b_dc_vz;   //!
   TBranch        *b_dc_vr;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_tl1_cz;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl1_r;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_dc_ccth;   //!
   TBranch        *b_dc_ccph;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!

   tester(TTree *tree=0);
   virtual ~tester();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     n_compare();
   virtual void     n_emp(Char_t runNum[10]);
   virtual void     n_carb(Char_t runNum[10]);
   virtual void     n_nh3(Char_t runNum[10]);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tester_cxx
tester::tester(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //*                                               //ONLY comment out when testing makeHistos2
      if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../root_51577_pass1.a16.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../root_51577_pass1.a16.root");
      }
      f->GetObject("h10",tree);

   }
   Init(tree);//*/
}

tester::~tester()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tester::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tester::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tester::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run_num", &run_num, &b_run_num);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
   fChain->SetBranchAddress("evntype", &evntype, &b_evntype);
   fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
   fChain->SetBranchAddress("l1bit", &l1bit, &b_l1bit);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("rf_time", &rf_time, &b_rf_time);
   fChain->SetBranchAddress("hlsc", &hlsc, &b_hlsc);
   fChain->SetBranchAddress("intt", &intt, &b_intt);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("sync", &sync, &b_sync);
   fChain->SetBranchAddress("mps", &mps, &b_mps);
   fChain->SetBranchAddress("qrt", &qrt, &b_qrt);
   fChain->SetBranchAddress("hel_flag", &hel_flag, &b_hel_flag);
   fChain->SetBranchAddress("hel_ofl", &hel_ofl, &b_hel_ofl);
   fChain->SetBranchAddress("hel_onl", &hel_onl, &b_hel_onl);
   fChain->SetBranchAddress("hel_r26", &hel_r26, &b_hel_r26);
   fChain->SetBranchAddress("clock_u", &clock_u, &b_clock_u);
   fChain->SetBranchAddress("fcup_u", &fcup_u, &b_fcup_u);
   fChain->SetBranchAddress("slm_u", &slm_u, &b_slm_u);
   fChain->SetBranchAddress("clock_g", &clock_g, &b_clock_g);
   fChain->SetBranchAddress("fcup_g", &fcup_g, &b_fcup_g);
   fChain->SetBranchAddress("slm_g", &slm_g, &b_slm_g);
   fChain->SetBranchAddress("raster_x", &raster_x, &b_raster_x);
   fChain->SetBranchAddress("raster_y", &raster_y, &b_raster_y);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("b", b, &b_b);
   fChain->SetBranchAddress("cx", cx, &b_cx);
   fChain->SetBranchAddress("cy", cy, &b_cy);
   fChain->SetBranchAddress("cz", cz, &b_cz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("tb_st", tb_st, &b_tb_st);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
   fChain->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
   fChain->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
   fChain->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
   fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
   fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
   fChain->SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
   fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
   fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
   fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
   fChain->SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
   fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
   fChain->SetBranchAddress("dc_ccth", dc_ccth, &b_dc_ccth);
   fChain->SetBranchAddress("dc_ccph", dc_ccph, &b_dc_ccph);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
   fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
   fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
   fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
   fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
   fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
   fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("edep", edep, &b_edep);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);
   fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
   fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
   fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
   Notify();
}

Bool_t tester::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tester::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tester::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tester_cxx
