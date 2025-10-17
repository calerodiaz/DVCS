// *********************************************************************************************************
// Writes the .txt datafile of Hall-A into a .root file   																								 *
// There are 45 kinematic set for the unpolarized data   																									 *
// Every TTree entry corresponds to one of the kimenatic sets																							 *
// Hall-A data reference: https://hallaweb.jlab.org/experiment/DVCS/documents/results/Frederic_thesis.pdf  *
// *********************************************************************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "TList.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

//_________________________________________________________________________________________________________
void genytree(Double_t var = 0.05) {


	Double_t M = 0.938272; //Mass of the proton in GeV
	const Int_t NumOfPts = 24;

	// File to save the TTree
	TFile fout("dvcs_HallA_E12-06-114.root","recreate");

	TCanvas *c1;

	struct kin_t {
      Double_t k;
      Double_t QQ;
			Double_t xB;
			Double_t t;
   };
	kin_t kin;
	Double_t phi[NumOfPts];
	Double_t F[NumOfPts], stat_F[NumOfPts];
	Double_t sys_high_F[NumOfPts], sys_low_F[NumOfPts], sys_corr_F[NumOfPts];


	TTree *t3 = new TTree("unpolarized","dvcs");
	t3->Branch("kinematics",&kin.k,"k/D:QQ:xB:t");
	//t3->Branch("kinematics",&kin.QQ,"QQ/D:xB:t");
	t3->Branch("phi",phi,"phi[24]/D");
	t3->Branch("F",F,"F[24]/D");
	t3->Branch("stat_F",stat_F,"stat_F[24]/D");
	t3->Branch("sys_high_F",sys_high_F,"sys_high_F[24]/D");
	t3->Branch("sys_low_F",sys_low_F,"sys_low_F[24]/D");
	t3->Branch("sys_corr_F",sys_corr_F,"sys_corr_F[24]/D");

	ifstream in;
	Double_t  f_xB, f_QQ, f_t, f_phi, f_F, f_stat_F, f_sys_high_F, f_sys_low_F, f_e1, f_e2, f_e3, f_e4;
	in.open("outdvcs.txt");


//	Float_t x,y,z;
Int_t nlines = 0, ipt = 0, iset = 0;
Double_t gg, t_min;

Int_t ik = 0;
Double_t Eb[] = {7.383, 8.521, 10.591, 4.487, 8.851, 8.847, 10.992, 8.521, 10.591};

while (1) {
	 in >> f_xB >> f_QQ >> f_t >> f_phi >> f_F >> f_stat_F >> f_sys_high_F >> f_sys_low_F >> f_e1 >> f_e2 >> f_e3 >> f_e4;
	 if (!in.good()) break;

	 phi[ipt] = f_phi;
	 F[ipt] = f_F;
	 stat_F[ipt] = f_stat_F;
	 sys_high_F[ipt] = f_sys_high_F;
	 sys_low_F[ipt] = f_sys_low_F;
	 sys_corr_F[ipt] = 0.03 * f_F;

	 ipt++;

	 if ( nlines ==  ( NumOfPts * ( iset + 1 ) - 1 ) )  {

		kin.k =  Eb[ik];
		kin.QQ = f_QQ;
		kin.xB = f_xB;
		// Get t from t' = t - tmin
		gg = 4. * M * M * f_xB * f_xB / f_QQ;
		t_min = -f_QQ * ( 1. - sqrt( 1. + gg ) + gg / 2. ) / ( f_xB * ( 1. - sqrt( 1. + gg ) + gg / ( 2. * f_xB ) ) );
		kin.t = f_t + t_min;

		iset++;

		t3 ->Fill();
		printf("set=%d, ipt=%d, line=%d, k=%8f, xB=%8f, QQ=%8f, t=%8f, t'=%8f, tmin=%8f, phi=%8f, F=%8f, stat_F=%8f, sys_high_F=%8f, sys_low_F=%8f\n", iset, ipt, nlines, Eb[ik], f_xB, f_QQ, f_t + t_min, f_t, t_min, f_phi, f_F, f_stat_F, f_sys_high_F, f_sys_low_F);

		if(iset % 5 == 0) ik++;

		ipt = 0;
	 }

	 nlines++;
}
printf(" found %d points\n",nlines);

in.close();

	t3->Print();
	t3->Show(0);
	// fout.cd();
	t3->Write();

}
