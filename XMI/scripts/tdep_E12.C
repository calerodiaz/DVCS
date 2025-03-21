using namespace std;
using namespace ROOT;
using namespace ROOT::RDF;
//___________________________________________________________________________
auto km15Graphs(RNode df){

   vector<TGraphErrors*>gkm15(4);

   for (Int_t ipar = 0; ipar < 4; ipar++) {
            gkm15[ipar] = new TGraphErrors();
            gkm15[ipar]->SetName("true");
            gkm15[ipar]->SetLineColor(kGray+3);
            gkm15[ipar]->SetLineWidth(1);
            gkm15[ipar]->SetMarkerColor(kGray+3);
            gkm15[ipar]->SetMarkerSize(0.001);
   }
   auto fillkm15Graph = [&gkm15](double t, double ReH_km15, double ReE_km15, double ReHtilde_km15, double dvcs_km15){
      gkm15[0]->SetPoint(gkm15[0]->GetN(), -t, ReH_km15);
      gkm15[0]->SetPointError( gkm15[0]->GetN()-1, 0, 0 );
      gkm15[1]->SetPoint(gkm15[1]->GetN(), -t, ReE_km15);
      gkm15[1]->SetPointError( gkm15[1]->GetN()-1, 0, 0 );
      gkm15[2]->SetPoint(gkm15[2]->GetN(), -t, ReHtilde_km15);
      gkm15[2]->SetPointError( gkm15[2]->GetN()-1, 0, 0 );
      gkm15[3]->SetPoint(gkm15[3]->GetN(), -t, dvcs_km15);
      gkm15[3]->SetPointError( gkm15[3]->GetN()-1, 0, 0 );
   };
   df.Foreach(fillkm15Graph, {"t", "ReH_km15", "ReE_km15", "ReHt_km15", "dvcs_km15"});

   return gkm15;
}
//______________________________________________________________________________
auto xmiGraphs(RNode df){

   vector<TGraphMultiErrors*>gr(4);

   // Get kinematics
   auto QQ = df.Take<double>("QQ");
   auto xB = df.Take<double>("xB");
   auto t = df.Take<double>("t");
   // Get xmi results and errors (stats and syst)
   auto xmi_ReH = df.Take<double>("xmi_ReH");
   auto exmi_ReH = df.Take<double>("exmi_ReH");
   auto sys_ReH = df.Take<double>("sys_ReH");
   auto xmi_ReE = df.Take<double>("xmi_ReE");
   auto exmi_ReE = df.Take<double>("exmi_ReE");
   auto sys_ReE = df.Take<double>("sys_ReE");
   auto xmi_ReHt = df.Take<double>("xmi_ReHt");
   auto exmi_ReHt = df.Take<double>("exmi_ReHt");
   auto sys_ReHt = df.Take<double>("sys_ReHt");
   auto xmi_dvcs = df.Take<double>("xmi_dvcs");
   auto exmi_dvcs = df.Take<double>("exmi_dvcs");
   auto sys_dvcs = df.Take<double>("sys_dvcs");

   int npoints = 5;
   Double_t *eReH = new Double_t [2];   
   Double_t *eReE = new Double_t [2];  
   Double_t *eReHt = new Double_t [2]; 
   Double_t *edvcs = new Double_t [2]; 

   for (int ipar = 0; ipar<1; ipar++){
      gr[ipar] = new TGraphMultiErrors(npoints,2);
      gr[ipar]->SetMarkerStyle(20);
      gr[ipar]->SetMarkerColor(kAzure-7);
      gr[ipar]->SetLineColor(kAzure-7);
      gr[ipar]->GetAttLine(0)->SetLineColor(kAzure-7);
      gr[ipar]->GetAttLine(1)->SetLineColor(kAzure-7);
      //gr[ipar]->GetAttLine(1)->SetLineWidth(0);
      gr[ipar]->GetAttFill(1)->SetFillStyle(0);
      //gr[ipar]->GetAttFill(1)->SetFillColor(17);
      //gr[ipar]->GetAttFill(1)->SetFillStyle(1001);
      //gr[ipar]->GetAttFill(1)->SetFillColorAlpha(kBlack, 0.2);

      gr[ipar]->SetFillStyle(0);
      //gr[ipar]->SetFillColor(17);
      //gr[ipar]->SetFillColorAlpha(kBlack, 0.2);
   } 

   for (int i = 0; i<npoints; i++){
         
      eReH[0] = exmi_ReH->at(i);        
      eReH[1] = sys_ReH->at(i); 
      eReE[0] = exmi_ReE->at(i);        
      eReE[1] = sys_ReE->at(i); 
      eReHt[0] = exmi_ReHt->at(i);        
      eReHt[1] = sys_ReHt->at(i); 
      edvcs[0] = exmi_dvcs->at(i);        
      edvcs[1] = sys_dvcs->at(i);

      gr[0]->SetPoint(i, -1*t->at(i), xmi_ReH->at(i));
      gr[0]->SetPointError(i, 2, 0.0001, 0.0001, eReH, eReH);

      gr[1]->SetPoint(i, -1*t->at(i), xmi_ReE->at(i));
      gr[1]->SetPointError(i, 2, 0.0001, 0.0001, eReE, eReE);

      gr[2]->SetPoint(i, -1*t->at(i), xmi_ReHt->at(i));
      gr[2]->SetPointError(i, 2, 0.0001, 0.0001, eReHt, eReHt);

      gr[3]->SetPoint(i, -1*t->at(i), xmi_dvcs->at(i));
      gr[3]->SetPointError(i, 2, 0.0001, 0.0001, edvcs, edvcs);
   }      

   return gr;
}  
//______________________________________________________________________________
Int_t xBdep2(const Int_t initial_set = 1, const Int_t last_set = 85)
{
   //ROOT::EnableImplicitMT();
   TH2::AddDirectory(false); // No dot write histogram automaticaly to the root file

   cout<<"Reading csv into dfs...";
   // fit results df
   auto df0_hA = FromCSV("/media/lily/Data/GPDs/HallA/realdata/maps/results_HallA_sets_1to85_data_5pct.csv");
   cout<<"done!"<<endl;
   TMultiGraph *mgr_results[4];
   Int_t fMarker[] = {22, 20, 20, 20, 20, 21, 23, 22, 34};
   Color_t col[] = {kBlack, kGreen+2, kBlue, kOrange+2, kMagenta+1, kRed-2, kAzure+9, kGray+2, kYellow+1};
   TString parsname[] = {"#displaystyle #mathfrak{Re}#mathcal{H}", "#displaystyle #mathfrak{Re}#mathcal{E}", "#displaystyle #mathfrak{Re}#mathcal{#widetilde{H}}",
                     "#displaystyle d^{4}#sigma^{DVCS}_{UU} [nb/GeV^{4}]"};
   TString kFormulation = "BKM10";
   TString kModel = "KM15";

   for (Int_t ipar = 0; ipar < 4; ipar++){
      mgr_results[ipar] = new TMultiGraph();
      mgr_results[ipar] ->SetTitle(Form("; #displaystyle -t[GeV^{2}]; %s", parsname[ipar].Data()));
   }   

   auto df_hA = df0_hA.Filter("exmi_ReH>0");
   // Filter desired sets
   auto setrange = [&initial_set, &last_set](Long64_t set) { return set >= initial_set && set <= last_set; };
   //auto df1_hA = df_hA.Filter(setrange, {"set"});
   auto df_hA_xB1_E00 = df_hA.Filter({"set <= 20"});
   auto df_hA_xB1_E12 = df_hA.Filter({"set > 20 && set <= 35"});
   auto df_hA_xB1_E07 = df_hA.Filter({"set > 65 && set <= 85"});
   auto df_hA_xB2_E12 = df_hA.Filter({"set > 35 && set <= 55"});
   auto df_hA_xB3_E12 = df_hA.Filter({"set > 55 && set <= 65"});

   // mean xB
   auto xB1_E00 = df_hA_xB1_E00.Mean("xB");
   auto xB1_E12 = df_hA_xB1_E12.Mean("xB");
   auto xB1_E07 = df_hA_xB1_E07.Mean("xB");
   auto xB2_E12 = df_hA_xB2_E12.Mean("xB");
   auto xB3_E12 = df_hA_xB3_E12.Mean("xB");

   vector<TGraphMultiErrors*>gE00_xB(4);
   vector<TGraphMultiErrors*>gE07_xB(4);
   vector<TGraphMultiErrors*>gE12_xB(4);
   // ----------------------- ReH ---------------------------------------
   // Mean ReH
   auto ReH_hA_xB1_E00 = df_hA_xB1_E00.Mean("xmi_ReH");
   auto ReH_hA_xB1_E12 = df_hA_xB1_E12.Mean("xmi_ReH");
   auto ReH_hA_xB1_E07 = df_hA_xB1_E07.Mean("xmi_ReH");
   auto ReH_hA_xB2_E12 = df_hA_xB2_E12.Mean("xmi_ReH");
   auto ReH_hA_xB3_E12 = df_hA_xB3_E12.Mean("xmi_ReH");
   // Mean sigma ReH
   auto eReH_hA_xB1_E00 = df_hA_xB1_E00.Mean("exmi_ReH");
   auto eReH_hA_xB1_E12 = df_hA_xB1_E12.Mean("exmi_ReH");
   auto eReH_hA_xB1_E07 = df_hA_xB1_E07.Mean("exmi_ReH");
   auto eReH_hA_xB2_E12 = df_hA_xB2_E12.Mean("exmi_ReH");
   auto eReH_hA_xB3_E12 = df_hA_xB3_E12.Mean("exmi_ReH");
   // Mean syst ReH
   auto sys_ReH_hA_xB1_E00 = df_hA_xB1_E00.Mean("sys_ReH");
   auto sys_ReH_hA_xB1_E12 = df_hA_xB1_E12.Mean("sys_ReH");
   auto sys_ReH_hA_xB1_E07 = df_hA_xB1_E07.Mean("sys_ReH");
   auto sys_ReH_hA_xB2_E12 = df_hA_xB2_E12.Mean("sys_ReH");
   auto sys_ReH_hA_xB3_E12 = df_hA_xB3_E12.Mean("sys_ReH");

   Double_t ReH_xmi_E00[1] = {*ReH_hA_xB1_E00};
   Double_t eReH_xmi_E00[1] = {*eReH_hA_xB1_E00};
   Double_t sys_ReH_xmi_E00[1] = {*sys_ReH_hA_xB1_E00};
   Double_t zero[1] = {0.001};
   Double_t xB_E00[1] = {*xB1_E00-0.01};

   Double_t ReH_xmi_E07[1] = {*ReH_hA_xB1_E07};
   Double_t eReH_xmi_E07[1] = {*eReH_hA_xB1_E07};
   Double_t sys_ReH_xmi_E07[1] = {*sys_ReH_hA_xB1_E07};
   Double_t xB_E07[1] = {*xB1_E07+0.02};

   Double_t ReH_xmi_E12[3] = {*ReH_hA_xB1_E12, *ReH_hA_xB2_E12, *ReH_hA_xB3_E12};
   Double_t eReH_xmi_E12[3] = {*eReH_hA_xB1_E12, *eReH_hA_xB2_E12, *eReH_hA_xB3_E12};
   Double_t sys_ReH_xmi_E12[3] = {*sys_ReH_hA_xB1_E12, *sys_ReH_hA_xB2_E12, *sys_ReH_hA_xB3_E12};
   Double_t xB_E12[3] = {*xB1_E12, *xB2_E12, *xB3_E12};
   Double_t zero3[3] = {0.001, 0.001, 0.001};
  
   gE00_xB[0] = new TGraphMultiErrors("gE00_xB_ReH", "gE00_xB_ReH", 1, xB_E00, ReH_xmi_E00, zero, zero, eReH_xmi_E00, eReH_xmi_E00);
   gE07_xB[0] = new TGraphMultiErrors("gE07_xB_ReH", "gE07_xB_ReH", 1, xB_E07, ReH_xmi_E07, zero, zero, eReH_xmi_E07, eReH_xmi_E07);
   gE12_xB[0] = new TGraphMultiErrors("gE12_xB_ReH", "gE12_xB_ReH", 3, xB_E12, ReH_xmi_E12, zero3, zero3, eReH_xmi_E12, eReH_xmi_E12);
   gE00_xB[0]->AddYError(1, sys_ReH_xmi_E00, sys_ReH_xmi_E00);
   gE07_xB[0]->AddYError(1, sys_ReH_xmi_E07, sys_ReH_xmi_E07);
   gE12_xB[0]->AddYError(3, sys_ReH_xmi_E12, sys_ReH_xmi_E12);

   // ----------------------- ReE ---------------------------------------
   // Mean ReE
   auto ReE_hA_xB1_E00 = df_hA_xB1_E00.Mean("xmi_ReE");
   auto ReE_hA_xB1_E12 = df_hA_xB1_E12.Mean("xmi_ReE");
   auto ReE_hA_xB1_E07 = df_hA_xB1_E07.Mean("xmi_ReE");
   auto ReE_hA_xB2_E12 = df_hA_xB2_E12.Mean("xmi_ReE");
   auto ReE_hA_xB3_E12 = df_hA_xB3_E12.Mean("xmi_ReE");
   // Mean sigma ReE
   auto eReE_hA_xB1_E00 = df_hA_xB1_E00.Mean("exmi_ReE");
   auto eReE_hA_xB1_E12 = df_hA_xB1_E12.Mean("exmi_ReE");
   auto eReE_hA_xB1_E07 = df_hA_xB1_E07.Mean("exmi_ReE");
   auto eReE_hA_xB2_E12 = df_hA_xB2_E12.Mean("exmi_ReE");
   auto eReE_hA_xB3_E12 = df_hA_xB3_E12.Mean("exmi_ReE");
   // Mean syst ReE
   auto sys_ReE_hA_xB1_E00 = df_hA_xB1_E00.Mean("sys_ReE");
   auto sys_ReE_hA_xB1_E12 = df_hA_xB1_E12.Mean("sys_ReE");
   auto sys_ReE_hA_xB1_E07 = df_hA_xB1_E07.Mean("sys_ReE");
   auto sys_ReE_hA_xB2_E12 = df_hA_xB2_E12.Mean("sys_ReE");
   auto sys_ReE_hA_xB3_E12 = df_hA_xB3_E12.Mean("sys_ReE");

   Double_t ReE_xmi_E00[1] = {*ReE_hA_xB1_E00};
   Double_t eReE_xmi_E00[1] = {*eReE_hA_xB1_E00};
   Double_t sys_ReE_xmi_E00[1] = {*sys_ReE_hA_xB1_E00};
   //Double_t zero[1] = {0.001};
   //Double_t xB_E00[1] = {*xB1_E00};

   Double_t ReE_xmi_E07[1] = {*ReE_hA_xB1_E07};
   Double_t eReE_xmi_E07[1] = {*eReE_hA_xB1_E07};
   Double_t sys_ReE_xmi_E07[1] = {*sys_ReE_hA_xB1_E07};
   //Double_t xB_E07[1] = {*xB1_E07};

   Double_t ReE_xmi_E12[3] = {*ReE_hA_xB1_E12, *ReE_hA_xB2_E12, *ReE_hA_xB3_E12};
   Double_t eReE_xmi_E12[3] = {*eReE_hA_xB1_E12, *eReE_hA_xB2_E12, *eReE_hA_xB3_E12};
   Double_t sys_ReE_xmi_E12[3] = {*sys_ReE_hA_xB1_E12, *sys_ReE_hA_xB2_E12, *sys_ReE_hA_xB3_E12};
   //Double_t xB_E12[3] = {*xB1_E12, *xB2_E12, *xB3_E12};
   //Double_t zero3[3] = {0.001, 0.001, 0.001};
  
   gE00_xB[1] = new TGraphMultiErrors("gE00_xB_ReE", "gE00_xB_ReE", 1, xB_E00, ReE_xmi_E00, zero, zero, eReE_xmi_E00, eReE_xmi_E00);
   gE07_xB[1] = new TGraphMultiErrors("gE07_xB_ReE", "gE07_xB_ReE", 1, xB_E07, ReE_xmi_E07, zero, zero, eReE_xmi_E07, eReE_xmi_E07);
   gE12_xB[1] = new TGraphMultiErrors("gE12_xB_ReE", "gE12_xB_ReE", 3, xB_E12, ReE_xmi_E12, zero3, zero3, eReE_xmi_E12, eReE_xmi_E12);
   gE00_xB[1]->AddYError(1, sys_ReE_xmi_E00, sys_ReE_xmi_E00);
   gE07_xB[1]->AddYError(1, sys_ReE_xmi_E07, sys_ReE_xmi_E07);
   gE12_xB[1]->AddYError(3, sys_ReE_xmi_E12, sys_ReE_xmi_E12);

   // ----------------------- ReHt ---------------------------------------
   // Mean ReHt
   auto ReHt_hA_xB1_E00 = df_hA_xB1_E00.Mean("xmi_ReHt");
   auto ReHt_hA_xB1_E12 = df_hA_xB1_E12.Mean("xmi_ReHt");
   auto ReHt_hA_xB1_E07 = df_hA_xB1_E07.Mean("xmi_ReHt");
   auto ReHt_hA_xB2_E12 = df_hA_xB2_E12.Mean("xmi_ReHt");
   auto ReHt_hA_xB3_E12 = df_hA_xB3_E12.Mean("xmi_ReHt");
   // Mean sigma ReHt
   auto eReHt_hA_xB1_E00 = df_hA_xB1_E00.Mean("exmi_ReHt");
   auto eReHt_hA_xB1_E12 = df_hA_xB1_E12.Mean("exmi_ReHt");
   auto eReHt_hA_xB1_E07 = df_hA_xB1_E07.Mean("exmi_ReHt");
   auto eReHt_hA_xB2_E12 = df_hA_xB2_E12.Mean("exmi_ReHt");
   auto eReHt_hA_xB3_E12 = df_hA_xB3_E12.Mean("exmi_ReHt");
   // Mean syst ReHt
   auto sys_ReHt_hA_xB1_E00 = df_hA_xB1_E00.Mean("sys_ReHt");
   auto sys_ReHt_hA_xB1_E12 = df_hA_xB1_E12.Mean("sys_ReHt");
   auto sys_ReHt_hA_xB1_E07 = df_hA_xB1_E07.Mean("sys_ReHt");
   auto sys_ReHt_hA_xB2_E12 = df_hA_xB2_E12.Mean("sys_ReHt");
   auto sys_ReHt_hA_xB3_E12 = df_hA_xB3_E12.Mean("sys_ReHt");

   Double_t ReHt_xmi_E00[1] = {*ReHt_hA_xB1_E00};
   Double_t eReHt_xmi_E00[1] = {*eReHt_hA_xB1_E00};
   Double_t sys_ReHt_xmi_E00[1] = {*sys_ReHt_hA_xB1_E00};
   //Double_t zero[1] = {0.001};
   //Double_t xB_E00[1] = {*xB1_E00};

   Double_t ReHt_xmi_E07[1] = {*ReHt_hA_xB1_E07};
   Double_t eReHt_xmi_E07[1] = {*eReHt_hA_xB1_E07};
   Double_t sys_ReHt_xmi_E07[1] = {*sys_ReHt_hA_xB1_E07};
   //Double_t xB_E07[1] = {*xB1_E07};

   Double_t ReHt_xmi_E12[3] = {*ReHt_hA_xB1_E12, *ReHt_hA_xB2_E12, *ReHt_hA_xB3_E12};
   Double_t eReHt_xmi_E12[3] = {*eReHt_hA_xB1_E12, *eReHt_hA_xB2_E12, *eReHt_hA_xB3_E12};
   Double_t sys_ReHt_xmi_E12[3] = {*sys_ReHt_hA_xB1_E12, *sys_ReHt_hA_xB2_E12, *sys_ReHt_hA_xB3_E12};
   //Double_t xB_E12[3] = {*xB1_E12, *xB2_E12, *xB3_E12};
   //Double_t zero3[3] = {0.001, 0.001, 0.001};
  
   gE00_xB[2] = new TGraphMultiErrors("gE00_xB_ReHt", "gE00_xB_ReHt", 1, xB_E00, ReHt_xmi_E00, zero, zero, eReHt_xmi_E00, eReHt_xmi_E00);
   gE07_xB[2] = new TGraphMultiErrors("gE07_xB_ReHt", "gE07_xB_ReHt", 1, xB_E07, ReHt_xmi_E07, zero, zero, eReHt_xmi_E07, eReHt_xmi_E07);
   gE12_xB[2] = new TGraphMultiErrors("gE12_xB_ReHt", "gE12_xB_ReHt", 3, xB_E12, ReHt_xmi_E12, zero3, zero3, eReHt_xmi_E12, eReHt_xmi_E12);
   gE00_xB[2]->AddYError(1, sys_ReHt_xmi_E00, sys_ReHt_xmi_E00);
   gE07_xB[2]->AddYError(1, sys_ReHt_xmi_E07, sys_ReHt_xmi_E07);
   gE12_xB[2]->AddYError(3, sys_ReHt_xmi_E12, sys_ReHt_xmi_E12);

   // ----------------------- dvcs ---------------------------------------
   // Mean dvcs
   auto dvcs_hA_xB1_E00 = df_hA_xB1_E00.Mean("xmi_dvcs");
   auto dvcs_hA_xB1_E12 = df_hA_xB1_E12.Mean("xmi_dvcs");
   auto dvcs_hA_xB1_E07 = df_hA_xB1_E07.Mean("xmi_dvcs");
   auto dvcs_hA_xB2_E12 = df_hA_xB2_E12.Mean("xmi_dvcs");
   auto dvcs_hA_xB3_E12 = df_hA_xB3_E12.Mean("xmi_dvcs");
   // Mean sigma dvcs
   auto edvcs_hA_xB1_E00 = df_hA_xB1_E00.Mean("exmi_dvcs");
   auto edvcs_hA_xB1_E12 = df_hA_xB1_E12.Mean("exmi_dvcs");
   auto edvcs_hA_xB1_E07 = df_hA_xB1_E07.Mean("exmi_dvcs");
   auto edvcs_hA_xB2_E12 = df_hA_xB2_E12.Mean("exmi_dvcs");
   auto edvcs_hA_xB3_E12 = df_hA_xB3_E12.Mean("exmi_dvcs");
   // Mean syst dvcs
   auto sys_dvcs_hA_xB1_E00 = df_hA_xB1_E00.Mean("sys_dvcs");
   auto sys_dvcs_hA_xB1_E12 = df_hA_xB1_E12.Mean("sys_dvcs");
   auto sys_dvcs_hA_xB1_E07 = df_hA_xB1_E07.Mean("sys_dvcs");
   auto sys_dvcs_hA_xB2_E12 = df_hA_xB2_E12.Mean("sys_dvcs");
   auto sys_dvcs_hA_xB3_E12 = df_hA_xB3_E12.Mean("sys_dvcs");

   Double_t dvcs_xmi_E00[1] = {*dvcs_hA_xB1_E00};
   Double_t edvcs_xmi_E00[1] = {*edvcs_hA_xB1_E00};
   Double_t sys_dvcs_xmi_E00[1] = {*sys_dvcs_hA_xB1_E00};
   //Double_t zero[1] = {0.001};
   //Double_t xB_E00[1] = {*xB1_E00};

   Double_t dvcs_xmi_E07[1] = {*dvcs_hA_xB1_E07};
   Double_t edvcs_xmi_E07[1] = {*edvcs_hA_xB1_E07};
   Double_t sys_dvcs_xmi_E07[1] = {*sys_dvcs_hA_xB1_E07};
   //Double_t xB_E07[1] = {*xB1_E07};

   Double_t dvcs_xmi_E12[3] = {*dvcs_hA_xB1_E12, *dvcs_hA_xB2_E12, *dvcs_hA_xB3_E12};
   Double_t edvcs_xmi_E12[3] = {*edvcs_hA_xB1_E12, *edvcs_hA_xB2_E12, *edvcs_hA_xB3_E12};
   Double_t sys_dvcs_xmi_E12[3] = {*sys_dvcs_hA_xB1_E12, *sys_dvcs_hA_xB2_E12, *sys_dvcs_hA_xB3_E12};
   //Double_t xB_E12[3] = {*xB1_E12, *xB2_E12, *xB3_E12};
   //Double_t zero3[3] = {0.001, 0.001, 0.001};
  
   gE00_xB[3] = new TGraphMultiErrors("gE00_xB_dvcs", "gE00_xB_dvcs", 1, xB_E00, dvcs_xmi_E00, zero, zero, edvcs_xmi_E00, edvcs_xmi_E00);
   gE07_xB[3] = new TGraphMultiErrors("gE07_xB_dvcs", "gE07_xB_dvcs", 1, xB_E07, dvcs_xmi_E07, zero, zero, edvcs_xmi_E07, edvcs_xmi_E07);
   gE12_xB[3] = new TGraphMultiErrors("gE12_xB_dvcs", "gE12_xB_dvcs", 3, xB_E12, dvcs_xmi_E12, zero3, zero3, edvcs_xmi_E12, edvcs_xmi_E12);
   gE00_xB[3]->AddYError(1, sys_dvcs_xmi_E00, sys_dvcs_xmi_E00);
   gE07_xB[3]->AddYError(1, sys_dvcs_xmi_E07, sys_dvcs_xmi_E07);
   gE12_xB[3]->AddYError(3, sys_dvcs_xmi_E12, sys_dvcs_xmi_E12);

   
   // ----------------------- KM15 ---------------------------------------
   // Mean ReH_km15   
   auto ReH_km15_xB1_E12 = df_hA_xB1_E12.Mean("ReH_km15");
   auto ReH_km15_xB2_E12 = df_hA_xB2_E12.Mean("ReH_km15");
   auto ReH_km15_xB3_E12 = df_hA_xB3_E12.Mean("ReH_km15");
   // Mean ReE_km15   
   auto ReE_km15_xB1_E12 = df_hA_xB1_E12.Mean("ReE_km15");
   auto ReE_km15_xB2_E12 = df_hA_xB2_E12.Mean("ReE_km15");
   auto ReE_km15_xB3_E12 = df_hA_xB3_E12.Mean("ReE_km15");
  // Mean ReH_km15   
   auto ReHt_km15_xB1_E12 = df_hA_xB1_E12.Mean("ReHt_km15");
   auto ReHt_km15_xB2_E12 = df_hA_xB2_E12.Mean("ReHt_km15");
   auto ReHt_km15_xB3_E12 = df_hA_xB3_E12.Mean("ReHt_km15");

   Double_t ReH_km15_E12[3] = {*ReH_km15_xB1_E12, *ReH_km15_xB2_E12, *ReH_km15_xB3_E12};
   Double_t ReE_km15_E12[3] = {*ReE_km15_xB1_E12, *ReE_km15_xB2_E12, *ReE_km15_xB3_E12};
   Double_t ReHt_km15_E12[3] = {*ReHt_km15_xB1_E12, *ReHt_km15_xB2_E12, *ReHt_km15_xB3_E12};
   
   vector<TGraph*>gkm15(4);

   gkm15[0] = new TGraph(3, xB_E12, ReH_km15_E12);
   gkm15[1] = new TGraph(3, xB_E12, ReE_km15_E12);
   gkm15[2] = new TGraph(3, xB_E12, ReHt_km15_E12);

   for (Int_t ipar = 0; ipar < 3; ipar++) {          
            gkm15[ipar]->SetName("true");
            gkm15[ipar]->SetLineColor(1);
            gkm15[ipar]->SetLineWidth(1);
            gkm15[ipar]->SetMarkerColor(1);
            gkm15[ipar]->SetMarkerSize(0.001);
   }

   // --------------- E12 Paper results F. Georges et. al.---------------------------------------------
   vector<TGraphMultiErrors*>gE00_xB_FG(4);
   vector<TGraphMultiErrors*>gE12_xB_FG(4);

   // ReH
   Double_t xB_E00_FG[1] = {0.354-0.0055};
   Double_t ReH_E00_FG[1] = {-1.929};
   Double_t eReH_E00_FG[1] = {-1.6475 + 1.929};
   Double_t sys_ReH_E00_FG[1] = {-1.7207 + 1.929};

   Double_t xB_E12_FG[3] = {0.3649-0.003, 0.4897-0.005, 0.6098-0.005};
   Double_t ReH_E12_FG[3] = {-1.3990, -1.7977, -1.1036};
   Double_t eReH_E12_FG[3] = {-0.9202 + 1.399, -1.0822 + 1.7977,  -0.04454 + 1.1036};
   Double_t sys_ReH_E12_FG[3] = {-1.134 + 1.399, -1.239 + 1.7977, -0.760 + 1.1036};

   // ReE
   Double_t ReE_E00_FG[1] = {-4.367};
   Double_t eReE_E00_FG[1] = {-2.411+4.367};
   Double_t sys_ReE_E00_FG[1] = { -2.7683+4.367};

   Double_t ReE_E12_FG[3] = {0.4791, 1.0807, -1.412};
   Double_t eReE_E12_FG[3] = {2.9451 - 0.4791, 4.8902 - 1.0807, 1.5127 + 1.412};
   Double_t sys_ReE_E12_FG[3] = {1.737 - 0.4791, 2.747 - 1.0807, -0.1368 + 1.412};
  
   // ReHt
   Double_t ReHt_E00_FG[1] = {2.0387};
   Double_t eReHt_E00_FG[1] = {2.5622 - 2.0387};
   Double_t sys_ReHt_E00_FG[1] = {2.4616 - 2.0387};

   Double_t ReHt_E12_FG[3] = {0.5156, 2.0701, 0.4232};
   Double_t eReHt_E12_FG[3] = {1.4619 - 0.5156, 3.2446 - 2.0701, 1.7117 - 0.4232};
   Double_t sys_ReHt_E12_FG[3] = {1.2338 - 0.5156, 2.6272 - 2.0701, 1.2151 - 0.4232};
 

   gE00_xB_FG[0] = new TGraphMultiErrors("gE00_xB_ReH", "gE00_xB_ReH", 1, xB_E00_FG, ReH_E00_FG, zero, zero, eReH_E00_FG, eReH_E00_FG);
   gE12_xB_FG[0] = new TGraphMultiErrors("gE12_xB_ReH", "gE12_xB_ReH", 3, xB_E12_FG, ReH_E12_FG, zero3, zero3, eReH_E12_FG, eReH_E12_FG);
   gE00_xB_FG[0]->AddYError(1, sys_ReH_E00_FG, sys_ReH_E00_FG);
   gE12_xB_FG[0]->AddYError(3, sys_ReH_E12_FG, sys_ReH_E12_FG);

   gE00_xB_FG[1] = new TGraphMultiErrors("gE00_xB_ReE", "gE00_xB_ReE", 1, xB_E00_FG, ReE_E00_FG, zero, zero, eReE_E00_FG, eReE_E00_FG);
   gE12_xB_FG[1] = new TGraphMultiErrors("gE12_xB_ReE", "gE12_xB_ReE", 3, xB_E12_FG, ReE_E12_FG, zero3, zero3, eReE_E12_FG, eReE_E12_FG);
   gE00_xB_FG[1]->AddYError(1, sys_ReE_E00_FG, sys_ReE_E00_FG);
   gE12_xB_FG[1]->AddYError(3, sys_ReE_E12_FG, sys_ReE_E12_FG);

   gE00_xB_FG[2] = new TGraphMultiErrors("gE00_xB_ReHt", "gE00_xB_ReHt", 1, xB_E00_FG, ReHt_E00_FG, zero, zero, eReHt_E00_FG, eReHt_E00_FG);
   gE12_xB_FG[2] = new TGraphMultiErrors("gE12_xB_ReHt", "gE12_xB_ReHt", 3, xB_E12_FG, ReHt_E12_FG, zero3, zero3, eReHt_E12_FG, eReHt_E12_FG);
   gE00_xB_FG[2]->AddYError(1, sys_ReHt_E00_FG, sys_ReHt_E00_FG);
   gE12_xB_FG[2]->AddYError(3, sys_ReHt_E12_FG, sys_ReHt_E12_FG);

    for (int ipar = 0; ipar<3; ipar++){

      gE00_xB_FG[ipar]->SetMarkerStyle(24);
      gE00_xB_FG[ipar]->SetMarkerColor(2);
      gE00_xB_FG[ipar]->SetLineColor(2);
      gE00_xB_FG[ipar]->SetMarkerSize(1.3);
      gE00_xB_FG[ipar]->GetAttLine(0)->SetLineColor(2);
      gE00_xB_FG[ipar]->GetAttLine(1)->SetLineColor(2);
      gE00_xB_FG[ipar]->GetAttFill(1)->SetFillStyle(0); 

      gE12_xB_FG[ipar]->SetMarkerStyle(20);
      gE12_xB_FG[ipar]->SetMarkerColor(2);
      gE12_xB_FG[ipar]->SetLineColor(2);
      gE12_xB_FG[ipar]->SetMarkerSize(1.3);
      gE12_xB_FG[ipar]->GetAttLine(0)->SetLineColor(2);
      gE12_xB_FG[ipar]->GetAttLine(1)->SetLineColor(2);
      gE12_xB_FG[ipar]->GetAttFill(1)->SetFillStyle(0);  
   
      gE00_xB[ipar]->SetMarkerStyle(24);
      gE00_xB[ipar]->SetMarkerColor(kAzure-7);
      gE00_xB[ipar]->SetLineColor(kAzure-7);
      gE00_xB[ipar]->SetMarkerSize(1.3);
      gE00_xB[ipar]->GetAttLine(0)->SetLineColor(kAzure-7);
      gE00_xB[ipar]->GetAttLine(1)->SetLineColor(kAzure-7);
      gE00_xB[ipar]->GetAttFill(1)->SetFillStyle(0); 

      gE07_xB[ipar]->SetMarkerStyle(22);
      gE07_xB[ipar]->SetMarkerColor(kAzure-7);
      gE07_xB[ipar]->SetLineColor(kAzure-7);
      gE07_xB[ipar]->SetMarkerSize(1.3);
      gE07_xB[ipar]->GetAttLine(0)->SetLineColor(kAzure-7);
      gE07_xB[ipar]->GetAttLine(1)->SetLineColor(kAzure-7);
      gE07_xB[ipar]->GetAttFill(1)->SetFillStyle(0);  

      gE12_xB[ipar]->SetMarkerStyle(20);
      gE12_xB[ipar]->SetMarkerColor(kAzure-7);
      gE12_xB[ipar]->SetLineColor(kAzure-7);
      gE12_xB[ipar]->SetMarkerSize(1.3);
      gE12_xB[ipar]->GetAttLine(0)->SetLineColor(kAzure-7);
      gE12_xB[ipar]->GetAttLine(1)->SetLineColor(kAzure-7);
      gE12_xB[ipar]->GetAttFill(1)->SetFillStyle(0);  
   } 


    
   // TCanvas *cgr;
   //          //cgr = new TCanvas("CFFs_HallA", "CFFs_HallA",415,257,1513,987);
   //          cgr = new TCanvas("CFFs_HallA", "CFFs_HallA",76,64,2483,824); // horizontal
   //          //cgr = new TCanvas("xmi_sys_HallA_E00", "xmi_sys_HallA_E00",463,64,703,1371); // vertical
   //          cgr->Divide(3,1);   

   TLine *l0 = new TLine(0.35, 0, 0.62, 0);
          l0->SetLineStyle(7);

   for (Int_t ipar = 0; ipar < 3; ipar++){
      mgr_results[ipar]->Add(gkm15[ipar],"l");
      mgr_results[ipar]->Add(gE00_xB[ipar], "PS ; z ; 5 s=5");
      mgr_results[ipar]->Add(gE07_xB[ipar], "PS ; z ; 5 s=5");
      mgr_results[ipar]->Add(gE12_xB[ipar], "PS ; z ; 5 s=5");
      mgr_results[ipar]->Add(gE00_xB_FG[ipar], "PS ; z ; 5 s=5"); 
      mgr_results[ipar]->Add(gE12_xB_FG[ipar], "PS ; z ; 5 s=5");

      //mgr_results[ipar]->Add(gkm15[ipar],"L");
      // if(ipar==0){
      //    mgr_results[ipar]->Add(gkum15,"Pz");
      //    mgr_results[ipar]->Add(gkum15_sys, "[]");
      //    mgr_results[ipar]->Add(gDupre);
      //    mgr_results[ipar]->Add(gBG);
      //    mgr_results[ipar]->Add(gMou);
      // }
      

      mgr_results[ipar]->GetXaxis()->SetLimits(0.32,0.63); 
      //mgr_results[ipar]->GetXaxis()->SetNdivisions(20);
      if(ipar==0) mgr_results[ipar]->GetYaxis()->SetRangeUser(-5.9, 4.6);
      if(ipar==1) mgr_results[ipar]->GetYaxis()->SetRangeUser(-11, 11);
      if(ipar==2) mgr_results[ipar]->GetYaxis()->SetRangeUser(-3.1, 7.1);
      // if(ipar==3) mgr_results[ipar]->GetYaxis()->SetRangeUser(0, 0.039);
      mgr_results[ipar]->GetXaxis()->SetLabelSize(0.06);
      mgr_results[ipar]->GetXaxis()->SetTitleSize(0.07);
      //mgr_results[ipar]->GetXaxis()->SetTitleOffset(0.94);
      mgr_results[ipar]->GetYaxis()->SetLabelSize(0.06);
      mgr_results[ipar]->GetYaxis()->SetTitleSize(0.08);
      //mgr_results[ipar]->GetYaxis()->SetTitleOffset(0.94);
      //if(ipar==3) mgr_results[ipar]->GetYaxis()->SetTitleOffset(1.1);
      mgr_results[ipar]->GetXaxis()->SetAxisColor(kGray);
      mgr_results[ipar]->GetYaxis()->SetAxisColor(kGray);

      // cgr->cd(ipar+1);
      // gPad->SetFrameLineColor(kGray);
      // gStyle->SetLineWidth(1);
      // gStyle->SetLineColor(kGray);
      // gPad->SetFrameLineWidth(1);
      // mgr_results[ipar]->Draw("AP");     
      // l0->Draw("same");
   } 
  
   gStyle->SetLineScalePS(2);
   gStyle->SetLineColor(kGray);
   gStyle->SetEndErrorSize(5);

   TLegend *leg = new TLegend(0.73,0.73,0.88,0.78);   
   leg->AddEntry(gE00_xB[0], "#chi MI [15]", "fpl");
   leg->AddEntry(gE07_xB[0], "#chi MI [17]", "fpl");
   leg->AddEntry(gE12_xB[0], "#chi MI [22]", "fpl");
   leg->AddEntry(gE00_xB_FG[0], "F.Georges [15]", "fpl");    
   leg->AddEntry(gE12_xB_FG[0], "F.Georges [22]", "fpl");    
   // leg->AddEntry(gBG, "Boer/Guidal 14", "pl");
   // leg->AddEntry(gMou, "Moutarde 09", "pl");
   leg->AddEntry(gkm15[0], "KM15*", "l");
   leg->SetBorderSize(1);
   leg->SetTextSize(0.07);
   leg->Draw("same"); 

   TH1D *hDummy = new TH1D("hDummy", "hDummy", 50, 0, 50);
   hDummy->GetXaxis()->SetTitle("x");
   hDummy->GetYaxis()->SetTitle("y");

   TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 1200);
   canvas->cd();
   canvas->Draw();

   Double_t l = 0.1, r = 0.1, t = 0.1, b = 0.1;

   TPad *pad_13 = new TPad("pad_13", "pad_13", 0.67, 0.67, 1., 1.);
   pad_13->SetLeftMargin(0.);
   pad_13->SetRightMargin(r);
   pad_13->SetTopMargin(t);
   pad_13->SetBottomMargin(b);
   pad_13->SetFrameLineColor(kGray);
   pad_13->Draw();
   pad_13->cd();

   hDummy->Draw();

   canvas->cd();

   TPad *pad_12 = new TPad("pad_12", "pad_12", 0.34, 0.67, 0.67, 1.);
   pad_12->SetLeftMargin(0.);
   pad_12->SetRightMargin(0.);
   pad_12->SetTopMargin(t);
   pad_12->SetBottomMargin(b);
   pad_12->SetFrameLineColor(kGray);
   pad_12->Draw();
   pad_12->cd();

   hDummy->Draw();

   canvas->cd();

   TPad *pad_11 = new TPad("pad_11", "pad_11", 0.01, 0.67, 0.34, 1.);
   pad_11->SetLeftMargin(l);
   pad_11->SetRightMargin(0.);
   pad_11->SetTopMargin(t);
   pad_11->SetBottomMargin(b);
   pad_11->SetFrameLineColor(kGray);
   pad_11->Draw();
   pad_11->cd();

   hDummy->Draw();

   canvas->cd();

   TPad *pad_23 = new TPad("pad_23", "pad_23", 0.67, 0.34, 1., 0.67);
   pad_23->SetLeftMargin(0.);
   pad_23->SetRightMargin(r);
   pad_23->SetTopMargin(t);
   pad_23->SetBottomMargin(b);
   pad_23->SetFrameLineColor(kGray);
   pad_23->Draw();
   pad_23->cd();
   hDummy->Draw();

   canvas->cd();

   TPad *pad_22 = new TPad("pad_22", "pad_22", 0.34, 0.34, 0.67, 0.67);
   pad_22->SetLeftMargin(0.);
   pad_22->SetRightMargin(0.);
   pad_22->SetTopMargin(t);
   pad_22->SetBottomMargin(b);
   pad_22->SetFrameLineColor(kGray);
   pad_22->Draw();
   pad_22->cd();
   hDummy->Draw();

   canvas->cd();

   TPad *pad_21 = new TPad("pad_21", "pad_21", 0.01, 0.34, 0.34, 0.67);
   pad_21->SetLeftMargin(l);
   pad_21->SetRightMargin(0.);
   pad_21->SetTopMargin(t);
   pad_21->SetBottomMargin(b);
   pad_21->SetFrameLineColor(kGray);
   pad_21->Draw();
   pad_21->cd();
   hDummy->Draw();

   canvas->cd();

   TPad *pad_33 = new TPad("pad_33", "pad_33", 0.67, 0.01, 1., 0.34);
   pad_33->SetLeftMargin(0.);
   pad_33->SetRightMargin(r);
   pad_33->SetTopMargin(t);
   pad_33->SetBottomMargin(b); // this is reducing the pad's height
   pad_33->SetFrameLineColor(kGray);
   pad_33->Draw();
   pad_33->cd();
   hDummy->Draw();

   canvas->cd();

   TPad *pad_32 = new TPad("pad_32", "pad_32", 0.34, 0.01, 0.67, 0.34);
   pad_32->SetLeftMargin(0.);
   pad_32->SetRightMargin(0.);
   pad_32->SetTopMargin(t);
   pad_32->SetBottomMargin(b); // this is reducing the pad's height
   pad_32->SetFrameLineColor(kGray);
   pad_32->Draw();
   pad_32->cd();
   hDummy->Draw();

   canvas->cd();

   TPad *pad_31 = new TPad("pad_31", "pad_31", 0.01, 0.01, 0.34, 0.34);
   pad_31->SetLeftMargin(l);
   pad_31->SetRightMargin(0.);
   pad_31->SetTopMargin(t);
   pad_31->SetBottomMargin(b); // this is reducing the pad's height
   pad_31->SetFrameLineColor(kGray);
   pad_31->Draw();
   pad_31->cd();
   hDummy->Draw();

   canvas->cd();


   
   return 0;



}