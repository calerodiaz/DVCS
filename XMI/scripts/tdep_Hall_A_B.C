using namespace std;
using namespace ROOT;
using namespace ROOT::RDF;
//___________________________________________________________________________
auto km15Graphs(RNode df){

   vector<TGraphErrors*>gkm15_hA(4);

   for (Int_t ipar = 0; ipar < 4; ipar++) {
            gkm15_hA[ipar] = new TGraphErrors();
            gkm15_hA[ipar]->SetName("true");
            gkm15_hA[ipar]->SetLineColor(kGray+2);
            gkm15_hA[ipar]->SetLineWidth(2);
            gkm15_hA[ipar]->SetMarkerColor(kGray+2);
            gkm15_hA[ipar]->SetMarkerSize(0.001);
   }
   auto fillkm15Graph = [&gkm15_hA](double t, double ReH_km15, double ReE_km15, double ReHtilde_km15, double dvcs_km15){
      gkm15_hA[0]->SetPoint(gkm15_hA[0]->GetN(), -t, ReH_km15);
      gkm15_hA[0]->SetPointError( gkm15_hA[0]->GetN()-1, 0, 0 );
      gkm15_hA[1]->SetPoint(gkm15_hA[1]->GetN(), -t, ReE_km15);
      gkm15_hA[1]->SetPointError( gkm15_hA[1]->GetN()-1, 0, 0 );
      gkm15_hA[2]->SetPoint(gkm15_hA[2]->GetN(), -t, ReHtilde_km15);
      gkm15_hA[2]->SetPointError( gkm15_hA[2]->GetN()-1, 0, 0 );
      gkm15_hA[3]->SetPoint(gkm15_hA[3]->GetN(), -t, dvcs_km15);
      gkm15_hA[3]->SetPointError( gkm15_hA[3]->GetN()-1, 0, 0 );
   };
   df.Foreach(fillkm15Graph, {"t", "ReH_km15", "ReE_km15", "ReHt_km15", "dvcs_km15"});

   return gkm15_hA;
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

   int nsets = df.Count().GetValue();
   Double_t *eReH = new Double_t [2];   
   Double_t *eReE = new Double_t [2];  
   Double_t *eReHt = new Double_t [2]; 
   Double_t *edvcs = new Double_t [2]; 

   Int_t dark   = TColor::GetColorDark(kAzure-6);
   Int_t bright = TColor::GetColorBright(kAzure-6);

   for (int ipar = 0; ipar<4; ipar++){
      gr[ipar] = new TGraphMultiErrors(nsets,2);
      gr[ipar]->SetMarkerStyle(20);
      gr[ipar]->SetMarkerSize(1.3);
      gr[ipar]->SetMarkerColor(kAzure-7);
      gr[ipar]->SetLineColor(kAzure-7);
      gr[ipar]->GetAttLine(0)->SetLineColor(kAzure-7);
      gr[ipar]->GetAttLine(1)->SetLineColor(33);
      //gr[ipar]->GetAttLine(1)->SetLineWidth(0);
      gr[ipar]->GetAttFill(1)->SetFillStyle(0);
      //gr[ipar]->GetAttFill(1)->SetFillColor(17);
      //gr[ipar]->GetAttFill(1)->SetFillStyle(1001);
      //gr[ipar]->GetAttFill(1)->SetFillColorAlpha(kBlack, 0.2);

      gr[ipar]->SetFillStyle(0);
      //gr[ipar]->SetFillColor(17);
      //gr[ipar]->SetFillColorAlpha(kBlack, 0.2);
   } 

   for (int i = 0; i<nsets; i++){
         
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
Int_t tdep_Hall_A_B(const Int_t initial_set = 6, const Int_t last_set = 10)
{
   ROOT::EnableImplicitMT();
   TH2::AddDirectory(false); // No dot write histogram automaticaly to the root file

   cout<<"Reading csv into dfs...";
   // fit results df
   auto df_hA = FromCSV("/media/lily/Data/GPDs/HallA/realdata/maps/results_HallA_sets_1to85_data_5pct.csv");
   auto df_hB = FromCSV("/media/lily/Data/GPDs/HallB/realdata/maps/results_HallB_sets_1to110_data_5pct.csv");
   cout<<"done!"<<endl;
   TMultiGraph *mgr_hA[4];
   TMultiGraph *mgr_hB[4];
   Int_t fMarker[] = {22, 20, 20, 20, 20, 21, 23, 22, 34};
   Color_t col[] = {kBlack, kGreen+2, kBlue, kOrange+2, kMagenta+1, kRed-2, kAzure+9, kGray+2, kYellow+1};
   TString parsname[] = {"#displaystyle #mathfrak{Re}#mathcal{H}", "#displaystyle #mathfrak{Re}#mathcal{E}", "#displaystyle #mathfrak{Re}#mathcal{#widetilde{H}}",
                     "#displaystyle d^{4}#sigma^{DVCS}_{UU} [nb/GeV^{4}]"};
   TString kFormulation = "BKM10";
   TString kModel = "KM15";

   for (Int_t ipar = 0; ipar < 4; ipar++){
      mgr_hA[ipar] = new TMultiGraph();
      mgr_hA[ipar] ->SetTitle(Form("; ; %s", parsname[ipar].Data()));
      if (ipar == 3 ) mgr_hA[ipar] ->SetTitle(Form("; #displaystyle -t[GeV^{2}]; %s", parsname[ipar].Data()));
      mgr_hB[ipar] = new TMultiGraph();
      if (ipar == 3 ) mgr_hB[ipar] ->SetTitle(Form("; #displaystyle -t[GeV^{2}]; "));
   }   
   // Filter desired sets
   //auto setrange_hA = [&initial_set, &last_set](Long64_t set) { return set >= initial_set && set <= last_set; };
   //auto df1_hA = df_hA.Filter(setrange_hA, {"set"});
   auto df1_hA = df_hA.Filter("set >= 6 && set <= 10");
   auto df1_hB = df_hB.Filter("set >= 68 && set <= 72");
   
   // Get mean xB and QQ
   auto k_mean_hA = df1_hA.Mean("k");
   auto xB_mean_hA = df1_hA.Mean("xB");
   auto QQ_mean_hA = df1_hA.Mean("QQ");
   auto k_mean_hB = df1_hB.Mean("k");
   auto xB_mean_hB = df1_hB.Mean("xB");
   auto QQ_mean_hB = df1_hB.Mean("QQ");

   auto gkm15_hA = km15Graphs(df1_hA);
   auto gxmi_hA = xmiGraphs(df1_hA);
   auto gkm15_hB = km15Graphs(df1_hB);
   auto gxmi_hB = xmiGraphs(df1_hB);

   // ======================== REFERENCES RESULTS ===========================================
   // HALL-A
   // ----------------------------------------------------------------
   // Dupre_hA point at t = -0.371, xB = 0.334, QQ = 2.23
   Double_t Dupre_hA_ReH[1] = {-3.32};
   Double_t Dupre_hA_dReH[1] = {4.456-3.32};
   Double_t t_Dupre_hA[1] = {0.369 - 0.004}; // Shift added
   auto gDupre_hA = new TGraphErrors(1, t_Dupre_hA, Dupre_hA_ReH, 0, Dupre_hA_dReH);
   gDupre_hA->SetMarkerColor(kGreen+2);
   gDupre_hA->SetLineColor(kGreen+2);
   gDupre_hA->SetMarkerStyle(23);
   gDupre_hA->SetMarkerSize(1.3);
   // ----------------------------------------------------------------
   //Boer-Guidal 14: 2 points at xB = 0.36, Q2 = 2.3
   Double_t t_BG[2]   = {0.279 + 0.004, 0.328 + 0.004}; // Shift added
   Double_t BG_ReH[2]   = {-1.098, -3.039};
   Double_t et[2] = {0., 0.};
   Double_t BG_ReH_el[2] = {2.862-1.098, 3.496-3.039};
   Double_t BG_ReH_eh[2] = {1.098-0.576, 3.039-2.517};
   auto gBG = new TGraphAsymmErrors(2, t_BG, BG_ReH, et, et, BG_ReH_el, BG_ReH_eh);
   gBG->SetMarkerColor(kMagenta);
   gBG->SetLineColor(kMagenta);
   gBG->SetMarkerStyle(33);
   gBG->SetMarkerSize(1.3);
   // ----------------------------------------------------------------
   // Moutarde 09 - local: 4 points at xB = 0.36, Q2 = 2.3 
   Double_t ReH_Mou[4] = {0.725, -1.353, -2.217, -3.211};
   Double_t eReH_Mou[4] = {0.966-0.725, 1.353-1.068, 2.217-1.999, 3.211-2.983};
   Double_t t_Mou[4] = {0.172 - 0.004, 0.231 - 0.006, 0.280 - 0.006, 0.330 - 0.012}; // Added 0.004 shift for visualization
   auto gMou = new TGraphErrors(4, t_Mou, ReH_Mou, 0, eReH_Mou);
   gMou->SetMarkerColor(2);
   gMou->SetLineColor(2);
   gMou->SetMarkerStyle(22);
   gMou->SetMarkerSize(1.3);
   // ----------------------------------------------------------------
   //KM15 paper fir
   Double_t ReH_kum15[5] = {-1.181, -1.937, -2.818, -3.039, -2.330};
   Double_t t_kum15[5] = {0.175, 0.231, 0.278, 0.324, 0.371};
   Double_t et_kum15[5] = {0.0001, 0.0001, 0.0001, 0.0001, 0.0001};
   Double_t eReH_kum15[5] = {-0.755+1.181, -1.472+1.937, -2.157+2.818, -2.212+3.039, -1.102+2.330};
   Double_t esysReH_kum15[5] = {1.322+1.181, 0.669+1.937, -0.039+2.818, -0.157+3.039, 0.866+2.330};
   auto gkum15 = new TGraphErrors(5, t_kum15, ReH_kum15, et_kum15, eReH_kum15);
   auto gkum15_sys = new TGraphErrors(5, t_kum15, ReH_kum15, et_kum15, esysReH_kum15);
   // gkum15->AddYError(5, esysReH_kum15, esysReH_kum15);
   gkum15->SetMarkerStyle(21);
   gkum15->SetMarkerColor(1);
   gkum15->SetLineColor(1);
   gkum15_sys->SetMarkerStyle(21);
   gkum15_sys->SetMarkerColor(1);
   gkum15_sys->SetLineColor(1);
   gkum15_sys->SetMarkerSize(1.3);
   // ----------------------------------------------------------------
   // HALL-B
   // ----------------------------------------------------------------
   // HallB paper ReH results at xB = 0.304, QQ = 2.10
   Double_t ReH_e1[4] = {-0.3783, 0.1391, 0.1425, 0.1728};
   Double_t t_e1[4] = {0.203, 0.262, 0.341, 0.447};
   Double_t et_e1[4] = {0., 0., 0., 0.};
   Double_t eReH_low_e1[4] = {1.485 - 0.3783, 0.929+0.1391, 1.016+0.1425, 1.166+0.1728};
   Double_t eReH_high_e1[4] = {0.3783 + 1.359, -0.1391 + 1.439, -0.1425 + 0.979, -0.1728 + 0.675};

   auto ge1 = new TGraphAsymmErrors(4, t_e1, ReH_e1, et_e1, et_e1, eReH_low_e1, eReH_high_e1);
   ge1->SetMarkerColor(kOrange+2);
   ge1->SetLineColor(kOrange+2);
   ge1->SetMarkerStyle(21);
   ge1->SetMarkerSize(1.3);

   // Dupre results
   Double_t ReH_Dupre[5] = {-0.884, -0.2115, 0.1346, -0.1346, -0.5961};
   Double_t t_Dupre[5] = {0.163, 0.213, 0.2719, 0.351, 0.458};
   Double_t eReH_Dupre[5] = {2-0.884, 1.403-0.2115, 0.884+0.1346, 1.115-0.1346, 1.519-0.5961};
   Double_t ezero[5] = {0., 0., 0., 0., 0.};
   auto gDupre_hB = new TGraphErrors(5, t_Dupre, ReH_Dupre, ezero, eReH_Dupre);
   gDupre_hB->SetMarkerColor(kGreen+2);
   gDupre_hB->SetLineColor(kGreen+2);
   gDupre_hB->SetMarkerStyle(23);
   gDupre_hB->SetMarkerSize(1.3);

   // ===========================================================================================

    for (Int_t ipar = 0; ipar < 4; ipar++){
      mgr_hA[ipar]->Add(gxmi_hA[ipar], "PS ; z ; 5 s=30");
      mgr_hA[ipar]->Add(gkm15_hA[ipar],"L");
      mgr_hB[ipar]->Add(gxmi_hB[ipar], "PS ; z ; 5 s=40");
      mgr_hB[ipar]->Add(gkm15_hB[ipar],"L");

      if(ipar==0){
         mgr_hA[ipar]->Add(gkum15,"Pz");
         mgr_hA[ipar]->Add(gkum15_sys, "[]");
         mgr_hA[ipar]->Add(gDupre_hA);
         mgr_hA[ipar]->Add(gBG);
         mgr_hA[ipar]->Add(gMou);
         mgr_hB[ipar]->Add(gDupre_hB);
         mgr_hB[ipar]->Add(ge1);
      }      

      mgr_hA[ipar]->GetXaxis()->SetLimits(0.13,0.41); 
      mgr_hB[ipar]->GetXaxis()->SetLimits(0.11,0.48); 

      if(ipar==0) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==1) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==2) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==3) mgr_hA[ipar]->GetYaxis()->SetRangeUser(0.001, 0.039); 
      mgr_hA[ipar]->GetYaxis()->SetNdivisions(505);

      if(ipar==0) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==1) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==2) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-9, 11);
      if(ipar==3) mgr_hB[ipar]->GetYaxis()->SetRangeUser(0.001, 0.039);  
      mgr_hB[ipar]->GetYaxis()->SetNdivisions(505);
 

      mgr_hA[ipar]->GetXaxis()->SetLabelSize(0.08);
      mgr_hA[ipar]->GetXaxis()->SetTitleSize(0.09);
      //mgr_hA[ipar]->GetXaxis()->SetTitleOffset(0.94);
      mgr_hA[ipar]->GetYaxis()->SetLabelSize(0.08);
      mgr_hA[ipar]->GetYaxis()->SetTitleSize(0.09);
      mgr_hA[ipar]->GetXaxis()->SetAxisColor(kGray);
      mgr_hA[ipar]->GetYaxis()->SetAxisColor(kGray);

      mgr_hB[ipar]->GetXaxis()->SetLabelSize(0.08);
      mgr_hB[ipar]->GetXaxis()->SetTitleSize(0.09);
      mgr_hB[ipar]->GetYaxis()->SetLabelSize(0.08);
      mgr_hB[ipar]->GetYaxis()->SetTitleSize(0.09);
      mgr_hB[ipar]->GetXaxis()->SetAxisColor(kGray);
      mgr_hB[ipar]->GetYaxis()->SetAxisColor(kGray);

      // cgr->cd(ipar+1);
      // gPad->SetFrameLineColor(kGray);
      // gStyle->SetLineWidth(1);
      // gStyle->SetLineColor(kGray);
      // gPad->SetFrameLineWidth(1);
      // mgr_hA[ipar]->Draw("AP");       
   } 

   gStyle->SetLineScalePS(2);
   gStyle->SetLineColor(kGray);
   gStyle->SetEndErrorSize(5);

   TLegend *leg = new TLegend(0.73,0.6,0.88,0.78);   
   leg->SetHeader("Hall-A E00-110");
   leg->SetNColumns(2);
   leg->AddEntry(gxmi_hA[0], "#chi MI", "fpl");
   leg->AddEntry(gDupre_hA, "Dupre 17", "pl");
   leg->AddEntry(gkum15, "Kumericki/Muller 15", "pl");    
   leg->AddEntry(gBG, "Boer/Guidal 14", "pl");
   leg->AddEntry(gMou, "Moutarde 09", "pl");
   leg->AddEntry(gkm15_hA[0], "KM15*", "l");
   leg->SetBorderSize(1);
   leg->SetTextSize(0.07);

   TLegend *leg3 = new TLegend(0.73,0.6,0.88,0.78);   
   leg3->SetHeader("Hall-B e1-DVCS1");
   leg3->AddEntry(gxmi_hB[0], "#chi MI", "fpl");
   leg3->AddEntry(gDupre_hB, "Dupre 17", "pl");
   leg3->AddEntry(ge1, "H. S. Jo 15", "pl");
   leg3->AddEntry(gkm15_hA[0], "KM15*", "l");    
   leg3->SetBorderSize(1);
   leg3->SetTextSize(0.07);


   Double_t l = 0.15, r = 0.15, t = 0., b = 0.1;

   //gPad->SetFrameLineWidth(1);
   
   TCanvas *canvas = new TCanvas("canvas", "canvas", 1024, 1300);
   canvas->cd();
   canvas->Draw();   

   TPad *pad_11 = new TPad("pad_11", "pad_11", 0.5, 0.725, 1., 0.95);
   pad_11->SetLeftMargin(0.);
   pad_11->SetRightMargin(r);
   pad_11->SetTopMargin(t);
   pad_11->SetBottomMargin(0.);
   pad_11->SetFrameLineColor(kGray);
   pad_11->Draw();
   pad_11->cd();

   mgr_hB[0]->Draw("AP");  
   leg3->Draw("same");
   canvas->cd();

   TPad *pad_1 = new TPad("pad_1", "pad_01", 0., 0.725, 0.5, 0.95);
   pad_1->SetLeftMargin(l);
   pad_1->SetRightMargin(0.);
   pad_1->SetTopMargin(t);
   pad_1->SetBottomMargin(0.);
   pad_1->SetFrameLineColor(kGray);
   pad_1->Draw();
   pad_1->cd();

   mgr_hA[0]->Draw("AP");  
   leg->Draw("same");
   
   canvas->cd();

   TPad *pad_22 = new TPad("pad_22", "pad_22", 0.5, 0.5, 1., 0.725);
   pad_22->SetLeftMargin(0.);
   pad_22->SetRightMargin(r);
   pad_22->SetTopMargin(0.);
   pad_22->SetBottomMargin(0.);
   pad_22->SetFrameLineColor(kGray);
   pad_22->Draw();
   pad_22->cd();

   mgr_hB[1]->Draw("AP");  

   canvas->cd();

   TPad *pad_2 = new TPad("pad_2", "pad_2", 0., 0.5, 0.5, 0.725);
   pad_2->SetLeftMargin(l);
   pad_2->SetRightMargin(0.);
   pad_2->SetTopMargin(0.);
   pad_2->SetBottomMargin(0.);
   pad_2->SetFrameLineColor(kGray);
   pad_2->Draw();
   pad_2->cd();

   mgr_hA[1]->Draw("AP");  

   canvas->cd();

   TPad *pad_33 = new TPad("pad_33", "pad_33", 0.5, 0.275, 1., 0.5);
   pad_33->SetLeftMargin(0.);
   pad_33->SetRightMargin(r);
   pad_33->SetTopMargin(0.);
   pad_33->SetBottomMargin(0.);
   pad_33->SetFrameLineColor(kGray);
   pad_33->Draw();
   pad_33->cd();

   mgr_hB[2]->Draw("AP");  

   canvas->cd();
    
   TPad *pad_3 = new TPad("pad_3", "pad_3", 0., 0.275, 0.5, 0.5);
   pad_3->SetLeftMargin(l);
   pad_3->SetRightMargin(0.);
   pad_3->SetTopMargin(0.);
   pad_3->SetBottomMargin(0.);
   pad_3->SetFrameLineColor(kGray);
   pad_3->Draw();
   pad_3->cd();

   mgr_hA[2]->Draw("AP");  

   canvas->cd();

   TPad *pad_44 = new TPad("pad_44", "pad_44", 0.5, 0.05, 1., 0.275);
   pad_44->SetLeftMargin(0.);
   pad_44->SetRightMargin(r);
   pad_44->SetTopMargin(0.);
   pad_44->SetBottomMargin(0.);
   pad_44->SetFrameLineColor(kGray);
   pad_44->Draw();
   pad_44->cd();

   mgr_hB[3]->Draw("AP");  

   canvas->cd();

   TPad *pad_4 = new TPad("pad_4", "pad_4", 0., 0.05, 0.5, 0.275);
   pad_4->SetLeftMargin(l);
   pad_4->SetRightMargin(0.);
   pad_4->SetTopMargin(0.);
   pad_4->SetBottomMargin(0.);
   pad_4->SetFrameLineColor(kGray);
   pad_4->Draw();
   pad_4->cd();

   mgr_hA[3]->Draw("AP");  

   canvas->cd();

   TPad *pad_kin_hA = new TPad("pad_kin_hA", "pad_kin_hA",0.06,0.925,0.5,0.99);
      pad_kin_hA->Draw();
      pad_kin_hA->cd();
      pad_kin_hA->SetFillStyle(0);
      pad_kin_hA->SetBorderSize(1);   
   auto tex_kinA = new TLatex(0.5,0.5,Form("k = %.2f GeV, Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *k_mean_hA, *QQ_mean_hA, *xB_mean_hA));
      tex_kinA->SetTextAlign(22);
      tex_kinA->SetTextSize(0.35);
      tex_kinA->Draw();

   canvas->cd();

   TPad *pad_kin_hB = new TPad("pad_kin_hB", "pad_kin_hB",0.45,0.925,0.99,0.99);
      pad_kin_hB->Draw();
      pad_kin_hB->cd();
      pad_kin_hB->SetFillStyle(0);
      pad_kin_hB->SetBorderSize(1);      
   auto tex_kinB = new TLatex(0.5,0.5,Form("k = %.2f GeV, Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *k_mean_hB, *QQ_mean_hB, *xB_mean_hB));
      tex_kinB->SetTextAlign(22);
      tex_kinB->SetTextSize(0.35);
      tex_kinB->Draw();

   //  TPad *pad_1 = new TPad("pad_1", "pad_1", 0., 0.75, 0.5, 1.);
   //  pad_1->SetLeftMargin(l);
   //  pad_1->SetTopMargin(t);
   //  pad_1->SetBottomMargin(0.);
   //  pad_1->Draw();
   //  pad_1->cd();

   //  hDummy->Draw();

   //  canvas->cd();

   //  double yp = 0.22;
   //  TPad *pad_2 = new TPad("pad_2", "pad_2", 0., 0.5, 0.5, 0.75);
   //  pad_2->SetLeftMargin(l);
   //  pad_2->SetRightMargin(r);
   //  pad_2->SetTopMargin(0.);
   //  pad_2->SetBottomMargin(0.);
   //  pad_2->Draw();
   //  pad_2->cd();
   //  hDummy->Draw();

   //  canvas->cd();

   //  TPad *pad_3 = new TPad("pad_3", "pad_3", 0., 0.25, 0.5, 0.5);
   //  pad_3->SetLeftMargin(l);
   //  pad_3->SetRightMargin(r);
   //  pad_3->SetTopMargin(0.);
   //  pad_3->SetBottomMargin(0.); // this is reducing the pad's height
   //  pad_3->Draw();
   //  pad_3->cd();
   //  hDummy->Draw();


   // TPad *pad_4 = new TPad("pad_4", "pad_4", 0., 0., 0.5, 0.25);
   //  pad_4->SetLeftMargin(l);
   //  pad_4->SetRightMargin(r);
   //  pad_4->SetTopMargin(0.);
   //  pad_4->SetBottomMargin(0.); // this is reducing the pad's height
   //  pad_4->Draw();
   //  pad_4->cd();
   //  hDummy->Draw();
    
   //  TCanvas *cgr;
   //          cgr = new TCanvas("CFFs_HallA", "CFFs_HallA",415,257,1513,987);
   //          //cgr = new TCanvas("xmi_sys_HallA_E00", "xmi_sys_HallA_E00",463,64,703,1371); // vertical
   //          cgr->Divide(2,2);
   // gStyle->SetEndErrorSize(5);
   // for (Int_t ipar = 0; ipar < 4; ipar++){
   //    mgr_hA[ipar]->Add(gxmi_hA[ipar], "PS ; z ; 5 s=30");
   //    mgr_hA[ipar]->Add(gkm15_hA[ipar],"L");
   //    mgr_hB[ipar]->Add(gxmi_hB[ipar], "PS ; z ; 5 s=30");
   //    mgr_hB[ipar]->Add(gkm15_hB[ipar],"L");

   //    if(ipar==0){
   //       mgr_hA[ipar]->Add(gkum15,"Pz");
   //       mgr_hA[ipar]->Add(gkum15_sys, "[]");
   //       mgr_hA[ipar]->Add(gDupre_hA);
   //       mgr_hA[ipar]->Add(gBG);
   //       mgr_hA[ipar]->Add(gMou);
   //       mgr_hB[ipar]->Add(gDupre_hB);
   //       mgr_hB[ipar]->Add(ge1);
   //    }      

   //    mgr_hA[ipar]->GetXaxis()->SetLimits(0.13,0.41); 
   //    mgr_hB[ipar]->GetXaxis()->SetLimits(0.11,0.48); 

   //    if(ipar==0) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==1) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==2) mgr_hA[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==3) {mgr_hA[ipar]->GetYaxis()->SetRangeUser(0.001, 0.039); mgr_hA[ipar]->GetYaxis()->SetNdivisions(505);}

   //    if(ipar==0) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==1) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==2) mgr_hB[ipar]->GetYaxis()->SetRangeUser(-11, 11);
   //    if(ipar==3) { mgr_hB[ipar]->GetYaxis()->SetRangeUser(0.001, 0.039);  mgr_hB[ipar]->GetYaxis()->SetNdivisions(505);}
 

   //    mgr_hA[ipar]->GetXaxis()->SetLabelSize(0.06);
   //    mgr_hA[ipar]->GetXaxis()->SetTitleSize(0.06);
   //    //mgr_hA[ipar]->GetXaxis()->SetTitleOffset(0.94);
   //    mgr_hA[ipar]->GetYaxis()->SetLabelSize(0.06);
   //    mgr_hA[ipar]->GetYaxis()->SetTitleSize(0.06);
   //    mgr_hA[ipar]->GetXaxis()->SetAxisColor(kGray);
   //    mgr_hA[ipar]->GetYaxis()->SetAxisColor(kGray);

   //    mgr_hB[ipar]->GetXaxis()->SetLabelSize(0.06);
   //    mgr_hB[ipar]->GetXaxis()->SetTitleSize(0.06);
   //    mgr_hB[ipar]->GetYaxis()->SetLabelSize(0.06);
   //    mgr_hB[ipar]->GetYaxis()->SetTitleSize(0.06);
   //    mgr_hB[ipar]->GetXaxis()->SetAxisColor(kGray);
   //    mgr_hB[ipar]->GetYaxis()->SetAxisColor(kGray);

   //    cgr->cd(ipar+1);
   //    gPad->SetFrameLineColor(kGray);
   //    gStyle->SetLineWidth(1);
   //    gStyle->SetLineColor(kGray);
   //    gPad->SetFrameLineWidth(1);
   //    mgr_hA[ipar]->Draw("AP");       
   // } 
   // cgr->cd(0);
   // // TPad *padtitle = new TPad("padtitle", "padtitle",0.1,0.92,0.9,0.99);
   // // padtitle->Draw();
   // // padtitle->cd();
   // // padtitle->SetFillStyle(0);
   // // padtitle->SetBorderSize(1);
   // // gPad->SetFrameLineWidth(1);
   // // gStyle->SetLineScalePS(2);
   // // auto tex = new TLatex(0.5,0.5,Form("Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *QQ_mean_hA, *xB_mean_hA));
   // // tex->SetTextAlign(22);
   // // tex->SetTextSize(0.8);
   // // tex->Draw();
   // // cgr->cd(2); 


   // TLegend *leg = new TLegend(0.73,0.73,0.88,0.78);   
   // leg->SetHeader("Hall-A E00-110");
   // leg->SetNColumns(2);
   // leg->AddEntry(gxmi_hA[0], "#chi MI", "fpl");
   // leg->AddEntry(gDupre_hA, "Dupre 17", "pl");
   // leg->AddEntry(gkum15, "Kumericki/Muller 15", "pl");    
   // leg->AddEntry(gBG, "Boer/Guidal 14", "pl");
   // leg->AddEntry(gMou, "Moutarde 09", "pl");
   // leg->AddEntry(gkm15_hA[0], "KM15*", "l");
   // //leg->AddEntry(gkm15_hA[0], "KM15*", "l");
   // leg->SetBorderSize(1);
   // leg->SetTextSize(0.07);
   // leg->Draw("same"); 

   // TLegend *leg3 = new TLegend(0.73,0.73,0.88,0.78);   
   // leg3->SetHeader("Hall-B e1-DVCS1");
   // leg3->AddEntry(gxmi_hB[0], "#chi MI", "fpl");
   // leg3->AddEntry(gDupre_hB, "Dupre 17", "pl");
   // leg3->AddEntry(ge1, "H. S. Jo 15", "pl");
   // leg3->AddEntry(gkm15_hA[0], "KM15*", "l");    
   // //leg3->AddEntry(gBG, "Boer/Guidal 14", "pl");
   // //leg3->AddEntry(gMou, "Moutarde 09", "pl");
   // //leg->AddEntry(gkm15_hA[0], "KM15*", "l");
   // leg3->SetBorderSize(1);
   // leg3->SetTextSize(0.07);
   // //leg->Draw("same"); 

   // cgr->cd(2);
   // TLegend *leg2 = new TLegend(0.73,0.73,0.88,0.78);   
   // leg2->AddEntry(gkm15_hA[0], "KM15*", "l");
   // leg2->SetBorderSize(1);
   // leg2->SetTextSize(0.07);
   // leg2->Draw("same"); 


   // // Test canvas pad divide
   // auto C = (TCanvas*) gROOT->FindObject("C");
   // if (C) delete C;
   // C = new TCanvas("C","canvas",1024,1300);
   // C->SetFillStyle(4000);
 
   // // Number of PADS
   // const Int_t Nx = 2;
   // const Int_t Ny = 4;
 
   // // Margins
   // Float_t lMargin = 0.12;
   // Float_t rMargin = 0.05;
   // Float_t bMargin = 0.05;
   // Float_t tMargin = 0.05;
 
   // // Canvas setup
   // CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
 
   // // Dummy histogram.
   // auto h = (TH1F*) gROOT->FindObject("histo");
   // if (h) delete h;
   // h = new TH1F("histo","",100,-5.0,5.0);
   // h->FillRandom("gaus",10000);
   // h->GetXaxis()->SetTitle("x axis");
   // h->GetYaxis()->SetTitle("y axis");
 
   // TPad *pad[Nx][Ny];
 
   // for (Int_t i = 0; i < Nx; i++) {
   //    for (Int_t j = 0; j < Ny; j++) {
   //       C->cd(0);
 
   //       // Get the pads previously created.
   //       pad[i][j] = (TPad*) C->FindObject(TString::Format("pad_%d_%d",i,j).Data());
   //       pad[i][j]->Draw();
   //       pad[i][j]->SetFillStyle(4000);
   //       pad[i][j]->SetFrameFillStyle(4000);
   //       pad[i][j]->cd();

   //       pad[i][j]->SetFrameLineColor(kGray);
   //       //pad->SetFrameLineWidth(1);
 
   //       // Size factors
   //       Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
   //       Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();
 
   //       if (i == 0){
   //          mgr_hA[3-j]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
   //          mgr_hA[3-j]->Draw("AP");  
   //          if (3-j == 0)  {leg->Draw("same");}

   //       }
   //       if (i == 1){
   //          mgr_hB[3-j]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
   //          mgr_hB[3-j]->Draw("AP");   
   //          if (3-j == 0)  {leg3->Draw("same");  }    

   //       }


         //TH1F *hFrame = (TH1F*) h->Clone(TString::Format("h_%d_%d",i,j).Data());
 
         // // y axis range
         // hFrame->SetMinimum(0.0001); // do not show 0
         // hFrame->SetMaximum(1.2*h->GetMaximum());
 
         // // Format for y axis
         // hFrame->GetYaxis()->SetLabelFont(43);
         // hFrame->GetYaxis()->SetLabelSize(16);
         // hFrame->GetYaxis()->SetLabelOffset(0.02);
         // hFrame->GetYaxis()->SetTitleFont(43);
         // hFrame->GetYaxis()->SetTitleSize(16);
         // hFrame->GetYaxis()->SetTitleOffset(2);
 
         // hFrame->GetYaxis()->CenterTitle();
         // hFrame->GetYaxis()->SetNdivisions(505);
 
         // // TICKS Y Axis
         // hFrame->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
 
         // // Format for x axis
         // hFrame->GetXaxis()->SetLabelFont(43);
         // hFrame->GetXaxis()->SetLabelSize(16);
         // hFrame->GetXaxis()->SetLabelOffset(0.02);
         // hFrame->GetXaxis()->SetTitleFont(43);
         // hFrame->GetXaxis()->SetTitleSize(16);
         // hFrame->GetXaxis()->SetTitleOffset(1);
         // hFrame->GetXaxis()->CenterTitle();
         // hFrame->GetXaxis()->SetNdivisions(505);
 
         // // TICKS X Axis
         // hFrame->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
 
         // // Draw cloned histogram with individual settings
         // hFrame->Draw();
 
         // TText text;
         // text.SetTextAlign(31);
         // text.SetTextFont(43);
         // text.SetTextSize(10);
         // text.DrawTextNDC(XtoPad(0.9), YtoPad(0.8), gPad->GetName());

         // auto tex = new TLatex(0.5,0.5,Form("Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *QQ_mean_hA, *xB_mean_hA));
         // tex->SetTextAlign(22);
         // tex->SetTextSize(0.4);
         // if(i==0 && j == 3)
         // tex->Draw();
         //cgr->cd(2); 
   //    }
   // }
   // C->cd();

   // TPad *pad_kin_hA = new TPad("pad_kin_hA", "pad_kin_hA",0.1,0.92,0.5,0.99);
   // pad_kin_hA->Draw();
   // pad_kin_hA->cd();
   // pad_kin_hA->SetFillStyle(0);
   // pad_kin_hA->SetBorderSize(1);
   // gPad->SetFrameLineWidth(1);
   // gStyle->SetLineScalePS(2);
   // auto tex_kinA = new TLatex(0.5,0.5,Form("k = %.2f GeV, Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *k_mean_hA, *QQ_mean_hA, *xB_mean_hA));
   // tex_kinA->SetTextAlign(22);
   // tex_kinA->SetTextSize(0.3);
   // tex_kinA->Draw();
   // C->cd();
   // TPad *pad_kin_hB = new TPad("pad_kin_hB", "pad_kin_hB",0.5,0.92,0.99,0.99);
   // pad_kin_hB->Draw();
   // pad_kin_hB->cd();
   // pad_kin_hB->SetFillStyle(0);
   // pad_kin_hB->SetBorderSize(1);
   // gPad->SetFrameLineWidth(1);
   // gStyle->SetLineScalePS(2);
   // auto tex_kinB = new TLatex(0.5,0.5,Form("k = %.2f GeV, Q^{2} = %.2f GeV^{2}, x_{B} = %.2f", *k_mean_hB, *QQ_mean_hB, *xB_mean_hB));
   // tex_kinB->SetTextAlign(22);
   // tex_kinB->SetTextSize(0.3);
   // tex_kinB->Draw();
  
   return 0;

}

