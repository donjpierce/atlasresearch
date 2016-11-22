{

  // use plain black on white colors
  Int_t icol=0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetStatColor(icol);
  //gStyle->SetFillColor(icol);
  gStyle->SetFrameFillColor(icol);
  gStyle->SetHistFillColor(icol);
  gStyle->SetTitleFillColor(icol);

  //   // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  //   gStyle->SetPadTopMargin(0.05);
  //   gStyle->SetPadRightMargin(0.05);
  //   gStyle->SetPadBottomMargin(0.16);
  //   gStyle->SetPadLeftMargin(0.12);

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.05;
  gStyle->SetTextFont(font);
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetTitleFont(font,"z");
  gStyle->SetLabelSize(0.9*tsize,"x");
  gStyle->SetLabelSize(0.9*tsize,"y");
  gStyle->SetLabelSize(0.9*tsize,"z");
  gStyle->SetTitleSize(0.9*tsize,"x");
  gStyle->SetTitleSize(0.9*tsize,"y");
  gStyle->SetTitleSize(0.9*tsize,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  // get rid of X error bars and y error bar caps
  //  gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  //  gStyle->SetOptTitle(0);
  //  gStyle->SetOptStat(1111);
  //  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(1111);
  //  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  ////////////////////////////////////////////////////////////

  // rainbow palette
  gStyle->SetPalette(1);

//   gStyle->SetOptStat(1111111);  // Set stat options
//   gStyle->SetStatY(0.9);  // Set y-position (fraction of pad size)
//   gStyle->SetStatX(0.9);  // Set x-position (fraction of pad size)
//   gStyle->SetStatW(0.4);  // Set width of stat-box (fraction of pad size)
//   gStyle->SetStatH(0.2);  // Set height of stat-box (fraction of pad size)

///////////////////////////////////////////////////////////////////////
  /* Colors:
     1 Black
     2 red
     3 light green
     4 blue
     5 yellow
     6 magenta
     7 cyan
     8 dark green
     9 dark blue
     10 clear
     11 grey
     12 black

     28 brown

     46 orange-red-brown
  */

  int colorfeb = 1;
  int colorcell = 2;
  int colortopocl =3;
  int colortopoclps = 4;
  int colortopoclpuc = 6;
  int colorMHT = 7;
  int colorl1 = 30;
  int colorrefwmu = 44;
  int coloroffrecal =41;

  float xminstatbox = 0.75;
  float xmaxstatbox = 0.9;
  //float deltastatboxy = 0.08;
  //float minstatboxy = 0.18;
  float deltastatboxy = 0.08;
  float minstatboxy = 0.26;
  int linewidth =3;

  int nmetbins=100;
  Double_t metmin=0.0;
  Double_t metmax=200.0;

  //int nmetbins=70;
  //Double_t metmin=0.0;
  //Double_t metmax=140.0;

  double zval;

  TString PlotCut("passrndm>0.0");
  //TString PlotCut("passrndm>0.0 && mettopoclpuc>0.1");
  //TString PlotCut("passrndm>0.0 && metl1>30.0");

  //char* mettotTitle = "Runs 282712 282992 283429 283780 284213 284285  HLT_noalg_zb_L1ZB";
  //char* mettotTitle = "Period J (19 Runs) HLT_noalg_zb_L1ZB";
  //char* mettotTitle = "Period J (19 Runs) HLT_j40_L1ZB"; 

  char* mettotTitle = "Period J Runs 283429 284285 HLT_noalg_zb_L1ZB";
  char* intplotTitle = "Period J Runs 283429 284285 HLT_noalg_zb_L1ZB";

  //char* mettotTitle = "Period J Runs 283429 284285 HLT_noalg_zb_L1ZB L1 MET > 30";
  //char* intplotTitle = "Period J Run 283429 284285 HLT_noalg_zb_L1ZB L1 MET > 30";

  //char* mettotTitle = "Period J Run 284285 HLT_noalg_zb_L1ZB TopoclPUC MET > 0.1";
  //char* intplotTitle = "Period J Run 284285 HLT_noalg_zb_L1ZB TopoclPUC MET > 0.1";

  //TFile *File1 = TFile::Open("Run284285.ZeroBias.ntuples.root");
  TFile *File1 = TFile::Open("Run283429R284285.ZeroBias.ntuples.root");

  TCanvas *METCELLVSL1 = new TCanvas("METCELLVSL1","METCELLVSL1");

    TH2F* metTotCellvsL1 = new TH2F("metTotCellvsL1"," Cell vs L1 MET", 100,0.,100,100,0.,100);
    metTotCellvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotCellvsL1->GetYaxis()->SetTitle("Cell MET [GeV]");
    tree->Draw("metcell:metl1>>metTotCellvsL1",PlotCut,"colz");

  TCanvas *METTOPOCLVSL1 = new TCanvas("METTOPOCLVSL1","METTOPOCLVSL1");

    TH2F* metTotTopoclvsL1 = new TH2F("metTotTopoclvsL1"," Topocl vs L1 MET", 100,0.,100,100,0.,100);
    metTotTopoclvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotTopoclvsL1->GetYaxis()->SetTitle("Topocl MET [GeV]");
    tree->Draw("mettopocl:metl1>>metTotTopoclvsL1",PlotCut,"colz");

  TCanvas *METPSVSL1 = new TCanvas("METPSVSL1","METPSVSL1");

    TH2F* metTotPSvsL1 = new TH2F("metTotPSvsL1"," TopoclPS vs L1 MET", 100,0.,100,100,0.,100);
    metTotPSvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotPSvsL1->GetYaxis()->SetTitle("TopoclPS MET [GeV]");
    tree->Draw("mettopoclps:metl1>>metTotPSvsL1",PlotCut,"colz");

  TCanvas *METPUCVSL1 = new TCanvas("METPUCVSL1","METPUCVSL1");

    TH2F* metTotPUCvsL1 = new TH2F("metTotPUCvsL1"," TopoclPUC vs L1 MET", 100,0.,100,100,0.,100);
    metTotPUCvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotPUCvsL1->GetYaxis()->SetTitle("TopoclPUC MET [GeV]");
    tree->Draw("mettopoclpuc:metl1>>metTotPUCvsL1",PlotCut,"colz");

  TCanvas *METMHTVSL1 = new TCanvas("METMHTVSL1","METMHTVSL1");

    TH2F* metTotMHTvsL1 = new TH2F("metTotMHTvsL1"," MHT vs L1 MET", 100,0.,100,100,0.,100);
    metTotMHTvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotMHTvsL1->GetYaxis()->SetTitle("MHT MET [GeV]");
    tree->Draw("metmht:metl1>>metTotMHTvsL1",PlotCut,"colz");

  TCanvas *METOFFRECALVSL1 = new TCanvas("METOFFRECALVSL1","METOFFRECALVSL1");

    TH2F* metTotOFFRECALvsL1 = new TH2F("metTotOFFRECALvsL1"," OFFRECAL vs L1 MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotOFFRECALvsL1->GetYaxis()->SetTitle("OFFRECAL MET [GeV]");
    tree->Draw("metoffrecal:metl1>>metTotOFFRECALvsL1",PlotCut,"colz");

  TCanvas *METTOPOCLVSCELL = new TCanvas("METTOPOCLVSCELL","METTOPOCLVSCELL");

    TH2F* metTotTOPOCLvsCELL = new TH2F("metTotTOPOCLvsCELL"," TOPOCL vs CELL MET", 100,0.,100,100,0.,100);
    metTotTOPOCLvsCELL->GetXaxis()->SetTitle("CELL MET [GeV]");
    metTotTOPOCLvsCELL->GetYaxis()->SetTitle("TOPOCL MET [GeV]");
    tree->Draw("mettopocl:metcell>>metTotTOPOCLvsCELL",PlotCut,"colz");   

  TCanvas *METPSVSCELL = new TCanvas("METPSVSCELL","METPSVSCELL");

    TH2F* metTotPSvsCELL = new TH2F("metTotPSvsCELL"," TopoclPS vs CELL MET", 100,0.,100,100,0.,100);
    metTotPSvsCELL->GetXaxis()->SetTitle("CELL MET [GeV]");
    metTotPSvsCELL->GetYaxis()->SetTitle("PS MET [GeV]");
    tree->Draw("mettopocl:metcell>>metTotPSvsCELL",PlotCut,"colz"); 

  TCanvas *METPUCVSCELL = new TCanvas("PUCVSCELL","PUCVSCELL");

    TH2F* metTotPUCvsCELL = new TH2F("metTotPUCvsCELL"," TopoclPUC vs CELL MET", 100,0.,100,100,0.,100);
    metTotPUCvsCELL->GetXaxis()->SetTitle("CELL MET [GeV]");
    metTotPUCvsCELL->GetYaxis()->SetTitle("PUC MET [GeV]");
    tree->Draw("mettopoclpuc:metcell>>metTotPUCvsCELL",PlotCut,"colz");  

  TCanvas *METMHTVSCELL = new TCanvas("METMHTVSCELL","METMHTVSCELL");

    TH2F* metTotMHTvsCELL = new TH2F("metTotMHTvsCELL"," MHT vs CELL MET", 100,0.,100,100,0.,100);
    metTotMHTvsCELL->GetXaxis()->SetTitle("CELL MET [GeV]");
    metTotMHTvsCELL->GetYaxis()->SetTitle("MHT MET [GeV]");
    tree->Draw("metmht:metcell>>metTotMHTvsCELL",PlotCut,"colz"); 

  TCanvas *METOFFRECALVSCELL = new TCanvas("METOFFRECALVSCELL","METOFFRECALVSCELL");

    TH2F* metTotOFFRECALvsCELL = new TH2F("metTotOFFRECALvsCELL"," OFFRECAL vs CELL MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsCELL->GetXaxis()->SetTitle("CELL MET [GeV]");
    metTotOFFRECALvsCELL->GetYaxis()->SetTitle("OFFRECAL MET [GeV]");
    tree->Draw("metoffrecal:metcell>>metTotOFFRECALvsCELL",PlotCut,"colz"); 

  TCanvas *METPSVSTOPOCL = new TCanvas("METPSVSTOPOCL","METPSVSTOPOCL");

    TH2F* metTotPSvsTOPOCL = new TH2F("metTotPSvsTOPOCL"," TopoclPS vs TOPOCL MET", 100,0.,100,100,0.,100);
    metTotPSvsTOPOCL->GetXaxis()->SetTitle("PS MET [GeV]");
    metTotPSvsTOPOCL->GetYaxis()->SetTitle("TOPOCL MET [GeV]");
    tree->Draw("mettopoclps:mettopocl>>metTotPSvsTOPOCL",PlotCut,"colz"); 

  TCanvas *METPUCVSTOPOCL = new TCanvas("PUCVSTOPOCL","PUCVSTOPOCL");

    TH2F* metTotPUCvsTOPOCL = new TH2F("metTotPUCvsTOPOCL"," TopoclPUC vs TOPOCL MET", 100,0.,100,100,0.,100);
    metTotPUCvsTOPOCL->GetXaxis()->SetTitle("PUC MET [GeV]");
    metTotPUCvsTOPOCL->GetYaxis()->SetTitle("TOPOCL MET [GeV]");
    tree->Draw("mettopoclpuc:mettopocl>>metTotPUCvsTOPOCL",PlotCut,"colz"); 

  TCanvas *METMHTVSTOPOCL = new TCanvas("METMHTVSTOPOCL","METMHTVSTOPOCL");

    TH2F* metTotMHTvsTOPOCL = new TH2F("metTotMHTvsTOPOCL"," MHT vs TOPOCL MET", 100,0.,100,100,0.,100);
    metTotMHTvsTOPOCL->GetXaxis()->SetTitle("MHT MET [GeV]");
    metTotMHTvsTOPOCL->GetYaxis()->SetTitle("TOPOCL MET [GeV]");
    tree->Draw("metmht:mettopocl>>metTotMHTvsTOPOCL",PlotCut,"colz"); 

  TCanvas *METOFFRECALVSTOPOCL = new TCanvas("METOFFRECALVSTOPOCL","METOFFRECALVSTOPOCL");

    TH2F* metTotOFFRECALvsTOPOCL = new TH2F("metTotOFFRECALvsTOPOCL"," OFFRECAL vs TOPOCL MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsTOPOCL->GetXaxis()->SetTitle("OFFRECAL MET [GeV]");
    metTotOFFRECALvsTOPOCL->GetYaxis()->SetTitle("TOPOCL MET [GeV]");
    tree->Draw("metoffrecal:mettopocl>>metTotOFFRECALvsTOPOCL",PlotCut,"colz"); 

  TCanvas *METPUCVSPS = new TCanvas("METPUCVSPS","METPUCVSPS");

    TH2F* metTotPUCvsPS = new TH2F("metTotPUCvsPS"," TopoclPUC vs TopoclPS MET", 100,0.,100,100,0.,100);
    metTotPUCvsPS->GetXaxis()->SetTitle("PUC MET [GeV]");
    metTotPUCvsPS->GetYaxis()->SetTitle("PS MET [GeV]");
    tree->Draw("mettopoclpuc:mettopoclps>>metTotPUCvsPS",PlotCut,"colz");

  TCanvas *METMHTVSPS = new TCanvas("METMHTVSPS","METMHTVSPS");

    TH2F* metTotMHTvsPS = new TH2F("metTotMHTvsPS"," MHT vs TopoclPS MET", 100,0.,100,100,0.,100);
    metTotMHTvsPS->GetXaxis()->SetTitle("MHT MET [GeV]");
    metTotMHTvsPS->GetYaxis()->SetTitle("PS MET [GeV]");
    tree->Draw("metmht:mettopoclps>>metTotMHTvsPS",PlotCut,"colz"); 

  TCanvas *METOFFRECALVSPS = new TCanvas("METOFFRECALVSPS","METOFFRECALVSPS");

    TH2F* metTotOFFRECALvsPS = new TH2F("metTotOFFRECALvsPS"," OFFRECAL vs TopoclPS MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsPS->GetXaxis()->SetTitle("OFFRECAL MET [GeV]");
    metTotOFFRECALvsPS->GetYaxis()->SetTitle("PS MET [GeV]");
    tree->Draw("metoffrecal:mettopoclps>>metTotOFFRECALvsPS",PlotCut,"colz"); 

   TCanvas *METMHTVSPUC = new TCanvas("METMHTVSPUC","METMHTVSPUC");

    TH2F* metTotMHTvsPUC = new TH2F("metTotMHTvsPUC"," MHT vs TopoclPUC MET", 100,0.,100,100,0.,100);
    metTotMHTvsPUC->GetXaxis()->SetTitle("MHT MET [GeV]");
    metTotMHTvsPUC->GetYaxis()->SetTitle("PUC MET [GeV]");
    tree->Draw("metmht:mettopoclpuc>>metTotMHTvsPUC",PlotCut,"colz"); 

  TCanvas *METOFFRECALVSPUC = new TCanvas("METOFFRECALVSPUC","METOFFRECALVSPUC");

    TH2F* metTotOFFRECALvsPUC = new TH2F("metTotOFFRECALvsPUC"," OFFRECAL vs TopoclPUC MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsPUC->GetXaxis()->SetTitle("OFFRECAL MET [GeV]");
    metTotOFFRECALvsPUC->GetYaxis()->SetTitle("PUC MET [GeV]");
    tree->Draw("metoffrecal:mettopoclpuc>>metTotOFFRECALvsPUC",PlotCut,"colz"); 

  TCanvas *METOFFRECALVSMHT = new TCanvas("METOFFRECALVSMHT","METOFFRECALVSMHT");

    TH2F* metTotOFFRECALvsMHT = new TH2F("metTotOFFRECALvsMHT"," OFFRECAL vs MHT MET", 100,0.,100,100,0.,100);
    metTotOFFRECALvsMHT->GetXaxis()->SetTitle("OFFRECAL MET [GeV]");
    metTotOFFRECALvsMHT->GetYaxis()->SetTitle("MHT MET [GeV]");
    tree->Draw("metoffrecal:metmht>>metTotOFFRECALvsMHT",PlotCut,"colz"); 

}