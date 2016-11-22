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

  /*
  int colorfeb = 1;
  int colorcell = 2;
  int colortopocl =3;
  int colortopoclps = 4;
  int colortopoclpuc = 6;
  int colorMHT = 7;
  int colorl1 = 30;
  int colorrefwmu = 44;
  int coloroffrecal =41;
  */
 
  int colorR1 = 1;
  int colorR2 = 2;
  int colorR3 = 3;
  int colorR4 = 4;
  int colorR5 = 6;
  int colorR6 = 44;

  float xminstatbox = 0.75;
  float xmaxstatbox = 0.9;
  float deltastatboxy = 0.08;
  //float minstatboxy = 0.18;
  float minstatboxy = 0.42;
  float minstatboxyR6 = minstatboxy;
  float maxstatboxyR6 = minstatboxy + deltastatboxy;
  float minstatboxyR5 = maxstatboxyR6;
  float maxstatboxyR5 = minstatboxyR5 + deltastatboxy;
  float minstatboxyR4 = maxstatboxyR5;
  float maxstatboxyR4 = minstatboxyR4 + deltastatboxy;
  float minstatboxyR3 = maxstatboxyR4;
  float maxstatboxyR3 = minstatboxyR3 + deltastatboxy;
  float minstatboxyR2 = maxstatboxyR3;
  float maxstatboxyR2 = minstatboxyR2 + deltastatboxy;
  float minstatboxyR1 = maxstatboxyR2;
  float maxstatboxyR1 = minstatboxyR1 + deltastatboxy;

  int linewidth =3;

  int nmetbins=100;
  Double_t metmin=0.0;
  Double_t metmax=200.0;

  char* mettotTitle = "Run XXX HLT_noalg_zb_L1ZB";
  char* mettotyaxisTitle = "Fraction of Events / 2 GeV";

  // First File
  TFile *FileR1 = TFile::Open("DperiodR276262R276689.root");

  TCanvas *METL1R1 = new TCanvas("METL1R1","METL1R1");

    TH1F *metTotL1R1=(TH1F*)metTotL1->Clone();
    metTotL1R1->SetName("metTotL1R1");
    metTotL1R1->SetDirectory(0);

    metTotL1R1->SetLineColor(colorR1);
    metTotL1R1->SetLineWidth(linewidth);
    metTotL1R1->SetTitle("L1 MET");
    metTotL1R1->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R1->SetNormFactor(1.0);
    metTotL1R1->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R1 = (TPaveStats*)metTotL1R1->GetListOfFunctions()->FindObject("stats");
    metTotL1R1->GetListOfFunctions()->Remove(statboxmetl1R1);
    metTotL1R1->SetStats(0);
    statboxmetl1R1->SetTextColor(colorR1);
    statboxmetl1R1->SetY1NDC(minstatboxyR1);
    statboxmetl1R1->SetY2NDC(maxstatboxyR1);
    statboxmetl1R1->SetX1NDC(xminstatbox);
    statboxmetl1R1->SetX2NDC(xmaxstatbox);
    statboxmetl1R1->Draw();

  TCanvas *METL2FEBR1 = new TCanvas("METL2FEBR1","METL2FEBR1");

    TH1F *metTotL2febR1=(TH1F*)metTotL2feb->Clone();
    metTotL2febR1->SetName("metTotL2febR1");
    metTotL2febR1->SetDirectory(0);

    metTotL2febR1->SetLineColor(colorR1);
    metTotL2febR1->SetLineWidth(linewidth);
    metTotL2febR1->SetTitle("FEB MET");
    metTotL2febR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR1->SetNormFactor(1.0);
    metTotL2febR1->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR1 = (TPaveStats*)metTotL2febR1->GetListOfFunctions()->FindObject("stats");
    metTotL2febR1->GetListOfFunctions()->Remove(statboxmetl2febR1);
    metTotL2febR1->SetStats(0);
    statboxmetl2febR1->SetTextColor(colorR1);
    statboxmetl2febR1->SetY1NDC(minstatboxyR1);
    statboxmetl2febR1->SetY2NDC(maxstatboxyR1);
    statboxmetl2febR1->SetX1NDC(xminstatbox);
    statboxmetl2febR1->SetX2NDC(xmaxstatbox);
    statboxmetl2febR1->Draw();

  TCanvas *METCELLR1 = new TCanvas("METCELLR1","METCELLR1");

    TH1F *metTotCellR1=(TH1F*)metTotCell->Clone();
    metTotCellR1->SetName("metTotCellR1");
    metTotCellR1->SetDirectory(0);

    metTotCellR1->SetLineColor(colorR1);
    metTotCellR1->SetLineWidth(linewidth);
    metTotCellR1->SetTitle("CELL MET");
    metTotCellR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR1->SetNormFactor(1.0);
    metTotCellR1->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR1 = (TPaveStats*)metTotCellR1->GetListOfFunctions()->FindObject("stats");
    metTotCellR1->GetListOfFunctions()->Remove(statboxmetcellR1);
    metTotCellR1->SetStats(0);
    statboxmetcellR1->SetTextColor(colorR1);
    statboxmetcellR1->SetY1NDC(minstatboxyR1);
    statboxmetcellR1->SetY2NDC(maxstatboxyR1);
    statboxmetcellR1->SetX1NDC(xminstatbox);
    statboxmetcellR1->SetX2NDC(xmaxstatbox);
    statboxmetcellR1->Draw();

  TCanvas *METTOPOCLR1 = new TCanvas("METTOPOCLR1","METTOPOCLR1");

    TH1F *metTotTopoclR1=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR1->SetName("metTotTopoclR1");
    metTotTopoclR1->SetDirectory(0);

    metTotTopoclR1->SetLineColor(colorR1);
    metTotTopoclR1->SetLineWidth(linewidth);
    metTotTopoclR1->SetTitle("TOPOCL MET");
    metTotTopoclR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR1->SetNormFactor(1.0);
    metTotTopoclR1->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR1 = (TPaveStats*)metTotTopoclR1->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR1->GetListOfFunctions()->Remove(statboxmettopoclR1);
    metTotTopoclR1->SetStats(0);
    statboxmettopoclR1->SetTextColor(colorR1);
    statboxmettopoclR1->SetY1NDC(minstatboxyR1);
    statboxmettopoclR1->SetY2NDC(maxstatboxyR1);
    statboxmettopoclR1->SetX1NDC(xminstatbox);
    statboxmettopoclR1->SetX2NDC(xmaxstatbox);
    statboxmettopoclR1->Draw();

  TCanvas *METTOPOCLPSR1 = new TCanvas("METTOPOCLPSR1","METTOPOCLPSR1");

    TH1F *metTotTopoclPSR1=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR1->SetName("metTotTopoclPSR1");
    metTotTopoclPSR1->SetDirectory(0);

    metTotTopoclPSR1->SetLineColor(colorR1);
    metTotTopoclPSR1->SetLineWidth(linewidth);
    metTotTopoclPSR1->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR1->SetNormFactor(1.0);
    metTotTopoclPSR1->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR1 = (TPaveStats*)metTotTopoclPSR1->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR1->GetListOfFunctions()->Remove(statboxmettopoclpsR1);
    metTotTopoclPSR1->SetStats(0);
    statboxmettopoclpsR1->SetTextColor(colorR1);
    statboxmettopoclpsR1->SetY1NDC(minstatboxyR1);
    statboxmettopoclpsR1->SetY2NDC(maxstatboxyR1);
    statboxmettopoclpsR1->SetX1NDC(xminstatbox);
    statboxmettopoclpsR1->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR1->Draw();

  TCanvas *METTOPOCLPUCR1 = new TCanvas("METTOPOCLPUCR1","METTOPOCLPUCR1");

    TH1F *metTotTopoclPUCR1=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR1->SetName("metTotTopoclPUCR1");
    metTotTopoclPUCR1->SetDirectory(0);

    metTotTopoclPUCR1->SetLineColor(colorR1);
    metTotTopoclPUCR1->SetLineWidth(linewidth);
    metTotTopoclPUCR1->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR1->SetNormFactor(1.0);
    metTotTopoclPUCR1->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR1 = (TPaveStats*)metTotTopoclPUCR1->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR1->GetListOfFunctions()->Remove(statboxmettopoclpucR1);
    metTotTopoclPUCR1->SetStats(0);
    statboxmettopoclpucR1->SetTextColor(colorR1);
    statboxmettopoclpucR1->SetY1NDC(minstatboxyR1);
    statboxmettopoclpucR1->SetY2NDC(maxstatboxyR1);
    statboxmettopoclpucR1->SetX1NDC(xminstatbox);
    statboxmettopoclpucR1->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR1->Draw();

  TCanvas *METMHTR1 = new TCanvas("METMHTR1","METMHTR1");

    TH1F *metTotMHTR1=(TH1F*)metTotMHT->Clone();
    metTotMHTR1->SetName("metTotMHTR1");
    metTotMHTR1->SetDirectory(0);

    metTotMHTR1->SetLineColor(colorR1);
    metTotMHTR1->SetLineWidth(linewidth);
    metTotMHTR1->SetTitle("MHT MET");
    metTotMHTR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR1->SetNormFactor(1.0);
    metTotMHTR1->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR1 = (TPaveStats*)metTotMHTR1->GetListOfFunctions()->FindObject("stats");
    metTotMHTR1->GetListOfFunctions()->Remove(statboxmetmhtR1);
    metTotMHTR1->SetStats(0);
    statboxmetmhtR1->SetTextColor(colorR1);
    statboxmetmhtR1->SetY1NDC(minstatboxyR1);
    statboxmetmhtR1->SetY2NDC(maxstatboxyR1);
    statboxmetmhtR1->SetX1NDC(xminstatbox);
    statboxmetmhtR1->SetX2NDC(xmaxstatbox);
    statboxmetmhtR1->Draw();

  TCanvas *METREFWMUR1 = new TCanvas("METREFWMUR1","METREFWMUR1");

    TH1F *metTotRefwmuR1=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR1->SetName("metTotRefwmuR1");
    metTotRefwmuR1->SetDirectory(0);

    metTotRefwmuR1->SetLineColor(colorR1);
    metTotRefwmuR1->SetLineWidth(linewidth);
    metTotRefwmuR1->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR1->SetNormFactor(1.0);
    metTotRefwmuR1->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR1 = (TPaveStats*)metTotRefwmuR1->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR1->GetListOfFunctions()->Remove(statboxmetrefwmuR1);
    metTotRefwmuR1->SetStats(0);
    statboxmetrefwmuR1->SetTextColor(colorR1);
    statboxmetrefwmuR1->SetY1NDC(minstatboxyR1);
    statboxmetrefwmuR1->SetY2NDC(maxstatboxyR1);
    statboxmetrefwmuR1->SetX1NDC(xminstatbox);
    statboxmetrefwmuR1->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR1->Draw();

  TCanvas *METOFFRECALR1 = new TCanvas("METOFFRECALR1","METOFFRECALR1");

    TH1F *metTotOffRecalR1=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR1->SetName("metTotOffRecalR1");
    metTotOffRecalR1->SetDirectory(0);

    metTotOffRecalR1->SetLineColor(colorR1);
    metTotOffRecalR1->SetLineWidth(linewidth);
    metTotOffRecalR1->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR1->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR1->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR1->SetNormFactor(1.0);
    metTotOffRecalR1->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR1 = (TPaveStats*)metTotOffRecalR1->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR1->GetListOfFunctions()->Remove(statboxmetoffrecalR1);
    metTotOffRecalR1->SetStats(0);
    statboxmetoffrecalR1->SetTextColor(colorR1);
    statboxmetoffrecalR1->SetY1NDC(minstatboxyR1);
    statboxmetoffrecalR1->SetY2NDC(maxstatboxyR1);
    statboxmetoffrecalR1->SetX1NDC(xminstatbox);
    statboxmetoffrecalR1->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR1->Draw();

  // Second File
  TFile *FileR2 = TFile::Open("Run279169.physics_ZeroBias.merge.AOD.r7562_p2521.root");

  TCanvas *METL1R2 = new TCanvas("METL1R2","METL1R2");

    TH1F *metTotL1R2=(TH1F*)metTotL1->Clone();
    metTotL1R2->SetName("metTotL1R2");
    metTotL1R2->SetDirectory(0);

    metTotL1R2->SetLineColor(colorR2);
    metTotL1R2->SetLineWidth(linewidth);
    metTotL1R2->SetTitle("L1 MET");
    metTotL1R2->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R2->SetNormFactor(1.0);
    metTotL1R2->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R2 = (TPaveStats*)metTotL1R2->GetListOfFunctions()->FindObject("stats");
    metTotL1R2->GetListOfFunctions()->Remove(statboxmetl1R2);
    metTotL1R2->SetStats(0);
    statboxmetl1R2->SetTextColor(colorR2);
    statboxmetl1R2->SetY1NDC(minstatboxyR2);
    statboxmetl1R2->SetY2NDC(maxstatboxyR2);
    statboxmetl1R2->SetX1NDC(xminstatbox);
    statboxmetl1R2->SetX2NDC(xmaxstatbox);
    statboxmetl1R2->Draw();

  TCanvas *METL2FEBR2 = new TCanvas("METL2FEBR2","METL2FEBR2");

    TH1F *metTotL2febR2=(TH1F*)metTotL2feb->Clone();
    metTotL2febR2->SetName("metTotL2febR2");
    metTotL2febR2->SetDirectory(0);

    metTotL2febR2->SetLineColor(colorR2);
    metTotL2febR2->SetLineWidth(linewidth);
    metTotL2febR2->SetTitle("FEB MET");
    metTotL2febR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR2->SetNormFactor(1.0);
    metTotL2febR2->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR2 = (TPaveStats*)metTotL2febR2->GetListOfFunctions()->FindObject("stats");
    metTotL2febR2->GetListOfFunctions()->Remove(statboxmetl2febR2);
    metTotL2febR2->SetStats(0);
    statboxmetl2febR2->SetTextColor(colorR2);
    statboxmetl2febR2->SetY1NDC(minstatboxyR2);
    statboxmetl2febR2->SetY2NDC(maxstatboxyR2);
    statboxmetl2febR2->SetX1NDC(xminstatbox);
    statboxmetl2febR2->SetX2NDC(xmaxstatbox);
    statboxmetl2febR2->Draw();

  TCanvas *METCELLR2 = new TCanvas("METCELLR2","METCELLR2");

    TH1F *metTotCellR2=(TH1F*)metTotCell->Clone();
    metTotCellR2->SetName("metTotCellR2");
    metTotCellR2->SetDirectory(0);

    metTotCellR2->SetLineColor(colorR2);
    metTotCellR2->SetLineWidth(linewidth);
    metTotCellR2->SetTitle("CELL MET");
    metTotCellR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR2->SetNormFactor(1.0);
    metTotCellR2->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR2 = (TPaveStats*)metTotCellR2->GetListOfFunctions()->FindObject("stats");
    metTotCellR2->GetListOfFunctions()->Remove(statboxmetcellR2);
    metTotCellR2->SetStats(0);
    statboxmetcellR2->SetTextColor(colorR2);
    statboxmetcellR2->SetY1NDC(minstatboxyR2);
    statboxmetcellR2->SetY2NDC(maxstatboxyR2);
    statboxmetcellR2->SetX1NDC(xminstatbox);
    statboxmetcellR2->SetX2NDC(xmaxstatbox);
    statboxmetcellR2->Draw();

  TCanvas *METTOPOCLR2 = new TCanvas("METTOPOCLR2","METTOPOCLR2");

    TH1F *metTotTopoclR2=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR2->SetName("metTotTopoclR2");
    metTotTopoclR2->SetDirectory(0);

    metTotTopoclR2->SetLineColor(colorR2);
    metTotTopoclR2->SetLineWidth(linewidth);
    metTotTopoclR2->SetTitle("TOPOCL MET");
    metTotTopoclR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR2->SetNormFactor(1.0);
    metTotTopoclR2->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR2 = (TPaveStats*)metTotTopoclR2->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR2->GetListOfFunctions()->Remove(statboxmettopoclR2);
    metTotTopoclR2->SetStats(0);
    statboxmettopoclR2->SetTextColor(colorR2);
    statboxmettopoclR2->SetY1NDC(minstatboxyR2);
    statboxmettopoclR2->SetY2NDC(maxstatboxyR2);
    statboxmettopoclR2->SetX1NDC(xminstatbox);
    statboxmettopoclR2->SetX2NDC(xmaxstatbox);
    statboxmettopoclR2->Draw();

  TCanvas *METTOPOCLPSR2 = new TCanvas("METTOPOCLPSR2","METTOPOCLPSR2");

    TH1F *metTotTopoclPSR2=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR2->SetName("metTotTopoclPSR2");
    metTotTopoclPSR2->SetDirectory(0);

    metTotTopoclPSR2->SetLineColor(colorR2);
    metTotTopoclPSR2->SetLineWidth(linewidth);
    metTotTopoclPSR2->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR2->SetNormFactor(1.0);
    metTotTopoclPSR2->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR2 = (TPaveStats*)metTotTopoclPSR2->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR2->GetListOfFunctions()->Remove(statboxmettopoclpsR2);
    metTotTopoclPSR2->SetStats(0);
    statboxmettopoclpsR2->SetTextColor(colorR2);
    statboxmettopoclpsR2->SetY1NDC(minstatboxyR2);
    statboxmettopoclpsR2->SetY2NDC(maxstatboxyR2);
    statboxmettopoclpsR2->SetX1NDC(xminstatbox);
    statboxmettopoclpsR2->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR2->Draw();

  TCanvas *METTOPOCLPUCR2 = new TCanvas("METTOPOCLPUCR2","METTOPOCLPUCR2");

    TH1F *metTotTopoclPUCR2=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR2->SetName("metTotTopoclPUCR2");
    metTotTopoclPUCR2->SetDirectory(0);

    metTotTopoclPUCR2->SetLineColor(colorR2);
    metTotTopoclPUCR2->SetLineWidth(linewidth);
    metTotTopoclPUCR2->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR2->SetNormFactor(1.0);
    metTotTopoclPUCR2->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR2 = (TPaveStats*)metTotTopoclPUCR2->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR2->GetListOfFunctions()->Remove(statboxmettopoclpucR2);
    metTotTopoclPUCR2->SetStats(0);
    statboxmettopoclpucR2->SetTextColor(colorR2);
    statboxmettopoclpucR2->SetY1NDC(minstatboxyR2);
    statboxmettopoclpucR2->SetY2NDC(maxstatboxyR2);
    statboxmettopoclpucR2->SetX1NDC(xminstatbox);
    statboxmettopoclpucR2->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR2->Draw();

  TCanvas *METMHTR2 = new TCanvas("METMHTR2","METMHTR2");

    TH1F *metTotMHTR2=(TH1F*)metTotMHT->Clone();
    metTotMHTR2->SetName("metTotMHTR2");
    metTotMHTR2->SetDirectory(0);

    metTotMHTR2->SetLineColor(colorR2);
    metTotMHTR2->SetLineWidth(linewidth);
    metTotMHTR2->SetTitle("MHT MET");
    metTotMHTR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR2->SetNormFactor(1.0);
    metTotMHTR2->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR2 = (TPaveStats*)metTotMHTR2->GetListOfFunctions()->FindObject("stats");
    metTotMHTR2->GetListOfFunctions()->Remove(statboxmetmhtR2);
    metTotMHTR2->SetStats(0);
    statboxmetmhtR2->SetTextColor(colorR2);
    statboxmetmhtR2->SetY1NDC(minstatboxyR2);
    statboxmetmhtR2->SetY2NDC(maxstatboxyR2);
    statboxmetmhtR2->SetX1NDC(xminstatbox);
    statboxmetmhtR2->SetX2NDC(xmaxstatbox);
    statboxmetmhtR2->Draw();

  TCanvas *METREFWMUR2 = new TCanvas("METREFWMUR2","METREFWMUR2");

    TH1F *metTotRefwmuR2=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR2->SetName("metTotRefwmuR2");
    metTotRefwmuR2->SetDirectory(0);

    metTotRefwmuR2->SetLineColor(colorR2);
    metTotRefwmuR2->SetLineWidth(linewidth);
    metTotRefwmuR2->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR2->SetNormFactor(1.0);
    metTotRefwmuR2->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR2 = (TPaveStats*)metTotRefwmuR2->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR2->GetListOfFunctions()->Remove(statboxmetrefwmuR2);
    metTotRefwmuR2->SetStats(0);
    statboxmetrefwmuR2->SetTextColor(colorR2);
    statboxmetrefwmuR2->SetY1NDC(minstatboxyR2);
    statboxmetrefwmuR2->SetY2NDC(maxstatboxyR2);
    statboxmetrefwmuR2->SetX1NDC(xminstatbox);
    statboxmetrefwmuR2->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR2->Draw();

  TCanvas *METOFFRECALR2 = new TCanvas("METOFFRECALR2","METOFFRECALR2");

    TH1F *metTotOffRecalR2=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR2->SetName("metTotOffRecalR2");
    metTotOffRecalR2->SetDirectory(0);

    metTotOffRecalR2->SetLineColor(colorR2);
    metTotOffRecalR2->SetLineWidth(linewidth);
    metTotOffRecalR2->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR2->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR2->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR2->SetNormFactor(1.0);
    metTotOffRecalR2->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR2 = (TPaveStats*)metTotOffRecalR2->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR2->GetListOfFunctions()->Remove(statboxmetoffrecalR2);
    metTotOffRecalR2->SetStats(0);
    statboxmetoffrecalR2->SetTextColor(colorR2);
    statboxmetoffrecalR2->SetY1NDC(minstatboxyR2);
    statboxmetoffrecalR2->SetY2NDC(maxstatboxyR2);
    statboxmetoffrecalR2->SetX1NDC(xminstatbox);
    statboxmetoffrecalR2->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR2->Draw();

  // Third File
  TFile *FileR3 = TFile::Open("Run280231.physics_ZeroBias.merge.AOD.r7562_p2521.root");

  TCanvas *METL1R3 = new TCanvas("METL1R3","METL1R3");

    TH1F *metTotL1R3=(TH1F*)metTotL1->Clone();
    metTotL1R3->SetName("metTotL1R3");
    metTotL1R3->SetDirectory(0);

    metTotL1R3->SetLineColor(colorR3);
    metTotL1R3->SetLineWidth(linewidth);
    metTotL1R3->SetTitle("L1 MET");
    metTotL1R3->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R3->SetNormFactor(1.0);
    metTotL1R3->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R3 = (TPaveStats*)metTotL1R3->GetListOfFunctions()->FindObject("stats");
    metTotL1R3->GetListOfFunctions()->Remove(statboxmetl1R3);
    metTotL1R3->SetStats(0);
    statboxmetl1R3->SetTextColor(colorR3);
    statboxmetl1R3->SetY1NDC(minstatboxyR3);
    statboxmetl1R3->SetY2NDC(maxstatboxyR3);
    statboxmetl1R3->SetX1NDC(xminstatbox);
    statboxmetl1R3->SetX2NDC(xmaxstatbox);
    statboxmetl1R3->Draw();

  TCanvas *METL2FEBR3 = new TCanvas("METL2FEBR3","METL2FEBR3");

    TH1F *metTotL2febR3=(TH1F*)metTotL2feb->Clone();
    metTotL2febR3->SetName("metTotL2febR3");
    metTotL2febR3->SetDirectory(0);

    metTotL2febR3->SetLineColor(colorR3);
    metTotL2febR3->SetLineWidth(linewidth);
    metTotL2febR3->SetTitle("FEB MET");
    metTotL2febR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR3->SetNormFactor(1.0);
    metTotL2febR3->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR3 = (TPaveStats*)metTotL2febR3->GetListOfFunctions()->FindObject("stats");
    metTotL2febR3->GetListOfFunctions()->Remove(statboxmetl2febR3);
    metTotL2febR3->SetStats(0);
    statboxmetl2febR3->SetTextColor(colorR3);
    statboxmetl2febR3->SetY1NDC(minstatboxyR3);
    statboxmetl2febR3->SetY2NDC(maxstatboxyR3);
    statboxmetl2febR3->SetX1NDC(xminstatbox);
    statboxmetl2febR3->SetX2NDC(xmaxstatbox);
    statboxmetl2febR3->Draw();

  TCanvas *METCELLR3 = new TCanvas("METCELLR3","METCELLR3");

    TH1F *metTotCellR3=(TH1F*)metTotCell->Clone();
    metTotCellR3->SetName("metTotCellR3");
    metTotCellR3->SetDirectory(0);

    metTotCellR3->SetLineColor(colorR3);
    metTotCellR3->SetLineWidth(linewidth);
    metTotCellR3->SetTitle("CELL MET");
    metTotCellR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR3->SetNormFactor(1.0);
    metTotCellR3->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR3 = (TPaveStats*)metTotCellR3->GetListOfFunctions()->FindObject("stats");
    metTotCellR3->GetListOfFunctions()->Remove(statboxmetcellR3);
    metTotCellR3->SetStats(0);
    statboxmetcellR3->SetTextColor(colorR3);
    statboxmetcellR3->SetY1NDC(minstatboxyR3);
    statboxmetcellR3->SetY2NDC(maxstatboxyR3);
    statboxmetcellR3->SetX1NDC(xminstatbox);
    statboxmetcellR3->SetX2NDC(xmaxstatbox);
    statboxmetcellR3->Draw();

  TCanvas *METTOPOCLR3 = new TCanvas("METTOPOCLR3","METTOPOCLR3");

    TH1F *metTotTopoclR3=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR3->SetName("metTotTopoclR3");
    metTotTopoclR3->SetDirectory(0);

    metTotTopoclR3->SetLineColor(colorR3);
    metTotTopoclR3->SetLineWidth(linewidth);
    metTotTopoclR3->SetTitle("TOPOCL MET");
    metTotTopoclR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR3->SetNormFactor(1.0);
    metTotTopoclR3->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR3 = (TPaveStats*)metTotTopoclR3->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR3->GetListOfFunctions()->Remove(statboxmettopoclR3);
    metTotTopoclR3->SetStats(0);
    statboxmettopoclR3->SetTextColor(colorR3);
    statboxmettopoclR3->SetY1NDC(minstatboxyR3);
    statboxmettopoclR3->SetY2NDC(maxstatboxyR3);
    statboxmettopoclR3->SetX1NDC(xminstatbox);
    statboxmettopoclR3->SetX2NDC(xmaxstatbox);
    statboxmettopoclR3->Draw();

  TCanvas *METTOPOCLPSR3 = new TCanvas("METTOPOCLPSR3","METTOPOCLPSR3");

    TH1F *metTotTopoclPSR3=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR3->SetName("metTotTopoclPSR3");
    metTotTopoclPSR3->SetDirectory(0);

    metTotTopoclPSR3->SetLineColor(colorR3);
    metTotTopoclPSR3->SetLineWidth(linewidth);
    metTotTopoclPSR3->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR3->SetNormFactor(1.0);
    metTotTopoclPSR3->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR3 = (TPaveStats*)metTotTopoclPSR3->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR3->GetListOfFunctions()->Remove(statboxmettopoclpsR3);
    metTotTopoclPSR3->SetStats(0);
    statboxmettopoclpsR3->SetTextColor(colorR3);
    statboxmettopoclpsR3->SetY1NDC(minstatboxyR3);
    statboxmettopoclpsR3->SetY2NDC(maxstatboxyR3);
    statboxmettopoclpsR3->SetX1NDC(xminstatbox);
    statboxmettopoclpsR3->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR3->Draw();

  TCanvas *METTOPOCLPUCR3 = new TCanvas("METTOPOCLPUCR3","METTOPOCLPUCR3");

    TH1F *metTotTopoclPUCR3=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR3->SetName("metTotTopoclPUCR3");
    metTotTopoclPUCR3->SetDirectory(0);

    metTotTopoclPUCR3->SetLineColor(colorR3);
    metTotTopoclPUCR3->SetLineWidth(linewidth);
    metTotTopoclPUCR3->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR3->SetNormFactor(1.0);
    metTotTopoclPUCR3->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR3 = (TPaveStats*)metTotTopoclPUCR3->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR3->GetListOfFunctions()->Remove(statboxmettopoclpucR3);
    metTotTopoclPUCR3->SetStats(0);
    statboxmettopoclpucR3->SetTextColor(colorR3);
    statboxmettopoclpucR3->SetY1NDC(minstatboxyR3);
    statboxmettopoclpucR3->SetY2NDC(maxstatboxyR3);
    statboxmettopoclpucR3->SetX1NDC(xminstatbox);
    statboxmettopoclpucR3->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR3->Draw();

  TCanvas *METMHTR3 = new TCanvas("METMHTR3","METMHTR3");

    TH1F *metTotMHTR3=(TH1F*)metTotMHT->Clone();
    metTotMHTR3->SetName("metTotMHTR3");
    metTotMHTR3->SetDirectory(0);

    metTotMHTR3->SetLineColor(colorR3);
    metTotMHTR3->SetLineWidth(linewidth);
    metTotMHTR3->SetTitle("MHT MET");
    metTotMHTR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR3->SetNormFactor(1.0);
    metTotMHTR3->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR3 = (TPaveStats*)metTotMHTR3->GetListOfFunctions()->FindObject("stats");
    metTotMHTR3->GetListOfFunctions()->Remove(statboxmetmhtR3);
    metTotMHTR3->SetStats(0);
    statboxmetmhtR3->SetTextColor(colorR3);
    statboxmetmhtR3->SetY1NDC(minstatboxyR3);
    statboxmetmhtR3->SetY2NDC(maxstatboxyR3);
    statboxmetmhtR3->SetX1NDC(xminstatbox);
    statboxmetmhtR3->SetX2NDC(xmaxstatbox);
    statboxmetmhtR3->Draw();

  TCanvas *METREFWMUR3 = new TCanvas("METREFWMUR3","METREFWMUR3");

    TH1F *metTotRefwmuR3=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR3->SetName("metTotRefwmuR3");
    metTotRefwmuR3->SetDirectory(0);

    metTotRefwmuR3->SetLineColor(colorR3);
    metTotRefwmuR3->SetLineWidth(linewidth);
    metTotRefwmuR3->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR3->SetNormFactor(1.0);
    metTotRefwmuR3->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR3 = (TPaveStats*)metTotRefwmuR3->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR3->GetListOfFunctions()->Remove(statboxmetrefwmuR3);
    metTotRefwmuR3->SetStats(0);
    statboxmetrefwmuR3->SetTextColor(colorR3);
    statboxmetrefwmuR3->SetY1NDC(minstatboxyR3);
    statboxmetrefwmuR3->SetY2NDC(maxstatboxyR3);
    statboxmetrefwmuR3->SetX1NDC(xminstatbox);
    statboxmetrefwmuR3->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR3->Draw();

  TCanvas *METOFFRECALR3 = new TCanvas("METOFFRECALR3","METOFFRECALR3");

    TH1F *metTotOffRecalR3=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR3->SetName("metTotOffRecalR3");
    metTotOffRecalR3->SetDirectory(0);

    metTotOffRecalR3->SetLineColor(colorR3);
    metTotOffRecalR3->SetLineWidth(linewidth);
    metTotOffRecalR3->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR3->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR3->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR3->SetNormFactor(1.0);
    metTotOffRecalR3->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR3 = (TPaveStats*)metTotOffRecalR3->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR3->GetListOfFunctions()->Remove(statboxmetoffrecalR3);
    metTotOffRecalR3->SetStats(0);
    statboxmetoffrecalR3->SetTextColor(colorR3);
    statboxmetoffrecalR3->SetY1NDC(minstatboxyR3);
    statboxmetoffrecalR3->SetY2NDC(maxstatboxyR3);
    statboxmetoffrecalR3->SetX1NDC(xminstatbox);
    statboxmetoffrecalR3->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR3->Draw();

  // Fourth File
  TFile *FileR4 = TFile::Open("Run280862.physics_ZeroBias.merge.AOD.r7562_p2521.root");

  TCanvas *METL1R4 = new TCanvas("METL1R4","METL1R4");

    TH1F *metTotL1R4=(TH1F*)metTotL1->Clone();
    metTotL1R4->SetName("metTotL1R4");
    metTotL1R4->SetDirectory(0);

    metTotL1R4->SetLineColor(colorR4);
    metTotL1R4->SetLineWidth(linewidth);
    metTotL1R4->SetTitle("L1 MET");
    metTotL1R4->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R4->SetNormFactor(1.0);
    metTotL1R4->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R4 = (TPaveStats*)metTotL1R4->GetListOfFunctions()->FindObject("stats");
    metTotL1R4->GetListOfFunctions()->Remove(statboxmetl1R4);
    metTotL1R4->SetStats(0);
    statboxmetl1R4->SetTextColor(colorR4);
    statboxmetl1R4->SetY1NDC(minstatboxyR4);
    statboxmetl1R4->SetY2NDC(maxstatboxyR4);
    statboxmetl1R4->SetX1NDC(xminstatbox);
    statboxmetl1R4->SetX2NDC(xmaxstatbox);
    statboxmetl1R4->Draw();

  TCanvas *METL2FEBR4 = new TCanvas("METL2FEBR4","METL2FEBR4");

    TH1F *metTotL2febR4=(TH1F*)metTotL2feb->Clone();
    metTotL2febR4->SetName("metTotL2febR4");
    metTotL2febR4->SetDirectory(0);

    metTotL2febR4->SetLineColor(colorR4);
    metTotL2febR4->SetLineWidth(linewidth);
    metTotL2febR4->SetTitle("FEB MET");
    metTotL2febR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR4->SetNormFactor(1.0);
    metTotL2febR4->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR4 = (TPaveStats*)metTotL2febR4->GetListOfFunctions()->FindObject("stats");
    metTotL2febR4->GetListOfFunctions()->Remove(statboxmetl2febR4);
    metTotL2febR4->SetStats(0);
    statboxmetl2febR4->SetTextColor(colorR4);
    statboxmetl2febR4->SetY1NDC(minstatboxyR4);
    statboxmetl2febR4->SetY2NDC(maxstatboxyR4);
    statboxmetl2febR4->SetX1NDC(xminstatbox);
    statboxmetl2febR4->SetX2NDC(xmaxstatbox);
    statboxmetl2febR4->Draw();

  TCanvas *METCELLR4 = new TCanvas("METCELLR4","METCELLR4");

    TH1F *metTotCellR4=(TH1F*)metTotCell->Clone();
    metTotCellR4->SetName("metTotCellR4");
    metTotCellR4->SetDirectory(0);

    metTotCellR4->SetLineColor(colorR4);
    metTotCellR4->SetLineWidth(linewidth);
    metTotCellR4->SetTitle("CELL MET");
    metTotCellR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR4->SetNormFactor(1.0);
    metTotCellR4->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR4 = (TPaveStats*)metTotCellR4->GetListOfFunctions()->FindObject("stats");
    metTotCellR4->GetListOfFunctions()->Remove(statboxmetcellR4);
    metTotCellR4->SetStats(0);
    statboxmetcellR4->SetTextColor(colorR4);
    statboxmetcellR4->SetY1NDC(minstatboxyR4);
    statboxmetcellR4->SetY2NDC(maxstatboxyR4);
    statboxmetcellR4->SetX1NDC(xminstatbox);
    statboxmetcellR4->SetX2NDC(xmaxstatbox);
    statboxmetcellR4->Draw();

  TCanvas *METTOPOCLR4 = new TCanvas("METTOPOCLR4","METTOPOCLR4");

    TH1F *metTotTopoclR4=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR4->SetName("metTotTopoclR4");
    metTotTopoclR4->SetDirectory(0);

    metTotTopoclR4->SetLineColor(colorR4);
    metTotTopoclR4->SetLineWidth(linewidth);
    metTotTopoclR4->SetTitle("TOPOCL MET");
    metTotTopoclR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR4->SetNormFactor(1.0);
    metTotTopoclR4->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR4 = (TPaveStats*)metTotTopoclR4->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR4->GetListOfFunctions()->Remove(statboxmettopoclR4);
    metTotTopoclR4->SetStats(0);
    statboxmettopoclR4->SetTextColor(colorR4);
    statboxmettopoclR4->SetY1NDC(minstatboxyR4);
    statboxmettopoclR4->SetY2NDC(maxstatboxyR4);
    statboxmettopoclR4->SetX1NDC(xminstatbox);
    statboxmettopoclR4->SetX2NDC(xmaxstatbox);
    statboxmettopoclR4->Draw();

  TCanvas *METTOPOCLPSR4 = new TCanvas("METTOPOCLPSR4","METTOPOCLPSR4");

    TH1F *metTotTopoclPSR4=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR4->SetName("metTotTopoclPSR4");
    metTotTopoclPSR4->SetDirectory(0);

    metTotTopoclPSR4->SetLineColor(colorR4);
    metTotTopoclPSR4->SetLineWidth(linewidth);
    metTotTopoclPSR4->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR4->SetNormFactor(1.0);
    metTotTopoclPSR4->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR4 = (TPaveStats*)metTotTopoclPSR4->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR4->GetListOfFunctions()->Remove(statboxmettopoclpsR4);
    metTotTopoclPSR4->SetStats(0);
    statboxmettopoclpsR4->SetTextColor(colorR4);
    statboxmettopoclpsR4->SetY1NDC(minstatboxyR4);
    statboxmettopoclpsR4->SetY2NDC(maxstatboxyR4);
    statboxmettopoclpsR4->SetX1NDC(xminstatbox);
    statboxmettopoclpsR4->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR4->Draw();

  TCanvas *METTOPOCLPUCR4 = new TCanvas("METTOPOCLPUCR4","METTOPOCLPUCR4");

    TH1F *metTotTopoclPUCR4=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR4->SetName("metTotTopoclPUCR4");
    metTotTopoclPUCR4->SetDirectory(0);

    metTotTopoclPUCR4->SetLineColor(colorR4);
    metTotTopoclPUCR4->SetLineWidth(linewidth);
    metTotTopoclPUCR4->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR4->SetNormFactor(1.0);
    metTotTopoclPUCR4->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR4 = (TPaveStats*)metTotTopoclPUCR4->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR4->GetListOfFunctions()->Remove(statboxmettopoclpucR4);
    metTotTopoclPUCR4->SetStats(0);
    statboxmettopoclpucR4->SetTextColor(colorR4);
    statboxmettopoclpucR4->SetY1NDC(minstatboxyR4);
    statboxmettopoclpucR4->SetY2NDC(maxstatboxyR4);
    statboxmettopoclpucR4->SetX1NDC(xminstatbox);
    statboxmettopoclpucR4->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR4->Draw();

  TCanvas *METMHTR4 = new TCanvas("METMHTR4","METMHTR4");

    TH1F *metTotMHTR4=(TH1F*)metTotMHT->Clone();
    metTotMHTR4->SetName("metTotMHTR4");
    metTotMHTR4->SetDirectory(0);

    metTotMHTR4->SetLineColor(colorR4);
    metTotMHTR4->SetLineWidth(linewidth);
    metTotMHTR4->SetTitle("MHT MET");
    metTotMHTR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR4->SetNormFactor(1.0);
    metTotMHTR4->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR4 = (TPaveStats*)metTotMHTR4->GetListOfFunctions()->FindObject("stats");
    metTotMHTR4->GetListOfFunctions()->Remove(statboxmetmhtR4);
    metTotMHTR4->SetStats(0);
    statboxmetmhtR4->SetTextColor(colorR4);
    statboxmetmhtR4->SetY1NDC(minstatboxyR4);
    statboxmetmhtR4->SetY2NDC(maxstatboxyR4);
    statboxmetmhtR4->SetX1NDC(xminstatbox);
    statboxmetmhtR4->SetX2NDC(xmaxstatbox);
    statboxmetmhtR4->Draw();

  TCanvas *METREFWMUR4 = new TCanvas("METREFWMUR4","METREFWMUR4");

    TH1F *metTotRefwmuR4=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR4->SetName("metTotRefwmuR4");
    metTotRefwmuR4->SetDirectory(0);

    metTotRefwmuR4->SetLineColor(colorR4);
    metTotRefwmuR4->SetLineWidth(linewidth);
    metTotRefwmuR4->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR4->SetNormFactor(1.0);
    metTotRefwmuR4->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR4 = (TPaveStats*)metTotRefwmuR4->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR4->GetListOfFunctions()->Remove(statboxmetrefwmuR4);
    metTotRefwmuR4->SetStats(0);
    statboxmetrefwmuR4->SetTextColor(colorR4);
    statboxmetrefwmuR4->SetY1NDC(minstatboxyR4);
    statboxmetrefwmuR4->SetY2NDC(maxstatboxyR4);
    statboxmetrefwmuR4->SetX1NDC(xminstatbox);
    statboxmetrefwmuR4->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR4->Draw();

  TCanvas *METOFFRECALR4 = new TCanvas("METOFFRECALR4","METOFFRECALR4");

    TH1F *metTotOffRecalR4=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR4->SetName("metTotOffRecalR4");
    metTotOffRecalR4->SetDirectory(0);

    metTotOffRecalR4->SetLineColor(colorR4);
    metTotOffRecalR4->SetLineWidth(linewidth);
    metTotOffRecalR4->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR4->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR4->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR4->SetNormFactor(1.0);
    metTotOffRecalR4->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR4 = (TPaveStats*)metTotOffRecalR4->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR4->GetListOfFunctions()->Remove(statboxmetoffrecalR4);
    metTotOffRecalR4->SetStats(0);
    statboxmetoffrecalR4->SetTextColor(colorR4);
    statboxmetoffrecalR4->SetY1NDC(minstatboxyR4);
    statboxmetoffrecalR4->SetY2NDC(maxstatboxyR4);
    statboxmetoffrecalR4->SetX1NDC(xminstatbox);
    statboxmetoffrecalR4->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR4->Draw();

  // Fifth File

  //TFile *FileR5 = TFile::Open("Run281143.physics_ZeroBias.merge.AOD.r7562_p2521.root");
  TFile *FileR5 = TFile::Open("JperiodR283429R284285.root");


  TCanvas *METL1R5 = new TCanvas("METL1R5","METL1R5");

    TH1F *metTotL1R5=(TH1F*)metTotL1->Clone();
    metTotL1R5->SetName("metTotL1R5");
    metTotL1R5->SetDirectory(0);

    metTotL1R5->SetLineColor(colorR5);
    metTotL1R5->SetLineWidth(linewidth);
    metTotL1R5->SetTitle("L1 MET");
    metTotL1R5->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R5->SetNormFactor(1.0);
    metTotL1R5->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R5 = (TPaveStats*)metTotL1R5->GetListOfFunctions()->FindObject("stats");
    metTotL1R5->GetListOfFunctions()->Remove(statboxmetl1R5);
    metTotL1R5->SetStats(0);
    statboxmetl1R5->SetTextColor(colorR5);
    statboxmetl1R5->SetY1NDC(minstatboxyR5);
    statboxmetl1R5->SetY2NDC(maxstatboxyR5);
    statboxmetl1R5->SetX1NDC(xminstatbox);
    statboxmetl1R5->SetX2NDC(xmaxstatbox);
    statboxmetl1R5->Draw();

  TCanvas *METL2FEBR5 = new TCanvas("METL2FEBR5","METL2FEBR5");

    TH1F *metTotL2febR5=(TH1F*)metTotL2feb->Clone();
    metTotL2febR5->SetName("metTotL2febR5");
    metTotL2febR5->SetDirectory(0);

    metTotL2febR5->SetLineColor(colorR5);
    metTotL2febR5->SetLineWidth(linewidth);
    metTotL2febR5->SetTitle("FEB MET");
    metTotL2febR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR5->SetNormFactor(1.0);
    metTotL2febR5->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR5 = (TPaveStats*)metTotL2febR5->GetListOfFunctions()->FindObject("stats");
    metTotL2febR5->GetListOfFunctions()->Remove(statboxmetl2febR5);
    metTotL2febR5->SetStats(0);
    statboxmetl2febR5->SetTextColor(colorR5);
    statboxmetl2febR5->SetY1NDC(minstatboxyR5);
    statboxmetl2febR5->SetY2NDC(maxstatboxyR5);
    statboxmetl2febR5->SetX1NDC(xminstatbox);
    statboxmetl2febR5->SetX2NDC(xmaxstatbox);
    statboxmetl2febR5->Draw();

  TCanvas *METCELLR5 = new TCanvas("METCELLR5","METCELLR5");

    TH1F *metTotCellR5=(TH1F*)metTotCell->Clone();
    metTotCellR5->SetName("metTotCellR5");
    metTotCellR5->SetDirectory(0);

    metTotCellR5->SetLineColor(colorR5);
    metTotCellR5->SetLineWidth(linewidth);
    metTotCellR5->SetTitle("CELL MET");
    metTotCellR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR5->SetNormFactor(1.0);
    metTotCellR5->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR5 = (TPaveStats*)metTotCellR5->GetListOfFunctions()->FindObject("stats");
    metTotCellR5->GetListOfFunctions()->Remove(statboxmetcellR5);
    metTotCellR5->SetStats(0);
    statboxmetcellR5->SetTextColor(colorR5);
    statboxmetcellR5->SetY1NDC(minstatboxyR5);
    statboxmetcellR5->SetY2NDC(maxstatboxyR5);
    statboxmetcellR5->SetX1NDC(xminstatbox);
    statboxmetcellR5->SetX2NDC(xmaxstatbox);
    statboxmetcellR5->Draw();

  TCanvas *METTOPOCLR5 = new TCanvas("METTOPOCLR5","METTOPOCLR5");

    TH1F *metTotTopoclR5=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR5->SetName("metTotTopoclR5");
    metTotTopoclR5->SetDirectory(0);

    metTotTopoclR5->SetLineColor(colorR5);
    metTotTopoclR5->SetLineWidth(linewidth);
    metTotTopoclR5->SetTitle("TOPOCL MET");
    metTotTopoclR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR5->SetNormFactor(1.0);
    metTotTopoclR5->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR5 = (TPaveStats*)metTotTopoclR5->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR5->GetListOfFunctions()->Remove(statboxmettopoclR5);
    metTotTopoclR5->SetStats(0);
    statboxmettopoclR5->SetTextColor(colorR5);
    statboxmettopoclR5->SetY1NDC(minstatboxyR5);
    statboxmettopoclR5->SetY2NDC(maxstatboxyR5);
    statboxmettopoclR5->SetX1NDC(xminstatbox);
    statboxmettopoclR5->SetX2NDC(xmaxstatbox);
    statboxmettopoclR5->Draw();

  TCanvas *METTOPOCLPSR5 = new TCanvas("METTOPOCLPSR5","METTOPOCLPSR5");

    TH1F *metTotTopoclPSR5=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR5->SetName("metTotTopoclPSR5");
    metTotTopoclPSR5->SetDirectory(0);

    metTotTopoclPSR5->SetLineColor(colorR5);
    metTotTopoclPSR5->SetLineWidth(linewidth);
    metTotTopoclPSR5->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR5->SetNormFactor(1.0);
    metTotTopoclPSR5->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR5 = (TPaveStats*)metTotTopoclPSR5->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR5->GetListOfFunctions()->Remove(statboxmettopoclpsR5);
    metTotTopoclPSR5->SetStats(0);
    statboxmettopoclpsR5->SetTextColor(colorR5);
    statboxmettopoclpsR5->SetY1NDC(minstatboxyR5);
    statboxmettopoclpsR5->SetY2NDC(maxstatboxyR5);
    statboxmettopoclpsR5->SetX1NDC(xminstatbox);
    statboxmettopoclpsR5->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR5->Draw();

  TCanvas *METTOPOCLPUCR5 = new TCanvas("METTOPOCLPUCR5","METTOPOCLPUCR5");

    TH1F *metTotTopoclPUCR5=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR5->SetName("metTotTopoclPUCR5");
    metTotTopoclPUCR5->SetDirectory(0);

    metTotTopoclPUCR5->SetLineColor(colorR5);
    metTotTopoclPUCR5->SetLineWidth(linewidth);
    metTotTopoclPUCR5->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR5->SetNormFactor(1.0);
    metTotTopoclPUCR5->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR5 = (TPaveStats*)metTotTopoclPUCR5->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR5->GetListOfFunctions()->Remove(statboxmettopoclpucR5);
    metTotTopoclPUCR5->SetStats(0);
    statboxmettopoclpucR5->SetTextColor(colorR5);
    statboxmettopoclpucR5->SetY1NDC(minstatboxyR5);
    statboxmettopoclpucR5->SetY2NDC(maxstatboxyR5);
    statboxmettopoclpucR5->SetX1NDC(xminstatbox);
    statboxmettopoclpucR5->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR5->Draw();

  TCanvas *METMHTR5 = new TCanvas("METMHTR5","METMHTR5");

    TH1F *metTotMHTR5=(TH1F*)metTotMHT->Clone();
    metTotMHTR5->SetName("metTotMHTR5");
    metTotMHTR5->SetDirectory(0);

    metTotMHTR5->SetLineColor(colorR5);
    metTotMHTR5->SetLineWidth(linewidth);
    metTotMHTR5->SetTitle("MHT MET");
    metTotMHTR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR5->SetNormFactor(1.0);
    metTotMHTR5->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR5 = (TPaveStats*)metTotMHTR5->GetListOfFunctions()->FindObject("stats");
    metTotMHTR5->GetListOfFunctions()->Remove(statboxmetmhtR5);
    metTotMHTR5->SetStats(0);
    statboxmetmhtR5->SetTextColor(colorR5);
    statboxmetmhtR5->SetY1NDC(minstatboxyR5);
    statboxmetmhtR5->SetY2NDC(maxstatboxyR5);
    statboxmetmhtR5->SetX1NDC(xminstatbox);
    statboxmetmhtR5->SetX2NDC(xmaxstatbox);
    statboxmetmhtR5->Draw();

  TCanvas *METREFWMUR5 = new TCanvas("METREFWMUR5","METREFWMUR5");

    TH1F *metTotRefwmuR5=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR5->SetName("metTotRefwmuR5");
    metTotRefwmuR5->SetDirectory(0);

    metTotRefwmuR5->SetLineColor(colorR5);
    metTotRefwmuR5->SetLineWidth(linewidth);
    metTotRefwmuR5->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR5->SetNormFactor(1.0);
    metTotRefwmuR5->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR5 = (TPaveStats*)metTotRefwmuR5->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR5->GetListOfFunctions()->Remove(statboxmetrefwmuR5);
    metTotRefwmuR5->SetStats(0);
    statboxmetrefwmuR5->SetTextColor(colorR5);
    statboxmetrefwmuR5->SetY1NDC(minstatboxyR5);
    statboxmetrefwmuR5->SetY2NDC(maxstatboxyR5);
    statboxmetrefwmuR5->SetX1NDC(xminstatbox);
    statboxmetrefwmuR5->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR5->Draw();

  TCanvas *METOFFRECALR5 = new TCanvas("METOFFRECALR5","METOFFRECALR5");

    TH1F *metTotOffRecalR5=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR5->SetName("metTotOffRecalR5");
    metTotOffRecalR5->SetDirectory(0);

    metTotOffRecalR5->SetLineColor(colorR5);
    metTotOffRecalR5->SetLineWidth(linewidth);
    metTotOffRecalR5->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR5->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR5->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR5->SetNormFactor(1.0);
    metTotOffRecalR5->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR5 = (TPaveStats*)metTotOffRecalR5->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR5->GetListOfFunctions()->Remove(statboxmetoffrecalR5);
    metTotOffRecalR5->SetStats(0);
    statboxmetoffrecalR5->SetTextColor(colorR5);
    statboxmetoffrecalR5->SetY1NDC(minstatboxyR5);
    statboxmetoffrecalR5->SetY2NDC(maxstatboxyR5);
    statboxmetoffrecalR5->SetX1NDC(xminstatbox);
    statboxmetoffrecalR5->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR5->Draw();

  /*
  // Sixth File
  TFile *FileR6 = TFile::Open("JperiodR283429R284285.root");

  TCanvas *METL1R6 = new TCanvas("METL1R6","METL1R6");

    TH1F *metTotL1R6=(TH1F*)metTotL1->Clone();
    metTotL1R6->SetName("metTotL1R6");
    metTotL1R6->SetDirectory(0);

    metTotL1R6->SetLineColor(colorR6);
    metTotL1R6->SetLineWidth(linewidth);
    metTotL1R6->SetTitle("L1 MET");
    metTotL1R6->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1R6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL1R6->SetNormFactor(1.0);
    metTotL1R6->Draw();
    gPad->Update();
    TPaveStats *statboxmetl1R6 = (TPaveStats*)metTotL1R6->GetListOfFunctions()->FindObject("stats");
    metTotL1R6->GetListOfFunctions()->Remove(statboxmetl1R6);
    metTotL1R6->SetStats(0);
    statboxmetl1R6->SetTextColor(colorR6);
    statboxmetl1R6->SetY1NDC(minstatboxyR6);
    statboxmetl1R6->SetY2NDC(maxstatboxyR6);
    statboxmetl1R6->SetX1NDC(xminstatbox);
    statboxmetl1R6->SetX2NDC(xmaxstatbox);
    statboxmetl1R6->Draw();

  TCanvas *METL2FEBR6 = new TCanvas("METL2FEBR6","METL2FEBR6");

    TH1F *metTotL2febR6=(TH1F*)metTotL2feb->Clone();
    metTotL2febR6->SetName("metTotL2febR6");
    metTotL2febR6->SetDirectory(0);

    metTotL2febR6->SetLineColor(colorR6);
    metTotL2febR6->SetLineWidth(linewidth);
    metTotL2febR6->SetTitle("FEB MET");
    metTotL2febR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2febR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotL2febR6->SetNormFactor(1.0);
    metTotL2febR6->Draw();
    gPad->Update();
    TPaveStats *statboxmetl2febR6 = (TPaveStats*)metTotL2febR6->GetListOfFunctions()->FindObject("stats");
    metTotL2febR6->GetListOfFunctions()->Remove(statboxmetl2febR6);
    metTotL2febR6->SetStats(0);
    statboxmetl2febR6->SetTextColor(colorR6);
    statboxmetl2febR6->SetY1NDC(minstatboxyR6);
    statboxmetl2febR6->SetY2NDC(maxstatboxyR6);
    statboxmetl2febR6->SetX1NDC(xminstatbox);
    statboxmetl2febR6->SetX2NDC(xmaxstatbox);
    statboxmetl2febR6->Draw();

  TCanvas *METCELLR6 = new TCanvas("METCELLR6","METCELLR6");

    TH1F *metTotCellR6=(TH1F*)metTotCell->Clone();
    metTotCellR6->SetName("metTotCellR6");
    metTotCellR6->SetDirectory(0);

    metTotCellR6->SetLineColor(colorR6);
    metTotCellR6->SetLineWidth(linewidth);
    metTotCellR6->SetTitle("CELL MET");
    metTotCellR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotCellR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotCellR6->SetNormFactor(1.0);
    metTotCellR6->Draw();
    gPad->Update();
    TPaveStats *statboxmetcellR6 = (TPaveStats*)metTotCellR6->GetListOfFunctions()->FindObject("stats");
    metTotCellR6->GetListOfFunctions()->Remove(statboxmetcellR6);
    metTotCellR6->SetStats(0);
    statboxmetcellR6->SetTextColor(colorR6);
    statboxmetcellR6->SetY1NDC(minstatboxyR6);
    statboxmetcellR6->SetY2NDC(maxstatboxyR6);
    statboxmetcellR6->SetX1NDC(xminstatbox);
    statboxmetcellR6->SetX2NDC(xmaxstatbox);
    statboxmetcellR6->Draw();

  TCanvas *METTOPOCLR6 = new TCanvas("METTOPOCLR6","METTOPOCLR6");

    TH1F *metTotTopoclR6=(TH1F*)metTotTopocl->Clone();
    metTotTopoclR6->SetName("metTotTopoclR6");
    metTotTopoclR6->SetDirectory(0);

    metTotTopoclR6->SetLineColor(colorR6);
    metTotTopoclR6->SetLineWidth(linewidth);
    metTotTopoclR6->SetTitle("TOPOCL MET");
    metTotTopoclR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclR6->SetNormFactor(1.0);
    metTotTopoclR6->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclR6 = (TPaveStats*)metTotTopoclR6->GetListOfFunctions()->FindObject("stats");
    metTotTopoclR6->GetListOfFunctions()->Remove(statboxmettopoclR6);
    metTotTopoclR6->SetStats(0);
    statboxmettopoclR6->SetTextColor(colorR6);
    statboxmettopoclR6->SetY1NDC(minstatboxyR6);
    statboxmettopoclR6->SetY2NDC(maxstatboxyR6);
    statboxmettopoclR6->SetX1NDC(xminstatbox);
    statboxmettopoclR6->SetX2NDC(xmaxstatbox);
    statboxmettopoclR6->Draw();

  TCanvas *METTOPOCLPSR6 = new TCanvas("METTOPOCLPSR6","METTOPOCLPSR6");

    TH1F *metTotTopoclPSR6=(TH1F*)metTotTopoclPS->Clone();
    metTotTopoclPSR6->SetName("metTotTopoclPSR6");
    metTotTopoclPSR6->SetDirectory(0);

    metTotTopoclPSR6->SetLineColor(colorR6);
    metTotTopoclPSR6->SetLineWidth(linewidth);
    metTotTopoclPSR6->SetTitle("TOPOCLPS MET");
    metTotTopoclPSR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPSR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPSR6->SetNormFactor(1.0);
    metTotTopoclPSR6->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpsR6 = (TPaveStats*)metTotTopoclPSR6->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPSR6->GetListOfFunctions()->Remove(statboxmettopoclpsR6);
    metTotTopoclPSR6->SetStats(0);
    statboxmettopoclpsR6->SetTextColor(colorR6);
    statboxmettopoclpsR6->SetY1NDC(minstatboxyR6);
    statboxmettopoclpsR6->SetY2NDC(maxstatboxyR6);
    statboxmettopoclpsR6->SetX1NDC(xminstatbox);
    statboxmettopoclpsR6->SetX2NDC(xmaxstatbox);
    statboxmettopoclpsR6->Draw();

  TCanvas *METTOPOCLPUCR6 = new TCanvas("METTOPOCLPUCR6","METTOPOCLPUCR6");

    TH1F *metTotTopoclPUCR6=(TH1F*)metTotTopoclPUC->Clone();
    metTotTopoclPUCR6->SetName("metTotTopoclPUCR6");
    metTotTopoclPUCR6->SetDirectory(0);

    metTotTopoclPUCR6->SetLineColor(colorR6);
    metTotTopoclPUCR6->SetLineWidth(linewidth);
    metTotTopoclPUCR6->SetTitle("TOPOCLPUC MET");
    metTotTopoclPUCR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUCR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotTopoclPUCR6->SetNormFactor(1.0);
    metTotTopoclPUCR6->Draw();
    gPad->Update();
    TPaveStats *statboxmettopoclpucR6 = (TPaveStats*)metTotTopoclPUCR6->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUCR6->GetListOfFunctions()->Remove(statboxmettopoclpucR6);
    metTotTopoclPUCR6->SetStats(0);
    statboxmettopoclpucR6->SetTextColor(colorR6);
    statboxmettopoclpucR6->SetY1NDC(minstatboxyR6);
    statboxmettopoclpucR6->SetY2NDC(maxstatboxyR6);
    statboxmettopoclpucR6->SetX1NDC(xminstatbox);
    statboxmettopoclpucR6->SetX2NDC(xmaxstatbox);
    statboxmettopoclpucR6->Draw();

  TCanvas *METMHTR6 = new TCanvas("METMHTR6","METMHTR6");

    TH1F *metTotMHTR6=(TH1F*)metTotMHT->Clone();
    metTotMHTR6->SetName("metTotMHTR6");
    metTotMHTR6->SetDirectory(0);

    metTotMHTR6->SetLineColor(colorR6);
    metTotMHTR6->SetLineWidth(linewidth);
    metTotMHTR6->SetTitle("MHT MET");
    metTotMHTR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHTR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotMHTR6->SetNormFactor(1.0);
    metTotMHTR6->Draw();
    gPad->Update();
    TPaveStats *statboxmetmhtR6 = (TPaveStats*)metTotMHTR6->GetListOfFunctions()->FindObject("stats");
    metTotMHTR6->GetListOfFunctions()->Remove(statboxmetmhtR6);
    metTotMHTR6->SetStats(0);
    statboxmetmhtR6->SetTextColor(colorR6);
    statboxmetmhtR6->SetY1NDC(minstatboxyR6);
    statboxmetmhtR6->SetY2NDC(maxstatboxyR6);
    statboxmetmhtR6->SetX1NDC(xminstatbox);
    statboxmetmhtR6->SetX2NDC(xmaxstatbox);
    statboxmetmhtR6->Draw();

  TCanvas *METREFWMUR6 = new TCanvas("METREFWMUR6","METREFWMUR6");

    TH1F *metTotRefwmuR6=(TH1F*)metTotRefwmu->Clone();
    metTotRefwmuR6->SetName("metTotRefwmuR6");
    metTotRefwmuR6->SetDirectory(0);

    metTotRefwmuR6->SetLineColor(colorR6);
    metTotRefwmuR6->SetLineWidth(linewidth);
    metTotRefwmuR6->SetTitle("REFFINAL CONTAINER MET");
    metTotRefwmuR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmuR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotRefwmuR6->SetNormFactor(1.0);
    metTotRefwmuR6->Draw();
    gPad->Update();
    TPaveStats *statboxmetrefwmuR6 = (TPaveStats*)metTotRefwmuR6->GetListOfFunctions()->FindObject("stats");
    metTotRefwmuR6->GetListOfFunctions()->Remove(statboxmetrefwmuR6);
    metTotRefwmuR6->SetStats(0);
    statboxmetrefwmuR6->SetTextColor(colorR6);
    statboxmetrefwmuR6->SetY1NDC(minstatboxyR6);
    statboxmetrefwmuR6->SetY2NDC(maxstatboxyR6);
    statboxmetrefwmuR6->SetX1NDC(xminstatbox);
    statboxmetrefwmuR6->SetX2NDC(xmaxstatbox);
    statboxmetrefwmuR6->Draw();

  TCanvas *METOFFRECALR6 = new TCanvas("METOFFRECALR6","METOFFRECALR6");

    TH1F *metTotOffRecalR6=(TH1F*)metTotOffRecal->Clone();
    metTotOffRecalR6->SetName("metTotOffRecalR6");
    metTotOffRecalR6->SetDirectory(0);

    metTotOffRecalR6->SetLineColor(colorR6);
    metTotOffRecalR6->SetLineWidth(linewidth);
    metTotOffRecalR6->SetTitle("OFFLINE RECALIBRATED MET");
    metTotOffRecalR6->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecalR6->GetYaxis()->SetTitle(mettotyaxisTitle);
    metTotOffRecalR6->SetNormFactor(1.0);
    metTotOffRecalR6->Draw();
    gPad->Update();
    TPaveStats *statboxmetoffrecalR6 = (TPaveStats*)metTotOffRecalR6->GetListOfFunctions()->FindObject("stats");
    metTotOffRecalR6->GetListOfFunctions()->Remove(statboxmetoffrecalR6);
    metTotOffRecalR6->SetStats(0);
    statboxmetoffrecalR6->SetTextColor(colorR6);
    statboxmetoffrecalR6->SetY1NDC(minstatboxyR6);
    statboxmetoffrecalR6->SetY2NDC(maxstatboxyR6);
    statboxmetoffrecalR6->SetX1NDC(xminstatbox);
    statboxmetoffrecalR6->SetX2NDC(xmaxstatbox);
    statboxmetoffrecalR6->Draw();
    */

  TCanvas *METL1ALLRUNS = new TCanvas("METL1ALLRUNS","METL1ALLRUNS");

    metTotL1R5->Draw();
    statboxmetl1R5->Draw();

    metTotL1R1->Draw("same");
    statboxmetl1R1->Draw();
    metTotL1R2->Draw("same");
    statboxmetl1R2->Draw();
    metTotL1R3->Draw("same");
    statboxmetl1R3->Draw();
    metTotL1R4->Draw("same");
    statboxmetl1R4->Draw();
    //metTotL1R5->Draw("same");
    //statboxmetl1R5->Draw();
    //metTotL1R6->Draw("same");
    //statboxmetl1R6->Draw();

    TLegend* metl1allruns=new TLegend(0.55,0.6,0.75,0.9);
    metl1allruns->AddEntry("metTotL1R1","Run D","L");
    metl1allruns->AddEntry("metTotL1R2","Run E","L");
    metl1allruns->AddEntry("metTotL1R3","Run F","L");
    metl1allruns->AddEntry("metTotL1R4","Run G","L");
    //metl1allruns->AddEntry("metTotL1R5","Run H","L");
    metl1allruns->AddEntry("metTotL1R5","Run J","L");
    //metl1allruns->AddEntry("metTotL1R6","Run J","L");
    metl1allruns->Draw();

  TCanvas *METFEBALLRUNS = new TCanvas("METFEBALLRUNS","METFEBALLRUNS");

    metTotL2febR5->Draw();
    statboxmetl2febR5->Draw();

    metTotL2febR1->Draw("same" );
    statboxmetl2febR1->Draw();
    metTotL2febR2->Draw("same");
    statboxmetl2febR2->Draw();
    metTotL2febR3->Draw("same");
    statboxmetl2febR3->Draw();
    metTotL2febR4->Draw("same");
    statboxmetl2febR4->Draw();
    //metTotL2febR5->Draw("same");
    //statboxmetl2febR5->Draw();
    //metTotL2febR6->Draw("same");
    //statboxmetl2febR6->Draw();
    TLegend* metl2feballruns=new TLegend(0.55,0.6,0.75,0.9);
    metl2feballruns->AddEntry("metTotL2febR1","Run D","L");
    metl2feballruns->AddEntry("metTotL2febR2","Run E","L");
    metl2feballruns->AddEntry("metTotL2febR3","Run F","L");
    metl2feballruns->AddEntry("metTotL2febR4","Run G","L");
    //metl2feballruns->AddEntry("metTotL2febR5","Run H","L");
    metl2feballruns->AddEntry("metTotL2febR5","Run J","L");
    //metl2feballruns->AddEntry("metTotL2febR6","Run J","L");
    metl2feballruns->Draw();

  TCanvas *METCELLALLRUNS = new TCanvas("METCELLALLRUNS","METCELLALLRUNS");

    metTotCellR5->Draw();
    statboxmetcellR5->Draw();

    metTotCellR1->Draw("same");
    statboxmetcellR1->Draw();
    metTotCellR2->Draw("same");
    statboxmetcellR2->Draw();
    metTotCellR3->Draw("same");
    statboxmetcellR3->Draw();
    metTotCellR4->Draw("same");
    statboxmetcellR4->Draw();
    //metTotCellR5->Draw("same");
    //statboxmetcellR5->Draw();
    //metTotCellR6->Draw("same");
    //statboxmetcellR6->Draw();
    TLegend* metcellallruns=new TLegend(0.55,0.6,0.75,0.9);
    metcellallruns->AddEntry("metTotCellR1","Run D","L");
    metcellallruns->AddEntry("metTotCellR2","Run E","L");
    metcellallruns->AddEntry("metTotCellR3","Run F","L");
    metcellallruns->AddEntry("metTotCellR4","Run G","L");
    //metcellallruns->AddEntry("metTotCellR5","Run H","L");
    metcellallruns->AddEntry("metTotCellR5","Run J","L");
    //metcellallruns->AddEntry("metTotCellR6","Run J","L");

    metcellallruns->Draw();

  TCanvas *METTOPOCLALLRUNS = new TCanvas("METTOPOCLALLRUNS","METTOPOCLALLRUNS");

    metTotTopoclR5->Draw();
    statboxmettopoclR5->Draw();

    metTotTopoclR1->Draw("same" );
    statboxmettopoclR1->Draw();
    metTotTopoclR2->Draw("same");
    statboxmettopoclR2->Draw();
    metTotTopoclR3->Draw("same");
    statboxmettopoclR3->Draw();
    metTotTopoclR4->Draw("same");
    statboxmettopoclR4->Draw();
    //metTotTopoclR5->Draw("same");
    //statboxmettopoclR5->Draw();
    //metTotTopoclR6->Draw("same");
    //statboxmettopoclR6->Draw();
    TLegend* mettopoclallruns=new TLegend(0.55,0.6,0.75,0.9);
    mettopoclallruns->AddEntry("metTotTopoclR1","Run D","L");
    mettopoclallruns->AddEntry("metTotTopoclR2","Run E","L");
    mettopoclallruns->AddEntry("metTotTopoclR3","Run F","L");
    mettopoclallruns->AddEntry("metTotTopoclR4","Run G","L");
    //mettopoclallruns->AddEntry("metTotTopoclR5","Run H","L");
    mettopoclallruns->AddEntry("metTotTopoclR5","Run J","L");
    //mettopoclallruns->AddEntry("metTotTopoclR6","Run J","L");
    mettopoclallruns->Draw();

  TCanvas *METTOPOCLPSALLRUNS = new TCanvas("METTOPOCLPSALLRUNS","METTOPOCLPSALLRUNS");

    metTotTopoclPSR5->Draw();
    statboxmettopoclpsR5->Draw();

    metTotTopoclPSR1->Draw("same");
    statboxmettopoclpsR1->Draw();
    metTotTopoclPSR2->Draw("same");
    statboxmettopoclpsR2->Draw();
    metTotTopoclPSR3->Draw("same");
    statboxmettopoclpsR3->Draw();
    metTotTopoclPSR4->Draw("same");
    statboxmettopoclpsR4->Draw();
    //metTotTopoclPSR5->Draw("same");
    //statboxmettopoclpsR5->Draw();
    //metTotTopoclPSR6->Draw("same");
    //statboxmettopoclpsR6->Draw();
    TLegend* mettopoclpsallruns=new TLegend(0.55,0.6,0.75,0.9);
    mettopoclpsallruns->AddEntry("metTotTopoclPSR1","Run D","L");
    mettopoclpsallruns->AddEntry("metTotTopoclPSR2","Run E","L");
    mettopoclpsallruns->AddEntry("metTotTopoclPSR3","Run F","L");
    mettopoclpsallruns->AddEntry("metTotTopoclPSR4","Run G","L");
    //mettopoclpsallruns->AddEntry("metTotTopoclPSR5","Run H","L");
    mettopoclpsallruns->AddEntry("metTotTopoclPSR5","Run J","L");
    //mettopoclpsallruns->AddEntry("metTotTopoclPSR6","Run J","L");
    mettopoclpsallruns->Draw();

  TCanvas *METTOPOCLPUCALLRUNS = new TCanvas("METTOPOCLPUCALLRUNS","METTOPOCLPUCALLRUNS");

    metTotTopoclPUCR5->Draw();
    statboxmettopoclpucR5->Draw();

    metTotTopoclPUCR1->Draw("same");
    statboxmettopoclpucR1->Draw();
    metTotTopoclPUCR2->Draw("same");
    statboxmettopoclpucR2->Draw();
    metTotTopoclPUCR3->Draw("same");
    statboxmettopoclpucR3->Draw();
    metTotTopoclPUCR4->Draw("same");
    statboxmettopoclpucR4->Draw();
    //metTotTopoclPUCR5->Draw("same");
    //statboxmettopoclpucR5->Draw();
    //metTotTopoclPUCR6->Draw("same");
    //statboxmettopoclpucR6->Draw();
    TLegend* mettopoclpucallruns=new TLegend(0.55,0.6,0.75,0.9);
    mettopoclpucallruns->AddEntry("metTotTopoclPUCR1","Run D","L");
    mettopoclpucallruns->AddEntry("metTotTopoclPUCR2","Run E","L");
    mettopoclpucallruns->AddEntry("metTotTopoclPUCR3","Run F","L");
    mettopoclpucallruns->AddEntry("metTotTopoclPUCR4","Run G","L");
    //mettopoclpucallruns->AddEntry("metTotTopoclPUCR5","Run H","L");
    mettopoclpucallruns->AddEntry("metTotTopoclPUCR5","Run J","L");
    //mettopoclpucallruns->AddEntry("metTotTopoclPUCR6","Run J","L");
    mettopoclpucallruns->Draw();

  TCanvas *METMHTALLRUNS = new TCanvas("METMHTALLRUNS","METMHTALLRUNS");

    metTotMHTR5->Draw();
    statboxmetmhtR5->Draw();

    metTotMHTR1->Draw("same");
    statboxmetmhtR1->Draw();
    metTotMHTR2->Draw("same");
    statboxmetmhtR2->Draw();
    metTotMHTR3->Draw("same");
    statboxmetmhtR3->Draw();
    metTotMHTR4->Draw("same");
    statboxmetmhtR4->Draw();
    //metTotMHTR5->Draw("same");
    //statboxmetmhtR5->Draw();
    //metTotMHTR6->Draw("same");
    //statboxmetmhtR6->Draw();
    TLegend* metmhtallruns=new TLegend(0.55,0.6,0.75,0.9);
    metmhtallruns->AddEntry("metTotMHTR1","Run D","L");
    metmhtallruns->AddEntry("metTotMHTR2","Run E","L");
    metmhtallruns->AddEntry("metTotMHTR3","Run F","L");
    metmhtallruns->AddEntry("metTotMHTR4","Run G","L");
    //metmhtallruns->AddEntry("metTotMHTR5","Run H","L");
    metmhtallruns->AddEntry("metTotMHTR5","Run J","L");
    //metmhtallruns->AddEntry("metTotMHTR6","Run J","L");
    metmhtallruns->Draw();

  TCanvas *METREFWMUALLRUNS = new TCanvas("METREFWMUALLRUNS","METREFWMUALLRUNS");

    metTotRefwmuR5->Draw();
    statboxmetrefwmuR5->Draw();

    metTotRefwmuR1->Draw("same");
    statboxmetrefwmuR1->Draw();
    metTotRefwmuR2->Draw("same");
    statboxmetrefwmuR2->Draw();
    metTotRefwmuR3->Draw("same");
    statboxmetrefwmuR3->Draw();
    metTotRefwmuR4->Draw("same");
    statboxmetrefwmuR4->Draw();
    //metTotRefwmuR5->Draw("same");
    //statboxmetrefwmuR5->Draw();
    //metTotRefwmuR6->Draw("same");
    //statboxmetrefwmuR6->Draw();
    TLegend* metrefwmuallruns=new TLegend(0.55,0.6,0.75,0.9);
    metrefwmuallruns->AddEntry("metTotRefwmuR1","Run D","L");
    metrefwmuallruns->AddEntry("metTotRefwmuR2","Run E","L");
    metrefwmuallruns->AddEntry("metTotRefwmuR3","Run F","L");
    metrefwmuallruns->AddEntry("metTotRefwmuR4","Run G","L");
    //metrefwmuallruns->AddEntry("metTotRefwmuR5","Run H","L");
    metrefwmuallruns->AddEntry("metTotRefwmuR5","Run J","L");
    //metrefwmuallruns->AddEntry("metTotRefwmuR6","Run J","L");
    metrefwmuallruns->Draw();

  TCanvas *METOFFRECALALLRUNS = new TCanvas("METOFFRECALALLRUNS","METOFFRECALALLRUNS");

    metTotOffRecalR5->Draw();
    statboxmetoffrecalR5->Draw();

    metTotOffRecalR1->Draw("same");
    statboxmetoffrecalR1->Draw();
    metTotOffRecalR2->Draw("same");
    statboxmetoffrecalR2->Draw();
    metTotOffRecalR3->Draw("same");
    statboxmetoffrecalR3->Draw();
    metTotOffRecalR4->Draw("same");
    statboxmetoffrecalR4->Draw();
    //metTotOffRecalR5->Draw("same");
    //statboxmetoffrecalR5->Draw();
    //metTotOffRecalR6->Draw("same");
    //statboxmetoffrecalR6->Draw();
    TLegend* metoffrecalallruns=new TLegend(0.55,0.6,0.75,0.9);
    metoffrecalallruns->AddEntry("metTotOffRecalR1","Run D","L");
    metoffrecalallruns->AddEntry("metTotOffRecalR2","Run E","L");
    metoffrecalallruns->AddEntry("metTotOffRecalR3","Run F","L");
    metoffrecalallruns->AddEntry("metTotOffRecalR4","Run G","L");
    //metoffrecalallruns->AddEntry("metTotOffRecalR5","Run H","L");
    metoffrecalallruns->AddEntry("metTotOffRecalR5","Run J","L");
    //metoffrecalallruns->AddEntry("metTotOffRecalR6","Run J","L");
    metoffrecalallruns->Draw();

}


