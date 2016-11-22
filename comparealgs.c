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

  TCanvas *METL1 = new TCanvas("METL1","METL1");

    TH1F* metTotL1 = new TH1F("metTotL1"," L1 MET",  nmetbins, metmin, metmax);
    //metTotL1->SetDirectory(0);
    metTotL1->SetLineColor(colorl1);
    metTotL1->SetLineWidth(linewidth);
    metTotL1->SetTitle(mettotTitle);
    metTotL1->GetXaxis()->SetTitle("MET [GeV]");
    metTotL1->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metl1>>metTotL1",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetl1 = (TPaveStats*)metTotL1->GetListOfFunctions()->FindObject("stats");
    metTotL1->GetListOfFunctions()->Remove(statboxmetl1);
    metTotL1->SetStats(0);
    statboxmetl1->SetTextColor(colorl1);
    statboxmetl1->SetY1NDC(minstatboxy);
    statboxmetl1->SetY2NDC(minstatboxy+deltastatboxy);
    statboxmetl1->SetX1NDC(xminstatbox);
    statboxmetl1->SetX2NDC(xmaxstatbox);
    statboxmetl1->Draw();

  TCanvas *METL1INT = new TCanvas("METL1INT","METL1INT");
    
    TH1F* l1metint = new TH1F("l1metint",intplotTitle, nmetbins, metmin, metmax);
    //l1metint->SetDirectory(0);
    l1metint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    l1metint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    l1metint->SetLineColor(colorl1);
    l1metint->SetLineWidth(linewidth);
    double zvall1, sumitl1, ratiol1;
    sumitl1=metTotL1->GetBinContent(nmetbins+1);
    Double_t* totnuml1 = metTotL1->GetEntries();
    std::cout << "L1 threshold plot \n";
    double zerobiasl1thresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvall1=metTotL1->GetBinContent(j);
      sumitl1 += zvall1;
      std::cout << "number of events in bin   " << j << "  content   " << zvall1 << "\n";
      ratiol1 = sumitl1/totnuml1;
      l1metint->SetBinContent(j,ratiol1);
      std::cout << "L1 Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiol1 <<"\n";
      zerobiasl1thresh[j-1]=ratiol1;
    }
    zval=metTotL1->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotL1->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotL1->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    l1metint->Draw();
    gPad->Update();
    TPaveStats *statboxl1int = (TPaveStats*)l1metint->GetListOfFunctions()->FindObject("stats");
    l1metint->GetListOfFunctions()->Remove(statboxl1int);
    l1metint->SetStats(0);
    statboxl1int->SetTextColor(colorl1);
    statboxl1int->SetY1NDC(minstatboxy);
    statboxl1int->SetY2NDC(minstatboxy+deltastatboxy);
    statboxl1int->SetX1NDC(xminstatbox);
    statboxl1int->SetX2NDC(xmaxstatbox);
    statboxl1int->Draw();

  TCanvas *METL2FEB = new TCanvas("METL2FEB","METL2FEB");

    TH1F* metTotL2feb = new TH1F("metTotL2feb"," L2 FEB MET",  nmetbins, metmin, metmax);
    //metTotL2feb->SetDirectory(0);
    metTotL2feb->SetLineColor(colorfeb);
    metTotL2feb->SetLineWidth(linewidth);
    metTotL2feb->SetTitle(mettotTitle);
    metTotL2feb->GetXaxis()->SetTitle("MET [GeV]");
    metTotL2feb->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metl2feb>>metTotL2feb",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetl2feb = (TPaveStats*)metTotL2feb->GetListOfFunctions()->FindObject("stats");
    metTotL2feb->GetListOfFunctions()->Remove(statboxmetl2feb);
    metTotL2feb->SetStats(0);
    statboxmetl2feb->SetTextColor(colorfeb);
    statboxmetl2feb->SetY1NDC(minstatboxy+deltastatboxy);
    statboxmetl2feb->SetY2NDC(minstatboxy+2*deltastatboxy);
    statboxmetl2feb->SetX1NDC(xminstatbox);
    statboxmetl2feb->SetX2NDC(xmaxstatbox);
    statboxmetl2feb->Draw();

  TCanvas *METL2FEBINT = new TCanvas("METL2FEBINT","METL2FEBINT");

    TH1F* l2febmetint = new TH1F("l2febmetint",intplotTitle, nmetbins, metmin, metmax);
    //l2febmetint->SetDirectory(0);
    l2febmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    l2febmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    l2febmetint->SetLineColor(colorfeb);
    l2febmetint->SetLineWidth(linewidth);
    double zvall2feb, sumitl2feb, ratiol2feb;
    sumitl2feb=metTotL2feb->GetBinContent(nmetbins+1);
    Double_t* totnuml2feb = metTotL2feb->GetEntries();
    std::cout << "L2 FEB threshold plot \n";
    double zerobiasl2febthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvall2feb=metTotL2feb->GetBinContent(j);
      sumitl2feb += zvall2feb;
      std::cout << "number of events in bin   " << j << "  content   " << zvall2feb << "\n";
      ratiol2feb = sumitl2feb/totnuml2feb;
      l2febmetint->SetBinContent(j,ratiol2feb);
      std::cout << "L2 FEB Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiol2feb <<"\n";
      zerobiasl2febthresh[j-1]=ratiol2feb;
    }
    zval=metTotL2feb->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotL2feb->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotL2feb->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    l2febmetint->Draw();
    gPad->Update();
    TPaveStats *statboxl2febint = (TPaveStats*)l2febmetint->GetListOfFunctions()->FindObject("stats");
    l2febmetint->GetListOfFunctions()->Remove(statboxl2febint);
    l2febmetint->SetStats(0);
    statboxl2febint->SetTextColor(colorfeb);
    statboxl2febint->SetY1NDC(minstatboxy+deltastatboxy);
    statboxl2febint->SetY2NDC(minstatboxy+2*deltastatboxy);
    statboxl2febint->SetX1NDC(xminstatbox);
    statboxl2febint->SetX2NDC(xmaxstatbox);
    statboxl2febint->Draw();

  TCanvas *METCELL = new TCanvas("METCELL","METCELL");

    TH1F* metTotCell = new TH1F("metTotCell"," Cell MET",  nmetbins, metmin, metmax);
    //metTotCell->SetDirectory(0);
    metTotCell->SetLineColor(colorcell);
    metTotCell->SetLineWidth(linewidth);
    metTotCell->SetTitle(mettotTitle);
    metTotCell->GetXaxis()->SetTitle("MET [GeV]");
    metTotCell->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metcell>>metTotCell",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetcell = (TPaveStats*)metTotCell->GetListOfFunctions()->FindObject("stats");
    metTotCell->GetListOfFunctions()->Remove(statboxmetcell);
    metTotCell->SetStats(0);
    statboxmetcell->SetTextColor(colorcell);
    statboxmetcell->SetY1NDC(minstatboxy+2*deltastatboxy);
    statboxmetcell->SetY2NDC(minstatboxy+3*deltastatboxy);
    statboxmetcell->SetX1NDC(xminstatbox);
    statboxmetcell->SetX2NDC(xmaxstatbox);
    statboxmetcell->Draw();

  TCanvas *METCELLINT = new TCanvas("METCELLINT","METCELLINT");

    TH1F* cellmetint = new TH1F("cellmetint",intplotTitle, nmetbins, metmin, metmax);
    //cellmetint->SetDirectory(0);
    cellmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    cellmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    cellmetint->SetLineColor(colorcell);
    cellmetint->SetLineWidth(linewidth);
    double zvalcell, sumitcell, ratiocell;
    sumitcell=metTotCell->GetBinContent(nmetbins+1);
    Double_t* totnumcell = metTotCell->GetEntries();
    std::cout << "CELL threshold plot \n";
    double zerobiascellthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvalcell=metTotCell->GetBinContent(j);
      sumitcell += zvalcell;
      std::cout << "number of events in bin   " << j << "  content   " << zvalcell << "\n";
      ratiocell = sumitcell/totnumcell;
      cellmetint->SetBinContent(j,ratiocell);
      std::cout << "CELL Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiocell <<"\n";
      zerobiascellthresh[j-1]=ratiocell;
    }
    zval=metTotCell->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotCell->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotCell->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    cellmetint->Draw();
    gPad->Update();
    TPaveStats *statboxcellint = (TPaveStats*)cellmetint->GetListOfFunctions()->FindObject("stats");
    cellmetint->GetListOfFunctions()->Remove(statboxcellint);
    cellmetint->SetStats(0);
    statboxcellint->SetTextColor(colorcell);
    statboxcellint->SetY1NDC(minstatboxy+2*deltastatboxy);
    statboxcellint->SetY2NDC(minstatboxy+3*deltastatboxy);
    statboxcellint->SetX1NDC(xminstatbox);
    statboxcellint->SetX2NDC(xmaxstatbox);
    statboxcellint->Draw();

  TCanvas *METTOPOCL = new TCanvas("METTOPOCL","METTOPOCL");

    TH1F* metTotTopocl = new TH1F("metTotTopocl"," Topocl MET",  nmetbins, metmin, metmax);
    //metTotTopocl->SetDirectory(0);
    metTotTopocl->SetLineColor(colortopocl);
    metTotTopocl->SetLineWidth(linewidth);
    metTotTopocl->SetTitle(mettotTitle);
    metTotTopocl->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopocl->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("mettopocl>>metTotTopocl",PlotCut);
    gPad->Update();
    TPaveStats *statboxmettopocl = (TPaveStats*)metTotTopocl->GetListOfFunctions()->FindObject("stats");
    metTotTopocl->GetListOfFunctions()->Remove(statboxmettopocl);
    metTotTopocl->SetStats(0);
    statboxmettopocl->SetTextColor(colortopocl);
    statboxmettopocl->SetY1NDC(minstatboxy+3*deltastatboxy);
    statboxmettopocl->SetY2NDC(minstatboxy+4*deltastatboxy);
    statboxmettopocl->SetX1NDC(xminstatbox);
    statboxmettopocl->SetX2NDC(xmaxstatbox);
    statboxmettopocl->Draw();

  TCanvas *METTOPOCLINT = new TCanvas("METTOPOCLINT","METTOPOCLINT");

    TH1F* topoclmetint = new TH1F("topoclmetint",intplotTitle, nmetbins, metmin, metmax);
    //topoclmetint->SetDirectory(0);
    topoclmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    topoclmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    topoclmetint->SetLineColor(colortopocl);
    topoclmetint->SetLineWidth(linewidth);
    double zvaltopocl, sumittopocl, ratiotopocl;
    sumittopocl=metTotTopocl->GetBinContent(nmetbins+1);
    Double_t* totnumtopocl = metTotTopocl->GetEntries();
    std::cout << "TOPOCL threshold plot \n";
    double zerobiastopoclthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvaltopocl=metTotTopocl->GetBinContent(j);
      sumittopocl += zvaltopocl;
      std::cout << "number of events in bin   " << j << "  content   " << zvaltopocl << "\n";
      ratiotopocl = sumittopocl/totnumtopocl;
      topoclmetint->SetBinContent(j,ratiotopocl);
      std::cout << "TOPOCL Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiotopocl <<"\n";
      zerobiastopoclthresh[j-1]=ratiotopocl;
    }
    zval=metTotTopocl->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotTopocl->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotTopocl->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    topoclmetint->Draw();
    gPad->Update();
    TPaveStats *statboxtopoclint = (TPaveStats*)topoclmetint->GetListOfFunctions()->FindObject("stats");
    topoclmetint->GetListOfFunctions()->Remove(statboxtopoclint);
    topoclmetint->SetStats(0);
    statboxtopoclint->SetTextColor(colortopocl);
    statboxtopoclint->SetY1NDC(minstatboxy+3*deltastatboxy);
    statboxtopoclint->SetY2NDC(minstatboxy+4*deltastatboxy);
    statboxtopoclint->SetX1NDC(xminstatbox);
    statboxtopoclint->SetX2NDC(xmaxstatbox);
    statboxtopoclint->Draw();

  TCanvas *METTOPOCLPS = new TCanvas("METTOPOCLPS","METTOPOCLPS");

    TH1F* metTotTopoclPS = new TH1F("metTotTopoclPS"," TopoclPS MET",  nmetbins, metmin, metmax);
    //metTotTopoclPS->SetDirectory(0);
    metTotTopoclPS->SetLineColor(colortopoclps);
    metTotTopoclPS->SetLineWidth(linewidth);
    metTotTopoclPS->SetTitle("mettotTitle");
    metTotTopoclPS->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPS->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("mettopoclps>>metTotTopoclPS",PlotCut);
    gPad->Update();
    TPaveStats *statboxmettopoclps = (TPaveStats*)metTotTopoclPS->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPS->GetListOfFunctions()->Remove(statboxmettopoclps);
    metTotTopoclPS->SetStats(0);
    statboxmettopoclps->SetTextColor(colortopoclps);
    statboxmettopoclps->SetY1NDC(minstatboxy+4*deltastatboxy);
    statboxmettopoclps->SetY2NDC(minstatboxy+5*deltastatboxy);
    statboxmettopoclps->SetX1NDC(xminstatbox);
    statboxmettopoclps->SetX2NDC(xmaxstatbox);
    statboxmettopoclps->Draw();

  TCanvas *METTOPOCLPSINT = new TCanvas("METTOPOCLPSINT","METTOPOCLPSINT");

    TH1F* topoclpsmetint = new TH1F("topoclpsmetint",intplotTitle, nmetbins, metmin, metmax);
    //topoclpsmetint->SetDirectory(0);
    topoclpsmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    topoclpsmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    topoclpsmetint->SetLineColor(colortopoclps);
    topoclpsmetint->SetLineWidth(linewidth);
    double zvaltopoclps, sumittopoclps, ratiotopoclps;
    sumittopoclps=metTotTopoclPS->GetBinContent(nmetbins+1);
    Double_t* totnumtopoclps = metTotTopoclPS->GetEntries();
    std::cout << "TOPOCLPS threshold plot \n";
    double zerobiastopoclpsthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvaltopoclps=metTotTopoclPS->GetBinContent(j);
      sumittopoclps += zvaltopoclps;
      std::cout << "number of events in bin   " << j << "  content   " << zvaltopoclps << "\n";
      ratiotopoclps = sumittopoclps/totnumtopoclps;
      topoclpsmetint->SetBinContent(j,ratiotopoclps);
      std::cout << "TOPOCLPS Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiotopoclps <<"\n";
      zerobiastopoclpsthresh[j-1]=ratiotopoclps;
    }
    zval=metTotTopoclPS->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotTopoclPS->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotTopoclPS->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    topoclpsmetint->Draw();
    gPad->Update();
    TPaveStats *statboxtopoclpsint = (TPaveStats*)topoclpsmetint->GetListOfFunctions()->FindObject("stats");
    topoclpsmetint->GetListOfFunctions()->Remove(statboxtopoclpsint);
    topoclpsmetint->SetStats(0);
    statboxtopoclpsint->SetTextColor(colortopoclps);
    statboxtopoclpsint->SetY1NDC(minstatboxy+4*deltastatboxy);
    statboxtopoclpsint->SetY2NDC(minstatboxy+5*deltastatboxy);
    statboxtopoclpsint->SetX1NDC(xminstatbox);
    statboxtopoclpsint->SetX2NDC(xmaxstatbox);
    statboxtopoclpsint->Draw();

  TCanvas *METTOPOCLPUC = new TCanvas("METTOPOCLPUC","METTOPOCLPUC");

    TH1F* metTotTopoclPUC = new TH1F("metTotTopoclPUC"," TopoclPUC MET",  nmetbins, metmin, metmax);
    //metTotTopoclPUC->SetDirectory(0);
    metTotTopoclPUC->SetLineColor(colortopoclpuc);
    metTotTopoclPUC->SetLineWidth(linewidth);
    metTotTopoclPUC->SetTitle(mettotTitle);
    metTotTopoclPUC->GetXaxis()->SetTitle("MET [GeV]");
    metTotTopoclPUC->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("mettopoclpuc>>metTotTopoclPUC",PlotCut);
    gPad->Update();
    TPaveStats *statboxmettopoclpuc = (TPaveStats*)metTotTopoclPUC->GetListOfFunctions()->FindObject("stats");
    metTotTopoclPUC->GetListOfFunctions()->Remove(statboxmettopoclpuc);
    metTotTopoclPUC->SetStats(0);
    statboxmettopoclpuc->SetTextColor(colortopoclpuc);
    statboxmettopoclpuc->SetY1NDC(minstatboxy+5*deltastatboxy);
    statboxmettopoclpuc->SetY2NDC(minstatboxy+6*deltastatboxy);
    statboxmettopoclpuc->SetX1NDC(xminstatbox);
    statboxmettopoclpuc->SetX2NDC(xmaxstatbox);
    statboxmettopoclpuc->Draw();

  TCanvas *METTOPOCLPUCINT = new TCanvas("METTOPOCLPUCINT","METTOPOCLPUCINT");

    TH1F* topoclpucmetint = new TH1F("topoclpucmetint",intplotTitle, nmetbins, metmin, metmax);
    //topoclpucmetint->SetDirectory(0);
    topoclpucmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    topoclpucmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    topoclpucmetint->SetLineColor(colortopoclpuc);
    topoclpucmetint->SetLineWidth(linewidth);
    double zvaltopoclpuc, sumittopoclpuc, ratiotopoclpuc;
    sumittopoclpuc=metTotTopoclPUC->GetBinContent(nmetbins+1);
    Double_t* totnumtopoclpuc = metTotTopoclPUC->GetEntries();
    std::cout << "TOPOCLPUC threshold plot \n";
    double zerobiastopoclpucthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvaltopoclpuc=metTotTopoclPUC->GetBinContent(j);
      sumittopoclpuc += zvaltopoclpuc;
      std::cout << "number of events in bin   " << j << "  content   " << zvaltopoclpuc << "\n";
      ratiotopoclpuc = sumittopoclpuc/totnumtopoclpuc;
      topoclpucmetint->SetBinContent(j,ratiotopoclpuc);
      std::cout << "TOPOCLPUC Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiotopoclpuc <<"\n";
      zerobiastopoclpucthresh[j-1]=ratiotopoclpuc;
    }
    zval=metTotTopoclPUC->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotTopoclPUC->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotTopoclPUC->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    topoclpucmetint->Draw();
    gPad->Update();
    TPaveStats *statboxtopoclpucint = (TPaveStats*)topoclpucmetint->GetListOfFunctions()->FindObject("stats");
    topoclpucmetint->GetListOfFunctions()->Remove(statboxtopoclpucint);
    topoclpucmetint->SetStats(0);
    statboxtopoclpucint->SetTextColor(colortopoclpuc);
    statboxtopoclpucint->SetY1NDC(minstatboxy+5*deltastatboxy);
    statboxtopoclpucint->SetY2NDC(minstatboxy+6*deltastatboxy);
    statboxtopoclpucint->SetX1NDC(xminstatbox);
    statboxtopoclpucint->SetX2NDC(xmaxstatbox);
    statboxtopoclpucint->Draw();

  TCanvas *METMHT = new TCanvas("METMHT","METMHT");

    TH1F* metTotMHT = new TH1F("metTotMHT"," MHT MET",  nmetbins, metmin, metmax);
    //metTotMHT->SetDirectory(0);
    metTotMHT->SetLineColor(colorMHT);
    metTotMHT->SetLineWidth(linewidth);
    metTotMHT->SetTitle(mettotTitle);
    metTotMHT->GetXaxis()->SetTitle("MET [GeV]");
    metTotMHT->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metmht>>metTotMHT",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetmht = (TPaveStats*)metTotMHT->GetListOfFunctions()->FindObject("stats");
    metTotMHT->GetListOfFunctions()->Remove(statboxmetmht);
    metTotMHT->SetStats(0);
    statboxmetmht->SetTextColor(colorMHT);
    statboxmetmht->SetY1NDC(minstatboxy+6*deltastatboxy);
    statboxmetmht->SetY2NDC(minstatboxy+7*deltastatboxy);
    statboxmetmht->SetX1NDC(xminstatbox);
    statboxmetmht->SetX2NDC(xmaxstatbox);
    statboxmetmht->Draw();

  TCanvas *METMHTINT = new TCanvas("METMHTINT","METMHTINT");

    TH1F* mhtmetint = new TH1F("mhtmetint",intplotTitle, nmetbins, metmin, metmax);
    //mhtmetint->SetDirectory(0);
    mhtmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    mhtmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    mhtmetint->SetLineColor(colorMHT);
    mhtmetint->SetLineWidth(linewidth);
    double zvalmht, sumitmht, ratiomht;
    sumitmht=metTotMHT->GetBinContent(nmetbins+1);
    Double_t* totnummht = metTotMHT->GetEntries();
    std::cout << "MHT threshold plot \n";
    double zerobiasmhtthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvalmht=metTotMHT->GetBinContent(j);
      sumitmht += zvalmht;
      std::cout << "number of events in bin   " << j << "  content   " << zvalmht << "\n";
      ratiomht = sumitmht/totnummht;
      mhtmetint->SetBinContent(j,ratiomht);
      std::cout << "MHT Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiomht <<"\n";
      zerobiasmhtthresh[j-1]=ratiomht;
    }
    zval=metTotMHT->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotMHT->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotMHT->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    mhtmetint->Draw();
    gPad->Update();
    TPaveStats *statboxmhtint = (TPaveStats*)mhtmetint->GetListOfFunctions()->FindObject("stats");
    mhtmetint->GetListOfFunctions()->Remove(statboxmhtint);
    mhtmetint->SetStats(0);
    statboxmhtint->SetTextColor(colorMHT);
    statboxmhtint->SetY1NDC(minstatboxy+6*deltastatboxy);
    statboxmhtint->SetY2NDC(minstatboxy+7*deltastatboxy);
    statboxmhtint->SetX1NDC(xminstatbox);
    statboxmhtint->SetX2NDC(xmaxstatbox);
    statboxmhtint->Draw();

  TCanvas *METREFWMU = new TCanvas("METREFWMU","METREFWMU");

    TH1F* metTotRefwmu = new TH1F("metTotRefwmu"," RefFinal with muons MET",  nmetbins, metmin, metmax);
    //metTotRefwmu->SetDirectory(0);
    metTotRefwmu->SetLineColor(colorrefwmu);
    metTotRefwmu->SetLineWidth(linewidth);
    metTotRefwmu->SetTitle(mettotTitle);
    metTotRefwmu->GetXaxis()->SetTitle("MET [GeV]");
    metTotRefwmu->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metrefwmuon>>metTotRefwmu",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetrefwmu = (TPaveStats*)metTotRefwmu->GetListOfFunctions()->FindObject("stats");
    metTotRefwmu->GetListOfFunctions()->Remove(statboxmetrefwmu);
    metTotRefwmu->SetStats(0);
    statboxmetrefwmu->SetTextColor(colorrefwmu);
    statboxmetrefwmu->SetY1NDC(minstatboxy+8*deltastatboxy);
    statboxmetrefwmu->SetY2NDC(minstatboxy+9*deltastatboxy);
    statboxmetrefwmu->SetX1NDC(xminstatbox);
    statboxmetrefwmu->SetX2NDC(xmaxstatbox);
    statboxmetrefwmu->Draw();

  TCanvas *METREFWMUINT = new TCanvas("METREFWMUINT","METREFWMUINT");

    TH1F* refwmumetint = new TH1F("refwmumetint",intplotTitle, nmetbins, metmin, metmax);
    //refwmumetint->SetDirectory(0);
    refwmumetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    refwmumetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    refwmumetint->SetLineColor(colorrefwmu);
    refwmumetint->SetLineWidth(linewidth);
    double zvalrefwmu, sumitrefwmu, ratiorefwmu;
    sumitrefwmu=metTotRefwmu->GetBinContent(nmetbins+1);
    Double_t* totnumrefwmu = metTotRefwmu->GetEntries();
    std::cout << "REF Wmu threshold plot \n";
    double zerobiasrefwmuthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvalrefwmu=metTotRefwmu->GetBinContent(j);
      sumitrefwmu += zvalrefwmu;
      std::cout << "number of events in bin   " << j << "  content   " << zvalrefwmu << "\n";
      ratiorefwmu = sumitrefwmu/totnumrefwmu;
      refwmumetint->SetBinContent(j,ratiorefwmu);
      std::cout << "Ref wmu Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiorefwmu <<"\n";
      zerobiasrefwmuthresh[j-1]=ratiorefwmu;
    }
    zval=metTotRefwmu->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotRefwmu->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotRefwmu->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    refwmumetint->Draw();
    gPad->Update();
    TPaveStats *statboxrefwmuint = (TPaveStats*)refwmumetint->GetListOfFunctions()->FindObject("stats");
    refwmumetint->GetListOfFunctions()->Remove(statboxrefwmuint);
    refwmumetint->SetStats(0);
    statboxrefwmuint->SetTextColor(colorrefwmu);
    statboxrefwmuint->SetY1NDC(minstatboxy+8*deltastatboxy);
    statboxrefwmuint->SetY2NDC(minstatboxy+9*deltastatboxy);
    statboxrefwmuint->SetX1NDC(xminstatbox);
    statboxrefwmuint->SetX2NDC(xmaxstatbox);
    statboxrefwmuint->Draw();

  TCanvas *METOFFRECAL = new TCanvas("METOFFRECAL","METOFFRECAL");

    TH1F* metTotOffRecal = new TH1F("metTotOffRecal"," Offlien Recalibrated MET",  nmetbins, metmin, metmax);
    //metTotOffRecal->SetDirectory(0);
    metTotOffRecal->SetLineColor(coloroffrecal);
    metTotOffRecal->SetLineWidth(linewidth);
    metTotOffRecal->SetTitle(mettotTitle);
    metTotOffRecal->GetXaxis()->SetTitle("MET [GeV]");
    metTotOffRecal->GetYaxis()->SetTitle("Number of Events / 2 GeV ");
    tree->Draw("metoffrecal>>metTotOffRecal",PlotCut);
    gPad->Update();
    TPaveStats *statboxmetoffrecal = (TPaveStats*)metTotOffRecal->GetListOfFunctions()->FindObject("stats");
    metTotOffRecal->GetListOfFunctions()->Remove(statboxmetoffrecal);
    metTotOffRecal->SetStats(0);
    statboxmetoffrecal->SetTextColor(coloroffrecal);
    statboxmetoffrecal->SetY1NDC(minstatboxy+7*deltastatboxy);
    statboxmetoffrecal->SetY2NDC(minstatboxy+8*deltastatboxy);
    statboxmetoffrecal->SetX1NDC(xminstatbox);
    statboxmetoffrecal->SetX2NDC(xmaxstatbox);
    statboxmetoffrecal->Draw();

  TCanvas *METOFFRECALINT = new TCanvas("METOFFRECALINT","METOFFRECALINT");

    TH1F* offrecalmetint = new TH1F("offrecalmetint",intplotTitle, nmetbins, metmin, metmax);
    //offrecalmetint->SetDirectory(0);
    offrecalmetint->GetXaxis()->SetTitle("Algorithm MET [GeV]");
    offrecalmetint->GetYaxis()->SetTitle("Fraction of Events Above Threshold");
    offrecalmetint->SetLineColor(coloroffrecal);
    offrecalmetint->SetLineWidth(linewidth);
    double zvaloffrecal, sumitoffrecal, ratiooffrecal;
    sumitoffrecal=metTotOffRecal->GetBinContent(nmetbins+1);
    Double_t* totnumoffrecal = metTotOffRecal->GetEntries();
    std::cout << "Offline Recal threshold plot \n";
    double zerobiasoffrecalthresh[nmetbins];
    for (int j=nmetbins; j>0; j--) {
      zvaloffrecal=metTotOffRecal->GetBinContent(j);
      sumitoffrecal += zvaloffrecal;
      std::cout << "number of events in bin   " << j << "  content   " << zvaloffrecal << "\n";
      ratiooffrecal = sumitoffrecal/totnumoffrecal;
      offrecalmetint->SetBinContent(j,ratiooffrecal);
      std::cout << "Off Recal Zero Bias " << j << "MET [GeV] " <<  float(metmin)+float((j-1)*(metmax-metmin))/float(nmetbins)    << "  fraction above threshold " << ratiooffrecal <<"\n";
      zerobiasoffrecalthresh[j-1]=ratiooffrecal;
    }
    zval=metTotOffRecal->GetBinContent(0);
    std::cout << "number of underflow events  " << zval << "\n";
    zval=metTotOffRecal->GetBinContent(nmetbins+1);
    std::cout << "number of overflow events  " << zval << "\n";
    zval=metTotOffRecal->GetEntries();
    std::cout << "number of entries = " << zval << "\n";
    offrecalmetint->Draw();
    gPad->Update();
    TPaveStats *statboxoffrecalint = (TPaveStats*)offrecalmetint->GetListOfFunctions()->FindObject("stats");
    offrecalmetint->GetListOfFunctions()->Remove(statboxoffrecalint);
    offrecalmetint->SetStats(0);
    statboxoffrecalint->SetTextColor(coloroffrecal);
    statboxoffrecalint->SetY1NDC(minstatboxy+7*deltastatboxy);
    statboxoffrecalint->SetY2NDC(minstatboxy+8*deltastatboxy);
    statboxoffrecalint->SetX1NDC(xminstatbox);
    statboxoffrecalint->SetX2NDC(xmaxstatbox);
    statboxoffrecalint->Draw();

  TCanvas *METALL = new TCanvas("METALL","METALL");

    metTotTopoclPUC->Draw();
    statboxmettopoclpuc->Draw();

    metTotL1->Draw("same");
    statboxmetl1->Draw();
    metTotL2feb->Draw("same");
    statboxmetl2feb->Draw();
    metTotCell->Draw("same");
    statboxmetcell->Draw();
    metTotTopocl->Draw("same");
    statboxmettopocl->Draw();
    metTotTopoclPS->Draw("same");
    statboxmettopoclps->Draw();
    //metTotTopoclPUC->Draw("same");
    //statboxmettopoclpuc->Draw();
    metTotMHT->Draw("same");
    statboxmetmht->Draw();
    //metTotRefwmu->Draw("same");
    //statboxmetrefwmu->Draw();
    metTotOffRecal->Draw("same");
    statboxmetoffrecal->Draw();
    TLegend* metallleg=new TLegend(0.55,0.6,0.75,0.9);
    metallleg->AddEntry("metTotOffRecal","Off Recal","L");
    //metallleg->AddEntry("metTotRefwmu","Ref Final Wmu","L");
    metallleg->AddEntry("metTotMHT","MHT","L");
    metallleg->AddEntry("metTotTopoclPUC","Topocl PUC","L");
    metallleg->AddEntry("metTotTopoclPS","Topocl PS","L");
    metallleg->AddEntry("metTotTopocl","Topocl","L");
    metallleg->AddEntry("metTotCell","Cell","L");
    metallleg->AddEntry("metTotL2feb","L2 FEB","L");
    metallleg->AddEntry("metTotL1","L1 RecEnergyROI","L");
    metallleg->Draw();

  TCanvas *METALLINT = new TCanvas("METALLINT","METALLINT");

    l1metint->Draw();
    //statboxl1int->Draw();
    l2febmetint->Draw("same");
    //statboxl2febint->Draw();
    cellmetint->Draw("same");
    //statboxcellint->Draw();
    topoclmetint->Draw("same");
    //statboxtopoclint->Draw();
    topoclpsmetint->Draw("same");
    //statboxtopoclpsint->Draw();
    topoclpucmetint->Draw("same");
    //statboxtopoclpucint->Draw();
    mhtmetint->Draw("same");
    //statboxmhtint->Draw();
    //refwmumetint->Draw("same");
    //statboxrefwmuint->Draw();
    offrecalmetint->Draw("same");
    //statboxoffrecalint->Draw();

    TLegend* metallintleg=new TLegend(0.55,0.6,0.75,0.9);
    metallintleg->AddEntry("offrecalmetint","Off Recal","L");
    //metallintleg->AddEntry("refwmumetint","Ref Wmu","L");
    metallintleg->AddEntry("mhtmetint","MHT","L");
    metallintleg->AddEntry("topoclpucmetint","Topocl PUC","L");
    metallintleg->AddEntry("topoclpsmetint","Topocl PS","L");
    metallintleg->AddEntry("topoclmetint","Topocl","L");
    metallintleg->AddEntry("cellmetint","Cell","L");
    metallintleg->AddEntry("l2febmetint","L2 FEB","L");
    metallintleg->AddEntry("l1metint","L1 RecEnergyROI","L");
    metallintleg->Draw();

  TCanvas *METCELLVSL1 = new TCanvas("METCELLVSL1","METCELLVSL1");

    TH2F* metTotCellvsL1 = new TH2F("metTotCellvsL1"," Cell vs L1 MET", nmetbins, metmin, metmax, nmetbins, metmin, metmax);
    metTotCellvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotCellvsL1->GetYaxis()->SetTitle("Cell MET [GeV]");
    tree->Draw("metcell:metl1>>metTotCellvsL1",PlotCut,"colz");

  TCanvas *METTOPOCLVSL1 = new TCanvas("METTOPOCLVSL1","METTOPOCLVSL1");

    TH2F* metTotTopoclvsL1 = new TH2F("metTotTopoclvsL1"," Topocl vs L1 MET", nmetbins, metmin, metmax, nmetbins, metmin, metmax);
    metTotTopoclvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotTopoclvsL1->GetYaxis()->SetTitle("Topocl MET [GeV]");
    tree->Draw("mettopocl:metl1>>metTotTopoclvsL1",PlotCut,"colz");

  TCanvas *METPSVSL1 = new TCanvas("METPSVSL1","METPSVSL1");

    TH2F* metTotPSvsL1 = new TH2F("metTotPSvsL1"," TopoclPS vs L1 MET", nmetbins, metmin, metmax, nmetbins, metmin, metmax);
    metTotPSvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotPSvsL1->GetYaxis()->SetTitle("TopoclPS MET [GeV]");
    tree->Draw("mettopoclps:metl1>>metTotPSvsL1",PlotCut,"colz");

  TCanvas *METPUCVSL1 = new TCanvas("METPUCVSL1","METPUCVSL1");

    TH2F* metTotPUCvsL1 = new TH2F("metTotPUCvsL1"," TopoclPUC vs L1 MET", nmetbins, metmin, metmax, nmetbins, metmin, metmax);
    metTotPUCvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotPUCvsL1->GetYaxis()->SetTitle("TopoclPUC MET [GeV]");
    tree->Draw("mettopoclpuc:metl1>>metTotPUCvsL1",PlotCut,"colz");

  TCanvas *METMHTVSL1 = new TCanvas("METMHTVSL1","METMHTVSL1");

    TH2F* metTotMHTvsL1 = new TH2F("metTotMHTvsL1"," MHT vs L1 MET", nmetbins, metmin, metmax, nmetbins, metmin, metmax);
    metTotMHTvsL1->GetXaxis()->SetTitle("L1 MET [GeV]");
    metTotMHTvsL1->GetYaxis()->SetTitle("MHT MET [GeV]");
    tree->Draw("metmht:metl1>>metTotMHTvsL1",PlotCut,"colz");


}


