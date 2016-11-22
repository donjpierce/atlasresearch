{

	// use plain black on white colors
	Int_t icol = 0;
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
	gStyle->SetPaperSize(20, 26);
	//   gStyle->SetPadTopMargin(0.05);
	//   gStyle->SetPadRightMargin(0.05);
	//   gStyle->SetPadBottomMargin(0.16);
	//   gStyle->SetPadLeftMargin(0.12);

	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	// use large fonts
	//Int_t font=72;
	Int_t font = 42;
	Double_t tsize = 0.05;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font, "x");
	gStyle->SetLabelFont(font, "y");
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "x");
	gStyle->SetTitleFont(font, "y");
	gStyle->SetTitleFont(font, "z");
	gStyle->SetLabelSize(0.9*tsize, "x");
	gStyle->SetLabelSize(0.9*tsize, "y");
	gStyle->SetLabelSize(0.9*tsize, "z");
	gStyle->SetTitleSize(0.9*tsize, "x");
	gStyle->SetTitleSize(0.9*tsize, "y");
	gStyle->SetTitleSize(0.9*tsize, "z");

	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2.);
	gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes
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
	int colortopocl = 3;
	int colortopoclps = 4;
	int colortopoclpuc = 6;
	int colorMHT = 7;
	int colorl1 = 30;
	int colorrefwmu = 44;
	int coloroffrecal = 41;

	float xminstatbox = 0.75;
	float xmaxstatbox = 0.9;
	//float deltastatboxy = 0.08;
	//float minstatboxy = 0.18;
	float deltastatboxy = 0.08;
	float minstatboxy = 0.26;
	int linewidth = 3;

	int nmetbins = 100;
	Double_t metmin = 0.0;
	Double_t metmax = 200.0;

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

	//defining values to calculate a correlation-coefficient for each graph would go here
	
	//beginning of comparisons to MET L1

	TCanvas *METCellvsL1 = new TCanvas("METCellvsL1", "MET Cell vs L1");

	TH2F* vMETCellvsL1 = new TH2F("vMETCellvsL1", "MET Cell vs L1", 100, 0., 100, 100, 0., 100);
	vMETCellvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vMETCellvsL1->GetYaxis()->SetTitle("MET Cell [GeV]");
	tree->Draw("metcell:metl1>>vMETCellvsL1", PlotCut, "colz");

	TCanvas *TopoclvsL1 = new TCanvas("TopoclvsL1", "Topocl vs L1");

	TH2F* vTopoclvsL1 = new TH2F("vTopoclvsL1", "Topocl vs L1", 100, 0., 100, 100, 0., 100);
	vTopoclvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vTopoclvsL1->GetYaxis()->SetTitle("Topocl [GeV]");
	tree->Draw("mettopocl:metl1>>vTopoclvsL1", PlotCut, "colz");

	TCanvas *TopoclPSvsL1 = new TCanvas("TopoclPSvsL1", "TopoclPS vs L1");

	TH2F* vTopoclPSvsL1 = new TH2F("vTopoclPSvsL1", " Topocl PS vs L1 MET", 100, 0., 100, 100, 0., 100);
	vTopoclPSvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vTopoclPSvsL1->GetYaxis()->SetTitle("Topocl PS [GeV]");
	tree->Draw("mettopoclps:metl1>>vTopoclPSvsL1", PlotCut, "colz");

	TCanvas *TopoclPUCvsL1 = new TCanvas("TopoclPUCvsL1", "Topocl PUC vs L1");

	TH2F* vTopoclPUCvsL1 = new TH2F("vTopoclPUCvsL1", " Topocl PUC vs L1 MET", 100, 0., 100, 100, 0., 100);
	vTopoclPUCvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vTopoclPUCvsL1->GetYaxis()->SetTitle("Topocl PUC [GeV]");
	tree->Draw("mettopoclpuc:metl1>>vTopoclPUCvsL1", PlotCut, "colz");

	TCanvas *MHTvsL1 = new TCanvas("MHTvsL1", "MHT vs L1");

	TH2F* vMHTvsL1 = new TH2F("vMHTvsL1", " MHT vs L1", 100, 0., 100, 100, 0., 100);
	vMHTvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vMHTvsL1->GetYaxis()->SetTitle("MHT [GeV]");
	tree->Draw("metmht:metl1>>vMHTvsL1", PlotCut, "colz");

	TCanvas *OffRecalvsL1 = new TCanvas("OffRecalvsL1", "Offline Recal vs L1");

	TH2F* vOffRecalvsL1 = new TH2F("vOffRecalvsL1", " Offline Recal vs L1", 100, 0., 100, 100, 0., 100);
	vOffRecalvsL1->GetXaxis()->SetTitle("L1 [GeV]");
	vOffRecalvsL1->GetYaxis()->SetTitle("Offline Recal [GeV]");
	tree->Draw("metoffrecal:metl1>>vOffRecalvsL1", PlotCut, "colz");

	//end of comparisons to MET L1
	//beginning of comparisons to MET Cell

	TCanvas *TopoclvsMetCell= new TCanvas("TopoclvsMetCell", "Topocl vs MET Cell");

	TH2F* vTopoclvsMetCell = new TH2F("vTopoclvsMetCell", "Topocl vs MET Cell", 100, 0., 100, 100, 0., 100);
	vTopoclvsMetCell->GetXaxis()->SetTitle("MET Cell [GeV]");
	vTopoclvsMetCell->GetYaxis()->SetTitle("Topocl [GeV]");
	tree->Draw("mettopocl:metcell>>vTopoclvsMetCell", PlotCut, "colz");

	TCanvas *TopoclPSvsMETCell = new TCanvas("TopoclPSvsMETCell", "Topocl PS vs MET Cell");

	TH2F* vTopoclPSvsMETCell = new TH2F("vTopoclPSvsMETCell", "Topocl PS vs MET Cell", 100, 0., 100, 100, 0., 100);
	vTopoclPSvsMETCell->GetXaxis()->SetTitle("MET Cell [GeV]");
	vTopoclPSvsMETCell->GetYaxis()->SetTitle("Topocl PS [GeV]");
	tree->Draw("mettopoclps:metcell>>vTopoclPSvsMETCell", PlotCut, "colz");

	TCanvas *TopoclPUCvsMETCell= new TCanvas("TopoclPUCvsMETCell", "Topocl PUC vs MET Cell");

	TH2F* vTopoclPUCvsMETCell = new TH2F("vTopoclPUCvsMETCell", "Topocl PUC vs MET Cell", 100, 0., 100, 100, 0., 100);
	vTopoclPUCvsMETCell->GetXaxis()->SetTitle("MET Cell [GeV]");
	vTopoclPUCvsMETCell->GetYaxis()->SetTitle("Topocl PUC [GeV]");
	tree->Draw("mettopoclpuc:metcell>>vTopoclPUCvsMETCell", PlotCut, "colz");

	TCanvas *MHTvsMETCell = new TCanvas("MHTvsMETCell", "MHT vs Met Cell");

	TH2F* vMHTvsMETCell = new TH2F("vMHTvsMETCell", "MET MHT vs MET Cell", 100, 0., 100, 100, 0., 100);
	vMHTvsMETCell->GetXaxis()->SetTitle("MET Cell [GeV]");
	vMHTvsMETCell->GetYaxis()->SetTitle("MHT [GeV]");
	tree->Draw("metmht:metcell>>vMHTvsMETCell", PlotCut, "colz");

	TCanvas *OffRecalvsMetCell = new TCanvas("OffRecalvsMetCell", "Offline Recal vs MET Cell");

	TH2F* vOffRecalvsMetCell = new TH2F("vOffRecalvsMetCell", "Offline Recal vs MET Cell", 100, 0., 100, 100, 0., 100);
	vOffRecalvsMetCell->GetXaxis()->SetTitle("MET Cell [GeV]");
	vOffRecalvsMetCell->GetYaxis()->SetTitle("Offline Recal [GeV]");
	tree->Draw("metoffrecal:metcell>>vOffRecalvsMetCell", PlotCut, "colz");

	//end of comparisons to MET CELL
	//beginning of comparisons to METtopocl

	TCanvas *TopoclPSvsTopocl = new TCanvas("TopoclPSvsTopocl", "Topocl PS vs Topocl");

	TH2F* vTopoclPSvsTopocl = new TH2F("vTopoclPSvsTopocl", "Topocl PS vs Topocl", 100, 0., 100, 100, 0., 100);
	vTopoclPSvsTopocl->GetXaxis()->SetTitle("Topocl [GeV]");
	vTopoclPSvsTopocl->GetYaxis()->SetTitle("Topocl PS [GeV]");
	tree->Draw("mettopoclps:mettopocl>>vTopoclPSvsTopocl", PlotCut, "colz");

	TCanvas *TopoclPUCvsTopocl = new TCanvas("TopoclPUCvsTopocl", "Topocl PUC vs Topocl");

	TH2F* vTopoclPUCvsTopocl = new TH2F("vTopoclPUCvsTopocl", "Topocl PUC vs Topocl", 100, 0., 100, 100, 0., 100);
	vTopoclPUCvsTopocl->GetXaxis()->SetTitle("Topocl [GeV]");
	vTopoclPUCvsTopocl->GetYaxis()->SetTitle("Topocl PUC [GeV]");
	tree->Draw("mettopoclpuc:mettopocl>>vTopoclPUCvsTopocl", PlotCut, "colz");

	TCanvas *MHTvsTopocl = new TCanvas("MHTvsTopocl", "MHT vs Topocl");

	TH2F* vMHTvsTopocl = new TH2F("vMHTvsTopocl", "MHT vs Topocl", 100, 0., 100, 100, 0., 100);
	vMHTvsTopocl->GetXaxis()->SetTitle("Topocl [GeV]");
	vMHTvsTopocl->GetYaxis()->SetTitle("MHT [GeV]");
	tree->Draw("metmht:mettopocl>>vMHTvsTopocl", PlotCut, "colz");

	TCanvas *OffRecalvsTopocl = new TCanvas("OffRecalvsTopocl", "Offline Recal vs Topocl");

	TH2F* vOffRecalvsTopocl = new TH2F("vOffRecalvsTopocl", "Offline Recal vs Topocl", 100, 0., 100, 100, 0., 100);
	vOffRecalvsTopocl->GetXaxis()->SetTitle("Topocl [GeV]");
	vOffRecalvsTopocl->GetYaxis()->SetTitle("Offline Recal [GeV}");
	tree->Draw("metoffrecal:mettopocl>>vOffRecalvsTopocl", PlotCut, "colz");

	//end of comparisons to MET topocl
	//beginning of comparisons to MET topoclPS

	TCanvas *TopoclPUCvsTopoclPS = new TCanvas("TopoclPUCvsTopoclPS", "Topocl PUC vs Topocl PS");

	TH2F* vTopoclPUCvsTopoclPS = new TH2F("vTopoclPUCvsTopoclPS", "Topocl PUC vs Topocl PS", 100, 0., 100, 100, 0., 100);
	vTopoclPUCvsTopoclPS->GetXaxis()->SetTitle("Topocl PS [GeV]");
	vTopoclPUCvsTopoclPS->GetYaxis()->SetTitle("Topocl PUC [GeV}");
	tree->Draw("mettopoclpuc:mettopoclps>>vTopoclPUCvsTopoclPS", PlotCut, "colz");

	TCanvas *MHTvsTopoclPS = new TCanvas("MHTvsTopoclPS", "MHT vs Topocl PS");

	TH2F* vMHTvsTopoclPS = new TH2F("vMHTvsTopoclPS", "MHT vs Topocl PS", 100, 0., 100, 100, 0., 100);
	vMHTvsTopoclPS->GetXaxis()->SetTitle("Topocl PS [GeV]");
	vMHTvsTopoclPS->GetYaxis()->SetTitle("MHT [GeV]");
	tree->Draw("metmht:mettopoclps>>vMHTvsTopoclPS", PlotCut, "colz");

	TCanvas *OffRecalvsTopoclPS = new TCanvas("OffRecalvsTopoclPS", "Offline Recal vs Topocl PS");

	TH2F* vOffRecalvsTopoclPS = new TH2F("vOffRecalvsTopoclPS", "Offline Recal vs TopoclPS", 100, 0., 100, 100, 0., 100);
	vOffRecalvsTopoclPS->GetXaxis()->SetTitle("Topocl PS [GeV]");
	vOffRecalvsTopoclPS->GetYaxis()->SetTitle("Offline Recal [GeV]");
	tree->Draw("metoffrecal:mettopoclps>>vOffRecalvsTopoclPS", PlotCut, "colz");

	//end comparisons to Topocl PS
	//begining of comparions to Topocl PUC

	TCanvas *MHTvsTopoclPUC = new TCanvas("MHTvsTopoclPUC", "MHT vs Topocl PUC");

	TH2F* vMHTvsTopoclPUC = new TH2F("vMHTvsTopoclPUC", "MHT vs Topocl PUC", 100, 0., 100, 100, 0., 100);
	vMHTvsTopoclPUC->GetXaxis()->SetTitle("Topocl PUC [GeV]");
	vMHTvsTopoclPUC->GetYaxis()->SetTitle("MHT [GeV]");
	tree->Draw("metmht:mettopoclpuc>>vMHTvsTopoclPUC", PlotCut, "colz");

	TCanvas *OffRecalvsTopoclPUC = new TCanvas("OffRecalvsTopoclPUC", "Offline Recal vs Topocl PUC");

	TH2F* vOffRecalvsTopoclPUC = new TH2F("vOffRecalvsTopoclPUC", "Offline Recal vs TopoclPUC", 100, 0., 100, 100, 0., 100);
	vOffRecalvsTopoclPUC->GetXaxis()->SetTitle("Topocl PUC [GeV]");
	vOffRecalvsTopoclPUC->GetYaxis()->SetTitle("Offline Recal [GeV]");
	tree->Draw("metoffrecal:mettopoclpuc>>vOffRecalvsTopoclPUC", PlotCut, "colz");

	//beginning of last comparison

	TCanvas *MHTvsOffRecal = new TCanvas("MHTvsOffRecal", "MHT vs Offline Recal");

	TH2F* vMHTvsOffRecal = new TH2F("vMHTvsOffRecal", "MHT vs Offline Recal", 100, 0., 100, 100, 0., 100);
	vMHTvsOffRecal->GetXaxis()->SetTitle("MET Offline Recal [GeV]");
	vMHTvsOffRecal->GetYaxis()->SetTitle("MET MHT [GeV]");
	tree->Draw("metmht:metoffrecal>>vMHTvsOffRecal", PlotCut, "colz");




}