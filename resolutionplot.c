{
	TString PlotCut("passrndm>0.1");

	// Define the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Defining a Linear Fit Function
	TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
	linfit->SetParameters(0, 1.);
	linfit->SetParameters(1, 1.);
	linfit->SetParLimits(0, -50., 50.);
	linfit->SetParLimits(1, -50, 50.);
	linfit->SetParName(0, "slope");
	linfit->SetParName(1, "intercept");

	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	
	//Produce fitting graphs 2015
	TH2F *L115 = new TH2F ("L115","", 20, 0., 20.,100,0.,100.);
		tree->Draw("metl1:sqrt(setl1)>>L115",PlotCut);
	/*TH2F *CELL15 = new TH2F ("CELL15","CELL2015", 50, 0., 50.,100,0.,100.);
		tree->Draw("metcell:sqrt(setcell)>>CELL15", PlotCut);
	TH2F *TOPOCL15 = new TH2F ("TOPOCL15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCL15",PlotCut);
	TH2F *TOPOCLPS15 = new TH2F ("TOPOCLPS15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPS15", PlotCut);
	TH2F *TOPOCLPUC15 = new TH2F ("TOPOCLPUC15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUC15", PlotCut);
	TH2F *MHT15 = new TH2F ("MHT15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("metmht:sqrt(setmht)>>MHT15", PlotCut);
	TH2F *OFFRECAL15 = new TH2F ("OFFRECAL15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL15",PlotCut);
		*/
	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");

	//Fitting graphs 2016
	TH2F *L116 = new TH2F ("L116","", 20, 0., 20.,100,0.,100.);
		tree->Draw("metl1:sqrt(setl1)>>L116",PlotCut);
	/*TH2F *CELL16 = new TH2F ("CELL16","",50,0.,50.,100,0.,100.);
		tree->Draw("metcell:sqrt(setcell)>>CELL16",PlotCut);
	TH2F *TOPOCL16 = new TH2F ("TOPOCL16","",50,0.,50.,100,0.,100.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCL16",PlotCut);
	TH2F *TOPOCLPS16 = new TH2F ("TOPOCLPS16","",50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPS16", PlotCut);
	TH2F *TOPOCLPUC16 = new TH2F ("TOPOCLPUC16","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUC16", PlotCut);
	TH2F *MHT16 = new TH2F ("MHT16","", 50, 0., 50.,100,0.,100.);
		tree->Draw("metmht:sqrt(setmht)>>MHT16",PlotCut);
	TH2F *OFFRECAL16 = new TH2F ("OFFRECAL16","", 50, 0., 50.,100,0.,100.);
		tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL16",PlotCut);
		*/

//L1 Algorithm resolutions in 2015 and 2016 
	TCanvas *cL115 = new TCanvas("cL115", "L1 2015");
	L115->Draw();
	L115->FitSlicesY(func, 0, -1, 10, "L");
	L115_1->Draw();
	L115_1->Fit("linfit");
	L115_1->SetTitle("Resolution of L1 in 2015");
	L115_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L115_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L115_1->SetLineColor(2);
	gPad->Update();
		TPaveStats *l115 = (TPaveStats*)L115_1 ->FindObject("stats");
		l115->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resl115 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl115->AddEntry("L115_1", "2015 Data", "L");
		resl115->Draw();

	TCanvas *cL116 = new TCanvas("cL116", "L1 2016 ");
	L116->Draw();
	L116->FitSlicesY(func, 0, -1, 10, "L");
	L116_1->Draw();
	L116_1->Fit("linfit");
	L116_1->SetTitle("Resolution of L1 in 2016 ");
	L116_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L116_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L116_1->SetLineColor(4);
	gPad->Update();
		TPaveStats *l116 = (TPaveStats*)L116_1 ->FindObject("stats");
		l116->SetTextColor(4);
		gStyle->SetOptFit(11);
	
		TLegend* resl116 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl116->AddEntry("L116_1", "2016 Data", "L");
		resl116->Draw();
/*
//CELL Algorithm resoltuions in 2015 and 2016 
		//TCanvas *cCELL15 = new TCanvas("cCELL15", "CELL 2015 ");
		TCanvas *cCELL15 = new TCanvas("cCELL15", "CELL Resolution in 2015 and 2016");
		CELL15->Draw();
		CELL15->FitSlicesY(func, 0, -1, 10, "L");
		CELL15_1->Draw();
		CELL15_1->SetTitle("Resolution of CELL in 2015 ");
		CELL15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELL15_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELL15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sCELL15 = (TPaveStats*)CELL15_1->FindObject("stats");
		sCELL15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resCELL15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELL15->AddEntry("CELL15_1", "2015 Data", "L");
		resCELL15->Draw();

		TCanvas *cCELL16 = new TCanvas("cCELL16", "CELL 2016 ");
		CELL16->Draw();
		CELL16->FitSlicesY(func, 0, -1, 10, "L");
		CELL16_1->Draw();
		CELL16_1->SetTitle("Resolution of CELL in 2016 ");
		CELL16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELL16_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELL16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sCELL16 = (TPaveStats*)CELL16_1->FindObject("stats");
		sCELL16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resCELL16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELL16->AddEntry("CELL16_1", "2016 Data", "L");
		resCELL16->Draw();
		
		//TOPOCL Algorithm resoltuions in 2015 and 2016 
		TCanvas *cTOPOCL15 = new TCanvas("cTOPOCL15", "TOPOCL 2015 ");
		TOPOCL15->Draw();
		TOPOCL15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCL15_1->Draw();
		TOPOCL15_1->SetTitle("Resolution of TOPOCL in 2015 ");
		TOPOCL15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCL15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCL15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCL15 = (TPaveStats*)TOPOCL15_1->FindObject("stats");
		sTOPOCL15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCL15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCL15->AddEntry("TOPOCL15_1", "2015 Data", "L");
		resTOPOCL15->Draw();

		TCanvas *cTOPOCL16 = new TCanvas("cTOPOCL16", "TOPOCL 2016 ");
		TOPOCL16->Draw();
		TOPOCL16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCL16_1->Draw();
		TOPOCL16_1->SetTitle("Resolution of TOPOCL in 2016 ");
		TOPOCL16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCL16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCL16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCL16 = (TPaveStats*)TOPOCL16_1->FindObject("stats");
		sTOPOCL16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCL16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCL16->AddEntry("TOPOCL16_1", "2016 Data", "L");
		resTOPOCL16->Draw();
	
		//TOPOCLPS Algorithm resoltuions in 2015 and 2016 
		TCanvas *cTOPOCLPS15 = new TCanvas("cTOPOCLPS15", "TOPOCLPS 2015 ");
		TOPOCLPS15->Draw();
		TOPOCLPS15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPS15_1->Draw();
		TOPOCLPS15_1->SetTitle("Resolution of TOPOCLPS in 2015 ");
		TOPOCLPS15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPS15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPS15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPS15 = (TPaveStats*)TOPOCLPS15_1->FindObject("stats");
		sTOPOCLPS15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPS15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPS15->AddEntry("TOPOCLPS15_1", "2015 Data", "L");
		resTOPOCLPS15->Draw();

		TCanvas *cTOPOCLPS16 = new TCanvas("cTOPOCLPS16", "TOPOCLPS 2016 ");
		TOPOCLPS16->Draw();
		TOPOCLPS16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPS16_1->Draw();
		TOPOCLPS16_1->SetTitle("Resolution of TOPOCLPS in 2016 ");
		TOPOCLPS16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPS16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPS16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPS16 = (TPaveStats*)TOPOCLPS16_1->FindObject("stats");
		sTOPOCLPS16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPS16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPS16->AddEntry("TOPOCLPS16_1", "2016 Data", "L");
		resTOPOCLPS16->Draw();
		
		//TOPOCLPUC Algorithm resoltuions in 2015 and 2016 
		TCanvas *cTOPOCLPUC15 = new TCanvas("cTOPOCLPUC15", "TOPOCLPUC 2015 ");
		TOPOCLPUC15->Draw();
		TOPOCLPUC15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPUC15_1->Draw();
		TOPOCLPUC15_1->SetTitle("Resolution of TOPOCLPUC in 2015 ");
		TOPOCLPUC15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUC15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUC15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPUC15 = (TPaveStats*)TOPOCLPUC15_1->FindObject("stats");
		sTOPOCLPUC15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUC15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUC15->AddEntry("TOPOCLPUC15_1", "2015 Data", "L");
		resTOPOCLPUC15->Draw();

		TCanvas *cTOPOCLPS16 = new TCanvas("cTOPOCLPS16", "TOPOCLPS 2016 ");
		TOPOCLPUC16->Draw();
		TOPOCLPUC16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPUC16_1->Draw();
		TOPOCLPUC16_1->SetTitle("Resolution of TOPOCLPUC in 2016 ");
		TOPOCLPUC16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUC16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUC16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPUC16 = (TPaveStats*)TOPOCLPUC16_1->FindObject("stats");
		sTOPOCLPUC16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUC16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUC16->AddEntry("TOPOCLPUC16_1", "2016 Data", "L");
		resTOPOCLPUC16->Draw();

		
		//MHT Algorithm resoltuions in 2015 and 2016 
		TCanvas *cMHT15 = new TCanvas("cMHT15", "MHT 2015 ");
		MHT15->Draw();
		MHT15->FitSlicesY(func, 0, -1, 10, "L");
		MHT15_1->Draw();
		MHT15_1->SetTitle("Resolution of MHT in 2015 ");
		MHT15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHT15_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHT15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sMHT15 = (TPaveStats*)MHT15_1->FindObject("stats");
		sMHT15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resMHT15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHT15->AddEntry("MHT15_1", "2015 Data", "L");
		resMHT15->Draw();

		TCanvas *cMHT16 = new TCanvas("cMHT16", "MHT 2016 ");
		MHT16->Draw();
		MHT16->FitSlicesY(func, 0, -1, 10, "L");
		MHT16_1->Draw();
		MHT16_1->SetTitle("Resolution of MHT in 2016 ");
		MHT16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHT16_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHT16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sMHT16 = (TPaveStats*)MHT16_1->FindObject("stats");
		sMHT16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resMHT16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHT16->AddEntry("MHT16_1", "2016 Data", "L");
		resMHT16->Draw();

	//OFFRECAL Algorithm resoltuions in 2015 and 2016 
		TCanvas *cOFFRECAL15 = new TCanvas("cOFFRECAL15", "OFFRECAL 2015 ");
		OFFRECAL15->Draw();
		OFFRECAL15->FitSlicesY(func, 0, -1, 10, "L");
		OFFRECAL15_1->Draw();
		OFFRECAL15_1->SetTitle("Resolution of OFFRECAL in 2015 ");
		OFFRECAL15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		OFFRECAL15_1->GetYaxis()->SetTitle("#sigma of Fit for OFFRECAL [GeV]");
		OFFRECAL15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sOFFRECAL15 = (TPaveStats*)OFFRECAL15_1->FindObject("stats");
		sOFFRECAL15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resOFFRECAL15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resOFFRECAL15->AddEntry("OFFRECAL15_1", "2015 Data", "L");
		resOFFRECAL15->Draw();

		TCanvas *cOFFRECAL16 = new TCanvas("cOFFRECAL16", "OFFRECAL 2016 ");
		OFFRECAL16->Draw();
		OFFRECAL16->FitSlicesY(func, 0, -1, 10, "L");
		OFFRECAL16_1->Draw();
		OFFRECAL16_1->SetTitle("Resolution of OFFRECAL in 2016 ");
		OFFRECAL16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		OFFRECAL16_1->GetYaxis()->SetTitle("#sigma of Fit for OFFRECAL [GeV]");
		OFFRECAL16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sOFFRECAL16 = (TPaveStats*)OFFRECAL16_1->FindObject("stats");
		sOFFRECAL16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resOFFRECAL16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resOFFRECAL16->AddEntry("OFFRECAL16_1", "2016 Data", "L");
		resOFFRECAL16->Draw();
		*/

}