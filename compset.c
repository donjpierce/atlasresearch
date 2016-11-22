{
	TString PlotCut ("passrndm>0.1&&runnum<283200.");
	TString PlotCut1 ("passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TString PlotCut2 ("passrndm>0.1&&runnum>283800.");
	TString PlotCut3 ("passrndm>0.1&&runnum<299000.");
	TString PlotCut4 ("passrndm>0.1&&runnum>299000.&&runnum<301000.");
	TString PlotCut5 ("passrndm>0.1&&runnum>301000.");

	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	TCanvas *CompareL1byrun15 = new TCanvas("CompareL1byrun15","CompareL1byrun15");
	TH1F *L115 = new TH1F ("L115","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L115","passrndm>0.1");
	TH1F *L115R1 = new TH1F ("L115R1","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L115R1", PlotCut);
	TH1F *L115R2 = new TH1F ("L115R2","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L115R2", PlotCut1);
	TH1F *L115R3 = new TH1F ("L115R3","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L115R3", PlotCut2);
	
		L115R1->Draw();
		L115R1->SetNormFactor(1.0);
		L115R1->SetLineColor(2);
		L115R2->Draw("sames");
		L115R2->SetNormFactor(1.0);
		L115R2->SetLineColor(3);
		L115R3->Draw("sames");
		L115R3->SetNormFactor(1.0);
		L115R3->SetLineColor(4);
		TLegend* compl1byrun=new TLegend(0.5,0.7,0.7,0.88);
    	compl1byrun->AddEntry("L115R1","2015 group 1","L");
    	compl1byrun->AddEntry("L115R2","2016 group 2","L");
    	compl1byrun->AddEntry("L115R3","2016 group 3","L");
    	compl1byrun->Draw();

	TCanvas *CompareCELLbyrun15 = new TCanvas("CompareCELLbyrun15","CompareCELLbyrun15");
	TH1F *CELL15 = new TH1F ("CELL15","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL15","passrndm>0.1");
	TH1F *CELL15R1 = new TH1F ("CELL15R1","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL15R1", PlotCut);
	TH1F *CELL15R2 = new TH1F ("CELL15R2","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL15R2", PlotCut1);
	TH1F *CELL15R3 = new TH1F ("CELL15R3","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL15R3", PlotCut2);
	
		CELL15R1->Draw();
		CELL15R1->SetNormFactor(1.0);
		CELL15R1->SetLineColor(2);
		CELL15R2->Draw("sames");
		CELL15R2->SetNormFactor(1.0);
		CELL15R2->SetLineColor(3);
		CELL15R3->Draw("sames");
		CELL15R3->SetNormFactor(1.0);
		CELL15R3->SetLineColor(4);
		TLegend* compcellbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	compcellbyrun->AddEntry("CELL15R1","2015 group 1","L");
    	compcellbyrun->AddEntry("CELL15R2","2016 group 2","L");
    	compcellbyrun->AddEntry("CELL15R3","2016 group 3","L");
    	compcellbyrun->Draw();


	TCanvas *CompareTOPOCLbyrun15 = new TCanvas("CompareTOPOCLbyrun15","CompareTOPOCLbyrun15");
	TH1F *TOPOCL15 = new TH1F ("TOPOCL15","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL15","passrndm>0.1");
	TH1F *TOPOCL15R1 = new TH1F ("TOPOCL15R1","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL15R1", PlotCut);
	TH1F *TOPOCL15R2 = new TH1F ("TOPOCL15R2","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL15R2", PlotCut1);
	TH1F *TOPOCL15R3 = new TH1F ("TOPOCL15R3","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL15R3", PlotCut2);
	
		TOPOCL15R1->Draw();
		TOPOCL15R1->SetNormFactor(1.0);
		TOPOCL15R1->SetLineColor(2);
		TOPOCL15R2->Draw("sames");
		TOPOCL15R2->SetNormFactor(1.0);
		TOPOCL15R2->SetLineColor(3);
		TOPOCL15R3->Draw("sames");
		TOPOCL15R3->SetNormFactor(1.0);
		TOPOCL15R3->SetLineColor(4);
		TLegend* comptopoclbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	comptopoclbyrun->AddEntry("TOPOCL15R1","2015 group 1","L");
    	comptopoclbyrun->AddEntry("TOPOCL15R2","2016 group 2","L");
    	comptopoclbyrun->AddEntry("TOPOCL15R3","2016 group 3","L");
    	comptopoclbyrun->Draw();


	TCanvas *ComparePSbyrun15 = new TCanvas("ComparePSbyrun15","ComparePSbyrun15");
	TH1F *PS15 = new TH1F ("PS15","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS15","passrndm>0.1");
	TH1F *PS15R1 = new TH1F ("PS15R1","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS15R1", PlotCut);
	TH1F *PS15R2 = new TH1F ("PS15R2","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS15R2", PlotCut1);
	TH1F *PS15R3 = new TH1F ("PS15R3","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS15R3", PlotCut2);
	
		PS15R1->Draw();
		PS15R1->SetNormFactor(1.0);
		PS15R1->SetLineColor(2);
		PS15R2->Draw("sames");
		PS15R2->SetNormFactor(1.0);
		PS15R2->SetLineColor(3);
		PS15R3->Draw("sames");
		PS15R3->SetNormFactor(1.0);
		PS15R3->SetLineColor(4);
		TLegend* comppsbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	comppsbyrun->AddEntry("PS15R1","2015 group 1","L");
    	comppsbyrun->AddEntry("PS15R2","2016 group 2","L");
    	comppsbyrun->AddEntry("PS15R3","2016 group 3","L");
    	comppsbyrun->Draw();


	TCanvas *ComparePUCbyrun15 = new TCanvas("ComparePUCbyrun15","ComparePUCbyrun15");
	TH1F *PUC15 = new TH1F ("PUC15","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC15","passrndm>0.1&&mettopoclpuc>0.1");
	TH1F *PUC15R1 = new TH1F ("PUC15R1","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC15R1", "passrndm>0.1&&runnum<283200.&&mettopoclpuc>0.1");
	TH1F *PUC15R2 = new TH1F ("PUC15R2","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC15R2", "passrndm>0.1&&runnum>283200.&&runnum<283800.&&mettopoclpuc>0.1");
	TH1F *PUC15R3 = new TH1F ("PUC15R3","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC15R3", "passrndm>0.1&&runnum>283800.&&mettopoclpuc>0.1");
	
		PUC15R1->Draw();
		PUC15R1->SetNormFactor(1.0);
		PUC15R1->SetLineColor(2);
		PUC15R2->Draw("sames");
		PUC15R2->SetNormFactor(1.0);
		PUC15R2->SetLineColor(3);
		PUC15R3->Draw("sames");
		PUC15R3->SetNormFactor(1.0);
		PUC15R3->SetLineColor(4);
		TLegend* comppucbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	comppucbyrun->AddEntry("PUC15R1","2015 group 1","L");
    	comppucbyrun->AddEntry("PUC15R2","2016 group 2","L");
    	comppucbyrun->AddEntry("PUC15R3","2016 group 2","L");
    	comppucbyrun->Draw();


	TCanvas *CompareMHTbyrun15 = new TCanvas("CompareMHTbyrun15","CompareMHTbyrun15");
	TH1F *MHT15 = new TH1F ("MHT15","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT15","passrndm>0.1");
	TH1F *MHT15R1 = new TH1F ("MHT15R1","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT15R1", PlotCut);
	TH1F *MHT15R2 = new TH1F ("MHT15R2","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT15R2", PlotCut1);
	TH1F *MHT15R3 = new TH1F ("MHT15R3","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT15R3", PlotCut2);
	
		MHT15R1->Draw();
		MHT15R1->SetNormFactor(1.0);
		MHT15R1->SetLineColor(2);
		MHT15R2->Draw("sames");
		MHT15R2->SetNormFactor(1.0);
		MHT15R2->SetLineColor(3);
		MHT15R3->Draw("sames");
		MHT15R3->SetNormFactor(1.0);
		MHT15R3->SetLineColor(4);
		TLegend* compmhtbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	compmhtbyrun->AddEntry("MHT15R1","2015 group 1","L");
    	compmhtbyrun->AddEntry("MHT15R2","2016 group 2","L");
    	compmhtbyrun->AddEntry("MHT15R3","2016 group 2","L");
    	compmhtbyrun->Draw();


	TCanvas *CompareOFFRECALbyrun15 = new TCanvas("CompareOFFRECALbyrun15","CompareOFFRECALbyrun15");
	TH1F *OFFRECAL15 = new TH1F ("OFFRECAL15","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL15","passrndm>0.1");
	TH1F *OFFRECAL15R1 = new TH1F ("OFFRECAL15R1","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL15R1", PlotCut);
	TH1F *OFFRECAL15R2 = new TH1F ("OFFRECAL15R2","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL15R2", PlotCut1);
	TH1F *OFFRECAL15R3 = new TH1F ("OFFRECAL15R3","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL15R3", PlotCut2);
	
		OFFRECAL15R1->Draw();
		OFFRECAL15R1->SetNormFactor(1.0);
		OFFRECAL15R1->SetLineColor(2);
		OFFRECAL15R2->Draw("sames");
		OFFRECAL15R2->SetNormFactor(1.0);
		OFFRECAL15R2->SetLineColor(3);
		OFFRECAL15R3->Draw("sames");
		OFFRECAL15R3->SetNormFactor(1.0);
		OFFRECAL15R3->SetLineColor(4);
		TLegend* compoffrecalbyrun=new TLegend(0.5,0.7,0.7,0.88);
    	compoffrecalbyrun->AddEntry("OFFRECAL15R1","2015 group 1","L");
    	compoffrecalbyrun->AddEntry("OFFRECAL15R2","2016 group 2","L");
    	compoffrecalbyrun->AddEntry("OFFRECAL15R3","2016 group 2","L");
    	compoffrecalbyrun->Draw();

	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");
	TCanvas *CompareL1byrun16 = new TCanvas("CompareL1byrun16","CompareL1byrun16");
	TH1F *L116 = new TH1F ("L116","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L116","passrndm>0.1");
	TH1F *L116R1 = new TH1F ("L116R1","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L116R1", PlotCut3);
	TH1F *L116R2 = new TH1F ("L116R2","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L116R2", PlotCut4);
	TH1F *L116R3 = new TH1F ("L116R3","",30,0.,30.);
		tree->Draw("sqrt(setl1)>>L116R3", PlotCut5);
	
		L116R1->Draw();
		L116R1->SetNormFactor(1.0);
		L116R1->SetLineColor(2);
		L116R2->Draw("sames");
		L116R2->SetNormFactor(1.0);
		L116R2->SetLineColor(3);
		L116R3->Draw("sames");
		L116R3->SetNormFactor(1.0);
		L116R3->SetLineColor(4);
		TLegend* compl1byrun16=new TLegend(0.5,0.7,0.7,0.88);
    	compl1byrun16->AddEntry("L116R1","2015 group 1","L");
    	compl1byrun16->AddEntry("L116R2","2016 group 2","L");
    	compl1byrun16->AddEntry("L116R3","2016 group 3","L");
    	compl1byrun16->Draw();

	TCanvas *CompareCELLbyrun16 = new TCanvas("CompareCELLbyrun16","CompareCELLbyrun16");
	TH1F *CELL16 = new TH1F ("CELL16","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL16","passrndm>0.1");
	TH1F *CELL16R1 = new TH1F ("CELL16R1","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL16R1", PlotCut3);
	TH1F *CELL16R2 = new TH1F ("CELL16R2","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL16R2", PlotCut4);
	TH1F *CELL16R3 = new TH1F ("CELL16R3","",30,0.,30.);
		tree->Draw("sqrt(setcell)>>CELL16R3", PlotCut5);
	
		CELL16R1->Draw();
		CELL16R1->SetNormFactor(1.0);
		CELL16R1->SetLineColor(2);
		CELL16R2->Draw("sames");
		CELL16R2->SetNormFactor(1.0);
		CELL16R2->SetLineColor(3);
		CELL16R3->Draw("sames");
		CELL16R3->SetNormFactor(1.0);
		CELL16R3->SetLineColor(4);
		TLegend* compcellbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	compcellbyrun16->AddEntry("CELL16R1","2015 group 1","L");
    	compcellbyrun16->AddEntry("CELL16R2","2016 group 2","L");
    	compcellbyrun16->AddEntry("CELL16R3","2016 group 3","L");
    	compcellbyrun16->Draw();


	TCanvas *CompareTOPOCLbyrun16 = new TCanvas("CompareTOPOCLbyrun16","CompareTOPOCLbyrun16");
	TH1F *TOPOCL16 = new TH1F ("TOPOCL16","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL16","passrndm>0.1");
	TH1F *TOPOCL16R1 = new TH1F ("TOPOCL16R1","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL16R1", PlotCut3);
	TH1F *TOPOCL16R2 = new TH1F ("TOPOCL16R2","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL16R2", PlotCut4);
	TH1F *TOPOCL16R3 = new TH1F ("TOPOCL16R3","",30,0.,30.);
		tree->Draw("sqrt(settopocl)>>TOPOCL16R3", PlotCut5);
	
		TOPOCL16R1->Draw();
		TOPOCL16R1->SetNormFactor(1.0);
		TOPOCL16R1->SetLineColor(2);
		TOPOCL16R2->Draw("sames");
		TOPOCL16R2->SetNormFactor(1.0);
		TOPOCL16R2->SetLineColor(3);
		TOPOCL16R3->Draw("sames");
		TOPOCL16R3->SetNormFactor(1.0);
		TOPOCL16R3->SetLineColor(4);
		TLegend* comptopoclbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	comptopoclbyrun16->AddEntry("TOPOCL16R1","2015 group 1","L");
    	comptopoclbyrun16->AddEntry("TOPOCL16R2","2016 group 2","L");
    	comptopoclbyrun16->AddEntry("TOPOCL16R3","2016 group 3","L");
    	comptopoclbyrun16->Draw();


	TCanvas *ComparePSbyrun16 = new TCanvas("ComparePSbyrun16","ComparePSbyrun16");
	TH1F *PS16 = new TH1F ("PS16","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS16","passrndm>0.1");
	TH1F *PS16R1 = new TH1F ("PS16R1","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS16R1", PlotCut3);
	TH1F *PS16R2 = new TH1F ("PS16R2","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS16R2", PlotCut4);
	TH1F *PS16R3 = new TH1F ("PS16R3","",30,0.,30.);
		tree->Draw("sqrt(settopoclps)>>PS16R3", PlotCut5);
	
		PS16R1->Draw();
		PS16R1->SetNormFactor(1.0);
		PS16R1->SetLineColor(2);
		PS16R2->Draw("sames");
		PS16R2->SetNormFactor(1.0);
		PS16R2->SetLineColor(3);
		PS16R3->Draw("sames");
		PS16R3->SetNormFactor(1.0);
		PS16R3->SetLineColor(4);
		TLegend* comppsbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	comppsbyrun16->AddEntry("PS16R1","2015 group 1","L");
    	comppsbyrun16->AddEntry("PS16R2","2016 group 2","L");
    	comppsbyrun16->AddEntry("PS16R3","2016 group 3","L");
    	comppsbyrun16->Draw();


	TCanvas *ComparePUCbyrun16 = new TCanvas("ComparePUCbyrun16","ComparePUCbyrun16");
	TH1F *PUC16 = new TH1F ("PUC16","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC16","passrndm>0.1&&mettopoclpuc>0.1");
	TH1F *PUC16R1 = new TH1F ("PUC16R1","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC16R1", "passrndm>0.1&&runnum<299000.&&mettopoclpuc>0.1");
	TH1F *PUC16R2 = new TH1F ("PUC16R2","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC16R2", "passrndm>0.1&&runnum>299000.&&runnum<301000.&&mettopoclpuc>0.1");
	TH1F *PUC16R3 = new TH1F ("PUC16R3","",50,0.,50.);
		tree->Draw("sqrt(settopoclpuc)>>PUC16R3", "passrndm>0.1&&runnum>301000.&&mettopoclpuc>0.1");
	
		PUC16R1->Draw();
		PUC16R1->SetNormFactor(1.0);
		PUC16R1->SetLineColor(2);
		PUC16R2->Draw("sames");
		PUC16R2->SetNormFactor(1.0);
		PUC16R2->SetLineColor(3);
		PUC16R3->Draw("sames");
		PUC16R3->SetNormFactor(1.0);
		PUC16R3->SetLineColor(4);
		TLegend* comppucbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	comppucbyrun16->AddEntry("PUC16R1","2015 group 1","L");
    	comppucbyrun16->AddEntry("PUC16R2","2016 group 2","L");
    	comppucbyrun16->AddEntry("PUC16R3","2016 group 3","L");
    	comppucbyrun16->Draw();


	TCanvas *CompareMHTbyrun16 = new TCanvas("CompareMHTbyrun16","CompareMHTbyrun16");
	TH1F *MHT16 = new TH1F ("MHT16","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT16","passrndm>0.1");
	TH1F *MHT16R1 = new TH1F ("MHT16R1","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT16R1", PlotCut3);
	TH1F *MHT16R2 = new TH1F ("MHT16R2","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT16R2", PlotCut4);
	TH1F *MHT16R3 = new TH1F ("MHT16R3","",30,0.,30.);
		tree->Draw("sqrt(setmht)>>MHT16R3", PlotCut5);
	
		MHT16R1->Draw();
		MHT16R1->SetNormFactor(1.0);
		MHT16R1->SetLineColor(2);
		MHT16R2->Draw("sames");
		MHT16R2->SetNormFactor(1.0);
		MHT16R2->SetLineColor(3);
		MHT16R3->Draw("sames");
		MHT16R3->SetNormFactor(1.0);
		MHT16R3->SetLineColor(4);
		TLegend* compmhtbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	compmhtbyrun16->AddEntry("MHT16R1","2015 group 1","L");
    	compmhtbyrun16->AddEntry("MHT16R2","2016 group 2","L");
    	compmhtbyrun16->AddEntry("MHT16R3","2016 group 3","L");
    	compmhtbyrun16->Draw();


	TCanvas *CompareOFFRECALbyrun16 = new TCanvas("CompareOFFRECALbyrun16","CompareOFFRECALbyrun16");
	TH1F *OFFRECAL16 = new TH1F ("OFFRECAL16","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL16","passrndm>0.1");
	TH1F *OFFRECAL16R1 = new TH1F ("OFFRECAL16R1","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL16R1", PlotCut3);
	TH1F *OFFRECAL16R2 = new TH1F ("OFFRECAL16R2","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL16R2", PlotCut4);
	TH1F *OFFRECAL16R3 = new TH1F ("OFFRECAL16R3","",30,0.,30.);
		tree->Draw("sqrt(setoffrecal)>>OFFRECAL16R3", PlotCut5);
	
		OFFRECAL16R1->Draw();
		OFFRECAL16R1->SetNormFactor(1.0);
		OFFRECAL16R1->SetLineColor(2);
		OFFRECAL16R2->Draw("sames");
		OFFRECAL16R2->SetNormFactor(1.0);
		OFFRECAL16R2->SetLineColor(3);
		OFFRECAL16R3->Draw("sames");
		OFFRECAL16R3->SetNormFactor(1.0);
		OFFRECAL16R3->SetLineColor(4);
		TLegend* compoffrecalbyrun16=new TLegend(0.5,0.7,0.7,0.88);
    	compoffrecalbyrun16->AddEntry("OFFRECAL16R1","2015 group 1","L");
    	compoffrecalbyrun16->AddEntry("OFFRECAL16R2","2016 group 2","L");
    	compoffrecalbyrun16->AddEntry("OFFRECAL16R3","2016 group 3","L");
    	compoffrecalbyrun16->Draw();


	TCanvas *CompareL1byyear = new TCanvas("CompareL1byyear","CompareL1byyear");
		L115->Draw();
		L115->GetXaxis()->SetTitle("sqrt(L1 SumEt) [GeV]");
		L115->GetXaxis()->SetTitle("Fraction of Events");
		L115->SetNormFactor(1.0);
		L115->SetLineColor(2);
		L116->Draw("sames");
		L116->SetNormFactor(1.0);
		L116->SetLineColor(4);
		TLegend* compl1byyear=new TLegend(0.5,0.7,0.7,0.88);
    	compl1byyear->AddEntry("L115","2015","L");
    	compl1byyear->AddEntry("L116","2016","L");
    	compl1byyear->Draw();

	TCanvas *CompareCELLbyyear = new TCanvas("CompareCELLbyyear","CompareCELLbyyear");
		CELL15->Draw();
		CELL15->GetXaxis()->SetTitle("sqrt(CELL SumEt) [GeV]");
		CELL15->GetXaxis()->SetTitle("Fraction of Events");
		CELL15->SetNormFactor(1.0);
		CELL15->SetLineColor(2);
		CELL16->Draw("sames");
		CELL16->SetNormFactor(1.0);
		CELL16->SetLineColor(4);
		TLegend* compcellbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	compcellbyyear->AddEntry("CELL15","2015","L");
    	compcellbyyear->AddEntry("CELL16","2016","L");
    	compcellbyyear->Draw();
	
	TCanvas *CompareTOPOCLbyyear = new TCanvas("CompareTOPOCLbyyear","CompareTOPOCLbyyear");
		TOPOCL15->Draw();
		TOPOCL15->GetXaxis()->SetTitle("sqrt(TOPOCL SumEt) [GeV]");
		TOPOCL15->GetXaxis()->SetTitle("Fraction of Events");
		TOPOCL15->SetNormFactor(1.0);
		TOPOCL15->SetLineColor(2);
		TOPOCL16->Draw("sames");
		TOPOCL16->SetNormFactor(1.0);
		TOPOCL16->SetLineColor(4);
		TLegend* comptopoclbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	comptopoclbyyear->AddEntry("TOPOCL15","2015","L");
    	comptopoclbyyear->AddEntry("TOPOCL16","2016","L");
    	comptopoclbyyear->Draw();

	TCanvas *ComparePSbyyear = new TCanvas("ComparePSbyyear","ComparePSbyyear");
		PS15->Draw();
		PS15->GetXaxis()->SetTitle("sqrt(PS SumEt) [GeV]");
		PS15->GetXaxis()->SetTitle("Fraction of Events");
		PS15->SetNormFactor(1.0);
		PS15->SetLineColor(2);
		PS16->Draw("sames");
		PS16->SetNormFactor(1.0);
		PS16->SetLineColor(4);
		TLegend* comppsbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	comppsbyyear->AddEntry("PS15","2015","L");
    	comppsbyyear->AddEntry("PS16","2016","L");
    	comppsbyyear->Draw();

	TCanvas *ComparePUCbyyear = new TCanvas("ComparePUCbyyear","ComparePUCbyyear");
		PUC15->Draw();
		PUC15->GetXaxis()->SetTitle("sqrt(PUC SumEt) [GeV]");
		PUC15->GetXaxis()->SetTitle("Fraction of Events");
		PUC15->SetNormFactor(1.0);
		PUC15->SetLineColor(2);
		PUC16->Draw("sames");
		PUC16->SetNormFactor(1.0);
		PUC16->SetLineColor(4);
		TLegend* comppucbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	comppucbyyear->AddEntry("PUC15","2015","L");
    	comppucbyyear->AddEntry("PUC16","2016","L");
    	comppucbyyear->Draw();

	TCanvas *CompareMHTbyyear = new TCanvas("CompareMHTbyyear","CompareMHTbyyear");
		MHT15->Draw();
		MHT15->GetXaxis()->SetTitle("sqrt(MHT SumEt) [GeV]");
		MHT15->GetXaxis()->SetTitle("Fraction of Events");
		MHT15->SetNormFactor(1.0);
		MHT15->SetLineColor(2);
		MHT16->Draw("sames");
		MHT16->SetNormFactor(1.0);
		MHT16->SetLineColor(4);
		TLegend* compmhtbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	compmhtbyyear->AddEntry("MHT15","2015","L");
    	compmhtbyyear->AddEntry("MHT16","2016","L");
    	compmhtbyyear->Draw();

	TCanvas *CompareOFFRECALbyyear = new TCanvas("CompareOFFRECALbyyear","CompareOFFRECALbyyear");
		OFFRECAL15->Draw();
		OFFRECAL15->GetXaxis()->SetTitle("sqrt(OFFRECAL SumEt) [GeV]");
		OFFRECAL15->GetXaxis()->SetTitle("Fraction of Events");
		OFFRECAL15->SetNormFactor(1.0);
		OFFRECAL15->SetLineColor(2);
		OFFRECAL16->Draw("sames");
		OFFRECAL16->SetNormFactor(1.0);
		OFFRECAL16->SetLineColor(4);
		TLegend* compoffrecalbyyear=new TLegend(0.5,0.7,0.7,0.88);
    	compoffrecalbyyear->AddEntry("OFFRECAL15","2015","L");
    	compoffrecalbyyear->AddEntry("OFFRECAL16","2016","L");
    	compoffrecalbyyear->Draw();
}