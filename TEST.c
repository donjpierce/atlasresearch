{
	TString PlotCut ("passrndm>0.1&&runnum<283200.");
	TString PlotCut1 ("passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TString PlotCut2 ("passrndm>0.1&&runnum>283800.");
	TString PlotCut3 ("passrndm>0.1&&runnum<299000.");
	TString PlotCut4 ("passrndm>0.1&&runnum>299000.&&runnum<301000.");
	TString PlotCut5 ("passrndm>0.1&&runnum>301000.");

	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	TH2F *TOPOCL15R1 = new TH2F ("TOPOCL15R1","",40,15.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL15R1",PlotCut);
	TH2F *TOPOCL15R2 = new TH2F ("TOPOCL15R2","",40,0.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL15R2",PlotCut1);
	TH2F *TOPOCL15R3 = new TH2F ("TOPOCL15R3","",40,0.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL15R3",PlotCut2);

	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");
	TH2F *TOPOCL16R1 = new TH2F ("TOPOCL16R1","",40,0.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL16R1",PlotCut3);
	TH2F *TOPOCL16R2 = new TH2F ("TOPOCL16R2","",40,0.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL16R2",PlotCut4);
	TH2F *TOPOCL16R3 = new TH2F ("TOPOCL16R3","",40,0.,40.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopocl)>>TOPOCL16R3",PlotCut5);


	TF1 *func = new TF1 ("func","[0]*(x/[1])*(1/[1])*exp(-.5*(x/[1])*(x/[1]))");
		func->SetParameters(0,100000.);
    	func->SetParameters(1,5.);
    	func->SetParLimits(1,0.1,100.);


   	TCanvas *CompareTOPOCL = new TCanvas ("CompareTOPOCL","CompareTOPOCL");
	TOPOCL15R1->Draw();
	TOPOCL16R3->FitSlicesY(func,0,-1,10,"L");
	TOPOCL16R2->FitSlicesY(func,0,-1,10,"L");
	TOPOCL16R1->FitSlicesY(func,0,-1,10,"L");
	TOPOCL15R3->FitSlicesY(func,0,-1,10,"L");
	TOPOCL15R2->FitSlicesY(func,0,-1,10,"L");
	TOPOCL15R1->FitSlicesY(func,0,-1,10,"L");
	
	TOPOCL15R1_1->Draw();
	TOPOCL15R1_1->SetTitle("Compare TOPOCL vs. Sum Et Across Runs");
	TOPOCL15R1_1->GetXaxis()->SetTitle("Sum Et [GeV]");
	TOPOCL15R1_1->GetYaxis()->SetTitle("TOPOCL [GeV]");
	TOPOCL15R1_1->SetLineColor(2);
	pol1->SetLineColor(2);
	TOPOCL15R1_1->Fit("pol1");
	TOPOCL15R2_1->Draw("sames");
	TOPOCL15R2_1->SetLineColor(3);
	pol1->SetLineColor(3);
	TOPOCL15R2_1->Fit("pol1");
	TOPOCL15R3_1->Draw("sames");
	TOPOCL15R3_1->SetLineColor(4);
	pol1->SetLineColor(4);
	TOPOCL15R3_1->Fit("pol1");
	TOPOCL16R1_1->Draw("sames");
	TOPOCL16R1_1->SetLineColor(44);
	pol1->SetLineColor(44);
	TOPOCL16R1_1->Fit("pol1");
	TOPOCL16R2_1->Draw("sames");
	TOPOCL16R2_1->SetLineColor(6);
	pol1->SetLineColor(6);
	TOPOCL16R2_1->Fit("pol1");
	TOPOCL16R3_1->Draw("sames");
	TOPOCL16R3_1->SetLineColor(7);
	pol1->SetLineColor(7);
	TOPOCL16R3_1->Fit("pol1");
	gStyle->SetOptFit(11);

}