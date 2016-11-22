{
	//2015 Plot Cuts
	TString PlotCut("passrndm>0.1&&runnum<283200.");
	TString PlotCut1("passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TString PlotCut2("passrndm>0.1&&runnum>283800.");

	//Gather 2015 Histograms
	TFile *2015 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

	TH1F *SumEt15L1R1 = new TH1F("SumEt15L1", "", 30, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEtL1R1", PlotCut);
	TH1F *SumEt15L1R2 = new TH1F("SumEt15L1R2", "", 30, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEt15L1R2", PlotCut1);
	TH1F *SumEt15L1R3 = new TH1F("SumEt15L1R3", "", 30, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEt15L1R3", PlotCut2);

	TH1F *SumEt15CellR1 = new TH1F("SumEt15CellL1", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEt15CellR1", PlotCut);
	TH1F *SumEt15CellR2 = new TH1F("SumEt15CellR2", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEtCellR2", PlotCut1);
	TH1F *SumEt15CellR3 = new TH1F("SumEt15CellR3", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEt15CellR3", PlotCut2);

	//Declare the Rayleigh Distribution as func
	TF1 *func = new TF1("func", "(x/[1])*(1/[1])*exp(-.5*(x/[1])*(x/[1]))");
	//func->SetParameters(0,100000.);
	func->SetParameters(1, 5.);
	func->SetParLimits(1, 0.1, 100.);

	//Draw 2015 Graphics
	TCanvas *SumEt15L1Run = new TCanvas("SumEt15L1Run", "L1 Tridivisions in 2015");
		SumEt15L1R3->SetNormFactor(1.0);
		SumEt15L1R3->SetLineColor(4);
		SumEt15L1R3->Draw();
		SumEt15L1R2->SetNormFactor(1.0);
		SumEt15L1R2->SetLineColor(3);
		SumEt15L1R2->Draw("sames");
		SumEt15L1R1->SetNormFactor(1.0);
		SumEt15L1R1->SetLineColor(2);
		SumEt15L1R1->Draw("sames");

	TCanvas *SumEt15CellRun = new TCanvas("SumEt15cellRun", "cell Tridivisions in 2015");s
		SumEt15CellR1->SetNormFactor(1.0);
		SumEt15CellR2->SetNormFactor(1.0);
		SumEt15CellR3->SetNormFactor(1.0);
		SumEt15CellR1->SetLineColor(2);
		SumEt15CellR2->SetLineColor(3);
		SumEt15CellR3->SetLineColor(4);
		SumEt15CellR3->Draw();
		SumEt15CellR2->Draw("sames");
		SumEt15CellR1->Draw("sames");


}