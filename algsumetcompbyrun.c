{
	//2015 Plot Cuts
	TString PlotCut("passrndm>0.1&&runnum<283200.");
	TString PlotCut1("passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TString PlotCut2("passrndm>0.1&&runnum>283800.");
	//2016 Plot Cuts
	TString PlotCut3("passrndm>0.1&&runnum<299000.");
	TString PlotCut4("passrndm>0.1&&runnum>299000.&&runnum<301000.");
	TString PlotCut5("passrndm>0.1&&runnum>301000.");

	//Gather 2015 Graphics
	TFile *2015 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
		
	TH1F *SumEt15L1R1 = new TH1F("SumEt15L1", "", 10, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEtL1R1", PlotCut);
	TH1F *SumEt15L1R2 = new TH1F("SumEt15L1R2", "", 10, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEt15L1R2", PlotCut1);
	TH1F *SumEt15L1R3 = new TH1F("SumEt15L1R3", "", 10, 0., 10.);
	tree->Draw("sqrt(setl1)>>SumEt15L1R3", PlotCut2);

	TH1F *SumEt15CellR1 = new TH1F("SumEt15CellL1", "", 21, 4., 25.);
	tree->Draw("sqrt(setcell)>>SumEt15CellR1", PlotCut);
	TH1F *SumEt15CellR2 = new TH1F("SumEt15CellR2", "", 21, 4., 25.);
	tree->Draw("sqrt(setcell)>>SumEtCellR2", PlotCut1);
	TH1F *SumEt15CellR3 = new TH1F("SumEt15CellR3", "", 21, 4., 25.);
	tree->Draw("sqrt(setcell)>>SumEt15CellR3", PlotCut2);

	TH1F *SumEt15MHTR1 = new TH1F("SumEt15MHTL1", "", 10, 0., 20.);
	tree->Draw("sqrt(setmht)>>SumEt15MHTR1", PlotCut);
	TH1F *SumEt15MHTR2 = new TH1F("SumEt15MHTR2", "", 20, 0., 20.);
	tree->Draw("sqrt(setmht)>>SumEt15MHTR2", PlotCut1);
	TH1F *SumEt15MHTR3 = new TH1F("SumEt15MHTR3", "", 20, 0., 20.);
	tree->Draw("sqrt(setmht)>>SumEt15MHTR3", PlotCut2);

	TH1F *SumEt15topoclR1 = new TH1F("SumEt15topoclL1", "", 27, 0., 27.);
	tree->Draw("sqrt(settopocl)>>SumEt15topoclR1", PlotCut);
	TH1F *SumEt15topoclR2 = new TH1F("SumEt15topoclR2", "", 27, 0., 27.);
	tree->Draw("sqrt(settopocl)>>SumEt15topoclR2", PlotCut1);
	TH1F *SumEt15topoclR3 = new TH1F("SumEt15topoclR3", "", 27, 0., 27.);
	tree->Draw("sqrt(settopocl)>>SumEt15topoclR3", PlotCut2);

	TH1F *SumEt15topoclpsR1 = new TH1F("SumEt15topoclpsL1", "", 26, 4., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt15topoclpsR1", PlotCut);
	TH1F *SumEt15topoclpsR2 = new TH1F("SumEt15topoclpsR2", "", 26, 4., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt15topoclpsR2", PlotCut1);
	TH1F *SumEt15topoclpsR3 = new TH1F("SumEt15topoclpsR3", "", 26, 4., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt15topoclpsR3", PlotCut2);

	TH1F *SumEt15topoclpucR1 = new TH1F("SumEt15topoclpucL1", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt15topoclpucR1", "mettopoclpuc>0.1&&passrndm>0.1&&runnum<283200.");
	TH1F *SumEt15topoclpucR2 = new TH1F("SumEt15topoclpucR2", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt15topoclpucR2", "mettopoclpuc>0.1&&passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TH1F *SumEt15topoclpucR3 = new TH1F("SumEt15topoclpucR3", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt15topoclpucR3", "mettopoclpuc>0.1&&passrndm>0.1&&runnum>283800.");

	TH1F *SumEt15offrecalR1 = new TH1F("SumEt15offrecalL1", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt15offrecalR1", PlotCut);
	TH1F *SumEtoffrecalR2 = new TH1F("SumEt15topocloffrecalR2", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt15offrecalR2", PlotCut1);
	TH1F *SumEt15offrecalR3 = new TH1F("SumEt15offrecalR3", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt15offrecalR3", PlotCut2);

	//Gather 2016 Graphics
	TFile *2016 = TFile::Open("ZeroBias2016.13Runs.root");

	TH1F *SumEt16L1R1 = new TH1F("SumEt16L1", "", 30, 0., 30.);
	tree->Draw("sqrt(setl1)>>SumEt16L1R1", PlotCut);
	TH1F *SumEt16L1R2 = new TH1F("SumEt16L1R2", "", 30, 0., 30.);
	tree->Draw("sqrt(setl1)>>SumEt16L1R2", PlotCut1);
	TH1F *SumEtL1R3 = new TH1F("SumEt16L1R3", "", 30, 0., 30.);
	tree->Draw("sqrt(setl1)>>SumEt16L1R3", PlotCut2);

	TH1F *SumEtCellR1 = new TH1F("SumCellL1", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEtCellR1", PlotCut);
	TH1F *SumEtCellR2 = new TH1F("SumEtCellR2", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEtCellR2", PlotCut1);
	TH1F *SumEtCellR3 = new TH1F("SumEtCellR3", "", 30, 0., 30.);
	tree->Draw("sqrt(setcell)>>SumEtCellR3", PlotCut2);

	TH1F *SumEt16MHTR1 = new TH1F("SumEt16MHTL1", "", 30, 0., 30.);
	tree->Draw("sqrt(setmht)>>SumEt16MHTR1", PlotCut);
	TH1F *SumEt16MHTR2 = new TH1F("SumEt16MHTR2", "", 30, 0., 30.);
	tree->Draw("sqrt(setmht)>>SumEt16MHTR2", PlotCut1);
	TH1F *SumEt16MHTR3 = new TH1F("SumEt16MHTR3", "", 30, 0., 30.);
	tree->Draw("sqrt(setmht)>>SumEt16MHTR3", PlotCut2);

	TH1F *SumEt16topoclR1 = new TH1F("SumEt16topoclL1", "", 30, 0., 30.);
	tree->Draw("sqrt(settopocl)>>SumEt16topoclR1", PlotCut);
	TH1F *SumEt16topoclR2 = new TH1F("SumEt16topoclR2", "", 30, 0., 30.);
	tree->Draw("sqrt(settopocl)>>SumEt16topoclR2", PlotCut1);
	TH1F *SumEt16topoclR3 = new TH1F("SumEt16topoclR3", "", 30, 0., 30.);
	tree->Draw("sqrt(settopocl)>>SumEt16topoclR3", PlotCut2);

	TH1F *SumEt16topoclpsR1 = new TH1F("SumEt16topoclpsL1", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt16topoclpsR1", PlotCut);
	TH1F *SumEt16topoclpsR2 = new TH1F("SumEt16topoclpsR2", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt16topoclpsR2", PlotCut1);
	TH1F *SumEt16topoclpsR3 = new TH1F("SumEt16topoclpsR3", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclps)>>SumEt16topoclpsR3", PlotCut2);

	TH1F *SumEt16topoclpucR1 = new TH1F("SumEt16topoclpucL1", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt16topoclpucR1", PlotCut);
	TH1F *SumEt16topoclpucR2 = new TH1F("SumEt16topoclpucR2", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt16topoclpucR2", PlotCut1);
	TH1F *SumEt16topoclpucR3 = new TH1F("SumEt16topoclpucR3", "", 30, 0., 30.);
	tree->Draw("sqrt(settopoclpuc)>>SumEt16topoclpucR3", PlotCut2);

	TH1F *SumEt16offrecalR1 = new TH1F("SumEt16offrecalL1", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt16offrecalR1", PlotCut);
	TH1F *SumEttopocloffrecalR2 = new TH1F("SumEt16offrecalR2", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt16offrecalR2", PlotCut1);
	TH1F *SumEttopocloffrecalR3 = new TH1F("SumEt16offrecalR3", "", 30, 0., 30.);
	tree->Draw("sqrt(setoffrecal)>>SumEt16offrecalR3", PlotCut2);


	//Draw 2015 Graphics
	TCanvas *SumEt15L1Run = new TCanvas("SumEt15L1Run", "2015 L1 Runs 1 2 3");
		SumEt15L1R1->SetTitle("L1 in Three Subdivisions of 2015 Runs");
		SumEt15L1R1->GetXaxis()->SetTitle("Sqrt SumEt of L1");
		SumEt15L1R1->SetNormFactor(1.0);
		SumEt15L1R2->SetNormFactor(1.0);
		SumEt15L1R3->SetNormFactor(1.0);
		SumEt15L1R1->SetLineColor(2);
		SumEt15L1R2->SetLineColor(3);
		SumEt15L1R3->SetLineColor(4);
		SumEt15L1R1->Draw();
		SumEt15L1R2->Draw("sames");
		SumEt15L1R3->Draw("sames");

	TCanvas *SumEt15CellRun = new TCanvas("SumEt15CellRun", "2015 Cell Runs 1 2 3");
		SumEt15CellR1->SetTitle("Cell in Three Subdivisions of 2015 Runs");
		SumEt15CellR1->GetXaxis()->SetTitle("Sqrt SumEt of Cell");
		SumEt15CellR1->SetNormFactor(1.0);
		SumEt15CellR2->SetNormFactor(1.0);
		SumEt15CellR3->SetNormFactor(1.0);
		SumEt15CellR1->SetLineColor(2);
		SumEt15CellR2->SetLineColor(3);
		SumEt15CellR3->SetLineColor(4);
		SumEt15CellR1->Draw();
		SumEt15CellR2->Draw("sames");
		SumEt15CellR3->Draw("sames");

	TCanvas *SumEt15MHTRun = new TCanvas("SumEt15MHTRun", "2015 MHT Runs 1 2 3");
		SumEt15MHTR1->SetTitle("MHT in Three Subdivisions of 2015 Runs");
		SumEt15MHTR1->GetXaxis()->SetTitle("Sqrt SumEt of MHT");
		SumEt15MHTR1->SetNormFactor(1.0);
		SumEt15MHTR2->SetNormFactor(1.0);
		SumEt15MHTR3->SetNormFactor(1.0);
		SumEt15MHTR1->SetLineColor(2);
		SumEt15MHTR2->SetLineColor(3);
		SumEt15MHTR3->SetLineColor(4);
		SumEt15MHTR1->Draw();
		SumEt15MHTR2->Draw("sames");
		SumEt15MHTR3->Draw("sames");

	TCanvas *SumEt15topoclRun = new TCanvas("SumEt15topoclRun", "2015 topocl Runs 1 2 3");
		SumEt15topoclR1->SetTitle("topocl in Three Subdivisions of 2015 Runs");
		SumEt15topoclR1->GetXaxis()->SetTitle("Sqrt SumEt of topocl");
		SumEt15topoclR1->SetNormFactor(1.0);
		SumEt15topoclR2->SetNormFactor(1.0);
		SumEt15topoclR3->SetNormFactor(1.0);
		SumEt15topoclR1->SetLineColor(2);
		SumEt15topoclR2->SetLineColor(3);
		SumEt15topoclR3->SetLineColor(4);
		SumEt15topoclR1->Draw();
		SumEt15topoclR2->Draw("sames");
		SumEt15topoclR3->Draw("sames");

	TCanvas *SumEt15topoclpsRun = new TCanvas("SumEt15topoclpsRun", "2015 topoclps Runs 1 2 3");
		SumEt15topoclpsR1->SetTitle("topoclps in Three Subdivisions of 2015 Runs");
		SumEt15topoclpsR1->GetXaxis()->SetTitle("Sqrt SumEt of topoclps");
		SumEt15topoclpsR1->SetNormFactor(1.0);
		SumEt15topoclpsR2->SetNormFactor(1.0);
		SumEt15topoclpsR3->SetNormFactor(1.0);
		SumEt15topoclpsR1->SetLineColor(2);
		SumEt15topoclpsR2->SetLineColor(3);
		SumEt15topoclpsR3->SetLineColor(4);
		SumEt15topoclpsR1->Draw();
		SumEt15topoclpsR2->Draw("sames");
		SumEt15topoclpsR3->Draw("sames");

	TCanvas *SumEt15topoclpucRun = new TCanvas("SumEt15topoclpucRun", "2015 topoclpuc Runs 1 2 3");
		SumEt15topoclpucR1->SetTitle("topoclpuc in Three Subdivisions of 2015 Runs");
		SumEt15topoclpucR1->GetXaxis()->SetTitle("Sqrt SumEt of topoclpuc");
		SumEt15topoclpucR1->SetNormFactor(1.0);
		SumEt15topoclpucR2->SetNormFactor(1.0);
		SumEt15topoclpucR3->SetNormFactor(1.0);
		SumEt15topoclpucR1->SetLineColor(2);
		SumEt15topoclpucR2->SetLineColor(3);
		SumEt15topoclpucR3->SetLineColor(4);
		SumEt15topoclpucR1->Draw();
		SumEt15topoclpucR2->Draw("sames");
		SumEt15topoclpucR3->Draw("sames");

	TCanvas *SumEt15offrecalRun = new TCanvas("SumEt15offrecalRun", "2015 offrecal Runs 1 2 3");
		SumEt15offrecalR1->SetTitle("offrecal in Three Subdivisions of 2015 Runs");
		SumEt15offrecalR1->GetXaxis()->SetTitle("Sqrt SumEt of offrecal");
		SumEt15offrecalR1->SetNormFactor(1.0);
		SumEt15offrecalR2->SetNormFactor(1.0);
		SumEt15offrecalR3->SetNormFactor(1.0);
		SumEt15offrecalR1->SetLineColor(2);
		SumEt15offrecalR2->SetLineColor(3);
		SumEt15offrecalR3->SetLineColor(4);
		SumEt15offrecalR1->Draw();
		SumEt15offrecalR2->Draw("sames");
		SumEt15offrecalR3->Draw("sames");
}