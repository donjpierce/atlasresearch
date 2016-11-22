{
	//Analyze three bins in the 2015 runs
	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

	TH1F *L115R1 = new TH1F("L115R1", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R1", "passrndm>0.1&&sqrt(setl1)>0.1&&sqrt(setl1)<3.");
	TH1F *L115R2 = new TH1F("L115R2", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R2", "passrndm>0.1&&sqrt(setl1)>3.&&sqrt(setl1)<7.");
	TH1F *L115R3 = new TH1F("L115R3", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R3", "passrndm>0.1&&sqrt(setl1)>7.&&sqrt(setl1)<12.");

	TH1F *Cell15R1 = new TH1F("Cell15R1", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell15R1", "passrndm>0.1&&sqrt(setcell)>5.&&sqrt(setcell)<14.");
	TH1F *Cell15R2 = new TH1F("Cell15R2", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell15R2", "passrndm>0.1&&sqrt(setcell)>14.&&sqrt(setcell)<21.");
	TH1F *Cell15R3 = new TH1F("Cell15R3", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell15R3", "passrndm>0.1&&sqrt(setcell)>21.&&sqrt(setcell)<30.");

	TH1F *Topocl15R1 = new TH1F("Topocl15R1", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl15R1", "passrndm>0.1&&sqrt(settopocl)>5.&&sqrt(settopocl)<20.");
	TH1F *Topocl15R2 = new TH1F("Topocl15R2", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl15R2", "passrndm>0.1&&sqrt(settopocl)>20.&&sqrt(settopocl)<30.");
	TH1F *Topocl15R3 = new TH1F("Topocl15R3", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl15R3", "passrndm>0.1&&sqrt(settopocl)>30.&&sqrt(settopocl)<45.");

	//Analyze three bins in the 2016 runs
	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");

	TH1F *L116R1 = new TH1F("L116R1", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R1", "passrndm>0.1&&sqrt(setl1)>0.1&&sqrt(setl1)<3.");
	TH1F *L116R2 = new TH1F("L116R2", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R2", "passrndm>0.1&&sqrt(setl1)>3.&&sqrt(setl1)<7.");
	TH1F *L116R3 = new TH1F("L116R3", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R3", "passrndm>0.1&&sqrt(setl1)>7.&&sqrt(setl1)<12.");

	TH1F *Cell16R1 = new TH1F("Cell16R1", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell16R1", "passrndm>0.1&&sqrt(setcell)>5.&&sqrt(setcell)<14.");
	TH1F *Cell16R2 = new TH1F("Cell16R2", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell16R2", "passrndm>0.1&&sqrt(setcell)>14.&&sqrt(setcell)<21.");
	TH1F *Cell16R3 = new TH1F("Cell16R3", "", 30, 0., 30.);
	tree->Draw("metcell>>Cell16R3", "passrndm>0.1&&sqrt(setcell)>21.&&sqrt(setcell)<30.");

	TH1F *Topocl16R1 = new TH1F("Topocl16R1", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl16R1", "passrndm>0.1&&sqrt(settopocl)>5.&&sqrt(settopocl)<20.");
	TH1F *Topocl16R2 = new TH1F("Topocl16R2", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl16R2", "passrndm>0.1&&sqrt(settopocl)>20.&&sqrt(settopocl)<30.");
	TH1F *Topocl16R3 = new TH1F("Topocl16R3", "", 30, 0., 30.);
	tree->Draw("mettopocl>>Topocl16R3", "passrndm>0.1&&sqrt(settopocl)>30.&&sqrt(settopocl)<45.");

	//Declare the Rayleigh Distribution as func
	TF1 *func = new TF1("func", "(x/[1])*(1/[1])*exp(-.5*(x/[1])*(x/[1]))");
	//func->SetParameters(0,100000.);
	func->SetParameters(1, 5.);
	func->SetParLimits(1, 0.1, 100.);

	//Begin L1 Graphs
	TCanvas *CompareL1R1 = new TCanvas("CompareL1R1", "L1: 0<L1 2015<3 , 0<L1 2016<4");
	L115R1->Draw();
	L115R1->SetTitle("Comparing L1 in 2015 and 2016 when there is low sqrt Sum Et");
	L115R1->GetXaxis()->SetTitle("sqrt Sum Et");
	L115R1->GetYaxis()->SetTitle("Normalized Fraction of Events");
	L115R1->SetLineColor(2);
	L115R1->SetNormFactor(1.0);
	func->SetLineColor(2);
	L115R1->Fit("func", "L");
	L116R1->Draw("sames");
	L116R1->SetLineColor(3);
	L116R1->SetNormFactor(1.0);
	func->SetLineColor(3);
	L116R1->Fit("func", "L");
	TCanvas *CompareL1R2 = new TCanvas("CompareL1R2", "L1: 14<L1 2015<21, 4<L1 2016<8");
	L115R2->Draw();
	L115R2->SetTitle("Comparing L1 in 2015 and 2016 when there is medium sqrt(Sum Et)")
	L115R2->SetLineColor(6);
	L115R2->SetNormFactor(1.0);
	func->SetLineColor(6);
	L115R2->Fit("func", "L");
	L116R2->Draw("sames");
	L116R2->SetLineColor(4);
	L116R2->SetNormFactor(1.0);
	func->SetLineColor(4);
	L116R2->Fit("func", "L");
	TCanvas *CompareL1R3 = new TCanvas("CompareL1R3", "L1: 8<L1 2015<12, 20<L1 2016<30");
	L115R3->Draw();
	L115R3->SetTitle("Comparing L1 in 2015 and 2016 when there is high sqrt(Sum Et)")
	L115R3->SetLineColor(44);
	L115R3->SetNormFactor(1.0);
	func->SetLineColor(44);
	L115R3->Fit("func", "L");
	L116R3->Draw("sames");
	L116R3->SetLineColor(7);
	L116R3->SetNormFactor(1.0);
	func->SetLineColor(7);
	L116R3->Fit("func", "L");
	gStyle->SetOptFit(11);

	//Begin Cell Graphs
	TCanvas *CompareCellR1 = new TCanvas("CompareCellR1", "Cell: 5<Cell 2015<14 , 5<Cell 2016<12");
	Cell15R1->Draw();
	Cell15R1->SetTitle("Comparing Cell in 2015 and 2016 when there is low sqrt(Sum Et)")
	Cell15R1->SetLineColor(2);
	Cell15R1->SetNormFactor(1.0);
	func->SetLineColor(2);
	Cell15R1->Fit("func", "L");
	Cell16R1->Draw("sames");
	Cell16R1->SetLineColor(3);
	Cell16R1->SetNormFactor(1.0);
	func->SetLineColor(3);
	Cell16R1->Fit("func", "L");
	TCanvas *CompareCellR2 = new TCanvas("CompareCellR2", "Cell: 14<Cell 2015<21, 12<Cell 2016<20");
	Cell15R2->Draw();
	Cell15R2->SetTitle("Comparing Cell in 2015 and 2016 when there is medium sqrt(Sum Et)")
	Cell15R2->SetLineColor(6);
	Cell15R2->SetNormFactor(1.0);
	func->SetLineColor(6);
	Cell15R2->Fit("func", "L");
	Cell16R2->Draw("sames");
	Cell16R2->SetLineColor(4);
	Cell16R2->SetNormFactor(1.0);
	func->SetLineColor(4);
	Cell16R2->Fit("func", "L");
	TCanvas *CompareCellR3 = new TCanvas("CompareCellR3", "Cell: 21<Cell 2015<30, 20<Cell 2016<30");
	Cell15R3->Draw();
	Cell15R3->SetTitle("Comparing Cell in 2015 and 2016 when there is high sqrt(Sum Et)")
	Cell15R3->SetLineColor(44);
	Cell15R3->SetNormFactor(1.0);
	func->SetLineColor(44);
	Cell15R3->Fit("func", "L");
	Cell16R3->Draw("sames");
	Cell16R3->SetLineColor(7);
	Cell16R3->SetNormFactor(1.0);
	func->SetLineColor(7);
	Cell16R3->Fit("func", "L");
	gStyle->SetOptFit(11);

	//Begin Topocl Graphs
	TCanvas *CompareTopoclR1 = new TCanvas("CompareTopoclR1", "Topocl: 7<Topocl 2015<12 , 5<Topocl 2016<20");
	Topocl15R1->Draw();
	Topocl15R1->SetTitle("Comparing Topocl in 2015 and 2016 when there is low sqrt(Sum Et)")
	Topocl15R1->SetLineColor(2);
	Topocl15R1->SetNormFactor(1.0);
	func->SetLineColor(2);
	Topocl15R1->Fit("func", "L");
	Topocl16R1->Draw("sames");
	Topocl16R1->SetLineColor(3);
	Topocl16R1->SetNormFactor(1.0);
	func->SetLineColor(3);
	Topocl16R1->Fit("func", "L");
	TCanvas *CompareTopoclR2 = new TCanvas("CompareTopoclR2", "Topocl: 12<Topocl 2015<21, 20<Topocl 2016<30");
	Topocl15R2->Draw();
	Topocl15R2->SetTitle("Comparing Topocl in 2015 and 2016 when there is medium sqrt(Sum Et)")
	Topocl15R2->SetLineColor(6);
	Topocl15R2->SetNormFactor(1.0);
	func->SetLineColor(6);
	Topocl15R2->Fit("func", "L");
	Topocl16R2->Draw("sames");
	Topocl16R2->SetLineColor(4);
	Topocl16R2->SetNormFactor(1.0);
	func->SetLineColor(4);
	Topocl16R2->Fit("func", "L");
	TCanvas *CompareTopoclR3 = new TCanvas("CompareTopclR3", "Cell: 21<Topocl 2015<32, 30<Topocl 2016<45");
	Topocl15R3->Draw();
	Topocl15R3->SetTitle("Comparing Topocl in 2015 and 2016 when there is high sqrt(Sum Et)")
	Topocl15R3->SetLineColor(44);
	Topocl15R3->SetNormFactor(1.0);
	func->SetLineColor(44);
	Topocl15R3->Fit("func", "L");
	Topocl16R3->Draw("sames");
	Topocl16R3->SetLineColor(7);
	Topocl16R3->SetNormFactor(1.0);
	func->SetLineColor(7);
	Topocl16R3->Fit("func", "L");
	gStyle->SetOptFit(11);
}
