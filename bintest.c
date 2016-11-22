{
	
	//Analyze three bins in the 2015 runs
	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");


	TH1F *L115R1 = new TH1F("L115R1", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R1", "passrndm>0.1&&sqrt(setl1)>2.&&sqrt(setl1)<3.");
	TH1F *L115R2 = new TH1F("L115R2", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R2", "passrndm>0.1&&sqrt(setl1)>4.&&sqrt(setl1)<5.");
	TH1F *L115R3 = new TH1F("L115R3", "", 30, 0., 30.);
	tree->Draw("metl1>>L115R3", "passrndm>0.1&&sqrt(setl1)>9.&&sqrt(setl1)<10.");



	//Analyze three bins in the 2016 runs
	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");

	TH1F *L116R1 = new TH1F("L116R1", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R1", "passrndm>0.1&&sqrt(setl1)>1.&&sqrt(setl1)<2.");
	TH1F *L116R2 = new TH1F("L116R2", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R2", "passrndm>0.1&&sqrt(setl1)>4.&&sqrt(setl1)<5.");
	TH1F *L116R3 = new TH1F("L116R3", "", 30, 0., 30.);
	tree->Draw("metl1>>L116R3", "passrndm>0.1&&sqrt(setl1)>10.&&sqrt(setl1)<11.");
	


	//Declare the Rayleigh Distribution as func
	TF1 *func = new TF1("func", "(x/[1])*(1/[1])*exp(-.5*(x/[1])*(x/[1]))");
	//func->SetParameters(0,100000.);
	func->SetParameters(1, 5.);
	func->SetParLimits(1, 0.1, 100.);

	//Begin L1 Graphs
	TCanvas *CompareL1R1 = new TCanvas("CompareL1R1", "L1");
	L116R1->SetNormFactor(1.0);
	L116R1->SetLineColor(3);
	func->SetLineColor(3);
	L116R1->Draw();
	L116R1->Fit("func", "L");
	L115R1->SetNormFactor(1.0);
	L115R1->Draw("sames");
	L115R1->SetLineColor(2);
	func->SetLineColor(2);
	L115R1->Fit("func", "L");
	
	TCanvas *CompareL1R2 = new TCanvas("CompareL1R2", "L1");
	L115R2->Draw();
	L115R2->SetTitle("Comparing L1 at Medium sqrt(Sum Et)")
	L115R2->SetLineColor(6);
	L115R2->SetNormFactor(1.0);
	func->SetLineColor(6);
	L115R2->Fit("func", "L");
	L116R2->Draw("sames");
	L116R2->SetLineColor(4);
	L116R2->SetNormFactor(1.0);
	func->SetLineColor(4);
	L116R2->Fit("func", "L");
	
	TCanvas *CompareL1R3 = new TCanvas("CompareL1R3", "L1");
	L115R3->Draw();
	L115R3->SetTitle("Comparing L1 at High sqrt(Sum Et)")
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


	
}