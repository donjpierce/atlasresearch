{

//2015 (Period J) Bunch Crossing IDs
	TFile *2015J = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	
	TCanvas *15bcids = new TCanvas("15bcids", "2015 Period J Bunch Crossing IDs");
	TH1F *2015bcid = new TH1F("2015bcid", "2015 Period J Bunch Crossing IDs", 3500, 0., 3500.);
	tree->Draw("bcid>>2015bcid");
	/*2015bcid->SetTitle("2015 Period J Bunch Crossings");
	2015bcid->GetXaxis()->SetTitle("Bunch Crossing ID Number");
	2015bcid->GetYAxis()->SetTitle("Number of Events");*/

	/* Compare 2015 L1 vs sqrt(setl1)*/
	TCanvas *2015L1 = new TCanvas("2015L1", "2015L1");
	TH2F *2015L1vsSET = new TH2F("2015L1vsSET", "2015 L1 MET vs sqrtSET", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSET", "colz");


//2016 (13 Runs) Bunch Crossing IDs
	TFile *2016 = TFile::Open("ZeroBias2016.13Runs.root");

	TCanvas *16bcids = new TCanvas("16bcids", "2016 Bunch Crossing IDs");
	TH1F *2016bcid = new TH1F("2016bcid", "2016 Bunch Crossing IDs", 3500, 0., 3500.);
	tree->Draw("bcid>>2016bcid");
	/*2016bcid->SetTitle("2016 Bunch Crossings");
	2016bcid->GetXaxis()->SetTitle("Bunch Crossing ID Number");
	2016bcid->GetYAxis()->SetTitle("Number of Events");*/

	/* Compare 2016 L1 vs sqrt(setl1)*/
	TCanvas *2016L1 = new TCanvas("2016L1", "2016L1");
	TH2F *2016L1vsSET = new TH2F("2016L1vsSET", "2016 L1 MET vs sqrtSET", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSET");

	





}