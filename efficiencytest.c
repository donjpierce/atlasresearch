{
	TFile * 2015 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

	//based on line 596 of turnoncurvesE4L1cut0.c
	TCanvas *L1MET = new TCanvas("L1MET", "L1MET");
	TH1F* l1met = new TH1F("l1met", " L1 MET ", 100, 0., 100.);
	l1met->GetXaxis()->SetTitle("L1 MET [GeV]");
	l1met->GetYaxis()->SetTitle("Number of Events");
	tree->Draw("metl1>>l1met");
	gPad->Update();

	TCanvas *L1METL1CUT = new TCanvas("L1METL1CUT", "L1METL1CUT");
	TH1F* l1metl1cut = new TH1F("l1metl1cut", " L1 MET with L1 Cut", 100, 0., 100.);
	l1metl1cut->GetXaxis()->SetTitle("L1 MET [GeV]");
	l1metl1cut->GetYaxis()->SetTitle("Number of Recorded Events");
	tree->Draw("metl1>>l1metl1cut", "metl1>0.1");
	gPad->Update();

	TCanvas *L1METTHRESH = new TCanvas("L1METTHERSH", "L1METTHRESH");
	TH1F* l1metthresh = new TH1F("l1metthresh", "Fraction of events above Algorithm Threshold", 100, 0., 100.);
	l1metthresh->GetXaxis()->SetTitle("L1 MET [GeV]");
	l1metthresh->GetYaxis()->SetTitle("Fraction of events above threshold");
	Double_t* zvall1, sumitl1, ratiol1;
	sumitl1 = l1metl1cut->GetBinContent(101);
	Double_t* totnuml1 = l1metl1cut->GetEntries();
	Double_t* totnuml1wol1cut = l1met->GetEntries();
	std::cout << "L1 threshold plot \n";
	
	double ttbarl1thresh[100];
	for (int j = 100; j > 0; j--) {
		zvall1 = l1metl1cut->GetBinContent(j);
		sumitl1 += zvall1;
		ratiol1 = sumitl1 / totnuml1wol1cut;
		l1metthresh->SetBinContent(j, ratiol1);
		std::cout << "L1 Signal bin " << j << "MET [GeV] " << float(0) + float((j - 1)*(100)) / float(100) << "  fraction above threshold " << ratiol1 << "\n";
		ttbarl1thresh[j - 1] = ratiol1;
	l1metthresh->Draw();
	gPad->Update();

}