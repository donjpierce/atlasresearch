 {

//2015 (Period J) Bunch Crossing IDs
	TFile *2015J = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	
	TCanvas *15bcids = new TCanvas("15bcids", "2015 Period J Bunch Crossing IDs");
	TH1F *2015bcid = new TH1F("2015bcid", "2015 Period J Bunch Crossing IDs", 3500, 0., 3500.);
	tree->Draw("bcid>>2015bcid");

	/* Compare 2015 [algorithm] MET vs sqrt([algorithm] SET) where BCID>570 (bcid=1000 is the start of the 11th bunch group in 2015)*/
	TCanvas *2015L1 = new TCanvas("2015L1", "2015L1");
	gPad->SetLogz();
	TH2F *2015L1vsSET = new TH2F("2015L1vsSET", "2015 L1 MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSET", "bcid>570.", "colz");

	TCanvas * 2015CELL = new TCanvas("2015CELL", "2015CELL");
	gPad->SetLogz();
	TH2F * 2015CELLvsSET = new TH2F("2015CELLvsSET", "2015 CELL MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015CELLvsSET", "bcid>570.", "colz");

	TCanvas * 2015MHT = new TCanvas("2015MHT", "2015MHT");
	gPad->SetLogz();
	TH2F * 2015MHTvsSET = new TH2F("2015MHTvsSET", "2015 MHT MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSET", "bcid>570.", "colz");

	TCanvas * 2015topoclps = new TCanvas("2015topoclps", "2015topoclps");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSET = new TH2F("2015topoclpsvsSET", "2015 TOPOCL PS MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSET", "bcid>570.", "colz");

	TCanvas * 2015topoclpuc = new TCanvas("2015topoclpuc", "2015topoclpuc");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSET = new TH2F("2015topoclpucvsSET", "2015 TOPOCL PUC MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSET", "bcid>570.", "colz");

	TCanvas * 2015topocl = new TCanvas("2015topocl", "2015topocl");
	gPad->SetLogz();
	TH2F * 2015topoclvsSET = new TH2F("2015topoclvsSET", "2015 TOPOCL MET vs sqrtSET for BCID > 570", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSET", "bcid>570.", "colz");



	/*Topocl PUC Bunch Group Analysis*/
	TCanvas * 2015topoclpuca = new TCanvas("2015topoclpuca", "2015topoclpuca");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETa = new TH2F("2015topoclpucvsSETa", "2015 topoclpuc MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015topoclpucb = new TCanvas("2015topoclpucb", "2015topoclpucb");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETb = new TH2F("2015topoclpucvsSETb", "2015 topoclpuc MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015topoclpucc = new TCanvas("2015topoclpucc", "2015topoclpucc");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETc = new TH2F("2015topoclpucvsSETc", "2015 topoclpuc MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015topoclpucd = new TCanvas("2015topoclpucd", "2015topoclpucd");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETd = new TH2F("2015topoclpucvsSETd", "2015 topoclpuc MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015topoclpuce = new TCanvas("2015topoclpuce", "2015topoclpuce");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETe = new TH2F("2015topoclpucvsSETe", "2015 topoclpuc MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015topoclpucf = new TCanvas("2015topoclpucf", "2015topoclpucf");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETf = new TH2F("2015topoclpucvsSETf", "2015 topoclpuc MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015topoclpucg = new TCanvas("2015topoclpucg", "2015topoclpucg");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETg = new TH2F("2015topoclpucvsSETg", "2015 topoclpuc MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015topoclpuch = new TCanvas("2015topoclpuch", "2015topoclpuch");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETh = new TH2F("2015topoclpucvsSETh", "2015 topoclpuc MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015topoclpuci = new TCanvas("2015topoclpuci", "2015topoclpuci");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETi = new TH2F("2015topoclpucvsSETi", "2015 topoclpuc MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015topoclpucj = new TCanvas("2015topoclpucj", "2015topoclpucj");
	gPad->SetLogz();
	TH2F * 2015topoclpucvsSETj = new TH2F("2015topoclpucvsSETj", "2015 topoclpuc MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2015topoclpucvsSETj", "bcid>520.&&bcid<565.", "colz");


	/* Topocl PS Bunch Group Analysis*/
	TCanvas *2015topoclpsa = new TCanvas("2015topoclpsa", "2015topoclpsa");
	gPad->SetLogz();
	TH2F *2015topoclpsvsSETa = new TH2F("2015topoclpsvsSETa", "2015 topoclps MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015topoclpsb = new TCanvas("2015topoclpsb", "2015topoclpsb");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETb = new TH2F("2015topoclpsvsSETb", "2015 topoclps MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015topoclpsc = new TCanvas("2015topoclpsc", "2015topoclpsc");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETc = new TH2F("2015topoclpsvsSETc", "2015 topoclps MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015topoclpsd = new TCanvas("2015topoclpsd", "2015topoclpsd");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETd = new TH2F("2015topoclpsvsSETd", "2015 topoclps MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015topoclpse = new TCanvas("2015topoclpse", "2015topoclpse");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETe = new TH2F("2015topoclpsvsSETe", "2015 topoclps MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015topoclpsf = new TCanvas("2015topoclpsf", "2015topoclpsf");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETf = new TH2F("2015topoclpsvsSETf", "2015 topoclps MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015topoclpsg = new TCanvas("2015topoclpsg", "2015topoclpsg");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETg = new TH2F("2015topoclpsvsSETg", "2015 topoclps MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015topoclpsh = new TCanvas("2015topoclpsh", "2015topoclpsh");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETh = new TH2F("2015topoclpsvsSETh", "2015 topoclps MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015topoclpsi = new TCanvas("2015topoclpsi", "2015topoclpsi");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETi = new TH2F("2015topoclpsvsSETi", "2015 topoclps MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015topoclpsj = new TCanvas("2015topoclpsj", "2015topoclpsj");
	gPad->SetLogz();
	TH2F * 2015topoclpsvsSETj = new TH2F("2015topoclpsvsSETj", "2015 topoclps MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2015topoclpsvsSETj", "bcid>520.&&bcid<565.", "colz");
	
	/*CELL Bunch Group Analysis*/

	TCanvas * 2015cella = new TCanvas("2015cella", "2015cella");
	gPad->SetLogz();
	TH2F * 2015cellvsSETa = new TH2F("2015cellvsSETa", "2015 cell MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015cellb = new TCanvas("2015cellb", "2015cellb");
	gPad->SetLogz();
	TH2F * 2015cellvsSETb = new TH2F("2015cellvsSETb", "2015 cell MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015cellc = new TCanvas("2015cellc", "2015cellc");
	gPad->SetLogz();
	TH2F * 2015cellvsSETc = new TH2F("2015cellvsSETc", "2015 cell MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015celld = new TCanvas("2015celld", "2015celld");
	gPad->SetLogz();
	TH2F * 2015cellvsSETd = new TH2F("2015cellvsSETd", "2015 cell MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015celle = new TCanvas("2015celle", "2015celle");
	gPad->SetLogz();
	TH2F * 2015cellvsSETe = new TH2F("2015cellvsSETe", "2015 cell MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015cellf = new TCanvas("2015cellf", "2015cellf");
	gPad->SetLogz();
	TH2F * 2015cellvsSETf = new TH2F("2015cellvsSETf", "2015 cell MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015cellg = new TCanvas("2015cellg", "2015cellg");
	gPad->SetLogz();
	TH2F * 2015cellvsSETg = new TH2F("2015cellvsSETg", "2015 cell MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015cellh = new TCanvas("2015cellh", "2015cellh");
	gPad->SetLogz();
	TH2F * 2015cellvsSETh = new TH2F("2015cellvsSETh", "2015 cell MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015celli = new TCanvas("2015celli", "2015celli");
	gPad->SetLogz();
	TH2F * 2015cellvsSETi = new TH2F("2015cellvsSETi", "2015 cell MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015cellj = new TCanvas("2015cellj", "2015cellj");
	gPad->SetLogz();
	TH2F * 2015cellvsSETj = new TH2F("2015cellvsSETj", "2015 cell MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2015cellvsSETj", "bcid>520.&&bcid<565.", "colz");

	/*MHT Bunch Group Analysis*/

	TCanvas * 2015mhta = new TCanvas("2015mhta", "2015mhta");
	gPad->SetLogz();
	TH2F * 2015mhtvsSETa = new TH2F("2015mhtvsSETa", "2015 MHT MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015mhtvsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015MHTb = new TCanvas("2015MHTb", "2015MHTb");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETb = new TH2F("2015MHTvsSETb", "2015 MHT MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015MHTc = new TCanvas("2015MHTc", "2015MHTc");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETc = new TH2F("2015MHTvsSETc", "2015 MHT MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015MHTd = new TCanvas("2015MHTd", "2015MHTd");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETd = new TH2F("2015MHTvsSETd", "2015 MHT MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015MHTe = new TCanvas("2015MHTe", "2015MHTe");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETe = new TH2F("2015MHTvsSETe", "2015 MHT MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015MHTf = new TCanvas("2015MHTf", "2015MHTf");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETf = new TH2F("2015MHTvsSETf", "2015 MHT MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015MHTg = new TCanvas("2015MHTg", "2015MHTg");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETg = new TH2F("2015MHTvsSETg", "2015 MHT MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015MHTh = new TCanvas("2015MHTh", "2015MHTh");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETh = new TH2F("2015MHTvsSETh", "2015 MHT MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015MHTi = new TCanvas("2015MHTi", "2015MHTi");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETi = new TH2F("2015MHTvsSETi", "2015 MHT MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015MHTj = new TCanvas("2015MHTj", "2015MHTj");
	gPad->SetLogz();
	TH2F * 2015MHTvsSETj = new TH2F("2015MHTvsSETj", "2015 MHT MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2015MHTvsSETj", "bcid>520.&&bcid<565.", "colz");


	/*L1 Bunch Group Analysis*/

	TCanvas * 2015l1a = new TCanvas("2015l1a", "2015l1a");
	gPad->SetLogz();
	TH2F * 2015l1vsSETa = new TH2F("2015l1vsSETa", "2015 L1 MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015l1vsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015L1b = new TCanvas("2015L1b", "2015L1b");
	gPad->SetLogz();
	TH2F * 2015L1vsSETb = new TH2F("2015L1vsSETb", "2015 L1 MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015L1c = new TCanvas("2015L1c", "2015L1c");
	gPad->SetLogz();
	TH2F * 2015L1vsSETc = new TH2F("2015L1vsSETc", "2015 L1 MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015L1d = new TCanvas("2015L1d", "2015L1d");
	gPad->SetLogz();
	TH2F * 2015L1vsSETd = new TH2F("2015L1vsSETd", "2015 L1 MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015L1e = new TCanvas("2015L1e", "2015L1e");
	gPad->SetLogz();
	TH2F * 2015L1vsSETe = new TH2F("2015L1vsSETe", "2015 L1 MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015L1f = new TCanvas("2015L1f", "2015L1f");
	gPad->SetLogz();
	TH2F * 2015L1vsSETf = new TH2F("2015L1vsSETf", "2015 L1 MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015L1g = new TCanvas("2015L1g", "2015L1g");
	gPad->SetLogz();
	TH2F * 2015L1vsSETg = new TH2F("2015L1vsSETg", "2015 L1 MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015L1h = new TCanvas("2015L1h", "2015L1h");
	gPad->SetLogz();
	TH2F * 2015L1vsSETh = new TH2F("2015L1vsSETh", "2015 L1 MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015L1i = new TCanvas("2015L1i", "2015L1i");
	gPad->SetLogz();
	TH2F * 2015L1vsSETi = new TH2F("2015L1vsSETi", "2015 L1 MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015L1j = new TCanvas("2015L1j", "2015L1j");
	gPad->SetLogz();
	TH2F * 2015L1vsSETj = new TH2F("2015L1vsSETj", "2015 L1 MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2015L1vsSETj", "bcid>520.&&bcid<565.", "colz");

	/*Topocl Bunch Group Analysis*/
	TCanvas * 2015topocla = new TCanvas("2015topocla", "2015topocla");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETa = new TH2F("2015topoclvsSETa", "2015 TOPOCL MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETa", "bcid>60.&&bcid<104.", "colz");

	TCanvas * 2015topoclb = new TCanvas("2015topoclb", "2015topoclb");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETb = new TH2F("2015topoclvsSETb", "2015 topocl MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETb", "bcid>104.&&bcid<150.", "colz");

	TCanvas * 2015topoclc = new TCanvas("2015topoclc", "2015topoclc");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETc = new TH2F("2015topoclvsSETc", "2015 topocl MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETc", "bcid>150.&&bcid<200.", "colz");

	TCanvas * 2015topocld = new TCanvas("2015topocld", "2015topocld");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETd = new TH2F("2015topoclvsSETd", "2015 topocl MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETd", "bcid>200.&&bcid<266.", "colz");

	TCanvas * 2015topocle = new TCanvas("2015topocle", "2015topocle");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETe = new TH2F("2015topoclvsSETe", "2015 topocl MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETe", "bcid>266.&&bcid<316.", "colz");

	TCanvas * 2015topoclf = new TCanvas("2015topoclf", "2015topoclf");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETf = new TH2F("2015topoclvsSETf", "2015 topocl MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETf", "bcid>316.&&bcid<360.", "colz");

	TCanvas * 2015topoclg = new TCanvas("2015topoclg", "2015topoclg");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETg = new TH2F("2015topoclvsSETg", "2015 topocl MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETg", "bcid>360.&&bcid<410.", "colz");

	TCanvas * 2015topoclh = new TCanvas("2015topoclh", "2015topoclh");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETh = new TH2F("2015topoclvsSETh", "2015 topocl MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETh", "bcid>410.&&bcid<474.", "colz");

	TCanvas * 2015topocli = new TCanvas("2015topocli", "2015topocli");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETi = new TH2F("2015topoclvsSETi", "2015 topocl MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETi", "bcid>474.&&bcid<520.", "colz");

	TCanvas * 2015topoclj = new TCanvas("2015topoclj", "2015topoclj");
	gPad->SetLogz();
	TH2F * 2015topoclvsSETj = new TH2F("2015topoclvsSETj", "2015 topocl MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2015topoclvsSETj", "bcid>520.&&bcid<565.", "colz");


//2016 (13 Runs) Bunch Crossing IDs
	TFile *2016 = TFile::Open("ZeroBias2016.13Runs.root");

	TCanvas *16bcids = new TCanvas("16bcids", "2016 Bunch Crossing IDs");
	TH1F *2016bcid = new TH1F("2016bcid", "2016 Bunch Crossing IDs", 3500, 0., 3500.);
	tree->Draw("bcid>>2016bcid");

	
	
	/* Compare 2016 [algorithm] MET vs sqrt([algorithm] SET) where BCID>570 (bcid=1000 is the start of the 11th bunch group in 2016)*/
	TCanvas * 2016L1 = new TCanvas("2016L1", "2016L1");
	gPad->SetLogz();
	TH2F * 2016L1vsSET = new TH2F("2016L1vsSET", "2016 L1 MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSET", "bcid>1186.", "colz");

	TCanvas * 2016CELL = new TCanvas("2016CELL", "2016CELL");
	gPad->SetLogz();
	TH2F * 2016CELLvsSET = new TH2F("2016CELLvsSET", "2016 CELL MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016CELLvsSET", "bcid>1186.", "colz");

	TCanvas * 2016MHT = new TCanvas("2016MHT", "2016MHT");
	gPad->SetLogz();
	TH2F * 2016MHTvsSET = new TH2F("2016MHTvsSET", "2016 MHT MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016MHTvsSET", "bcid>1186.", "colz");

	TCanvas * 2016topoclpuc = new TCanvas("2016topoclpuc", "2016topoclpuc");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSET = new TH2F("2016topoclpucvsSET", "2016 topoclpuc MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSET", "bcid>1186.", "colz");

	TCanvas * 2016topoclps = new TCanvas("2016topoclps", "2016topoclps");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSET = new TH2F("2016topoclpsvsSET", "2016 topoclps MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSET", "bcid>1186.", "colz");

	TCanvas * 2016topocl = new TCanvas("2016topocl", "2016topocl");
	gPad->SetLogz();
	TH2F * 2016topoclvsSET = new TH2F("2016topoclvsSET", "2016 topocl MET vs sqrtSET for BCID > 1186", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSET", "bcid>1186.", "colz");


	/* Compare 2016 L1 MET vs sqrt(L1 SET) for all first 10 bunch groups*/
	TCanvas * 2016L1a = new TCanvas("2016L1a", "2016L1a");
	gPad->SetLogz();
	TH2F * 2016L1vsSETa = new TH2F("2016L1vsSETa", "2016 L1 MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016L1b = new TCanvas("2016L1b", "2016L1b");
	gPad->SetLogz();
	TH2F * 2016L1vsSETb = new TH2F("2016L1vsSETb", "2016 L1 MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016L1c = new TCanvas("2016L1c", "2016L1c");
	gPad->SetLogz();
	TH2F * 2016L1vsSETc = new TH2F("2016L1vsSETc", "2016 L1 MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016L1d = new TCanvas("2016L1d", "2016L1d");
	gPad->SetLogz();
	TH2F * 2016L1vsSETd = new TH2F("2016L1vsSETd", "2016 L1 MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016L1e = new TCanvas("2016L1e", "2016L1e");
	gPad->SetLogz();
	TH2F * 2016L1vsSETe = new TH2F("2016L1vsSETe", "2016 L1 MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016L1f = new TCanvas("2016L1f", "2016L1f");
	gPad->SetLogz();
	TH2F * 2016L1vsSETf = new TH2F("2016L1vsSETf", "2016 L1 MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016L1g = new TCanvas("2016L1g", "2016L1g");
	gPad->SetLogz();
	TH2F * 2016L1vsSETg = new TH2F("2016L1vsSETg", "2016 L1 MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016L1h = new TCanvas("2016L1h", "2016L1h");
	gPad->SetLogz();
	TH2F * 2016L1vsSETh = new TH2F("2016L1vsSETh", "2016 L1 MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETh", "bcid>890.&&bcid<1000.", "colz");

	TCanvas * 2016L1i = new TCanvas("2016L1i", "2016L1i");
	gPad->SetLogz();
	TH2F * 2016L1vsSETi = new TH2F("2016L1vsSETi", "2016 L1 MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETi", "bcid>1000.&&bcid<1100.", "colz");

	TCanvas * 2016L1j = new TCanvas("2016L1j", "2016L1j");
	gPad->SetLogz();
	TH2F * 2016L1vsSETj = new TH2F("2016L1vsSETj", "2016 L1 MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metl1:sqrt(setl1)>>2016L1vsSETj", "bcid>1100.&&bcid<1186.", "colz");

	/*CELL Bunch Group Analysis*/

	TCanvas * 2016cella = new TCanvas("2016cella", "2016cella");
	gPad->SetLogz();
	TH2F * 2016cellvsSETa = new TH2F("2016cellvsSETa", "2016 cell MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016cellb = new TCanvas("2016cellb", "2016cellb");
	gPad->SetLogz();
	TH2F * 2016cellvsSETb = new TH2F("2016cellvsSETb", "2016 cell MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016cellc = new TCanvas("2016cellc", "2016cellc");
	gPad->SetLogz();
	TH2F * 2016cellvsSETc = new TH2F("2016cellvsSETc", "2016 cell MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016celld = new TCanvas("2016celld", "2016celld");
	gPad->SetLogz();
	TH2F * 2016cellvsSETd = new TH2F("2016cellvsSETd", "2016 cell MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016celle = new TCanvas("2016celle", "2016celle");
	gPad->SetLogz();
	TH2F * 2016cellvsSETe = new TH2F("2016cellvsSETe", "2016 cell MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016cellf = new TCanvas("2016cellf", "2016cellf");
	gPad->SetLogz();
	TH2F * 2016cellvsSETf = new TH2F("2016cellvsSETf", "2016 cell MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016cellg = new TCanvas("2016cellg", "2016cellg");
	gPad->SetLogz();
	TH2F * 2016cellvsSETg = new TH2F("2016cellvsSETg", "2016 cell MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016cellh = new TCanvas("2016cellh", "2016cellh");
	gPad->SetLogz();
	TH2F * 2016cellvsSETh = new TH2F("2016cellvsSETh", "2016 cell MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETh", "bcid>890.&&bcid<1000." , "colz");

	TCanvas * 2016celli = new TCanvas("2016celli", "2016celli");
	gPad->SetLogz();
	TH2F * 2016cellvsSETi = new TH2F("2016cellvsSETi", "2016 cell MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETi", "bcid>1000.&&bcid<1100." , "colz");

	TCanvas * 2016cellj = new TCanvas("2016cellj", "2016cellj");
	gPad->SetLogz();
	TH2F * 2016cellvsSETj = new TH2F("2016cellvsSETj", "2016 cell MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metcell:sqrt(setcell)>>2016cellvsSETj", "bcid>1100.&&bcid<1186.", "colz");

	/*mht Bunch Group Analysis*/

	TCanvas * 2016mhta = new TCanvas("2016mhta", "2016mhta");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETa = new TH2F("2016mhtvsSETa", "2016 mht MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016mhtb = new TCanvas("2016mhtb", "2016mhtb");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETb = new TH2F("2016mhtvsSETb", "2016 mht MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016mhtc = new TCanvas("2016mhtc", "2016mhtc");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETc = new TH2F("2016mhtvsSETc", "2016 mht MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016mhtd = new TCanvas("2016mhtd", "2016mhtd");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETd = new TH2F("2016mhtvsSETd", "2016 mht MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016mhte = new TCanvas("2016mhte", "2016mhte");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETe = new TH2F("2016mhtvsSETe", "2016 mht MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016mhtf = new TCanvas("2016mhtf", "2016mhtf");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETf = new TH2F("2016mhtvsSETf", "2016 mht MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016mhtg = new TCanvas("2016mhtg", "2016mhtg");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETg = new TH2F("2016mhtvsSETg", "2016 mht MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016mhth = new TCanvas("2016mhth", "2016mhth");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETh = new TH2F("2016mhtvsSETh", "2016 mht MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETh", "bcid>890.&&bcid<1000.", "colz");

	TCanvas * 2016mhti = new TCanvas("2016mhti", "2016mhti");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETi = new TH2F("2016mhtvsSETi", "2016 mht MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETi", "bcid>1000.&&bcid<1100.", "colz");

	TCanvas * 2016mhtj = new TCanvas("2016mhtj", "2016mhtj");
	gPad->SetLogz();
	TH2F * 2016mhtvsSETj = new TH2F("2016mhtvsSETj", "2016 mht MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("metmht:sqrt(setmht)>>2016mhtvsSETj", "bcid>1100.&&bcid<1186.", "colz");
	
	/* Topocl PS Bunch Group Analysis*/
	TCanvas * 2016topoclpsa = new TCanvas("2016topoclpsa", "2016topoclpsa");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETa = new TH2F("2016topoclpsvsSETa", "2016 topoclps MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016topoclpsb = new TCanvas("2016topoclpsb", "2016topoclpsb");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETb = new TH2F("2016topoclpsvsSETb", "2016 topoclps MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016topoclpsc = new TCanvas("2016topoclpsc", "2016topoclpsc");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETc = new TH2F("2016topoclpsvsSETc", "2016 topoclps MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016topoclpsd = new TCanvas("2016topoclpsd", "2016topoclpsd");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETd = new TH2F("2016topoclpsvsSETd", "2016 topoclps MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016topoclpse = new TCanvas("2016topoclpse", "2016topoclpse");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETe = new TH2F("2016topoclpsvsSETe", "2016 topoclps MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016topoclpsf = new TCanvas("2016topoclpsf", "2016topoclpsf");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETf = new TH2F("2016topoclpsvsSETf", "2016 topoclps MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016topoclpsg = new TCanvas("2016topoclpsg", "2016topoclpsg");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETg = new TH2F("2016topoclpsvsSETg", "2016 topoclps MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016topoclpsh = new TCanvas("2016topoclpsh", "2016topoclpsh");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETh = new TH2F("2016topoclpsvsSETh", "2016 topoclps MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETh", "bcid>890.&&bcid<1000.", "colz");

	TCanvas * 2016topoclpsi = new TCanvas("2016topoclpsi", "2016topoclpsi");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETi = new TH2F("2016topoclpsvsSETi", "2016 topoclps MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETi", "bcid>1000.&&bcid<1100.", "colz");

	TCanvas * 2016topoclpsj = new TCanvas("2016topoclpsj", "2016topoclpsj");
	gPad->SetLogz();
	TH2F * 2016topoclpsvsSETj = new TH2F("2016topoclpsvsSETj", "2016 topoclps MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclps:sqrt(settopoclps)>>2016topoclpsvsSETj", "bcid>1100.&&bcid<1186.", "colz");

	/*topoclpuc Bunch Group Analysis*/

	TCanvas * 2016topoclpuca = new TCanvas("2016topoclpuca", "2016topoclpuca");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETa = new TH2F("2016topoclpucvsSETa", "2016 topoclpuc MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016topoclpucb = new TCanvas("2016topoclpucb", "2016topoclpucb");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETb = new TH2F("2016topoclpucvsSETb", "2016 topoclpuc MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016topoclpucc = new TCanvas("2016topoclpucc", "2016topoclpucc");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETc = new TH2F("2016topoclpucvsSETc", "2016 topoclpuc MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016topoclpucd = new TCanvas("2016topoclpucd", "2016topoclpucd");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETd = new TH2F("2016topoclpucvsSETd", "2016 topoclpuc MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016topoclpuce = new TCanvas("2016topoclpuce", "2016topoclpuce");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETe = new TH2F("2016topoclpucvsSETe", "2016 topoclpuc MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016topoclpucf = new TCanvas("2016topoclpucf", "2016topoclpucf");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETf = new TH2F("2016topoclpucvsSETf", "2016 topoclpuc MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016topoclpucg = new TCanvas("2016topoclpucg", "2016topoclpucg");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETg = new TH2F("2016topoclpucvsSETg", "2016 topoclpuc MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016topoclpuch = new TCanvas("2016topoclpuch", "2016topoclpuch");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETh = new TH2F("2016topoclpucvsSETh", "2016 topoclpuc MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETh", "bcid>890.&&bcid<1000.", "colz");

	TCanvas * 2016topoclpuci = new TCanvas("2016topoclpuci", "2016topoclpuci");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETi = new TH2F("2016topoclpucvsSETi", "2016 topoclpuc MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETi", "bcid>1000.&&bcid<1100.", "colz");

	TCanvas * 2016topoclpucj = new TCanvas("2016topoclpucj", "2016topoclpucj");
	gPad->SetLogz();
	TH2F * 2016topoclpucvsSETj = new TH2F("2016topoclpucvsSETj", "2016 topoclpuc MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>2016topoclpucvsSETj", "bcid>1100.&&bcid<1186.", "colz");


	/*topocl Bunch Group Analysis*/

	TCanvas * 2016topocla = new TCanvas("2016topocla", "2016topocla");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETa = new TH2F("2016topoclvsSETa", "2016 topocl MET vs sqrtSET BG = 1", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETa", "bcid>60.&&bcid<210.", "colz");

	TCanvas * 2016topoclb = new TCanvas("2016topoclb", "2016topoclb");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETb = new TH2F("2016topoclvsSETb", "2016 topocl MET vs sqrtSET BG = 2", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETb", "bcid>210.&&bcid<300.", "colz");

	TCanvas * 2016topoclc = new TCanvas("2016topoclc", "2016topoclc");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETc = new TH2F("2016topoclvsSETc", "2016 topocl MET vs sqrtSET BG = 3", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETc", "bcid>300.&&bcid<420.", "colz");

	TCanvas * 2016topocld = new TCanvas("2016topocld", "2016topocld");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETd = new TH2F("2016topoclvsSETd", "2016 topocl MET vs sqrtSET BG = 4", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETd", "bcid>420.&&bcid<520.", "colz");

	TCanvas * 2016topocle = new TCanvas("2016topocle", "2016topocle");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETe = new TH2F("2016topoclvsSETe", "2016 topocl MET vs sqrtSET BG = 5", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETe", "bcid>520.&&bcid<620.", "colz");

	TCanvas * 2016topoclf = new TCanvas("2016topoclf", "2016topoclf");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETf = new TH2F("2016topoclvsSETf", "2016 topocl MET vs sqrtSET BG = 6", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETf", "bcid>620.&&bcid<750.", "colz");

	TCanvas * 2016topoclg = new TCanvas("2016topoclg", "2016topoclg");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETg = new TH2F("2016topoclvsSETg", "2016 topocl MET vs sqrtSET BG = 7", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETg", "bcid>750.&&bcid<890.", "colz");

	TCanvas * 2016topoclh = new TCanvas("2016topoclh", "2016topoclh");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETh = new TH2F("2016topoclvsSETh", "2016 topocl MET vs sqrtSET BG = 8", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETh", "bcid>890.&&bcid<1000.", "colz");

	TCanvas * 2016topocli = new TCanvas("2016topocli", "2016topocli");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETi = new TH2F("2016topoclvsSETi", "2016 topocl MET vs sqrtSET BG = 9", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETi", "bcid>1000.&&bcid<1100.", "colz");

	TCanvas * 2016topoclj = new TCanvas("2016topoclj", "2016topoclj");
	gPad->SetLogz();
	TH2F * 2016topoclvsSETj = new TH2F("2016topoclvsSETj", "2016 topocl MET vs sqrtSET BG = 10", 100, 0., 100., 100, 0., 100.);
	tree->Draw("mettopocl:sqrt(settopocl)>>2016topoclvsSETj", "bcid>1100.&&bcid<1186.", "colz");





}