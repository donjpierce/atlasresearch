{
/*	//2015 Plot Cuts
	TString PlotCut("passrndm>0.1&&runnum<283200.");
	TString PlotCut1("passrndm>0.1&&runnum>283200.&&runnum<283800.");
	TString PlotCut2("passrndm>0.1&&runnum>283800.");
	//2016 Plot Cuts
	TString PlotCut3("passrndm>0.1&&runnum<299000.");
	TString PlotCut4("passrndm>0.1&&runnum>299000.&&runnum<301000.");
	TString PlotCut5("passrndm>0.1&&runnum>301000.");
*/

	//Comparing BCID across 2015 runs
	TFile *2015 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

	TCanvas *BCIDR1 = new TCanvas("BCIDR1", "BCIDR1");
	TH1F *2015bcid1 = new TH1F("2015bcid1", "2015 BCID Run=1", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid1", "passrndm>0.1&&runnum>282610&&runnum<282640.");

	TCanvas *BCIDR2 = new TCanvas("BCIDR2", "BCIDR2");
	TH1F *2015bcid2 = new TH1F("2015bcid2", "2015 BCID Run=2", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid2", "passrndm>0.1&&runnum>282695.&&runnum<282716.");

	TCanvas *BCIDR3 = new TCanvas("BCIDR3", "BCIDR3");
	TH1F *2015bcid3 = new TH1F("2015bcid3", "2015 BCID Run=3", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid3", "passrndm>0.1&&runnum>282750&&runnum<282800.");

	TCanvas *BCIDR4 = new TCanvas("BCIDR4", "BCIDR4");
	TH1F * 2015bcid4 = new TH1F("2015bcid4", "2015 BCID Run=4", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid4", "passrndm>0.1&&runnum>282980&&runnum<283010.");

	TCanvas *BCIDR5 = new TCanvas("BCIDR5", "BCIDR5");
	TH1F * 2015bcid5 = new TH1F("2015bcid5", "2015 BCID Run=5", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid5", "passrndm>0.1&&runnum>283050&&runnum<283090.");

	TCanvas *BCIDR6 = new TCanvas("BCIDR6", "BCIDR6");
	TH1F * 2015bcid6 = new TH1F("2015bcid6", "2015 BCID Run=6", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid6", "passrndm>0.1&&runnum>283140&&runnum<283170.");

	TCanvas *BCIDR7 = new TCanvas("BCIDR7", "BCIDR7");
	TH1F * 2015bcid7 = new TH1F("2015bcid7", "2015 BCID Run=7", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid7", "passrndm>0.1&&runnum>283240&&runnum<283300.");

	TCanvas *BCIDR8 = new TCanvas("BCIDR8", "BCIDR8");
	TH1F * 2015bcid8 = new TH1F("2015bcid8", "2015 BCID Run=8", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid8", "passrndm>0.1&&runnum>283410&&runnum<283440.");

	TCanvas *BCIDR9 = new TCanvas("BCIDR9", "BCIDR9");
	TH1F * 2015bcid9 = new TH1F("2015bcid9", "2015 BCID Run=9", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid9", "passrndm>0.1&&runnum>283550&&runnum<283650.");

	TCanvas *BCIDR10 = new TCanvas("BCIDR10", "BCIDR10");
	TH1F * 2015bcid10 = new TH1F("2015bcid10", "2015 BCID Run=10", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid10", "passrndm>0.1&&runnum>283760&&runnum<283800.");
	
	TCanvas *BCIDR11 = new TCanvas("BCIDR11", "BCIDR11");
	TH1F * 2015bcid11 = new TH1F("2015bcid11", "2015 BCID Run=11", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid11", "passrndm>0.1&&runnum>283980&&runnum<284020.");

	TCanvas *BCIDR12 = new TCanvas("BCIDR12", "BCIDR12");
	TH1F * 2015bcid12 = new TH1F("2015bcid12", "2015 BCID Run=12", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid12", "passrndm>0.1&&runnum>284150&&runnum<284170.");

	TCanvas *BCIDR13 = new TCanvas("BCIDR13", "BCIDR13");
	TH1F * 2015bcid13 = new TH1F("2015bcid13", "2015 BCID Run=13", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid13", "passrndm>0.1&&runnum>284200&&runnum<284220.");

	TCanvas *BCIDR14 = new TCanvas("BCIDR14", "BCIDR14");
	TH1F * 2015bcid14 = new TH1F("2015bcid14", "2015 BCID Run=14", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid14", "passrndm>0.1&&runnum>284280&&runnum<284300.");

	TCanvas *BCIDR15 = new TCanvas("BCIDR15", "BCIDR15");
	TH1F * 2015bcid15 = new TH1F("2015bcid15", "2015 BCID Run=15", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid15", "passrndm>0.1&&runnum>284400&&runnum<284450.");

	TCanvas *BCIDR16 = new TCanvas("BCIDR16", "BCIDR16");
	TH1F * 2015bcid16 = new TH1F("2015bcid16", "2015 BCID Run=16", 1000, 0., 3300.);
	tree->Draw("bcid>>2015bcid16", "passrndm>0.1&&runnum>284450&&runnum<284500.");

	//Comparing BCID Across 2015 Runs
	TFile * 2016 = TFile::Open("ZeroBias2016.13Runs.root");

	TCanvas *16BCIDR1 = new TCanvas("BCIDR1", "BCIDR1");
	TH1F * 2016bcid1 = new TH1F("2016bcid1", "2016 BCID Run=1", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid1", "passrndm>0.1&&runnum>297680&&runnum<297800.");

	TCanvas *16BCIDR2 = new TCanvas("BCIDR2", "BCIDR2");
	TH1F * 2016bcid2 = new TH1F("2016bcid2", "2016 BCID Run=2", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid2", "passrndm>0.1&&runnum>298500.&&runnum<298800.");

	TCanvas *16BCIDR3 = new TCanvas("BCIDR3", "BCIDR3");
	TH1F * 2016bcid3 = new TH1F("2016bcid3", "2016 BCID Run=3", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid3", "passrndm>0.1&&runnum>299000&&runnum<299200.");

	TCanvas *16BCIDR4 = new TCanvas("BCIDR4", "BCIDR4");
	TH1F * 2016bcid4 = new TH1F("2016bcid4", "2016 BCID Run=4", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid4", "passrndm>0.1&&runnum>299500&&runnum<299700.");

	TCanvas *16BCIDR5 = new TCanvas("BCIDR5", "BCIDR5");
	TH1F * 2016bcid5 = new TH1F("2016bcid5", "2016 BCID Run=5", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid5", "passrndm>0.1&&runnum>300800&&runnum<301000.");

	TCanvas *16BCIDR6 = new TCanvas("BCIDR6", "BCIDR6");
	TH1F * 2016bcid6 = new TH1F("2016bcid6", "2016 BCID Run=6", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid6", "passrndm>0.1&&runnum>302200&&runnum<302500.");

	TCanvas *16BCIDR7 = new TCanvas("BCIDR7", "BCIDR7");
	TH1F * 2016bcid7 = new TH1F("2016bcid7", "2016 BCID Run=7", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid7", "passrndm>0.1&&runnum>302800&&runnum<303100.");

	TCanvas *16BCIDR8 = new TCanvas("BCIDR8", "BCIDR8");
	TH1F * 2016bcid8 = new TH1F("2016bcid8", "2016 BCID Run=8", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid8", "passrndm>0.1&&runnum>303100&&runnum<303240.");

	TCanvas *16BCIDR9 = new TCanvas("BCIDR9", "BCIDR9");
	TH1F * 2016bcid9 = new TH1F("2016bcid9", "2016 BCID Run=9", 1000, 0., 3300.);
	tree->Draw("bcid>>2016bcid9", "passrndm>0.1&&runnum>303260&&runnum<303500.");

}