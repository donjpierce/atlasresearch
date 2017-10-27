//Int_t tailcalculator()
{
	#include <iostream>
	#include <fstream>
	#include <vector>
	#include "TF1.h"
	#include "TMath.h"
	#include "TH1.h"
	#include "TFile.h"
	#include "TTree.h"
	#include "TString.h"
	#include "TEfficiency.h"
	#include "TLegend.h"
	#include "TROOT.h"
	#include "TCanvas.h"
	#include "TSystem.h"
	#include "TH2F.h"
	#include "TPaveStats.h"
	#include "TStyle.h"

	//gROOT->ProcessLine("gROOT->SetBatch(kTRUE)"); // suppresses the drawing of graphs
	gROOT->ProcessLine("gROOT->Time();");

	// Define the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Defining a Linear Fit Function
	TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
	linfit->SetParameters(0, -80.);
	linfit->SetParameters(1, -80.);
	linfit->SetParLimits(0, -80., 80.);
	linfit->SetParLimits(1, -80., 80.);
	linfit->SetParName(0, "slope");
	linfit->SetParName(1, "intercept");

	// Defining a second-order fit Function
	TF1 *nfit = new TF1("nfit", "[0]*(x + [2])*(x + [2]) + [1]");
	nfit->SetParameters(0, -50.);
	nfit->SetParameters(1, -50.);
	nfit->SetParLimits(0, -50000., 50000.);
	nfit->SetParLimits(1, -50000., 50000.);
	nfit->SetParLimits(2, -100., 100.);
	nfit->SetParName(0, "amplitude");
	nfit->SetParName(1, "intercept");
	nfit->SetParName(2, "shift");

	//Fit Parameters
	//linear parametrs
	Double_t slope_ZeroBias[4];
	Double_t slope_Muon[4];
	Double_t intercept_ZeroBias[4];
	Double_t intercept_Muon[4];
	//second-order parametrs
	Double_t slope_ZeroBias[4];
	Double_t intercept_ZeroBias[4];
	Double_t shift_ZeroBias[4];

	TFile *zbFile = TFile::Open("../PhysicsMain.L1KFnoalgXEtriggers.2016.f731f758_m1659m1710.Run309759.48Runs-001.root");
	TTree *zbTree = (TTree*)zbFile->Get("tree");
	Int_t zbl1gt10, zbl1gt30, zbl1gt40, zbl1gt45;
	Float_t zbl1, zbint;
	zbTree->SetBranchAddress("passnoalgL1XE10", &zbl1gt10);
	zbTree->SetBranchAddress("passnoalgL1XE30", &zbl1gt30);
	zbTree->SetBranchAddress("passnoalgL1XE40", &zbl1gt40);
	zbTree->SetBranchAddress("passnoalgL1XE45", &zbl1gt45);
	zbTree->SetBranchAddress("metl1", &zbl1);
	zbTree->SetBranchAddress("actint", &zbint);

	TFile *muonFile = TFile::Open("../PhysicsMain.L1KFmuontriggers.2016.f731f758_m1659m1710.Run309759.48Runs-002.root");
	TTree* muonTree = (TTree*)muonFile->Get("tree");
	Int_t passmuonmed, passmuonvarmed, muonrecal, muonclean;
	Float_t ml1, muonint;
	muonTree->SetBranchAddress("passmu26med", &passmuonmed);
	muonTree->SetBranchAddress("passmu26varmed", &passmuonvarmed);
	muonTree->SetBranchAddress("metl1", &ml1);
	muonTree->SetBranchAddress("actint", &muonint);
	muonTree->SetBranchAddress("recalbroke", &muonrecal);
	muonTree->SetBranchAddress("passcleancuts", &muonclean);

	// choose with which file you're creating correlation plots
	// ZERO BIAS CORRELATIOON RUN SELECT
	TTree* runtree = (TTree*)zbFile->Get("tree");
	TString graphtitle = "2016 Prescaled (L1KFnoalgXEtriggers...48Runs-001) L1 > 50GeV";
	//TString runcut = (zbl1gt10 > 0.1 || zbl1gt30 > 0.1 || zbl1gt40 > 0.1 || zbl1gt45 > 0.1) && zbl1 > 50. ;
	// MUON CORRELATION RUN SELECT
	//TTree* runtree = (TTree*)muonFile->Get("tree");
	//TString graphtitle = "2016 Muons (L1KFmuontriggers...48Runs-002) for L1 > 50GeV, 40GeV < transversemass < 100GeV, and actint > 35.";
	//TString runcut = (passmuonmed > 0.1 || passmuonvarmed > 0.1) && ml1 > 50. && muonclean > 0.1 && muonrecal < 0.1;

	// initialize zerobias and muon cuts for resolution graphs
	TString zbPlotCut("(passnoalgL1XE10>0.5||passnoalgL1XE30>0.5||passnoalgL1XE40>0.5||passnoalgL1XE45>0.5)");
	TString muonsPlotCut("passmu26med>0.5||passmu26varmed>0.5");

	//Produce fitting graphs for zerobias events
	//TH2F *L1zb = new TH2F ("L1zb","", 60, 0., 60.,100,0.,100.);
		//zbTree->Draw("metl1:sqrt(setl1)>>L1zb","passrndm>0.5&&metl1>30.");
	TH2F *l1zb = new TH2F ("l1zb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("metl1:sqrt(setl1)>>l1zb", zbPlotCut);
	TH2F *l1kfzb = new TH2F("l1kfzb", "", 100, 0., 100., 1000, 0., 1000.);
		zbTree->Draw("metl1kf:sqrt(setl1kf)>>l1kfzb", zbPlotCut);
	TH2F *mhtl1zb = new TH2F ("mhtl1zb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("metl1mht:sqrt(setl1mht)>>mhtl1zb", zbPlotCut);
	TH2F *mhtl1kfzb = new TH2F("mhtl1kfzb", "", 100, 0., 100., 1000, 0., 1000.);
		zbTree->Draw("metl1mhtkf:sqrt(setl1mhtkf)>>mhtl1kfzb", zbPlotCut);

//L1 Algorithm resoltuions in ZeroBias
		TCanvas *cl1zb = new TCanvas("cl1zb", "L1 2016 ");
		l1zb->Draw();
		l1zb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *l1zb_1 = (TH1D*)gDirectory->Get("l1zb_1");
		l1zb_1->Draw();
		l1zb_1->Fit(linfit);
		slope_ZeroBias[0] = linfit->GetParameter(0);
		intercept_ZeroBias[0] = linfit->GetParameter(1);
		l1zb_1->SetTitle("Resolution of L1 in ZeroBias 2016 ");
		l1zb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		l1zb_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		l1zb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sl1zb = (TPaveStats*)l1zb_1->FindObject("stats");
		sl1zb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resl1zb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl1zb->AddEntry("l1zb_1", "Zero Bias Data", "L");
		resl1zb->Draw();

	//l1kf Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cl1kfzb = new TCanvas("cl1kfzb", "l1kf 2016 ");
		l1kfzb->Draw();
		l1kfzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *l1kfzb_1 = (TH1D*)gDirectory->Get("l1kfzb_1");
		l1kfzb_1->Draw();
		l1kfzb_1->Fit(linfit);
		slope_ZeroBias[1] = linfit->GetParameter(0);
		intercept_ZeroBias[1] = linfit->GetParameter(1);
		l1kfzb_1->SetTitle("Resolution of L1KF in ZeroBias 2016 ");
		l1kfzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		l1kfzb_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		l1kfzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sl1kfzb = (TPaveStats*)l1kfzb_1->FindObject("stats");
		sl1kfzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resl1kfzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl1kfzb->AddEntry("l1kfzb_1", "Zero Bias Data", "L");
		resl1kfzb->Draw();

		//mhtl1 Algorithm resoltuions in ZeroBias
		TCanvas *cmhtl1zb = new TCanvas("cmhtl1zb", "mhtl1 2016 ");
		mhtl1zb->Draw();
		mhtl1zb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *mhtl1zb_1 = (TH1D*)gDirectory->Get("mhtl1zb_1");
		mhtl1zb_1->Draw();
		mhtl1zb_1->Fit(linfit);
		slope_ZeroBias[2] = linfit->GetParameter(0);
		intercept_ZeroBias[2] = linfit->GetParameter(1);
		mhtl1zb_1->SetTitle("Resolution of MHT L1 in ZeroBias 2016 ");
		mhtl1zb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		mhtl1zb_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		mhtl1zb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *smhtl1zb = (TPaveStats*)mhtl1zb_1->FindObject("stats");
		smhtl1zb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resmhtl1zb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resmhtl1zb->AddEntry("mhtl1zb_1", "Zero Bias Data", "L");
		resmhtl1zb->Draw();

		//xMHT Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cmhtl1kfzb = new TCanvas("cmhtl1kfzb", "xMHT 2016 ");
		mhtl1kfzb->Draw();
		mhtl1kfzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *mhtl1kfzb_1 = (TH1D*)gDirectory->Get("mhtl1kfzb_1");
		mhtl1kfzb_1->Draw("");
		mhtl1kfzb_1->Fit(linfit);
		slope_ZeroBias[3] = linfit->GetParameter(0);
		intercept_ZeroBias[3] = linfit->GetParameter(1);
		mhtl1kfzb_1->SetTitle("Resolution of MHTKF L1 in ZeroBias 2016 ");
		mhtl1kfzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		mhtl1kfzb_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		mhtl1kfzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *smhtl1kfzb = (TPaveStats*)mhtl1kfzb_1->FindObject("stats");
		smhtl1kfzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resmhtl1kfzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resmhtl1kfzb->AddEntry("mhtl1kfzb_1", "Zero Bias Data", "L");
		resmhtl1kfzb->Draw();


		//___Calculate Tail Events Based on Resolutions___
		TString metalgName[4] = {"metl1", "metl1kf", "metl1mht", "metl1mhtkf"};
		TString setalgName[4] = {"setl1", "setl1kf", "setl1mht", "setl1mhtkf"};

		// create arrays for MET and SET branches
		Float_t met[4]; Float_t set[4];
		for (Int_t i = 0; i < 4; i++)
		{
			runtree->SetBranchAddress(metalgName[i], &met[i]);
			runtree->SetBranchAddress(setalgName[i], &set[i]);
		}

		// initialize variables for calculating transverse mass
		Double_t transversemass;
		Float_t metoff, metoffw, mexoff, mexoffw, meyoff, meyoffw;
		muonTree->SetBranchAddress("metoffrecal", &metoff);
		muonTree->SetBranchAddress("metoffrecalmuon", &metoffw);
		muonTree->SetBranchAddress("mexoffrecal", &mexoff);
		muonTree->SetBranchAddress("mexoffrecalmuon", &mexoffw);
		muonTree->SetBranchAddress("meyoffrecal", &meyoff);
		muonTree->SetBranchAddress("meyoffrecalmuon", &meyoffw);

		// create graphs which I will later populate with TailMET vs. MET of different algorithm pairs
		// correlationgraphs will be populated with the FULL dataset
		TH2F *correlationgraph[12];
		char *histname = new char[12];
		Int_t bins = 1000;
		Double_t min = 0.;
		Double_t max = 1000.;
		for (Int_t i = 0; i < 12; i++)
		{
			sprintf(histname, "histo%d", i+1);
			correlationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create oddcorrelationgraphs to be populated with the odd-numbered entries in the dataset
		TH2F *oddcorrelationgraph[12];
		for (int i = 0; i < 12; i++)
		{
			sprintf(histname, "oddhisto%d", i+1);
			oddcorrelationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create evencorrelationgraphs to be populated with the even-numbered entries in the dataset
		TH2F *evencorrelationgraph[12];
		for (int i = 0; i < 12; i++)
		{
			sprintf(histname, "evenhisto%d", i+1);
			evencorrelationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create a list whose entries correspond to algorithms and are the number of events in that algorithm's tail
		Double_t tailagreement[12];
		for (i = 0; i < 12; i++)
		{
			tailagreement[i] = 0;
		}

	int n = 0; // this variable will determine whether an event is even-numbered or odd-numbered-
	Long64_t nentries = runtree->GetEntries();
		for (Int_t i = 0; i < nentries; i++)
		{
			n = ( 1 - n ); // this logic changes n to be either 0 or 1
			runtree->GetEntry(i);
			if (i % 100000 == 0)
			{
				cout << "hey there good lookin'";
			}

		// tranvserse mass based on metoffrecal
		//transversemass = sqrt(2*metoff*metoffw*(1+((mexoff*mexoffw+meyoff*meyoffw) / (metoff*metoffw))));

		if ( /*(passmuonmed > 0.1 || passmuonvarmed > 0.1) && ml1 > 50. && muonint > 35. && muonclean > 0.1 && muonrecal < 0.1 && 40. < transversemass && transversemass < 100. *//*(zbl1gt10 > 0.1 || zbl1gt30 > 0.1 || zbl1gt40 > 0.1 || zbl1gt45 > 0.1) && zbl1 > 50.*/ )
		{
			Double_t sigma[4];
			Double_t metdist[4]; // metdist will be the distance of the event's MET from the median
			Double_t x[4]; // x = bulkmet and y = tailmet will be calculated for each algorithm
			Double_t y[4];

			// the following loop populates the sigma and metdist arrays
			for (Int_t j = 0; j < 4; j++)
			{
				if (sqrt(set[j]) >= 10.) // throw out events whose SET values are too low
				{
					// compute sigma of this event for all algorithms
					sigma[j] = slope_ZeroBias[j]*sqrt(set[j]) + intercept_ZeroBias[j];
					//compute metdist for all algorithms
					metdist[j] = TMath::Abs( met[j] - (sigma[j]*TMath::Sqrt(TMath::PiOver2())));
				}
			}

			// the following logic populates correlationgraphs with (x = met, y = tailmet) tuples
			// only if they exist for a given event in the tree
			Int_t h = 0; // this variable counts each correlationgraph
			for (Int_t l = 0; l < 4; l++)
			{
				if (metdist[l] >= 3*sigma[l]) // if the event is in the tail of alg A
				{
					y[l] = met[l]; // save to y = tail met
					for (Int_t m = 0; m < 4; m++)
					{
						if (l == m)	continue;
							if (metdist[m] >= 3*sigma[m]) // if the event is also in the tail of Alg B
							{
								tailagreement[h]++;
							}
								x[m] = met[m]; // save to x = met
								correlationgraph[h]->Fill(x[m], y[l]); // and populate the appropraite correlationgraph
								if (n == 0)
								{
									oddcorrelationgraph[h]->Fill(x[m], y[l]); // populate with odd-numbered entry
								}
								if (n == 1)
								{
									evencorrelationgraph[h]->Fill(x[m], y[l]); // populate with even-numbered entry
								}
							h++;
					}
				}
			}
		}
	}

//======================================================================================================================================//

		TString xaxisNames[4] = {"L1 MET [GeV]", "L1KF MET [GeV]", "MHT L1 MET [GeV]", "MHTKF L1 MET [GeV]"};
		TString yaxisNames[4] = {"L1 Tail MET [GeV]", "L1KF Tail MET [GeV]", "MHT L1 Tail MET [GeV]", "MHTKF L1 Tail MET [GeV]"};

		ofstream correlationcoefficients; // prepare log file of correlation coefficients
		correlationcoefficients.open("correlationvalues.txt"); // open log file
		correlationcoefficients << "FILE:" << " " << graphtitle << "\n\n";
		correlationcoefficients << "Graph" << "\t" << "Correlation" << " " << "±" << " " << "Approx. Uncertainty" << " " << "\t" << "Correlation Graph" << "\t\t\t\t\t" << "Tail Fractions" "\n"; // write title of table

		Double_t tailfractions[12];
		Double_t correlationentries[12];
		TCanvas *mycanv[12];
		char *canvname = new char[12];
		Double_t r[12]; // correlation coefficients
		Double_t oddvalue[12]; // correlation values from oddcorrelationgraphs
		Double_t evenvalue[12]; // correlation values from evencorrelationgraphs
		Double_t c[12]; // final confidence interval of original correlation coefficients (r values)
		int k = 0; // this variable counts correlationgraphs
		for (int q = 0; q < 4; q++)
		{
			for (int l = 0; l < 4; l++)
			{
				if (q == l)	continue;
					canvname = Form("canv%d",k+1);
					mycanv[k] = new TCanvas(canvname, "");
					correlationgraph[k]->Draw("colz"); // add "colz" in function if desired
					correlationgraph[k]->GetYaxis()->SetTitle(yaxisNames[q]);
					correlationgraph[k]->GetXaxis()->SetTitle(xaxisNames[l]);
					correlationgraph[k]->SetTitle(graphtitle);
					mycanv[k]->SetLogz();
					correlationentries[k] = correlationgraph[k]->GetEntries();
					tailfractions[k] = tailagreement[k]/correlationentries[k];
					r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of each graph
					oddvalue[k] = oddcorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of odd graphs
					evenvalue[k] = evencorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of even graphs
					c[k] = 0.5*(oddvalue[k] - evenvalue[k]);
					mycanv[k]->Print(Form("%d.png", k+1));
					correlationcoefficients << k+1 << "\t" << r[k] << " " << "±" << " " << c[k]	<< ',' << "\t\t" << yaxisNames[q] << " " << "vs." << " " << xaxisNames[l] << "\t\t" << tailfractions[k] << "\n";
					k++;
			}
		}

		correlationcoefficients.close();

		return 0;
}
