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
	Double_t slope_ZeroBias[6];
	Double_t slope_Muon[6];
	Double_t intercept_ZeroBias[6];
	Double_t intercept_Muon[6];
	//second-order parametrs
	Double_t slope_ZeroBias[6];
	Double_t intercept_ZeroBias[6];
	Double_t shift_ZeroBias[6];

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
	TString graphtitle = "2016 Prescaled (L1KFnoalgXEtriggers...48Runs-001) for L1 > 50GeV";
	//bool runcut = (zbl1gt10 > 0.1 || zbl1gt30 > 0.1 || zbl1gt40 > 0.1 || zbl1gt45 > 0.1) && zbl1 > 50. ;
	// MUON CORRELATION RUN SELECT
	//TTree* runtree = (TTree*)muonFile->Get("tree");
	//TString graphtitle = "2016 Muons (L1KFmuontriggers...48Runs-002) for 40 < REF mass < 100, L1 > 50GeV";
	//bool runcut = (passmuonmed > 0.1 || passmuonvarmed > 0.1) && ml1 > 50. && muonint > 35. && 20. < transversemass && transversemass < 100.;

	// initialize zerobias and muon cuts for resolution graphs
	TString zbPlotCut("(passnoalgL1XE10>0.5||passnoalgL1XE30>0.5||passnoalgL1XE40>0.5||passnoalgL1XE45>0.5)");
	TString muonsPlotCut("passmu26med>0.5||passmu26varmed>0.5");

	//Produce fitting graphs for zerobias events
	//TH2F *L1zb = new TH2F ("L1zb","", 60, 0., 60.,100,0.,100.);
		//zbTree->Draw("metl1:sqrt(setl1)>>L1zb","passrndm>0.5&&metl1>30.");
	TH2F *CELLzb = new TH2F ("CELLzb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("metcell:sqrt(setcell)>>CELLzb", zbPlotCut);
	TH2F *MHTzb = new TH2F("MHTzb", "", 100, 0., 100., 1000, 0., 1000.);
		zbTree->Draw("metmht:sqrt(setmht)>>MHTzb", zbPlotCut);
	TH2F *TOPOCLzb = new TH2F ("TOPOCLzb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLzb", zbPlotCut);
	TH2F *TopoclEMzb = new TH2F("TopoclEMzb", "", 100, 0., 100., 1000, 0., 1000.);
		zbTree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMzb", zbPlotCut);
	TH2F *TOPOCLPSzb = new TH2F ("TOPOCLPSzb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSzb", zbPlotCut);
	TH2F *TOPOCLPUCzb = new TH2F ("TOPOCLPUCzb","", 100, 0., 100.,1000,0.,1000.);
		zbTree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCzb", zbPlotCut);

/*
	//Produce fitting graphs for muon events
	TH2F *L1muon = new TH2F ("L1muon","", 100, 0., 100.,1000,0.,1000.);
		muonTree->Draw("metl1:sqrt(setl1)>>L1muon",muonsPlotCut);
	TH2F *CELLmuon = new TH2F ("CELLmuon","",100,0.,100.,1000,0.,1000.);
		muonTree->Draw("metcell:sqrt(setcell)>>CELLmuon",muonsPlotCut);
	TH2F *MHTmuon = new TH2F("MHTmuon", "", 100, 0., 100., 1000, 0., 1000.);
		muonTree->Draw("metmht:sqrt(setmht)>>MHTmuon", muonsPlotCut);
	TH2F *TOPOCLmuon = new TH2F ("TOPOCLmuon","",100,0.,100.,1000,0.,1000.);
		muonTree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLmuon",muonsPlotCut);
	TH2F *TopoclEMmuon = new TH2F("TopoclEMmuon", "", 100, 0., 100., 1000, 0., 1000.);
		muonTree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMmuon", muonsPlotCut);
	TH2F *TOPOCLPSmuon = new TH2F ("TOPOCLPSmuon","", 100, 0., 100., 1000, 0., 1000.);
		muonTree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSmuon", muonsPlotCut);
	TH2F *TOPOCLPUCmuon = new TH2F ("TOPOCLPUCmuon","", 100, 0., 100., 1000, 0., 1000.);
		muonTree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCmuon", muonsPlotCut);
*/

/*
//L1 Algorithm resolutions in ZeroBias and Muons
	TCanvas *cL115 = new TCanvas("cL115", "L1 2015");
	L115->Draw();
	L115->FitSlicesY(func, 0, -1, 10, "L");
	L115_1->Draw();
	L115_1->Fit("nfit");
	slope_ZeroBias[0] = nfit->GetParameter(0);
	intercept_ZeroBias[0] = nfit->GetParameter(1);
	shift_ZeroBias[0] = nfit->GetParameter(2);
	L115_1->SetTitle("Resolution of L1 in ZeroBias 2016");
	L115_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L115_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L115_1->SetLineColor(2);
	gPad->Update();
		TPaveStats *l115 = (TPaveStats*)L115_1 ->FindObject("stats");
		l115->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resl115 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl115->AddEntry("L115_1", "Zero", "L");
		resl115->Draw();

	TCanvas *cL116 = new TCanvas("cL116", "L1 2016 ");
	L116->Draw();
	L116->FitSlicesY(func, 0, -1, 10, "L");
	L116_1->Draw();
	L116_1->Fit("nfit");
	slope_Muon[0] = nfit->GetParameter(0);
	intercept_Muon[0] = nfit->GetParameter(1);
	shift_Muon[0] = nfit->GetParameter(2);
	L116_1->SetTitle("Resolution of L1 in Muons (L1XE45..Runs9B) 2016 ");
	L116_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L116_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L116_1->SetLineColor(4);
	gPad->Update();
		TPaveStats *l116 = (TPaveStats*)L116_1 ->FindObject("stats");
		l116->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resl116 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl116->AddEntry("L116_1", "Muon Data", "L");
		resl116->Draw();
*/


//CELL Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cCELLzb = new TCanvas("cCELLzb", "CELL 2015 ");
		CELLzb->Draw();
		CELLzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *CELLzb_1 = (TH1D*)gDirectory->Get("CELLzb_1");
		CELLzb_1->Draw();
		CELLzb_1->Fit(linfit);
		slope_ZeroBias[0] = linfit->GetParameter(0);
		intercept_ZeroBias[0] = linfit->GetParameter(1);
		CELLzb_1->SetTitle("Resolution of CELL in ZeroBias 2016 ");
		CELLzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELLzb_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELLzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sCELLzb = (TPaveStats*)CELLzb_1->FindObject("stats");
		sCELLzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resCELLzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELLzb->AddEntry("CELLzb_1", "Zero Bias Data", "L");
		resCELLzb->Draw();

/*
		TCanvas *cCELLmuon = new TCanvas("cCELLmuon", "CELL 2016 ");
		CELLmuon->Draw();
		CELLmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *CELLmuon_1 = (TH1D*)gDirectory->Get("CELLmuon_1");
		CELLmuon_1->Draw();
		CELLmuon_1->Fit(linfit);
		slope_Muon[0] = linfit->GetParameter(0);
		intercept_Muon[0] = linfit->GetParameter(1);
		CELLmuon_1->SetTitle("Resolution of CELL in Muons (L1XE45..Runs9B) 2016 ");
		CELLmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELLmuon_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELLmuon_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sCELLmuon = (TPaveStats*)CELLmuon_1->FindObject("stats");
		sCELLmuon->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resCELLmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELLmuon->AddEntry("CELLmuon_1", "Muon Data", "L");
		resCELLmuon->Draw();
*/

		//MHT Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cMHTzb = new TCanvas("cMHTzb", "MHT 2015 ");
		MHTzb->Draw();
		MHTzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *MHTzb_1 = (TH1D*)gDirectory->Get("MHTzb_1");
		MHTzb_1->Draw();
		MHTzb_1->Fit("linfit");
		slope_ZeroBias[1] = linfit->GetParameter(0);
		intercept_ZeroBias[1] = linfit->GetParameter(1);
		MHTzb_1->SetTitle("Resolution of MHT in ZeroBias 2016 ");
		MHTzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHTzb_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHTzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sMHTzb = (TPaveStats*)MHTzb_1->FindObject("stats");
		sMHTzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resMHTzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHTzb->AddEntry("MHTzb_1", "Zero Bias Data", "L");
		resMHTzb->Draw();
/*
		TCanvas *cMHTmuon = new TCanvas("cMHTmuon", "MHT 2016 ");
		MHTmuon->Draw();
		MHTmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *MHTmuon_1 = (TH1D*)gDirectory->Get("MHTmuon_1");
		MHTmuon_1->Draw();
		MHTmuon_1->Fit("linfit");
		slope_Muon[1] = linfit->GetParameter(0);
		intercept_Muon[1] = linfit->GetParameter(1);
		MHTmuon_1->SetTitle("Resolution of MHT in Muons (L1XE45..Runs9B) 2016 ");
		MHTmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHTmuon_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHTmuon_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sMHTmuon = (TPaveStats*)MHTmuon_1->FindObject("stats");
		sMHTmuon->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resMHTmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHTmuon->AddEntry("MHTmuon_1", "Muon Data", "L");
		resMHTmuon->Draw();
*/

		//TOPOCL Algorithm resoltuions in ZeroBias and Muon
		TCanvas *cTOPOCLzb = new TCanvas("cTOPOCLzb", "TOPOCL 2015 ");
		TOPOCLzb->Draw();
		TOPOCLzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLzb_1 = (TH1D*)gDirectory->Get("TOPOCLzb_1");
		TOPOCLzb_1->Draw();
		TOPOCLzb_1->Fit(linfit);
		slope_ZeroBias[2] = linfit->GetParameter(0);
		intercept_ZeroBias[2] = linfit->GetParameter(1);
		TOPOCLzb_1->SetTitle("Resolution of TOPOCL in ZeroBias 2016 ");
		TOPOCLzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLzb_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCLzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLzb = (TPaveStats*)TOPOCLzb_1->FindObject("stats");
		sTOPOCLzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLzb->AddEntry("TOPOCLzb_1", "Zero Bias Data", "L");
		resTOPOCLzb->Draw();

/*
		TCanvas *cTOPOCLmuon = new TCanvas("cTOPOCLmuon", "TOPOCL 2016 ");
		TOPOCLmuon->Draw();
		TOPOCLmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLmuon_1 = (TH1D*)gDirectory->Get("TOPOCLmuon_1");
		TOPOCLmuon_1->Draw();
		TOPOCLmuon_1->Fit(linfit);
		slope_Muon[2] = linfit->GetParameter(0);
		intercept_Muon[2] = linfit->GetParameter(1);
		TOPOCLmuon_1->SetTitle("Resolution of TOPOCL in Muons (L1XE45..Runs9B) 2016 ");
		TOPOCLmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLmuon_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCLmuon_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLmuon = (TPaveStats*)TOPOCLmuon_1->FindObject("stats");
		sTOPOCLmuon->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLmuon->AddEntry("TOPOCLmuon_1", "Muon Data", "L");
		resTOPOCLmuon->Draw();
		*/

		//TOPOCLPS Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cTOPOCLPSzb = new TCanvas("cTOPOCLPSzb", "TOPOCLPS 2015 ");
		TOPOCLPSzb->Draw();
		TOPOCLPSzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPSzb_1 = (TH1D*)gDirectory->Get("TOPOCLPSzb_1");
		TOPOCLPSzb_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPSzb_1->Draw();
		TOPOCLPSzb_1->Fit(linfit);
		slope_ZeroBias[3] = linfit->GetParameter(0);
		intercept_ZeroBias[3] = linfit->GetParameter(1);
		TOPOCLPSzb_1->SetTitle("Resolution of TOPOCLPS in ZeroBias 2016 ");
		TOPOCLPSzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPSzb_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPSzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPSzb = (TPaveStats*)TOPOCLPSzb_1->FindObject("stats");
		sTOPOCLPSzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPSzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPSzb->AddEntry("TOPOCLPSzb_1", "Zero Bias Data", "L");
		resTOPOCLPSzb->Draw();

/*
		TCanvas *cTOPOCLPSmuon = new TCanvas("cTOPOCLPSmuon", "TOPOCLPS 2016 ");
		TOPOCLPSmuon->Draw();
		TOPOCLPSmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPSmuon_1 = (TH1D*)gDirectory->Get("TOPOCLPSmuon_1");
		TOPOCLPSmuon_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPSmuon_1->Draw();
		TOPOCLPSmuon_1->Fit(linfit);
		slope_Muon[3] = linfit->GetParameter(0);
		intercept_Muon[3] = linfit->GetParameter(1);
		TOPOCLPSmuon_1->SetTitle("Resolution of TOPOCLPS in Muons (L1XE45..Runs9B) 2016 ");
		TOPOCLPSmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPSmuon_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPSmuon_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPSmuon = (TPaveStats*)TOPOCLPSmuon_1->FindObject("stats");
		sTOPOCLPSmuon->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPSmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPSmuon->AddEntry("TOPOCLPSmuon_1", "Muon Data", "L");
		resTOPOCLPSmuon->Draw();
*/

		//TOPOCLPUC Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cTOPOCLPUCzb = new TCanvas("cTOPOCLPUCzb", "TOPOCLPUC 2015 ");
		TOPOCLPUCzb->Draw();
		TOPOCLPUCzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPUCzb_1 = (TH1D*)gDirectory->Get("TOPOCLPUCzb_1");
		TOPOCLPUCzb_1->Draw("");
		TOPOCLPUCzb_1->Fit(linfit);
		slope_ZeroBias[4] = linfit->GetParameter(0);
		intercept_ZeroBias[4] = linfit->GetParameter(1);
		TOPOCLPUCzb_1->SetTitle("Resolution of TOPOCLPUC in ZeroBias 2016 ");
		TOPOCLPUCzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUCzb_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUCzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPUCzb = (TPaveStats*)TOPOCLPUCzb_1->FindObject("stats");
		sTOPOCLPUCzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUCzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUCzb->AddEntry("TOPOCLPUCzb_1", "Zero Bias Data", "L");
		resTOPOCLPUCzb->Draw();

/*
		TCanvas *cTOPOCLPUCmuon = new TCanvas("cTOPOCLPUCmuon", "TOPOCLPUC 2016 ");
		TOPOCLPUCmuon->Draw();
		TOPOCLPUCmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPUCmuon_1 = (TH1D*)gDirectory->Get("TOPOCLPUCmuon_1");
		TOPOCLPUCmuon_1->Draw();
		TOPOCLPUCmuon_1->Fit(linfit);
		slope_Muon[4] = linfit->GetParameter(0);
		intercept_Muon[4] = linfit->GetParameter(1);
		TOPOCLPUCmuon_1->SetTitle("Resolution of TOPOCLPUC in Muons (L1XE45..Runs9B) 2016 ");
		TOPOCLPUCmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUCmuon_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUCmuon_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPUCmuon = (TPaveStats*)TOPOCLPUCmuon_1->FindObject("stats");
		sTOPOCLPUCmuon->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUCmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUCmuon->AddEntry("TOPOCLPUCmuon_1", "Muon Data", "L");
		resTOPOCLPUCmuon->Draw();
*/

		//TopoclEM Algorithm resolutions in ZeroBias2016 and Muons2016
		TCanvas *cTopoclEMzb = new TCanvas("cTopoclEMzb", "TopoclEM ZeroBias2016");
		TopoclEMzb->Draw();
		TopoclEMzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TopoclEMzb_1 = (TH1D*)gDirectory->Get("TopoclEMzb_1");
		TopoclEMzb_1->Draw();
		TopoclEMzb_1->Fit(linfit);
		slope_ZeroBias[5] = linfit->GetParameter(0);
		intercept_ZeroBias[5] = linfit->GetParameter(1);
		TopoclEMzb_1->SetTitle("Resolution of TopoclEM in ZeroBias2016");
		TopoclEMzb_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TopoclEMzb_1->GetYaxis()->SetTitle("#sigma of Fit for TopoclEM [GeV]");
		TopoclEMzb_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *topoclemzb = (TPaveStats*)TopoclEMzb_1 ->FindObject("stats");
		topoclemzb->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* restopoclemzb = new TLegend(0.37, 0.7, 0.55, 0.88);
		restopoclemzb->AddEntry("TopoclEMzb_1", "Zero Bias Data", "L");
		restopoclemzb->Draw();

/*
		TCanvas *cTopoclEMmuon = new TCanvas("cTopoclEMmuon", "TopoclEMMuon 2016 ");
		TopoclEMmuon->Draw();
		TopoclEMmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TopoclEMmuon_1 = (TH1D*)gDirectory->Get("TopoclEMmuon_1");
		TopoclEMmuon_1->Draw();
		TopoclEMmuon_1->Fit(linfit);
		slope_Muon[5] = linfit->GetParameter(0);
		intercept_Muon[5] = linfit->GetParameter(1);
		TopoclEMmuon_1->SetTitle("Resolution of TopoclEM in Muons2016 ");
		TopoclEMmuon_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TopoclEMmuon_1->GetYaxis()->SetTitle("#sigma of Fit for TopoclEM [GeV]");
		TopoclEMmuon_1->SetLineColor(4);
		gPad->Update();
			TPaveStats *topoclemmuon = (TPaveStats*)TopoclEMmuon_1 ->FindObject("stats");
			topoclemmuon->SetTextColor(4);
			gStyle->SetOptFit(11);

			TLegend* restopoclemmuon = new TLegend(0.37, 0.7, 0.55, 0.88);
			restopoclemmuon->AddEntry("TopoclEMmuon_1", "Muon Data", "L");
			restopoclemmuon->Draw();
*/

		//___Calculate Tail Events Based on Resolutions___
		TString metalgName[6] = {"metcell", "metmht", "mettopocl", "mettopoclps", "mettopoclpuc", "mettopoclem"};
		TString setalgName[6] = {"setcell", "setmht", "settopocl", "settopoclps", "settopoclpuc", "settopoclem"};

		// create arrays for MET and SET branches
		Float_t met[6]; Float_t set[6];
		for (Int_t i = 0; i < 6; i++)
		{
			runtree->SetBranchAddress(metalgName[i], &met[i]);
			runtree->SetBranchAddress(setalgName[i], &set[i]);
		}

			// initialize variables for calculating transverse mass
			Float_t transversemass;

			//Float_t metref, metrefw, mexref, mexrefw, meyref, meyrefw;
			//muonTree->SetBranchAddress("metrefmuon", &metref);
			//muonTree->SetBranchAddress("metrefwmuon", &metrefw);
			//muonTree->SetBranchAddress("mexrefmuon", &mexref);
			//muonTree->SetBranchAddress("mexrefwmuon", &mexrefw);
			//muonTree->SetBranchAddress("meyrefmuon", &meyref);
			//muonTree->SetBranchAddress("meyrefwmuon", &meyrefw);

			Float_t metoff, metoffmuon, mexoff, mexoffmuon, meyoff, meyoffmuon;
			muonTree->SetBranchAddress("metoffrecal", &metoff);
			muonTree->SetBranchAddress("metoffrecalmuon", &metoffmuon);
			muonTree->SetBranchAddress("mexoffrecal", &mexoff);
			muonTree->SetBranchAddress("mexoffrecalmuon", &mexoffmuon);
			muonTree->SetBranchAddress("meyoffrecal", &meyoff);
			muonTree->SetBranchAddress("meyoffrecalmuon", &meyoffmuon);

		// create histograms which will later be populated with TailMET vs. MET of different algorithm pairs
		TH2F *correlationgraph[30];
		char *histname = new char[30];
		Int_t bins = 1000;
		Double_t min = 0.;
		Double_t max = 1000.;
		for (Int_t i = 0; i < 30; i++)
		{
			sprintf(histname, "histo%d", i+1);
			correlationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create oddcorrelationgraphs to be populated with the odd-numbered entries in the dataset
		TH2F *oddcorrelationgraph[30];
		for (int i = 0; i < 30; i++)
		{
			sprintf(histname, "oddhisto%d", i+1);
			oddcorrelationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create evencorrelationgraphs to be populated with the even-numbered entries in the dataset
		TH2F *evencorrelationgraph[30];
		for (int i = 0; i < 30; i++)
		{
			sprintf(histname, "evenhisto%d", i+1);
			evencorrelationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		// create a list whose entries correspond to algorithms and are the number of events in that algorithm's tail
	Double_t tailagreement[30];
		for (Int_t i = 0; i < 30; i++)
		{
			tailagreement[i] = 0;
		}

	int n = 0; // this variable will determine whether an event is even-numbered or odd-numbered
	Long64_t remainingentries;
	Long64_t nentries = runtree->GetEntries();
	for (Int_t i = 0; i < nentries; i++)
	{
		n = ( 1 - n ); // this logic changes n to be either 0 or 1
		runtree->GetEntry(i);

			// transverse mass based on metrefmuon:
			//transversemass = sqrt(2*metref*metrefw*(1+((mexref*mexrefw+meyref*meyrefw) / metref*metrefw))));
			// transverse mass based on metoffrecal:
			//transversemass = sqrt(2*metoff*metoffmuon*(1+((mexoff*mexoffmuon+meyoff*meyoffmuon) / (metoff*metoffmuon))));

		if ( /*(passmuonmed > 0.1 || passmuonvarmed > 0.1) && ml1 > 50. && 40. < transversemass && transversemass < 100. && muonclean > 0.1 && muonrecal < 0.1 && muonint > 35.*/ (zbl1gt10 > 0.1 || zbl1gt30 > 0.1 || zbl1gt40 > 0.1 || zbl1gt45 > 0.1) && zbl1 > 50.)
		{
			remainingentries++;
			Double_t sigma[6];
			Double_t metdist[6]; // metdist will be the distance of the event's MET from the median
			Double_t x[6]; // x = bulkmet and y = tailmet will be calculated for each algorithm
			Double_t y[6];

			// the following loop populates the sigma and metdist arrays
			for (Int_t j = 0; j < 6; j++)
			{
				if (sqrt(set[j]) >= 4.0) // throw out events whose SET values are too low
				{
					// compute sigma of this event for all algorithms
					sigma[j] = slope_ZeroBias[j]*sqrt(set[j]) + intercept_ZeroBias[j];
					//compute metdist for all algorithms
					metdist[j] = TMath::Abs( met[j] - (sigma[j]*TMath::Sqrt(TMath::PiOver2())));
				}
			}

/*
			if(metdist[0] >= 3*sigma[0]) // if the event is in the tail of Cell
			{
				y[0] = met[0];
				x[1] = met[1];
				correlationgraph[0]->Fill(x[1],y[0]);
			}
*/

			// the following logic populates correlationgraphs with (x = met, y = tailmet) tuples
			// only if they exist for a given event in the tree
			Int_t h = 0; // this variable counts each TH2F correlationgraph
			for (Int_t l = 0; l < 5; l++)
			{
				if (metdist[l] >= 3*sigma[l]) // if the event is in the tail of alg A
				{
					y[l] = met[l]; // save to y = tail met
					for (Int_t m = l+1; m < 6; m++)
					{
						if (metdist[m] >= 3*sigma[m]) // if the event is also in the tail of Alg B
						{
							tailagreement[h]++;
						}

							x[m] = met[m]; // save to x = met
							correlationgraph[h]->Fill(x[m], y[l]); // and populate the appropraite correlationgraph
							if (n == 1)
							{
								oddcorrelationgraph[h]->Fill(x[m], y[l]); // populate with odd-numbered entry
							}
							if (n == 0)
							{
								evencorrelationgraph[h]->Fill(x[m], y[l]); // populate with even-numbered entry
							}
						h++;
					}
				}
				else
				{
					h++;
				}
			 }

/*
			h = 15; // this variable counts each correlationgraph
			for (Int_t l = 0; l < 5; l++)
			{
				if (metdist[l] <= 3*sigma[l]) // if the event is in the bulk of A
				{
					x[l] = met[l]; // save to x = met
					for (Int_t m = l+1; m < 6; m++)
					{
						if (metdist[m] >= 3*sigma[m]) // if it's in the tail of B
						{
							y[m] = met[m]; // save to y = tail met
							correlationgraph[h]->Fill(y[l], x[m]); // and populate the appropraite correlationgraph
							if (n == 1)
							{
								oddcorrelationgraph[h]->Fill(x[l], y[m]); // populate with odd-numbered entry
							}
							if (n == 0)
							{
								evencorrelationgraph[h]->Fill(x[l], y[m]); // populate with even-numbered entry
							}
							h++;
						}
					}
				}
				else
				{
					h++;
				}
			 }
			 */
		 }
	  }


//======================================================================================================================================//

		TString xaxisNames[6] = {"Cell MET [GeV]", "MHT MET [GeV]", "Topocl MET [GeV]", "TopoclPS MET [GeV]", "TopoclPUC MET [GeV]", "TopoclEM MET [GeV]"};
		TString yaxisNames[6] = {"Cell Tail MET [GeV]", "MHT Tail MET [GeV]", "Topocl Tail MET [GeV]", "TopoclPS Tail MET [GeV]", "TopoclPUC Tail MET [GeV]", "TopoclEM Tail MET [GeV]"};

		ofstream correlationcoefficients; // prepare log file of correlation coefficients
		correlationcoefficients.open("correlationvalues.txt"); // open log file
		correlationcoefficients << "FILE:" << " " << graphtitle << "\n\n";
		correlationcoefficients << "Graph" << "\t" << "Correlation" << " " << "±" << " " << "Approx. Uncertainty" << " " << "\t" << "Correlation Graph" << "\t\t\t\t\t" << "Tail Fractions" "\n"; // write title of table



		TCanvas *mycanv[30];
		char *canvname = new char[30];
		Double_t r[30]; // correlation coefficients
		Double_t tailfractions[30]; // fraction of events in the tail of Alg A that are also in the tail of Alg B
		Int_t correlationentries;
		Double_t oddvalue[30]; // correlation values from oddcorrelationgraphs
		Double_t evenvalue[30]; // correlation values from evencorrelationgraphs
		Double_t c[30]; // final confidence interval of original correlation coefficients (r values)
		int k = 0; // this variable counts correlationgraphs
		for (int q = 0; q < 5; q++)
		{
			for (int l = q+1; l < 6; l++)
			{
				canvname = Form("canv%d",k+1);
				mycanv[k] = new TCanvas(canvname, "");
				correlationgraph[k]->Draw("colz"); // add "colz" in function if desired
				correlationgraph[k]->GetYaxis()->SetTitle(yaxisNames[q]);
				correlationgraph[k]->GetXaxis()->SetTitle(xaxisNames[l]);
				correlationgraph[k]->SetTitle(graphtitle);
				mycanv[k]->SetLogz();
				r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of each graph
				oddvalue[k] = oddcorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of odd graphs
				evenvalue[k] = evencorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of even graphs
				c[k] = 0.5*(oddvalue[k] - evenvalue[k]);
				correlationentries = correlationgraph[k]->GetEntries();
				tailfractions[k] = tailagreement[k]/correlationentries;
				mycanv[k]->Print(Form("%d.png", k+1));
				correlationcoefficients << k+1 << "\t" << r[k] << " " << "±" << " " << c[k]	<< ',' << "\t\t" << yaxisNames[q] << " " << "vs." << " " << xaxisNames[l] << "\t\t" << tailfractions[k] << "\n";
				k++;
			}
		}
		int k = 15;
		for (int u = 0; u < 5; u++)
		{
			for (int t = u+1; t < 6; t++)
			{
				canvname = Form("canv%d",k+1);
				mycanv[k] = new TCanvas(canvname, "");
				correlationgraph[k]->Draw("colz");  // add "colz" in function if desired
				correlationgraph[k]->GetYaxis()->SetTitle(xaxisNames[u]);
				correlationgraph[k]->GetXaxis()->SetTitle(yaxisNames[t]);
				correlationgraph[k]->SetTitle(graphtitle);
				mycanv[k]->SetLogz();
				r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of each graph
				oddvalue[k] = oddcorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of odd graphs
				evenvalue[k] = evencorrelationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of even graphs
				c[k] = 0.5*(oddvalue[k] - evenvalue[k]);
				correlationentries = correlationgraph[k]->GetEntries();
				tailfractions[k] = tailagreement[k]/correlationentries;
				mycanv[k]->Print(Form("%d.png", k+1));
				correlationcoefficients << k+1 << "\t" << r[k] << " " << "±" << " " << c[k]	<< ',' << "\t\t" << xaxisNames[u] << " " << "vs." << " " << yaxisNames[t] << "\t\t" << tailfractions[k] << "\n";
				k++;
			}
		}

		correlationcoefficients.close();

		return 0;
}
