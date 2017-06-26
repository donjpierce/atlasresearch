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


Int_t tailcalculator()
{
	gROOT->ProcessLine("gROOT->SetBatch(kTRUE)"); // suppresses the drawing of graphs
	TString PlotCut("passrndm>0.5"); // for 2015
	TString PlotCutmuons("passmu26med>0.5&&metl1>50."); // for 2016 muons

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

	// Defining a non-linear fit Function
	TF1 *nfit = new TF1("nfit", "[0]*(x + [2])*(x + [2]) + [1]");
	nfit->SetParameters(0, -50.);
	nfit->SetParameters(1, -50.);
	nfit->SetParLimits(0, -50000., 50000.);
	nfit->SetParLimits(1, -50000., 50000.);
	nfit->SetParLimits(2, -100., 100.);
	nfit->SetParName(0, "slope");
	nfit->SetParName(1, "intercept");
	nfit->SetParName(2, "shift");

	TFile *File1 = TFile::Open("../ZeroBias2016R307195R311481Runs56.root");
	TTree* tree = (TTree*)File1->Get("tree");
	//Produce fitting graphs 2015
	//TH2F *L1zb = new TH2F ("L1zb","", 60, 0., 60.,100,0.,100.);
		//tree->Draw("metl1:sqrt(setl1)>>L1zb","passrndm>0.5&&metl1>30.");
	TH2F *CELLzb = new TH2F ("CELLzb","", 50, 0., 50.,100,0.,100.);
		tree->Draw("metcell:sqrt(setcell)>>CELLzb", PlotCut);
	TH2F *MHTzb = new TH2F("MHTzb", "", 50, 0., 50., 100, 0., 100.);
		tree->Draw("metmht:sqrt(setmht)>>MHTzb", "passrndm>0.1&&metmht>0.1");
	TH2F *TOPOCLzb = new TH2F ("TOPOCLzb","", 70, 0., 70.,100,0.,100.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLzb",PlotCut);
	TH2F *TopoclEMzb = new TH2F("TopoclEMzb", "", 100, 0., 100., 100, 0., 100.);
		tree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMzb", PlotCut);
	TH2F *TOPOCLPSzb = new TH2F ("TOPOCLPSzb","", 70, 0., 70.,100,0.,100.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSzb", PlotCut);
	TH2F *TOPOCLPUCzb = new TH2F ("TOPOCLPUCzb","", 70, 0., 70.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCzb", "passrndm>0.1&&mettopoclpuc>0.1");
	//TH2F *OFFRECAL15 = new TH2F ("OFFRECAL15","", 50, 0., 50.,100,0.,100.);
		//tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL15",PlotCut);

	TFile *File2 = TFile::Open("../PhysicsMain2016.Muons.noalgL1XE45R3073065R311481Runs9B.root");
	//Fitting graphs 2016
	TH2F *L1muon = new TH2F ("L1muon","", 60, 0., 60.,1000,0.,1000.);
		tree->Draw("metl1:sqrt(setl1)>>L1muon",PlotCutmuons);
	TH2F *CELLmuon = new TH2F ("CELLmuon","",60,0.,60.,1000,0.,1000.);
		tree->Draw("metcell:sqrt(setcell)>>CELLmuon",PlotCutmuons);
	TH2F *MHTmuon = new TH2F("MHTmuon", "", 60, 0., 60., 1000, 0., 1000.);
		tree->Draw("metmht:sqrt(setmht)>>MHTmuon", PlotCutmuons);
	TH2F *TOPOCLmuon = new TH2F ("TOPOCLmuon","",60,0.,60.,1000,0.,1000.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLmuon",PlotCutmuons);
	TH2F *TopoclEMmuon = new TH2F("TopoclEMmuon", "", 100, 0., 100., 1000, 0., 1000.);
		tree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMmuon", PlotCutmuons);
	TH2F *TOPOCLPSmuon = new TH2F ("TOPOCLPSmuon","", 60, 0., 60., 1000, 0., 1000.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSmuon", PlotCutmuons);
	TH2F *TOPOCLPUCmuon = new TH2F ("TOPOCLPUCmuon","", 60, 0., 60., 1000, 0., 1000.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCmuon", PlotCutmuons);
	//TH2F *OFFRECAL16 = new TH2F ("OFFRECAL16","", 50, 0., 50.,100,0.,100.);
		//tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL16",PlotCut);

		//Fit Parameters
		Double_t shift_ZeroBias[6];
		Double_t shift_Muon[6];
		Double_t slope_ZeroBias[6];
		Double_t slope_Muon[6];
		Double_t intercept_ZeroBias[6];
		Double_t intercept_Muon[6];

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
		CELLzb_1->Fit("linfit");
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

		TCanvas *cCELLmuon = new TCanvas("cCELLmuon", "CELL 2016 ");
		CELLmuon->Draw();
		CELLmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *CELLmuon_1 = (TH1D*)gDirectory->Get("CELLmuon_1");
		CELLmuon_1->Draw();
		CELLmuon_1->Fit("linfit");
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

		//TOPOCL Algorithm resoltuions in ZeroBias and Muon
		TCanvas *cTOPOCLzb = new TCanvas("cTOPOCLzb", "TOPOCL 2015 ");
		TOPOCLzb->Draw();
		TOPOCLzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLzb_1 = (TH1D*)gDirectory->Get("TOPOCLzb_1");
		TOPOCLzb_1->Draw();
		TOPOCLzb_1->Fit("linfit");
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

		TCanvas *cTOPOCLmuon = new TCanvas("cTOPOCLmuon", "TOPOCL 2016 ");
		TOPOCLmuon->Draw();
		TOPOCLmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLmuon_1 = (TH1D*)gDirectory->Get("TOPOCLmuon_1");
		TOPOCLmuon_1->Draw();
		TOPOCLmuon_1->Fit("linfit");
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

		//TOPOCLPS Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cTOPOCLPSzb = new TCanvas("cTOPOCLPSzb", "TOPOCLPS 2015 ");
		TOPOCLPSzb->Draw();
		TOPOCLPSzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPSzb_1 = (TH1D*)gDirectory->Get("TOPOCLPSzb_1");
		TOPOCLPSzb_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPSzb_1->Draw();
		TOPOCLPSzb_1->Fit("linfit");
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

		TCanvas *cTOPOCLPSmuon = new TCanvas("cTOPOCLPSmuon", "TOPOCLPS 2016 ");
		TOPOCLPSmuon->Draw();
		TOPOCLPSmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPSmuon_1 = (TH1D*)gDirectory->Get("TOPOCLPSmuon_1");
		TOPOCLPSmuon_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPSmuon_1->Draw();
		TOPOCLPSmuon_1->Fit("linfit");
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

		//TOPOCLPUC Algorithm resoltuions in ZeroBias and Muons
		TCanvas *cTOPOCLPUCzb = new TCanvas("cTOPOCLPUCzb", "TOPOCLPUC 2015 ");
		TOPOCLPUCzb->Draw();
		TOPOCLPUCzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPUCzb_1 = (TH1D*)gDirectory->Get("TOPOCLPUCzb_1");
		TOPOCLPUCzb_1->Draw("");
		TOPOCLPUCzb_1->Fit("linfit");
		slope_ZeroBias[4] = linfit->GetParameter(0);
		intercept_ZeroBias[4] = linfit->GetParameter(1);
		//shift_ZeroBias[5] = nfit->GetParameter(2);
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

		TCanvas *cTOPOCLPUCmuon = new TCanvas("cTOPOCLPUCmuon", "TOPOCLPUC 2016 ");
		TOPOCLPUCmuon->Draw();
		TOPOCLPUCmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TOPOCLPUCmuon_1 = (TH1D*)gDirectory->Get("TOPOCLPUCmuon_1");
		TOPOCLPUCmuon_1->Draw();
		TOPOCLPUCmuon_1->Fit("linfit");
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


		//TopoclEM Algorithm resolutions in ZeroBias2016 and Muons2016
		TCanvas *cTopoclEMzb = new TCanvas("cTopoclEMzb", "TopoclEM ZeroBias2016");
		TopoclEMzb->Draw();
		TopoclEMzb->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TopoclEMzb_1 = (TH1D*)gDirectory->Get("TopoclEMzb_1");
		TopoclEMzb_1->Draw();
		TopoclEMzb_1->Fit("linfit");
		slope_ZeroBias[5] = linfit->GetParameter(0);
		intercept_ZeroBias[5] = linfit->GetParameter(1);
		//shift_ZeroBias[0] = nfit->GetParameter(2);
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

		TCanvas *cTopoclEMmuon = new TCanvas("cTopoclEMmuon", "TopoclEMMuon 2016 ");
		TopoclEMmuon->Draw();
		TopoclEMmuon->FitSlicesY(func, 0, -1, 10, "L");
		TH1D *TopoclEMmuon_1 = (TH1D*)gDirectory->Get("TopoclEMmuon_1");
		TopoclEMmuon_1->Draw();
		TopoclEMmuon_1->Fit("nfit");
		slope_Muon[5] = nfit->GetParameter(0);
		intercept_Muon[5] = nfit->GetParameter(1);
		shift_Muon[5] = nfit->GetParameter(2);
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


		//___Calculate Tail Events Based on Resolutions___
		Int_t passmu26Flag, passmuvarmedFlag;
		Float_t l1cut = 50.0; Float_t l1met;
/*
		TFile *zerobiasfile = TFile::Open("../ZeroBias2016R307195R311481Runs56.root");
		TTree* zeroBiasTree = NULL;
	    zerobiasfile->GetObject("tree",zeroBiasTree);
		TString zerobiasGraphTitle = "2016 ZeroBias (Runs56) NO L1 CUT";
*/
		TFile *muonFile = TFile::Open("../PhysicsMain2016.Muons.noalgL1XE45R3073065R311481Runs9B.root");
		TTree* muonTree = NULL;
		//EXPLICTLY SELECT THE TTREE CALLED "tree" FROM MUONFILE; STORE IT IN "muonTree"
		muonFile->GetObject("tree",muonTree);
		TString muonGraphTitle = "2016 Muons (L1XE45...Runs9B) L1 > 40GeV";

		muonTree->SetBranchAddress("passmu26med", &passmu26Flag);
		muonTree->SetBranchAddress("passmu26varmed", & passmuvarmedFlag);
		muonTree->SetBranchAddress("metl1", &l1met);

		TString metalgName[6] = {"metcell", "metmht", "mettopocl", "mettopoclps", "mettopoclpuc", "mettopoclem"};
		TString setalgName[6] = {"setcell", "setmht", "settopocl", "settopoclps", "settopoclpuc", "settopoclem"};

		// create arrays for MET and SET branches
		Float_t met[6];
		for (Int_t i = 0; i < 6; i++)
		{
			muonTree->SetBranchAddress(metalgName[i], &met[i]);
		}

		Float_t set[6];
		for (Int_t i = 0; i < 6; i++)
		{
			muonTree->SetBranchAddress(setalgName[i], &set[i]);
		}



		// create graphs which I will later populate with TailMET vs. MET of different algorithm pairs
		TH2F *correlationgraph[30];
		char *histname = new char[30];
		Int_t bins = 900;
		Double_t min = 0.;
		Double_t max = 900.;
		for (Int_t i = 0; i < 30; i++)
		{
			sprintf(histname, "histo%d", i+1);
			correlationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

	Long64_t nentries = muonTree->GetEntries();
	for (Int_t i = 0; i < nentries; i++)
	{
		muonTree->GetEntry(i);

		if ( ( passmu26Flag > 0.5 || passmuvarmedFlag > 0.5 ) && l1met > l1cut) // throw out events which don't pass either muon trigger
		{
			Double_t sigma[6];
			Double_t metdist[6]; // metdist will be the distance of the event's MET from the median
			Double_t x[6]; // x = bulkmet and y = tailmet will be calculated for each algorithm
			Double_t y[6];

			// the following loop populates the sigma and metdist arrays
			for (Int_t j = 0; j < 6; j++)
			{
				if (sqrt(set[j]) < 4.0) // throw out events whose SET values are too low
				{
					continue;
				}

				// commented-out regions allow for nonlinear calculations
				//if (j < 5)
				//{
					// compute sigma and metdist for l1, cell, mht, topocl, and topoclps
					sigma[j] = slope_ZeroBias[j]*sqrt(set[j]) + intercept_ZeroBias[j];
					metdist[j] = abs( met[j] - (sigma[j]*sqrt(TMath::PiOver2())));
				//}
				//else
				//{
					// compute sigma and metdist for topoclpuc whose fit is nonlinear
					//sigma[j] = slope_ZeroBias[j]*(sqrt(set[j]) + shift_ZeroBias[j])*(sqrt(set[j]) + shift_ZeroBias[j]) + intercept_ZeroBias[j];
					//metdist[j] = abs( met[j] - (sigma[j]*sqrt(1.57079633)) ); // 1.5707963 = pi/2
				//}

			}

			// the following logic populates correlationgraphs with (x = met, y = tailmet) touples only...
			// if they exist for a given event in the tree
			Int_t h = 0; // this variable counts each TH2F correlationgraph
			for (Int_t l = 0; l < 5; l++)
			{
				if (metdist[l] < 3*sigma[l]) // if the event is in the bulk
				{
					x[l] = met[l]; // save to x = met

					for (Int_t m = l+1; m < 6; m++)
					{
						if (metdist[m] > 3*sigma[m]) // if the event is in the tail of alg[m] (does not equal alg[l])
						{
							y[m] = met[m]; // save to y = tailmet
							correlationgraph[h]->Fill(x[l], y[m]); // and populate the appropraite correlationgraph
							h++;
						}
						else
						{
							h++;
						}
					}
				}
				else
				{
					y[l] = met[l]; // save to y = tailmet

					for (Int_t m = l+1; m < 6; m++) // for each remaining alg
					{
							x[m] = met[m]; // save to x = met
							correlationgraph[h]->Fill(x[l], y[m]); // and populate the appropraite correlationgraph
							h++;
					}
				}
			}

			h = 15; // this variable counts each correlationgraph
			for (Int_t l = 0; l < 5; l++)
			{
				if (metdist[l] > 3*sigma[l]) // if the event is in the tail of alg[l]
				{
					y[l] = met[l]; // save to y = tailmet

					for (Int_t m = l+1; m < 6; m++)
					{
						x[m] = met[m]; // save to x = met
						correlationgraph[h]->Fill(y[l], x[m]); // and populate the appropraite correlationgraph
						h++;
					}
				}
				else
				{
					x[l] = met[l]; // save to x = met

					for (Int_t m = l+1; m < 6; m++)
					{
						if (metdist[m] > 3*sigma[m]) // if the event is in the tail of alg[m]
						{
							y[m] = met[m]; // save to y = tailmet
							correlationgraph[h]->Fill(y[l], x[m]);
							h++;
						}
						else
						{
							h++;
						}
					}
				}
			}
		}
	}

		TString xaxisNames[6] = {"Cell MET [GeV]", "MHT MET [GeV]", "Topocl MET [GeV]", "TopoclPS MET [GeV]", "TopoclPUC MET [GeV]", "TopoclEM MET [GeV]"};
		TString yaxisNames[6] = {"Cell Tail MET [GeV]", "MHT Tail MET [GeV]", "Topocl Tail MET [GeV]", "TopoclPS Tail MET [GeV]", "TopoclPUC Tail MET [GeV]", "TopoclEM Tail MET [GeV]"};

		ofstream correlationcoefficients; // prepare log file of correlation coefficients
		correlationcoefficients.open("correlationvalues.txt"); // open log file
		correlationcoefficients << "Graph" << "\t" << "Correlation" << " " << "±" << " " << "90\% Confidence Range" << " " << "\t" << "FILE:" << "\t" << muonGraphTitle << "\n"; // write title of table

		TCanvas *mycanv[30];
		char *canvname = new char[30];
		Double_t r[30]; // correlation coefficients
		Double_t z[30]; // Fisher-transformed correlation variables
		Double_t zsigma[30]; // Fisher standard errors
		Double_t zlowerlimit[30]; // lower limit of Fisher variable
		Double_t zupperlimit[30]; // upper limit of Fisher variables
		Double_t c[30]; // confidence interval of original correlation coefficients
		Int_t entries[30]; // records number of entries in each correlationgraph
		Int_t k = 0; // this variable counts correlationgraphs
		for (Int_t q = 0; q < 5; q++)
		{
			for (Int_t l = q+1; l < 6; l++)
			{
				canvname = Form("canv%d",k+1);
				mycanv[k] = new TCanvas(canvname, "");
				correlationgraph[k]->Draw("colz"); // add "colz" in function if desired
				correlationgraph[k]->GetYaxis()->SetTitle(yaxisNames[q]);
				correlationgraph[k]->GetXaxis()->SetTitle(xaxisNames[l]);
				correlationgraph[k]->SetTitle(muonGraphTitle);
				entries[k] = correlationgraph[k]->GetEntries();
				mycanv[k]->SetLogz();
				r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of each graph
				z[k] = TMath::ATanH(r[k]); // Fisher-transformation of correlation factor
				zsigma[k] = 1/(sqrt(entries[k] - 3)); // sigma of gaussian-distributed Fisher variables
				zlowerlimit[k] = z[k] - 1.645*zsigma[k]/sqrt(entries[k]); // lower limit of Fisher variable (CALCULATED WITH 90% CONFIDENCE)
				zupperlimit[k] = z[k] + 1.645*zsigma[k]/sqrt(entries[k]); // upper limit of Fisher variable (CALCULATED WITH 90% CONFIDENCE)
				c[k] = ( exp(2*zlowerlimit[k]) - 1 ) / ( exp(2*zupperlimit[k]) + 1 ); // back-transformation of confidence interval
				mycanv[k]->Print(Form("%d.png", k+1));
				correlationcoefficients << k+1 << "\t" << r[k] << " " << "±" << " " << c[k]	<< ',' << "\t" << yaxisNames[q] << " " << "vs." << " " << xaxisNames[l] << "\n";
				k++;
			}
		}
		k = 15;
		for (Int_t u = 0; u < 5; u++)
		{
			for (Int_t t = u+1; t < 6; t++)
			{
				canvname = Form("canv%d",k+1);
				mycanv[k] = new TCanvas(canvname, "");
				correlationgraph[k]->Draw("colz");  // add "colz" in function if desired
				correlationgraph[k]->GetYaxis()->SetTitle(xaxisNames[u]);
				correlationgraph[k]->GetXaxis()->SetTitle(yaxisNames[t]);
				correlationgraph[k]->SetTitle(muonGraphTitle);
				entries[k] = correlationgraph[k]->GetEntries();
				mycanv[k]->SetLogz();
				r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2); // record correlation factor of each graph
				z[k] = TMath::ATanH(r[k]); // Fisher-transformation of correlation factor
				zsigma[k] = 1/(sqrt(entries[k] - 3)); // sigma of gaussian-distributed Fisher variables
				zlowerlimit[k] = z[k] - 1.645*zsigma[k]/sqrt(entries[k]); // lower limit of Fisher variable (CALCULATED WITH 90% CONFIDENCE)
				zupperlimit[k] = z[k] + 1.645*zsigma[k]/sqrt(entries[k]); // upper limit of Fisher variable (CALCULATED WITH 90% CONFIDENCE)
				c[k] = ( exp(2*zlowerlimit[k]) - 1 ) / ( exp(2*zupperlimit[k]) + 1 ); // back-transformation of confidence interval
				mycanv[k]->Print(Form("%d.png", k+1));
				correlationcoefficients << k+1 << "\t" << r[k] << " " << "±" << " " << c[k]	<< ',' << "\t" << xaxisNames[u] << " " << "vs." << " " << yaxisNames[t] << "\n";
				k++;
			}
		}

		correlationcoefficients.close();

		return 0;

}
