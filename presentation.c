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

	gROOT->ProcessLine("gROOT->SetBatch(kTRUE)"); // suppresses the drawing of graphs
	gROOT->ProcessLine("gROOT->Time();");
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
	TTree* zbTree = (TTree*)File1->Get("tree");
	//Produce fitting graphs 2015
	//TH2F *L1zb = new TH2F ("L1zb","", 60, 0., 60.,100,0.,100.);
		//zbTree->Draw("metl1:sqrt(setl1)>>L1zb","passrndm>0.5&&metl1>30.");
	TH2F *CELLzb = new TH2F ("CELLzb","", 50, 0., 50.,100,0.,100.);
		zbTree->Draw("metcell:sqrt(setcell)>>CELLzb", PlotCut);
	TH2F *MHTzb = new TH2F("MHTzb", "", 50, 0., 50., 100, 0., 100.);
		zbTree->Draw("metmht:sqrt(setmht)>>MHTzb", "passrndm>0.1&&metmht>0.1");
	TH2F *TOPOCLzb = new TH2F ("TOPOCLzb","", 70, 0., 70.,100,0.,100.);
		zbTree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLzb",PlotCut);
	TH2F *TopoclEMzb = new TH2F("TopoclEMzb", "", 100, 0., 100., 100, 0., 100.);
		zbTree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMzb", PlotCut);
	TH2F *TOPOCLPSzb = new TH2F ("TOPOCLPSzb","", 70, 0., 70.,100,0.,100.);
		zbTree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSzb", PlotCut);
	TH2F *TOPOCLPUCzb = new TH2F ("TOPOCLPUCzb","", 70, 0., 70.,100,0.,100.);
		zbTree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCzb", "passrndm>0.1&&mettopoclpuc>0.1");
	//TH2F *OFFRECAL15 = new TH2F ("OFFRECAL15","", 50, 0., 50.,100,0.,100.);
		//zbTree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL15",PlotCut);

/*
	TFile *muonFile = TFile::Open("../PhysicsMain2016.Muons.noalgL1XE45R3073065R311481Runs9B.root");
	//EXPLICTLY SELECT THE TTREE CALLED "tree" FROM MUONFILE; STORE IT IN "muonTree"
	TTree* muonTree = (TTree*)muonFile->Get("tree");
	//Fitting graphs 2016
	TH2F *L1muon = new TH2F ("L1muon","", 60, 0., 60.,1000,0.,1000.);
		muonTree->Draw("metl1:sqrt(setl1)>>L1muon",PlotCutmuons);
	TH2F *CELLmuon = new TH2F ("CELLmuon","",60,0.,60.,1000,0.,1000.);
		muonTree->Draw("metcell:sqrt(setcell)>>CELLmuon",PlotCutmuons);
	TH2F *MHTmuon = new TH2F("MHTmuon", "", 60, 0., 60., 1000, 0., 1000.);
		muonTree->Draw("metmht:sqrt(setmht)>>MHTmuon", PlotCutmuons);
	TH2F *TOPOCLmuon = new TH2F ("TOPOCLmuon","",60,0.,60.,1000,0.,1000.);
		muonTree->Draw("mettopocl:sqrt(settopocl)>>TOPOCLmuon",PlotCutmuons);
	TH2F *TopoclEMmuon = new TH2F("TopoclEMmuon", "", 100, 0., 100., 1000, 0., 1000.);
		muonTree->Draw("mettopoclem:sqrt(settopoclem)>>TopoclEMmuon", PlotCutmuons);
	TH2F *TOPOCLPSmuon = new TH2F ("TOPOCLPSmuon","", 60, 0., 60., 1000, 0., 1000.);
		muonTree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPSmuon", PlotCutmuons);
	TH2F *TOPOCLPUCmuon = new TH2F ("TOPOCLPUCmuon","", 60, 0., 60., 1000, 0., 1000.);
		muonTree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUCmuon", PlotCutmuons);
	//TH2F *OFFRECAL16 = new TH2F ("OFFRECAL16","", 50, 0., 50.,100,0.,100.);
		//muonTree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL16",PlotCut);
*/
		//Fit Parameters
		Double_t shift_ZeroBias[6];
		Double_t shift_Muon[6];
		Double_t slope_ZeroBias[6];
		Double_t slope_Muon[6];
		Double_t intercept_ZeroBias[6];
		Double_t intercept_Muon[6];









}
