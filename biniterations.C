#include <iostream>
#include <fstream>
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
#include "TF1.h"


Int_t biniterations( const TString& metAlgName = "metcell" , const TString& setAlgName = "setcell",
 const TString& muonFilePath = "../ZeroBias2016R307195R311481Runs56.root" )
{

    gROOT->ProcessLine("gROOT->SetBatch(kTRUE)");

	const TFile *muonFile = TFile::Open( muonFilePath ,"READ" );
    //explicitly get the ttree called "tree" in muonFile, and store its address in myMuonTree
    const TTree* myMuonTree = (TTree*)muonFile->Get("tree");

	Float_t setalg , metalg;
	Int_t passrndm, numEntries, passmuflag,passmuvarmed;
	Int_t tail = 1;

	//set branch addresses
	myMuonTree->SetBranchAddress(setAlgName, &setalg);
	myMuonTree->SetBranchAddress(metAlgName, &metalg);
	myMuonTree->SetBranchAddress("passrndm",&passrndm);
	myMuonTree->SetBranchAddress("passmu26med",&passmuflag);
	myMuonTree->SetBranchAddress("passmu26varmed",&passmuvarmed);

	//==================================================================================================================================================//
	//Initialize scatter plot of alg met vs. sqrt alg set

	//**NOTE: CHANGE RANGE FOR SIGNAL EVENTS***
	TH1F *algSetHist = new TH1F("algSetHist", setAlgName , 100, 0., 100.);

	Long64_t muonNentries = myMuonTree->GetEntries();
	for (Int_t i = 0; i < muonNentries; i++)
	{
		myMuonTree->GetEntry(i);
		//are you supposed to cut on passrndm for muon events?
		if ( ( (passmuflag > 0.5) || (passmuvarmed > 0.5) ) && (metalg > tail))
		{
			algSetHist->Fill(sqrt(setalg));
		}
	}

	//If the bin content is nonzero, save the number of that bin in the binnum array
	Int_t nHists = 0;
 	vector<Int_t> binArray(100);
	for (Int_t j = 0; j < 100; j++)
	{
		numEntries = algSetHist->GetBinContent(j);
		if (numEntries > 0)
		{
			binArray[j] = algSetHist->GetBin(j);
			nHists++;
		}
		else
		{
			binArray[j] = 0;
		}
	}
	std::cout << "nHists: " << nHists << std::endl;
	//Produce an array of histograms for the algorithm MET,
	//where the length of the array is the number of non-zero elements in binnum
	TH1F* histArray[100];
	TString histname;
	Int_t nhistbins = 300;
	Float_t xmin = 0.0;
	Float_t xmax = 300.0;
	for (Int_t l = 1; l < nHists; l++)
	{
		histname = Form("histo%d",l);
		histArray[l] = new TH1F(histname, "", nhistbins, xmin, xmax);
	}

	// Parse through the tree again, and fill each MET histogram only with events in their respective sqrt(set) bins

	for (Int_t i = 0; i < muonNentries; i++)
	{
		myMuonTree->GetEntry(i);
		for (Int_t p = 1; p < sizeof(binArray) ; p++)
		{
            //don't use passrndm on muon events. must use the random muon triggers!
			if ((passrndm > 0.5) && (binArray[p] < sqrt(setalg)) && (sqrt(setalg) < binArray[p] + 1))
			{
				histArray[p]->Fill(metalg);
			}
		}
	}

	for (Int_t n = 1; n < nHists; n++)
	{
		//histArray[n]->SetTitle(metalg "for bin %d of sq rt SET");
		histArray[n]->GetYaxis()->SetTitle("Number of Events");
		histArray[n]->GetXaxis()->SetTitle("MET [GeV]");
	}

	// Initialize the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Plot events vs. MET/SIGMA
	Double_t sigmaarray[100];
	TCanvas *canvasArray[100];
	TString canvname;
    gROOT->ProcessLine(".> biniterations.log");
	for (Int_t m = 1; m < nHists; m++)
	{
		canvname = Form("canv%d",m);
		canvasArray[m] = new TCanvas(canvname, "");
		canvasArray[m]->SetLogy();
		histArray[m]->Fit("func", "L");
		sigmaarray[m] = func->GetParameter(1);
		canvasArray[m]->Print(Form("./Pictures/%s.png", histArray[m]->GetName()));
	}
    gROOT->ProcessLine(".>");

	ofstream sigmafile;
	sigmafile.open("sigmaarray.txt");
	sigmafile << "sigma" << "\n";
	for (Int_t p = 1; p < nHists; p++)
	{
		sigmafile << sigmaarray[p] << "\n";
	}
	sigmafile.close();
return 0;
}
