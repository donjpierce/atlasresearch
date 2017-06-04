void biniterations(TString& metAlgName = "mettopocl" , TString& setAlgName = "settopocl")
{

	TFile * file = TFile::Open("C:/root_v5.34.36/myData/ZeroBias2016R307195R311481Runs56.root");

	Float_t setalg , metalg;
	Int_t passrndm, numEntries;
	Int_t tail = 1;

	tree->SetBranchAddress(setAlgName, &setalg);
	tree->SetBranchAddress(metAlgName, &metalg);
	tree->SetBranchAddress("passrndm",&passrndm);



	//==================================================================================================================================================//
	//Initialize scatter plot of alg met vs. sqrt alg set
	//**NOTE: CHANGE RANGE FOR SIGNAL EVENTS***
	TH1F *algSetHist = new TH1F("algSetHist", Form("sqrt of %s",setAlgName), 100, 0., 100.);

	Long64_t nentries = tree->GetEntries();
	for (Int_t i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);
		if ((passrndm > 0.1) && (metalg > tail))
		{
			algSetHist->Fill(sqrt(setalg));
		}
	}

	//If the bin content is nonzero, save the number of that bin in the binnum array
	Int_t nHists = 0;
	vector<Int_t> binArray;
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
			binnum[j] = 0;
		}
	}
	for (Int_t k = 0; k < 100; k++)
	{
		if (binnum[k] > 0)
		{
			binArray.push_back(binnum[k]);
		}
	}

	//Produce an array of histograms for the algorithm MET,
	//where the length of the array is the number of non-zero elements in binnum
	TH1F *myhist[nHists];
	TString *histname;
	Int_t nhistbins = 300;
	Float_t xmin = 0.0
	Float_t xmax = 300.0;
	for (Int_t l = 1; l < nHists; l++)
	{
		histname = Form("histo%d",l);
		myhist[l] = new TH1F(histname, "", nhistbins, xmin, xmax);
	}

	// Parse through the tree again, and fill each MET histogram only with events in their respective sqrt(set) bins

	for (Int_t i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);
		for (Int_t p = 1; p < binArray.size(); p++)
		{
			if ((passrndm > 0.5) && (binArray[p] < sqrt(setalg)) && (sqrt(setalg) < binArray[p] + 1))
			{
				myhist[p]->Fill(metalg);
			}
		}
	}

	for (Int_t n = 1; n < nHists; n++)
	{
		//myhist[n]->SetTitle(metalg "for bin %d of sq rt SET");
		myhist[n]->GetYaxis()->SetTitle("Number of Events");
		myhist[n]->GetXaxis()->SetTitle("MET [GeV]");
	}

	// Initialize the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Plot events vs. MET/SIGMA
	Double_t sigmaarray[nHists];
	TCanvas *mycanv[nHists];
	Char_t *canvname = new char[nHists];
	for (Int_t m = 1; m < nHists; m++)
	{
		canvname = Form("canv%d",m);
		mycanv[m] = new TCanvas(canvname, "");
		mycanv[m]->SetLogy();
		myhist[m]->Fit("func", "L");
		sigmaarray[m] = func->GetParameter(1);
		mycanv[m]->Print(Form("../Pictures/%s.png", myhist[m]->GetName()));
	}

	ofstream sigmafile;
	sigmafile.open("sigmaarray.txt");
	sigmafile << "sigma" << "\n";
	for (Int_t p = 1; p < nHists; p++)
	{
		sigmafile << sigmaarray[p] << "\n";
	}
	sigmafile.close();
}
