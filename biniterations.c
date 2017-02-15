{
	#include <vector>
	
	TFile * 2015 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

	Float_t setalg;
	TBranch *b_setalg = new TBranch();
	tree->SetBranchAddress("settopoclps", &setalg, &b_setalg);

	Float_t metalg;
	TBranch *b_metalg = new TBranch();
	tree->SetBranchAddress("mettopoclps", &metalg, &b_metalg);

	int tail = 121;

	//==================================================================================================================================================//
	//Initialize scatter plot of alg met vs. sqrt alg set ***NOTE: CHANGE RANGE FOR SIGNAL EVENTS***	
	TH1F *algset = new TH1F("algset", "sqrt of Algorithm SET", 100, 0., 100.);
	TH1F *algmet = new TH1F("algmet", "Algorith MET", 200, 0., 200.);

	bool pass;
	Long64_t nentries = tree->GetEntries();
	for (int i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);

		// Only pass the events that are in the specified tail of the algorithm
		if (metalg > tail)
		{
			algset->Fill(sqrt(setalg));
		}
	}

	//Get the bin content of each bin of sq.rt(SET)
	//If the bin content is nonzero, save the number of that bin in the binnum array
	int bincontents[100];
	int binnum[100];
	for (int j = 0; j < 100; j++)
	{
		int binpop = algset->GetBinContent(j);
		bincontents[j] = binpop;

		if (binpop != 0)
		{
			binnum[j] = algset->GetBin(j);
		}
		else
		{
			binnum[j] = 0;
		}
	}

	//Count the non-zero elements of binnum for the number of histograms we will generate
	
	int counter = 0;
	for (int k = 1; k < 100; k++)
	{
		if (binnum[k] > 0)
		{
			counter++;
		}
	}
	int nhist = counter;
	

	vector<Handle_t> binarray;
	for (int k = 1; k < 100; k++)
	{
		if (binnum[k] > 0)
		{
			binarray.push_back(binnum[k]);
		}
	}

	//Produce an array of histograms for the algorithm MET,
	//where the length of the array is the number of non-zero elements in binnum 
	TH1F *myhist[nhist];
	char *histname = new char[nhist];
	int nhistbins = 300;
	float xmin = 0., xmax = 300.;
	for (int l = 1; l < nhist; l++)
	{
		sprintf(histname, "histo%d", l);
		myhist[l] = new TH1F(histname, "", nhistbins, xmin, xmax);
	}

	// Parse through the tree again, and fill each MET histogram only with events in their respective sqrt(set) bins
	Long64_t nentries = tree->GetEntries();
	for (int i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);		
		for (int p = 1; p < nhist; p++)
		{
				if (binarray[p] < sqrt(setalg) && sqrt(setalg) < binarray[p] + 1)
				{
					myhist[p]->Fill(metalg);
				}
			
		}
	}

	// Fit each histgram to a rayleigh distribution and save their sigmas to an array
	// Initialize the Rayleigh Distribution
	TF1 *raydist = new TF1("raydist", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	raydist->SetParameters(0, 100000.);
	raydist->SetParameters(1, 1.);
	raydist->SetParLimits(0, 0.1, 10000000.);
	raydist->SetParLimits(1, 0.1, 10000000.);


	//TString a = "bin";
	//TString b = ".png";
	TCanvas *mycanv[nhist];
	char *canvname = new char[nhist];
	for (int m = 1; m < nhist; m++)
	{
		sprintf(canvname, "canv%d", m);
		mycanv[m] = new TCanvas(canvname, "");
		myhist[m]->Fit("raydist", "L");
		//TString filename = a + binarray[m] + b;
		myhist[m]->Draw();
	}



}