{
	#include <vector>

	TFile *file = TFile::Open("../ZeroBias2016R307195R311481Runs56.root");

	Float_t setalg;
	TBranch *b_setalg = new TBranch();
	tree->SetBranchAddress("settopoclpuc", &setalg, &b_setalg);

	Float_t metalg;
	TBranch *b_metalg = new TBranch();
	tree->SetBranchAddress("mettopoclpuc", &metalg, &b_metalg);

	int tail = 1;

	//==================================================================================================================================================//
	//Initialize scatter plot of alg met vs. sqrt alg set
	//***NOTE: CHANGE RANGE FOR SIGNAL EVENTS***
	TH1F *algset = new TH1F("algset", "sqrt of Algorithm SET", 100, 0., 100.);

	bool pass;
	Long64_t nentries = tree->GetEntries();
	for (int i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);

		// Only pass the events that are in the specified tail of the algorithm
		if ("passrndm" > 0.1)
		{
			if (metalg > tail)
			{
				algset->Fill(sqrt(setalg));
			}
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
			cout << binnum[j];
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

	for (int i = 0; i < nentries; i++)
	{
		if ("passrndm" > 0.1 && "mettopoclpuc" > 0.1)
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
	}

	for (int n = 1; n < nhist; n++)
	{
		//myhist[n]->SetTitle(metalg "for bin %d of sq rt SET");
		myhist[n]->GetYaxis()->SetTitle("Number of Events");
		myhist[n]->GetXaxis()->SetTitle("Cell MET [GeV]");
		myhist[n]->SetTitle(Form("Algorithm MET in bin %i of #sqrt{SumEt}", n));
	}

	// Initialize the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Plot events vs. MET/SIGMA
	double sigmaarray[nhist];
	TCanvas *mycanv[nhist];
	char *canvname = new char[nhist];
	for (int m = 1; m < nhist; m++)
	{
		sprintf(canvname, "canv%d", m);
		mycanv[m] = new TCanvas(canvname, "");
		mycanv[m]->SetLogy();
		myhist[m]->Fit("func", "L");
		sigmaarray[m] = func->GetParameter(1);
		mycanv[m]->Print(Form("%s.png", myhist[m]->GetName()));
	}

	ofstream sigmafile;
	sigmafile.open("sigmaarray.txt");
	sigmafile << "sigma" << "\n";
	for (int p = 1; p < nhist; p++)
	{
		sigmafile << sigmaarray[p] << "\n";
	}
	sigmafile.close();

}
