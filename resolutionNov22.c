{
		/*
		DP resolutionNov22.c
		This program will graph the resolution of various algorithms.
		Resolution is defined as the sigma of the Rayleigh-distributed algorithm's SUM Et for all sqrt(SUM Et)
		*/
	
		//FILE NAME (commented choice, 2015 or 2016)
		TString filename = "ZeroBias2015.p2634.PeriodJ.root";
		//TString filename = "ZeroBias2016.13Runs";
	
		//Bunch Group IDs
		TString bunch = "bcid";
		TString bunch_title = "Bunch Group ID";
	
		//Passed alg parameters
		TString passed_hist_name = "sqrt(setl1)";
		TString passed_hist_title = "sqrt(SUM Et L1)";

		//2D Histogram parameters
		int nxbins = 100;
		Double_t metxMin = 0.;
		Double_t metxMax = 100.;
		int nybins = 100;
		Double_t metyMin = 0.;
		Double_t metyMax = 100.;

		//Histogram Axis Labels
		TString yaxis "sigma of fit for algorithm sqrt(SUM Et) [sqrt(GeV)]";
		TString xaxis = "Algorithm sqrt(Sumet) [sqrt(GeV)]";

		//Assign Full Hist Params
		TString passed_hist_full_param = (passed_hist_title + " HIST;" + xaxis + ";" + yaxis);

	//===============================================================================================================================//
		//Open ROOT File
		TFile *file = TFile::Open(filename);

		//containers for an entry of set and a sigma of a fit
		Float_t passed_hist_sqrtset, raysigma;

		//Get tree address of variables
		tree->SetBranchAddress(passed_hist_name, &passed_hist_sqrtset);

		//Initialize histogram object
		TH2F *passed_hist = new TH2F(passed_hist_name, passed_hist_full_param, nxbins, metxMin, metxMax, nybins, metyMin, metyMax);

		//Initialize the Rayleigh Distribution:
		TF1 *raydist = new TF1("raydist", "[0]*(x/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
		raydist->SetParameters(0, 100000.);
		raydist->SetParameters(1, 1.);
		raydist->SetParLimits(0, 0.1, 10000000.);
		raydist->SetParLimits(1, 0.1, 10000000.);
		
		//Fit the distribution to the histogram
		passed_hist->Fit("raydist", "L");
		//Acquire the sigma of the fit
		Double_t sigma = raydist->GetParameter(1);
	
		//Fill the histogram with entries
		int nentries = tree->GetEntries();
		for (Long64_t i = 0; i < nentries; i++)
		{
			tree->GetEntry(i);
			//Only consider events in this loop that have been stored above a given Bunch Crossing ID (BCID)
			//For 2015, the 11th BCID begins at 2261. For 2016, the 11th BCID begins at 1514.
			if (bunch > 2261.)
			{
				passed_hist->Fill(passed_hist_sqrtset);
			}

			//Fill histogram with sqrt(set) bins in a one-to-one correspondence with a sigma of a fit of set
			

		}
	
		//Plot passed histogram
		passed_hist->Draw();
		
}