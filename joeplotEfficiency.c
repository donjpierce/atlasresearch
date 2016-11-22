{
	/*
	Joseph Corrado plotEffiency.c
	Plots effiencies of various ATLAS algorithms versus the MET OFFLINE RECAL algorithm
	*/

	//ROOT FileName
	TString fileName = "ZeroBias2015.p2634.PeriodJ.root";

	//Passed hist parameters
	TString passed_hist_name = "metcell";
	TString passed_hist_title = "MET CELL";

	//Reference hist parameters 
	TString reference_hist_name = "metoffrecal";
	TString reference_hist_title = "METOFFRECAL";

	//Histogram parameters
	int nbins = 250;
	Double_t metMin = 0.0;
	Double_t metMax = 100.0;

	//Scatter Plot Parameters
	int nscatterBinsx = 250;
	int nscatterBinsy = 250;
	int scatterxmin = -5.0;
	int scatterxmax = 160.0;
	int scatterymin = -5.0;
	int scatterymax = 100.0;

	//Histogram Axis labels 
	TString xaxis = "MET [GeV]";
	TString yaxis = "Events";

	//Calculate Full Hist Params
	TString passed_hist_full_param = (passed_hist_title + " HIST; " + xaxis + ";" + yaxis);
	TString reference_hist_full_param = (reference_hist_title + " HIST;" + xaxis + ";" + yaxis);

//========================================================================================================================
	//TEfficiency Labels
	TString teff_name, teff_title;
	teff_name = passed_hist_name + " vs " + reference_hist_name;
	teff_title = passed_hist_title + " / " + reference_hist_title + " Efficiency" + ";" + xaxis + ";" + " Efficiency";

	//Scatter Plot Labels 
	scatterPlot_title = passed_hist_title + " vs " + reference_hist_title + " Scatter Plot" + ";" + reference_hist_title + ";" + passed_hist_title;
	TString scatterxLabel, scatteryLabel;
	scatterxLabel = reference_hist_title + " [GeV]";
	scatteryLabel = passed_hist_title + " [GeV]";

	//Open ROOT file 
	TFile * 2016Data = TFile::Open(fileName, "READ");

	//containers for 1st entry of met 
	Float_t passed_hist_met, reference_hist_met;

	//Get addresses of variables
	tree->SetBranchAddress(passed_hist_name, &passed_hist_met);
	tree->SetBranchAddress(reference_hist_name, &reference_hist_met);

	//Initialize TH1F objects
	TH1F *reference_hist = new TH1F(reference_hist_name, reference_hist_full_param, nbins, metMin, metMax);
	TH1F *passed_hist = new TH1F(passed_hist_name, passed_hist_full_param, nbins, metMin, metMax);

	//Initialize TEfficiency Object 
	TEfficiency *teff = new TEfficiency(teff_name, teff_title, nbins, metMin, metMax);

	//2D scatter plot of passed hist versus reference hist
	TH2F *scatterPlot = new TH2F(teff_name, scatterPlot_title, nscatterBinsx, scatterxmin, scatterxmax,nscatterBinsy , scatterymin,scatterymax);

	//Fill hist and effiency with entries
	bool pass;
	//# of entries
	int nentries = tree->GetEntries();
	for (Long64_t i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);
		//If the passed_hist_met is above a certain value, add the reference_hist_met from the corresponding entry 
		if (passed_hist_met > 80.0)
		{
			pass = true;
		}
		else
		{
			pass = false;
		}

		//Fill passed hist; increments bin corresponding to met by 1
		passed_hist->Fill(passed_hist_met);

		//Fill reference hist 
		reference_hist->Fill(reference_hist_met);

		//Fill 2D scatter Plot
		scatterPlot->Fill(reference_hist_met, passed_hist_met);

		//Fill TEfficiency 
		teff->Fill(pass, reference_hist_met);
	}

	//Divide canvas into pads
	TCanvas *plotCanvas = new TCanvas();
	plotCanvas->Divide(2,2);

	//Plot Reference Histogram on first pad
	plotCanvas->cd(1);
	reference_hist->Draw();

	//Plot passed histogram on second pad
	plotCanvas->cd(2);
	passed_hist->Draw();

	//Plot tefficiency curve on third pad
	plotCanvas->cd(3);
	teff->Draw();

	//Plot 2D scatter of histograms on fourth pad
	plotCanvas->cd(4);
	scatterPlot->SetMarkerColor(kGreen);
	scatterPlot->Draw("SURF2");

/*FOR LATER (Drawing n graphs):
	TGraph *g[10];
	Double_t x[10] = { 0,1,2,3,4,5,6,7,8,9 };
	Double_t y[10] = { 1,2,3,4,5,6,4,3,2,1 };
	TMultiGraph *mg = new TMultiGraph();
	for (int i = 0; i<10; i++) {
		g[i] = new TGraph(10, x, y);
		y[i] = y[i] + 1;
		mg->Add(g[i]);
	}
	mg->Draw("AL");
===================================================================================================*/
}
	
