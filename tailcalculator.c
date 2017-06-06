{
	#include <vector>
	TString PlotCut("passrndm>0.1");

	// Define the Rayleigh Distribution
	TF1 *func = new TF1("func", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	func->SetParameters(0, 100000.);
	func->SetParameters(1, 1.);
	func->SetParLimits(0, 0.1, 10000000.);
	func->SetParLimits(1, 0.1, 10000000.);

	// Defining a Linear Fit Function
	TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
	linfit->SetParameters(0, -50.);
	linfit->SetParameters(1, -50.);
	linfit->SetParLimits(0, -50., 50.);
	linfit->SetParLimits(1, -50, 50.);
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

	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	//Produce fitting graphs 2015
	TH2F *L115 = new TH2F ("L115","", 20, 0., 20.,100,0.,100.);
		tree->Draw("metl1:sqrt(setl1)>>L115",PlotCut);
	TH2F *CELL15 = new TH2F ("CELL15","CELL2015", 50, 0., 50.,100,0.,100.);
		tree->Draw("metcell:sqrt(setcell)>>CELL15", PlotCut);
	TH2F *MHT15 = new TH2F("MHT15", "", 50, 0., 50., 100, 0., 100.);
		tree->Draw("metmht:sqrt(setmht)>>MHT15", "passrndm>0.1&&metmht>0.1");
	TH2F *TOPOCL15 = new TH2F ("TOPOCL15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCL15",PlotCut);
	TH2F *TOPOCLPS15 = new TH2F ("TOPOCLPS15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPS15", PlotCut);
	TH2F *TOPOCLPUC15 = new TH2F ("TOPOCLPUC15","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUC15", "passrndm>0.1&&mettopoclpuc>0.1");
	//TH2F *OFFRECAL15 = new TH2F ("OFFRECAL15","", 50, 0., 50.,100,0.,100.);
		//tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL15",PlotCut);

	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");
	//Fitting graphs 2016
	TH2F *L116 = new TH2F ("L116","", 20, 0., 20.,100,0.,100.);
		tree->Draw("metl1:sqrt(setl1)>>L116",PlotCut);
	TH2F *CELL16 = new TH2F ("CELL16","",50,0.,50.,100,0.,100.);
		tree->Draw("metcell:sqrt(setcell)>>CELL16",PlotCut);
	TH2F *MHT16 = new TH2F("MHT16", "", 50, 0., 50., 100, 0., 100.);
		tree->Draw("metmht:sqrt(setmht)>>MHT16", "passrndm>0.1&&metmht>0.1");
	TH2F *TOPOCL16 = new TH2F ("TOPOCL16","",50,0.,50.,100,0.,100.);
		tree->Draw("mettopocl:sqrt(settopocl)>>TOPOCL16",PlotCut);
	TH2F *TOPOCLPS16 = new TH2F ("TOPOCLPS16","",50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclps:sqrt(settopoclps)>>TOPOCLPS16", PlotCut);
	TH2F *TOPOCLPUC16 = new TH2F ("TOPOCLPUC16","", 50, 0., 50.,100,0.,100.);
		tree->Draw("mettopoclpuc:sqrt(settopoclpuc)>>TOPOCLPUC16", "passrndm>0.1&&mettopoclpuc>0.1");
	//TH2F *OFFRECAL16 = new TH2F ("OFFRECAL16","", 50, 0., 50.,100,0.,100.);
		//tree->Draw("metoffrecal:sqrt(setoffrecal)>>OFFRECAL16",PlotCut);

		//Linear Fit Parameters
		Double_t slope15[6];
		Double_t slope16[6];
		Double_t intercept15[6];
		Double_t intercept16[6];

//L1 Algorithm resolutions in 2015 and 2016
	TCanvas *cL115 = new TCanvas("cL115", "L1 2015");
	L115->Draw();
	L115->FitSlicesY(func, 0, -1, 10, "L");
	L115_1->Draw();
	L115_1->Fit("linfit");
	slope15[0] = linfit->GetParameter(0);
	intercept15[0] = linfit->GetParameter(1);
	L115_1->SetTitle("Resolution of L1 in 2015");
	L115_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L115_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L115_1->SetLineColor(2);
	gPad->Update();
		TPaveStats *l115 = (TPaveStats*)L115_1 ->FindObject("stats");
		l115->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resl115 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl115->AddEntry("L115_1", "2015 Data", "L");
		resl115->Draw();

	TCanvas *cL116 = new TCanvas("cL116", "L1 2016 ");
	L116->Draw();
	L116->FitSlicesY(func, 0, -1, 10, "L");
	L116_1->Draw();
	L116_1->Fit("linfit");
	slope16[0] = linfit->GetParameter(0);
	intercept16[0] = linfit->GetParameter(1);
	L116_1->SetTitle("Resolution of L1 in 2016 ");
	L116_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
	L116_1->GetYaxis()->SetTitle("#sigma of Fit for L1 [GeV]");
	L116_1->SetLineColor(4);
	gPad->Update();
		TPaveStats *l116 = (TPaveStats*)L116_1 ->FindObject("stats");
		l116->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resl116 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resl116->AddEntry("L116_1", "2016 Data", "L");
		resl116->Draw();


//CELL Algorithm resoltuions in 2015 and 2016
		TCanvas *cCELL15 = new TCanvas("cCELL15", "CELL 2015 ");
		CELL15->Draw();
		CELL15->FitSlicesY(func, 0, -1, 10, "L");
		CELL15_1->Draw();
		CELL15_1->Fit("linfit");
		slope15[1] = linfit->GetParameter(0);
		intercept15[1] = linfit->GetParameter(1);
		CELL15_1->SetTitle("Resolution of CELL in 2015 ");
		CELL15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELL15_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELL15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sCELL15 = (TPaveStats*)CELL15_1->FindObject("stats");
		sCELL15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resCELL15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELL15->AddEntry("CELL15_1", "2015 Data", "L");
		resCELL15->Draw();

		TCanvas *cCELL16 = new TCanvas("cCELL16", "CELL 2016 ");
		CELL16->Draw();
		CELL16->FitSlicesY(func, 0, -1, 10, "L");
		CELL16_1->Draw();
		CELL16_1->Fit("linfit");
		slope16[1] = linfit->GetParameter(0);
		intercept16[1] = linfit->GetParameter(1);
		CELL16_1->SetTitle("Resolution of CELL in 2016 ");
		CELL16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		CELL16_1->GetYaxis()->SetTitle("#sigma of Fit for CELL [GeV]");
		CELL16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sCELL16 = (TPaveStats*)CELL16_1->FindObject("stats");
		sCELL16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resCELL16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resCELL16->AddEntry("CELL16_1", "2016 Data", "L");
		resCELL16->Draw();

		//MHT Algorithm resoltuions in 2015 and 2016
		TCanvas *cMHT15 = new TCanvas("cMHT15", "MHT 2015 ");
		MHT15->Draw();
		MHT15->FitSlicesY(func, 0, -1, 10, "L");
		MHT15_1->Draw();
		MHT15_1->Fit("linfit");
		slope15[2] = linfit->GetParameter(0);
		intercept15[2] = linfit->GetParameter(1);
		MHT15_1->SetTitle("Resolution of MHT in 2015 ");
		MHT15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHT15_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHT15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sMHT15 = (TPaveStats*)MHT15_1->FindObject("stats");
		sMHT15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resMHT15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHT15->AddEntry("MHT15_1", "2015 Data", "L");
		resMHT15->Draw();

		TCanvas *cMHT16 = new TCanvas("cMHT16", "MHT 2016 ");
		MHT16->Draw();
		MHT16->FitSlicesY(func, 0, -1, 10, "L");
		MHT16_1->Draw();
		MHT16_1->Fit("linfit");
		slope16[2] = linfit->GetParameter(0);
		intercept16[2] = linfit->GetParameter(1);
		MHT16_1->SetTitle("Resolution of MHT in 2016 ");
		MHT16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		MHT16_1->GetYaxis()->SetTitle("#sigma of Fit for MHT [GeV]");
		MHT16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sMHT16 = (TPaveStats*)MHT16_1->FindObject("stats");
		sMHT16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resMHT16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resMHT16->AddEntry("MHT16_1", "2016 Data", "L");
		resMHT16->Draw();

		//TOPOCL Algorithm resoltuions in 2015 and 2016
		TCanvas *cTOPOCL15 = new TCanvas("cTOPOCL15", "TOPOCL 2015 ");
		TOPOCL15->Draw();
		TOPOCL15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCL15_1->Draw();
		TOPOCL15_1->Fit("linfit");
		slope15[3] = linfit->GetParameter(0);
		intercept15[3] = linfit->GetParameter(1);
		TOPOCL15_1->SetTitle("Resolution of TOPOCL in 2015 ");
		TOPOCL15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCL15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCL15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCL15 = (TPaveStats*)TOPOCL15_1->FindObject("stats");
		sTOPOCL15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCL15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCL15->AddEntry("TOPOCL15_1", "2015 Data", "L");
		resTOPOCL15->Draw();

		TCanvas *cTOPOCL16 = new TCanvas("cTOPOCL16", "TOPOCL 2016 ");
		TOPOCL16->Draw();
		TOPOCL16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCL16_1->Draw();
		TOPOCL16_1->Fit("linfit");
		slope16[3] = linfit->GetParameter(0);
		intercept16[3] = linfit->GetParameter(1);
		TOPOCL16_1->SetTitle("Resolution of TOPOCL in 2016 ");
		TOPOCL16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCL16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCL [GeV]");
		TOPOCL16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCL16 = (TPaveStats*)TOPOCL16_1->FindObject("stats");
		sTOPOCL16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCL16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCL16->AddEntry("TOPOCL16_1", "2016 Data", "L");
		resTOPOCL16->Draw();

		//TOPOCLPS Algorithm resoltuions in 2015 and 2016
		TCanvas *cTOPOCLPS15 = new TCanvas("cTOPOCLPS15", "TOPOCLPS 2015 ");
		TOPOCLPS15->Draw();
		TOPOCLPS15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPS15_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPS15_1->Draw();
		TOPOCLPS15_1->Fit("linfit");
		slope15[4] = linfit->GetParameter(0);
		intercept15[4] = linfit->GetParameter(1);
		TOPOCLPS15_1->SetTitle("Resolution of TOPOCLPS in 2015 ");
		TOPOCLPS15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPS15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPS15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPS15 = (TPaveStats*)TOPOCLPS15_1->FindObject("stats");
		sTOPOCLPS15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPS15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPS15->AddEntry("TOPOCLPS15_1", "2015 Data", "L");
		resTOPOCLPS15->Draw();

		TCanvas *cTOPOCLPS16 = new TCanvas("cTOPOCLPS16", "TOPOCLPS 2016 ");
		TOPOCLPS16->Draw();
		TOPOCLPS16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPS16_1->GetYaxis()->SetRange(0, 50.);
		TOPOCLPS16_1->Draw();
		TOPOCLPS16_1->Fit("linfit");
		slope16[4] = linfit->GetParameter(0);
		intercept16[4] = linfit->GetParameter(1);
		TOPOCLPS16_1->SetTitle("Resolution of TOPOCLPS in 2016 ");
		TOPOCLPS16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPS16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPS [GeV]");
		TOPOCLPS16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPS16 = (TPaveStats*)TOPOCLPS16_1->FindObject("stats");
		sTOPOCLPS16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPS16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPS16->AddEntry("TOPOCLPS16_1", "2016 Data", "L");
		resTOPOCLPS16->Draw();

		//TOPOCLPUC Algorithm resoltuions in 2015 and 2016
		TCanvas *cTOPOCLPUC15 = new TCanvas("cTOPOCLPUC15", "TOPOCLPUC 2015 ");
		TOPOCLPUC15->Draw();
		TOPOCLPUC15->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPUC15_1->Draw("");
		TOPOCLPUC15_1->Fit("nfit", "", "", 10., 100.);
		slope15[5] = nfit->GetParameter(0);
		intercept15[5] = nfit->GetParameter(1);
		Double_t shift15 = nfit->GetParameter(2);
		TOPOCLPUC15_1->SetTitle("Resolution of TOPOCLPUC in 2015 ");
		TOPOCLPUC15_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUC15_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUC15_1->SetLineColor(2);
		gPad->Update();
		TPaveStats *sTOPOCLPUC15 = (TPaveStats*)TOPOCLPUC15_1->FindObject("stats");
		sTOPOCLPUC15->SetTextColor(2);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUC15 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUC15->AddEntry("TOPOCLPUC15_1", "2015 Data", "L");
		resTOPOCLPUC15->Draw();

		TCanvas *cTOPOCLPUC16 = new TCanvas("cTOPOCLPUC16", "TOPOCLPUC 2016 ");
		TOPOCLPUC16->Draw();
		TOPOCLPUC16->FitSlicesY(func, 0, -1, 10, "L");
		TOPOCLPUC16_1->Draw();
		TOPOCLPUC16_1->Fit("nfit", "", "", 10., 100.);
		slope16[5] = nfit->GetParameter(0);
		intercept16[5] = nfit->GetParameter(1);
		Double_t shift16 = nfit->GetParameter(2);
		TOPOCLPUC16_1->SetTitle("Resolution of TOPOCLPUC in 2016 ");
		TOPOCLPUC16_1->GetXaxis()->SetTitle("#sqrt{SumEt} #left[#sqrt{GeV} #right]");
		TOPOCLPUC16_1->GetYaxis()->SetTitle("#sigma of Fit for TOPOCLPUC [GeV]");
		TOPOCLPUC16_1->SetLineColor(4);
		gPad->Update();
		TPaveStats *sTOPOCLPUC16 = (TPaveStats*)TOPOCLPUC16_1->FindObject("stats");
		sTOPOCLPUC16->SetTextColor(4);
		gStyle->SetOptFit(11);

		TLegend* resTOPOCLPUC16 = new TLegend(0.37, 0.7, 0.55, 0.88);
		resTOPOCLPUC16->AddEntry("TOPOCLPUC16_1", "2016 Data", "L");
		resTOPOCLPUC16->Draw();

		// Calculate Tail Events Based on Resolutions
		Handle_t setalg[6] = {"setl1", "setcell", "setmht", "settopocl", "settopoclps", "settopoclpuc"};
		Handle_t metalg[6] = {"metl1", "metcell", "metmht", "mettopocl", "mettopoclps", "mettopoclpuc"};

		TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");

		TH2F *correlationgraph[30];
		char *histname = new char[30];
		int bins = 200;
		Double_t min = 0.;
		Double_t max = 200.;
		for (int i = 0; i < 30; i++)
		{
			sprintf(histname, "histo%d", i);
			correlationgraph[i] = new TH2F(histname, "", bins, min, max, bins, min, max);
		}

		bool pass;
		Long64_t nentries = tree->GetEntries();
		for (int i = 0; i < nentries; i++)
		{
			tree->GetEntry(i);

			Double_t sigma[6];
			Double_t metdist[6]; // metdist will be the distance of the event's MET from the median
			Double_t x[6]; // x = bulkmet and y = tailmet will be calculated for each algorithm
			Double_t y[6];

			// the following loop populates the sigma and metdist arrays
			for (int j = 0; j < 6; j++)
			{
				if (j < 5)
				{
					// compute sigma and metdist for l1, cell, mht, topocl, and topoclps
					sigma[j] = slope15[j]*sqrt(setalg[j]) + intercept15[j];
					metdist[j] = abs( metalg[j] - (sigma[j]*sqrt(1.57079633)) ); // 1.5707963 = pi/2
				}
				else
				{
					// compute sigma and metdist for topoclpuc whose fit is nonlinear
					sigma[j] = slope15[j]*(x + shift15)*(x + shift15) + intercept15[j];
					metdist[j] = abs( metalg[j] - (sigma[j]*sqrt(1.57079633)) );
				}
			}

			// the following logic populates correlationgraphs with (x = bulkmet, y = tailmet) touples only...
			// if they exist for a given event in the tree

			int z = 0; // this variable counts each correlationgraph
			for (int l = 0; l < 5; l++)
			{
				if (metdist[l] < 3*sigma[l]) // if the event is in the bulk of alg[l]
				{
					x[l] = metalg[l]; // save to x = bulkmet

					for (int m = l+1; m < 6; m++)
					{
						if (metdist[m] > 3*sigma[m]) // if the event is in the tail of alg[m] DNE alg[l]
						{
							y[m] = metalg[m]; // save to y = tailmet
							correlationgraph[z]->Fill(x[l], y[m]); // and populate the appropraite correlationgraph
							z++;
						}
						else
						{
							z++;
						}
					}
				}
				else
				{
					y[l] = metalg[l]; // event is in the tail of alg[l]

					for (int m = l+1; m < 6; m++)
					{
						if (metdist[m] < 3*sigma[m])
						{
							// IF the event is in the BULK of alg[]
							x[m] = metalg[m];
							correlationgraph[z]->Fill(y[l], x[m]);
							z++;
						}
						else
						{
							z++;
						}
					}
				}
			}


			int k = 0;
				for (int l = 0; l < 5; l++)
				{
					for (int m = l+1; m < 6; m++)
					{
						correlationgraph[k]->Fill(x[l], y[m]);
						k ++;
					}
				}

				for (int l = 0; l < 5; l++)
				{
					for (int m = l+1; m < 6; m++)
					{
						correlationgraph[k]->Fill(y[l], x[m]);
						k ++;
					}
				}

				correlationgraph[0]->Fill(x[0], y[1]); // l1 bulk cell tail
				correlationgraph[1]->Fill(x[0], y[2]); // l1 bulk mht tail
				correlationgraph[2]->Fill(x[0], y[3]); // l1 bulk topocl tail
				correlationgraph[3]->Fill(x[0], y[4]); // l1 bulk topoclps tail
				correlationgraph[4]->Fill(x[0], y[5]); // l1 bulk topoclpuc tail

				correlationgraph[5]->Fill(x[1], y[2]); // cell bulk mht tail
				correlationgraph[6]->Fill(x[1], y[3]); // cell bulk topocl tail
				correlationgraph[7]->Fill(x[1], y[4]); // cell bulk topoclps tail
				correlationgraph[8]->Fill(x[1], y[5]); // cell bulk topoclpuc tail

				correlationgraph[9]->Fill(x[2], y[3]); // mht bulk topocl tail
				correlationgraph[10]->Fill(x[2], y[4]); // mht bulk topoclps tail
				correlationgraph[11]->Fill(x[2], y[5]); // mht bulk topoclpuc tail

				correlationgraph[12]->Fill(x[3], y[4]); // topocl bulk topoclps tail
				correlationgraph[13]->Fill(x[3], y[5]); // topocl bulk topoclpuc tail

				correlationgraph[14]->Fill(x[4], y[5]); // topoclps bulk topoclpuc tail
				//=====================================================================//

				correlationgraph[15]->Fill(y[0], x[1]); // l1 tail cell bulk
				correlationgraph[16]->Fill(y[0], x[2]); // l1 tail mht bulk
				correlationgraph[17]->Fill(y[0], x[3]); // l1 tail topocl bulk
				correlationgraph[18]->Fill(y[0], x[4]); // l1 tail topoclps bulk
				correlationgraph[19]->Fill(y[0], x[5]); // l1 tail topoclpuc bulk

				correlationgraph[20]->Fill(y[1], x[2]); // cell tail mht bulk
				correlationgraph[21]->Fill(y[1], x[3]); // cell tail topocl bulk
				correlationgraph[22]->Fill(y[1], x[4]); // cell tail topoclps bulk
				correlationgraph[23]->Fill(y[1], x[5]); // cell tail topoclpuc bulk

				correlationgraph[24]->Fill(y[2], x[3]); // mht tail topocl bulk
				correlationgraph[25]->Fill(y[2], x[4]); // mht tail topoclps bulk
				correlationgraph[26]->Fill(y[2], x[5]); // mht tail topoclpuc bulk

				correlationgraph[27]->Fill(y[3], x[4]); // topocl tail topoclps bulk
				correlationgraph[28]->Fill(y[3], x[5]); // topocl tail topoclpuc bulk

				correlationgraph[29]->Fill(y[4], x[5]); // topoclps tail topoclpuc bulk


		}

			// record correlation coefficients for each graph
			Double_t r[30];
			for (int k = 0; k < 30; k++)
			{
				r[k] = correlationgraph[k]->GetCorrelationFactor(1, 2);
			}

		}

		TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");




}
