  {   

    TFile *File1 = TFile::Open("ZeroBias2016.13Runs.root");

    TH1F *plot = new TH1F ("plot","",100,0.,100.);
    tree->Draw("metl1>>plot","passrndm>0.1&&bcid>410.&&bcid<620.");

	//Normal Distribution:
   //TF1 *func = new TF1 ("func","([0]/([2]*sqrt(2*pi)))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))");
	//Rayleigh Distribution:
     TF1 *func = new TF1 ("func","[0]*(x/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
        func->SetParameters(0,100000.);
        func->SetParameters(1,1.);
        //func->SetParameters(2,7.);
        func->SetParLimits(0,0.1,10000000.);
        func->SetParLimits(1,0.1,10000000.);
	plot->Fit("func","L");

	Double_t p1 = func->GetParameter(1);
	cout << endl << sqrt(3.141592653 / 2)*p1 << endl;
	Double_t mean = plot->GetMean();
	cout << mean << endl;
	Double_t p0 = func->GetParameter(0);
	cout << p0 << endl;
	Double_t entries = plot->GetEntries();
	cout << entries << endl;

	TH2F *resolution = new TH2F("resolution", "L1 resolution BG = 4", 100, 0., 100.);
	tree->Draw("metl1:p1>>resolution", "passrndm>0.1&&bcid>410.&&bcid<620.");//HOW CAN I GET p1 INTO THIS DRAW FUNCTION?


    //TCanvas *c0 = new TCanvas("c0,c0");
    //plot->Fit("gaus");

 /* TH2F *plot1 = new TH2F ("plot1","",25,0.,25.,100,0.,100.);
        tree->Draw("metl1:npv>>plot1","metl1<100.");
        func->SetParameters(0,100000.);
        func->SetParameters(1,1.);
        //func->SetParameters(2,7.);
        func->SetParLimits(0,0.1,10000000.);
        func->SetParLimits(1,0.1,10000000.);

        
        plot1->FitSlicesY(func,0,-1,10,"L");
    plot1_1->Draw();*/

    
 

    
}

