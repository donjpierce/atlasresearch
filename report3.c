//mean comparison for L1 compared across runs only  
{
	  TFile *File1 = TFile::Open("Run2016R00300800_myOutput.root");
		TH1F *plot = new TH1F("plot", "", 25, 0., 25.);
		 tree->Draw("metl1>>plot", "metl1<100.&&npv>0.9&&npv<1.1&&passrndmj40>0.1");

	  //Gaussian Fit
	  //TF1 *func = new TF1 ("func","([0]/([2]*sqrt(6.28)))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))");
	  //Rayleigh Fit
	  TF1 *func = new TF1("func", "[0]*(x/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
	  func->SetParameters(0, 100000.);
	  func->SetParameters(1, 1.);
	  func->SetParLimits(0, 0.1, 10000000.);
	  func->SetParLimits(1, 0.1, 10000000.);

	  Double_t p1 = func->GetParameter(1);
	  cout << endl << sqrt(3.141592653 / 2)*p1 << endl;
	  Double_t mean = plot->GetMean();
	  cout << mean << endl;

	  Double_t p0 = func->GetParameter(0);
	  cout << p0 << endl;
	  Double_t entries = plot->GetEntries();
	  cout << entries << endl;

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