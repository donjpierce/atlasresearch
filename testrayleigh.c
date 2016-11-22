int main()
{
	TFile *File1 = TFile::Open("Run2016R00300800_myOutput.root");
	

	TF1 *func = new TF1("func", "[0]*(x/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
		func->SetParameters(0, 100000.);
		func->SetParameters(1, 1.);
		//func->SetParameters(2,7.);
		func->SetParLimits(0, 0.1, 10000000.);
		func->SetParLimits(1, 0.1, 10000000.);
		plot->Fit("func", "L");

		int a = 0.1;
		int b = 1.1;

	while (a < 10, b < 11)
	{
		tree->Draw("metl1>>plot", "metl1<100.&&npv<a&&npv>b");

		Double_t p1 = func->GetParameter(1);
		Double_t m = sqrt(3.141592653 / 2)*p1;
		return m
		Double_t mean = plot->GetMean();
		return mean;
		Double_t p0 = func->GetParameter(0);
		return p0;
		Double_t entries = plot->GetEntries();
		return entries;

		a++;
		b++;
	}

}