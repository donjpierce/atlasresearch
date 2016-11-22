{
	TFile *File1 = TFile::Open("ZeroBias2015.p2634.PeriodJ.root");
	TCanvas *CompareL1 = new TCanvas ("CompareL1","CompareL1");
	TH1F *L115 = new TH1F ("L115","",30,0.,30.);
		tree->Draw("metl1>>L115","passrndm>0.1&&sqrt(setl1)>4.&&sqrt(setl1)<5.");
	TH1F *CELL15 = new TH1F ("CELL15","",50,0.,50.);
		tree->Draw("metcell>>CELL15","passrndm>0.1&&sqrt(setcell)>14.&&sqrt(setcell)<15.");
	TH1F *TOPOCL15 = new TH1F ("TOPOCL15","",100,0.,100.);
		tree->Draw("mettopocl>>TOPOCL15","passrndm>0.1&&sqrt(settopocl)>19.&&sqrt(settopocl)<20.");
	TH1F *PS15 = new TH1F ("PS15","",100,0.,100.);
		tree->Draw("mettopoclps>>PS15","passrndm>0.1&&sqrt(settopoclps)>17.&&sqrt(settopoclps)<18.");
	TH1F *PUC15 = new TH1F ("PUC15","",100,0.,100.);
		tree->Draw("mettopoclpuc>>PUC15","passrndm>0.1&&mettopoclpuc>0.1&&sqrt(settopoclpuc)>25.&&sqrt(settopoclpuc)<26.");
	TH1F *MHT15 = new TH1F ("MHT15","",100,0.,100.);
		tree->Draw("metmht>>MHT15","passrndm>0.1&&sqrt(setmht)>9.&&sqrt(setmht)<10.");
	TH1F *OFFRECAL15 = new TH1F ("OFFRECAL15","",100,0.,100.);
		tree->Draw("metoffrecal>>OFFRECAL15","passrndm>0.1&&sqrt(setoffrecal)>20.&&sqrt(setoffrecal)<21.");

	TFile *File2 = TFile::Open("ZeroBias2016.13Runs.root");
	TH1F *L116 = new TH1F ("L116","",30,0.,30.);
		tree->Draw("metl1>>L116","passrndm>0.1&&sqrt(setl1)>4.&&sqrt(setl1)<5.");
	TH1F *CELL16 = new TH1F ("CELL16","",50,0.,50.);
		tree->Draw("metcell>>CELL16","passrndm>0.1&&sqrt(setcell)>14.&&sqrt(setcell)<15.");
	TH1F *TOPOCL16 = new TH1F ("TOPOCL16","",100,0.,100.);
		tree->Draw("mettopocl>>TOPOCL16","passrndm>0.1&&sqrt(settopocl)>19.&&sqrt(settopocl)<20.");
	TH1F *PS16 = new TH1F ("PS16","",100,0.,100.);
		tree->Draw("mettopoclps>>PS16","passrndm>0.1&&sqrt(settopoclps)>17.&&sqrt(settopoclps)<18.");
	TH1F *PUC16 = new TH1F ("PUC16","",100,0.,100.);
		tree->Draw("mettopoclpuc>>PUC16","passrndm>0.1&&mettopoclpuc>0.1&&sqrt(settopoclpuc)>25.&&sqrt(settopoclpuc)<26.");
	TH1F *MHT16 = new TH1F ("MHT16","",100,0.,100.);
		tree->Draw("metmht>>MHT16","passrndm>0.1&&sqrt(setmht)>9.&&sqrt(setmht)<10.");
	TH1F *OFFRECAL16 = new TH1F ("OFFRECAL16","",100,0.,100.);
		tree->Draw("metoffrecal>>OFFRECAL16","passrndm>0.1&&sqrt(setoffrecal)>20.&&sqrt(setoffrecal)<21.");

	TF1 *func = new TF1 ("func","(x/[1])*(1/[1])*exp(-.5*(x/[1])*(x/[1]))");
	//func->SetParameters(0,100000.);
    func->SetParameters(1,5.);
    func->SetParLimits(1,0.1,100.);

    L115->Draw();
    L115->GetXaxis()->SetTitle("L1 [GeV]");
    L115->GetYaxis()->SetTitle("Fraction of Events per GeV");
    L115->SetNormFactor(1.0);
    L115->SetLineColor(2);
    L115->Fit("func","L");
    	gPad->Update();
		TPaveStats *l115 = (TPaveStats*)L115 ->FindObject("stats");
		l115->SetTextColor(2);
		l115->SetX1NDC(.78);
		l115->SetX2NDC(.98);
		l115->SetY1NDC(.95);
		l115->SetY2NDC(.78);
    L116->Draw("sames");
    L116->SetNormFactor(1.0);
    L116->SetLineColor(4);
    func->SetLineColor(4);
    L116->Fit("func","L");
    	gPad->Update();
		TPaveStats *l116 = (TPaveStats*)L116 ->FindObject("stats");
		l116->SetTextColor(4);
		l116->SetX1NDC(.78);
		l116->SetX2NDC(.98);
		l116->SetY1NDC(.6);
		l116->SetY2NDC(.76);
	gStyle->SetOptFit(11);

    	TLegend* compl1=new TLegend(0.5,0.7,0.7,0.88);
    	compl1->AddEntry("L115","2015 Data","L");
    	compl1->AddEntry("L116","2016 Data","L");
    	compl1->Draw();

    TCanvas *CompareCELL = new TCanvas ("CompareCELL","CompareCELL");
    CELL15->Draw();
    CELL15->GetXaxis()->SetTitle("CELL [GeV]");
    CELL15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    CELL15->SetNormFactor(1.0);
    CELL15->SetLineColor(2);
    CELL15->Fit("func","L");
        gPad->Update();
		TPaveStats *cell15 = (TPaveStats*)CELL15 ->FindObject("stats");
		cell15->SetTextColor(2);
		cell15->SetX1NDC(.78);
		cell15->SetX2NDC(.98);
		cell15->SetY1NDC(.95);
		cell15->SetY2NDC(.78);
    CELL16->Draw("sames");
    CELL16->SetNormFactor(1.0);
    CELL16->SetLineColor(4);
    func->SetLineColor(4);
    CELL16->Fit("func","L");
        gPad->Update();
		TPaveStats *cell16 = (TPaveStats*)CELL16 ->FindObject("stats");
		cell16->SetTextColor(4);
		cell16->SetX1NDC(.78);
		cell16->SetX2NDC(.98);
		cell16->SetY1NDC(.6);
		cell16->SetY2NDC(.76);

    	TLegend* compcell=new TLegend(0.5,0.7,0.7,0.88);
    	compcell->AddEntry("CELL15","2015 Data","L");
    	compcell->AddEntry("CELL16","2016 Data","L");
    	compcell->Draw();

    TCanvas *CompareTOPOCL = new TCanvas ("CompareTOPOCL","CompareTOPOCL");
    TOPOCL15->Draw();
    CELL15->GetXaxis()->SetTitle("TOPOCL [GeV]");
    CELL15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    TOPOCL15->SetNormFactor(1.0);
    TOPOCL15->SetLineColor(2);
    TOPOCL15->Fit("func","L");
        gPad->Update();
		TPaveStats *topocl15 = (TPaveStats*)TOPOCL15 ->FindObject("stats");
		topocl15->SetTextColor(2);
		topocl15->SetX1NDC(.78);
		topocl15->SetX2NDC(.98);
		topocl15->SetY1NDC(.95);
		topocl15->SetY2NDC(.78);
    TOPOCL16->Draw("sames");
    TOPOCL16->SetNormFactor(1.0);
    TOPOCL16->SetLineColor(4);
    func->SetLineColor(4);
    TOPOCL16->Fit("func","L");
        gPad->Update();
		TPaveStats *topocl16 = (TPaveStats*)TOPOCL16 ->FindObject("stats");
		topocl16->SetTextColor(4);
		topocl16->SetX1NDC(.78);
		topocl16->SetX2NDC(.98);
		topocl16->SetY1NDC(.6);
		topocl16->SetY2NDC(.76);

    	TLegend* comptopocl=new TLegend(0.5,0.7,0.7,0.88);
    	comptopocl->AddEntry("TOPOCL15","2015 Data","L");
    	comptopocl->AddEntry("TOPOCL16","2016 Data","L");
    	comptopocl->Draw();

    TCanvas *ComparePS = new TCanvas ("ComparePS","ComparePS");
    PS15->Draw();
    PS15->GetXaxis()->SetTitle("PS [GeV]");
    PS15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    PS15->SetNormFactor(1.0);
    PS15->SetLineColor(2);
    PS15->Fit("func","L");
        gPad->Update();
		TPaveStats *ps15 = (TPaveStats*)PS15 ->FindObject("stats");
		ps15->SetTextColor(2);
		ps15->SetX1NDC(.78);
		ps15->SetX2NDC(.98);
		ps15->SetY1NDC(.95);
		ps15->SetY2NDC(.78);
    PS16->Draw("sames");
    PS16->SetNormFactor(1.0);
    PS16->SetLineColor(4);
    func->SetLineColor(4);
    PS16->Fit("func","L");
        gPad->Update();
		TPaveStats *ps16 = (TPaveStats*)PS16 ->FindObject("stats");
		ps16->SetTextColor(4);
		ps16->SetX1NDC(.78);
		ps16->SetX2NDC(.98);
		ps16->SetY1NDC(.6);
		ps16->SetY2NDC(.76);

    	TLegend* compps=new TLegend(0.5,0.7,0.7,0.88);
    	compps->AddEntry("PS15","2015 Data","L");
    	compps->AddEntry("PS16","2016 Data","L");
    	compps->Draw();

    TCanvas *ComparePUC = new TCanvas ("ComparePUC","ComparePUC");
    PUC15->Draw();
    PUC15->GetXaxis()->SetTitle("PUC [GeV]");
    PUC15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    PUC15->SetNormFactor(1.0);
    PUC15->SetLineColor(2);
    PUC15->Fit("func","L");
        gPad->Update();
        TPaveStats *puc15 = (TPaveStats*)PUC15 ->FindObject("stats");
        puc15->SetTextColor(2);
		puc15->SetX1NDC(.78);
		puc15->SetX2NDC(.98);
		puc15->SetY1NDC(.95);
		puc15->SetY2NDC(.78);
		puc15->SetTextColor(2);
    PUC16->Draw("sames");
    PUC16->SetNormFactor(1.0);
    PUC16->SetLineColor(4);
    func->SetLineColor(4);
    PUC16->Fit("func","L");
        gPad->Update();
		TPaveStats *puc16 = (TPaveStats*)PUC16 ->FindObject("stats");
		puc16->SetTextColor(4);
		puc16->SetX1NDC(.78);
		puc16->SetX2NDC(.98);
		puc16->SetY1NDC(.6);
		puc16->SetY2NDC(.76);

    	TLegend* comppuc=new TLegend(0.5,0.7,0.7,0.88);
    	comppuc->AddEntry("PUC15","2015 Data","L");
    	comppuc->AddEntry("PUC16","2016 Data","L");
    	comppuc->Draw();

    TCanvas *CompareMHT = new TCanvas ("CompareMHT","CompareMHT");
    MHT15->Draw();
    MHT15->GetXaxis()->SetTitle("MHT [GeV]");
    MHT15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    MHT15->SetNormFactor(1.0);
    MHT15->SetLineColor(2);
    MHT15->Fit("func","L");
        gPad->Update();
		TPaveStats *mht15 = (TPaveStats*)MHT15 ->FindObject("stats");
		mht15->SetTextColor(2);
		mht15->SetX1NDC(.78);
		mht15->SetX2NDC(.98);
		mht15->SetY1NDC(.95);
		mht15->SetY2NDC(.78);
    MHT16->Draw("sames");
    MHT16->SetNormFactor(1.0);
    MHT16->SetLineColor(4);
    func->SetLineColor(4);
    MHT16->Fit("func","L");
        gPad->Update();
		TPaveStats *mht16 = (TPaveStats*)MHT16 ->FindObject("stats");
		mht16->SetTextColor(4);
		mht16->SetX1NDC(.78);
		mht16->SetX2NDC(.98);
		mht16->SetY1NDC(.6);
		mht16->SetY2NDC(.76);

    	TLegend* compmht=new TLegend(0.5,0.7,0.7,0.88);
    	compmht->AddEntry("MHT15","2015 Data","L");
    	compmht->AddEntry("MHT16","2016 Data","L");
    	compmht->Draw();

   	TCanvas *CompareOFFRECAL = new TCanvas ("CompareOFFRECAL","CompareOFFRECAL");
    OFFRECAL15->Draw();
    OFFRECAL15->GetXaxis()->SetTitle("OFFRECAL [GeV]");
    OFFRECAL15->GetYaxis()->SetTitle("Fraction of Events");
    func->SetLineColor(2);
    OFFRECAL15->SetNormFactor(1.0);
    OFFRECAL15->SetLineColor(2);
    OFFRECAL15->Fit("func","L");
        gPad->Update();
		TPaveStats *offrecal15 = (TPaveStats*)OFFRECAL15 ->FindObject("stats");
		offrecal15->SetTextColor(2);
		offrecal15->SetX1NDC(.78);
		offrecal15->SetX2NDC(.98);
		offrecal15->SetY1NDC(.95);
		offrecal15->SetY2NDC(.78);
    OFFRECAL16->Draw("sames");
    OFFRECAL16->SetNormFactor(1.0);
    OFFRECAL16->SetLineColor(4);
    func->SetLineColor(4);
    OFFRECAL16->Fit("func","L");
        gPad->Update();
		TPaveStats *offrecal16 = (TPaveStats*)OFFRECAL16 ->FindObject("stats");
		offrecal16->SetTextColor(4);
		offrecal16->SetX1NDC(.78);
		offrecal16->SetX2NDC(.98);
		offrecal16->SetY1NDC(.6);
		offrecal16->SetY2NDC(.76);

    	TLegend* compoffrecal=new TLegend(0.5,0.7,0.7,0.88);
    	compoffrecal->AddEntry("OFFRECAL15","2015 Data","L");
    	compoffrecal->AddEntry("OFFRECAL16","2016 Data","L");
    	compoffrecal->Draw();


}