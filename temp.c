{

	Int_t k = 0;
	for (int q = 0; q < 4; q++)
	{
		for (int l = 0; l < 4; l++)
		{
			if (q == l)	continue;
			correlationgraph[k]->SetTitle("2016 Muons (L1KFmuontriggers...48Runs-002) for L1 > 50GeV");
			k++;
		}
	}


}

	/*

	TFile *muonFile = TFile::Open("../PhysicsMain.L1KFmuontriggers.2016.f731f758_m1659m1710.Run309759.48Runs-002.root");
	TTree* muonTree = (TTree*)muonFile->Get("tree");
	Int_t passmuonmed, passmuonvarmed, muonrecal, muonclean;
	Float_t ml1, muonint;
	muonTree->SetBranchAddress("passmu26med", &passmuonmed);
	muonTree->SetBranchAddress("passmu26varmed", &passmuonvarmed);
	muonTree->SetBranchAddress("metl1", &ml1);
	muonTree->SetBranchAddress("actint", &muonint);
	muonTree->SetBranchAddress("recalbroke", &muonrecal);
	muonTree->SetBranchAddress("passcleancuts", &muonclean);

	TFile *zbFile = TFile::Open("../PhysicsMain.L1KFnoalgXEtriggers.2016.f731f758_m1659m1710.Run309759.48Runs-001.root");
	TTree *zbTree = (TTree*)zbFile->Get("tree");
	Int_t zbl1gt10, zbl1gt30, zbl1gt40, zbl1gt45;
	Float_t zbl1, zbint;
	zbTree->SetBranchAddress("passnoalgL1XE10", &zbl1gt10);
	zbTree->SetBranchAddress("passnoalgL1XE30", &zbl1gt30);
	zbTree->SetBranchAddress("passnoalgL1XE40", &zbl1gt40);
	zbTree->SetBranchAddress("passnoalgL1XE45", &zbl1gt45);
	zbTree->SetBranchAddress("actint", &zbint);

	//___Calculate Tail Events Based on Resolutions___
	TString metalgName[4] = {"metl1", "metl1kf", "metl1mht", "metl1mhtkf"};
	TString setalgName[4] = {"setl1", "setl1kf", "setl1mht", "setl1mhtkf"};

	// create arrays for MET and SET branches
	Float_t met[4]; Float_t set[4];
	for (Int_t i = 0; i < 4; i++)
	{
		zbTree->SetBranchAddress(metalgName[i], &met[i]);
		zbTree->SetBranchAddress(setalgName[i], &set[i]);
	}

		TH2F *correlationgraph = new TH2F("correlationgraph", "",100, 0., 100., 100, 0., 100.);

		TTree* runtree = zbTree;

		Int_t n = 0;
		Long64_t nentries = runtree->GetEntries();
		for (Int_t i = 0; i < nentries/10; i++)
		{
			runtree->GetEntry(i);
			//std::cout << zbl1gt10 << zbl1gt30 << zbl1gt40 << zbl1gt45  << std::endl;

			Double_t x[4]; // x = bulkmet and y = tailmet will be calculated for each algorithm
			Double_t y[4];

			if ( (zbl1gt10 > 0.1 || zbl1gt30 > 0.1 || zbl1gt40 > 0.1 || zbl1gt45 > 0.1) && met[0] > 50.)
			{
				x[0] = met[0];
				y[0] = met[1];
				correlationgraph->Fill(x[0], y[0]); // and populate the appropraite correlationgraph
			}
		}

}
