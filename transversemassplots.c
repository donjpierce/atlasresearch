{
  TFile *muonFile = TFile::Open("PhysicsMain.L1KFmuontriggers.2016.f731f758_m1659m1710.Run309759.48Runs-002.root");
  TTree* muonTree = (TTree*)muonFile->Get("tree");

  Float_t reftransversemass, offrecaltransversemass;
  Float_t metref, metrefw, mexref, mexrefw, meyref, meyrefw;
  muonTree->SetBranchAddress("metrefmuon", &metref);
  muonTree->SetBranchAddress("metrefwmuon", &metrefw);
  muonTree->SetBranchAddress("mexrefmuon", &mexref);
  muonTree->SetBranchAddress("mexrefwmuon", &mexrefw);
  muonTree->SetBranchAddress("meyrefmuon", &meyref);
  muonTree->SetBranchAddress("meyrefwmuon", &meyrefw);

  Float_t metoff, mexoff, meyoff, metoffmu, mexoffmu, meyoffmu;
  muonTree->SetBranchAddress("metoffrecal", &metoff);
  muonTree->SetBranchAddress("mexoffrecal", &mexoff);
  muonTree->SetBranchAddress("meyoffrecal", &meyoff);
  muonTree->SetBranchAddress("metoffrecalmuon", &metoffmu);
  muonTree->SetBranchAddress("mexoffrecalmuon", &mexoffmu);
  muonTree->SetBranchAddress("meyoffrecalmuon", &meyoffmu);

  Int_t passmu, passmuvar;
  muonTree->SetBranchAddress("passmu26med", &passmu);
  muonTree->SetBranchAddress("passmu26varmed", &passmuvar);

  Float_t l1;
  muonTree->SetBranchAddress("metl1", &l1);

  TH1F *refmtplot = new TH1F("refmtplot", "", 200, 0., 200.);
  TH1F *offrecalmtplot = new TH1F("offrecalmtplot", "", 200, 0., 200.);

  Long64_t nentries = muonTree->GetEntries();
  for(int i = 0; i < nentries; i++)
  {
    muonTree->GetEntry(i);

    reftransversemass = sqrt(2*metref*metrefw*(1+((mexref*mexrefw+meyref*meyrefw) / (metref*metrefw))));
    offrecaltransversemass = sqrt(2*metoff*metoffmu*(1+((mexoff*mexoffmu+meyoff*meyoffmu) / (metoff*metoffmu))));

      if ((passmu > 0.1 || passmuvar > 0.1) && l1 > 50. && 40. < offrecaltransversemass && offrecaltransversemass < 100.)
      {
        refmtplot->Fill(reftransversemass);
      }
      if ((passmu > 0.1 || passmuvar > 0.1) && l1 > 50. && 20. < reftransversemass && reftransversemass < 100.)
      {
        offrecalmtplot->Fill(offrecaltransversemass);
      }

  }

  refmtplot->SetTitle("Transverse Mass based on METREF (L1 > 50GeV, 40 < OFFRECAL mass < 100)");
  refmtplot->GetYaxis()->SetTitle("No. of Events");
  refmtplot->GetXaxis()->SetTitle("Mass [GeV]");
  TCanvas *reftransmass = new TCanvas ("reftransmass", "");
  refmtplot->Draw();

  offrecalmtplot->SetTitle("Transverse Mass based on METOFFRECAL (L1 > 50GeV, 20 < REF mass < 100)");
  offrecalmtplot->GetYaxis()->SetTitle("No. of Events");
  offrecalmtplot->GetXaxis()->SetTitle("Mass [GeV]");
  TCanvas *offrecaltransmass = new TCanvas("offrecaltransmass", "");
  offrecalmtplot->Draw();

}
