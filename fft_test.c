#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TMath.h>

void fft_test () {

    TCanvas *canv = new TCanvas("canv", "FFT");
    TH1::AddDirectory(kFALSE);

    TF1 *fsin = new TF1("fsin", "sin(x) + sin(2*x) + sin(0.5*x) + 1", 0, 4*TMath::Pi());
    fsin->Draw();

    Int_t n = 25;
    TH1D *hsin = 0;
    // funchist = new TH1D("hsin", "hsin", TMath::Sin(), lowval, hival);
    hsin = new TH1D("hsin", "hsin", n + 1, 0, 4*TMath::Pi());
    Double_t x;

    for (Int_t i = 0; i <= n; i++) {
        x = (Double_t(i) / n) * (4*TMath::Pi());
        hsin->SetBinContent(i+1, fsin->Eval(x));
    }

    hsin->Draw("same");

    TCanvas *fft_canv = new TCanvas("fft_canv", "fft canvas");
    TH1 *hm = 0;
    TVirtualFFT::SetTransform(0);
    hm = hsin->FFT(hm, "MAG");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->Draw();


}
