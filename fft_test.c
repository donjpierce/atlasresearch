#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TMath.h>

TH1* histogram () {

    TCanvas *canv = new TCanvas("canv", "FFT");
    TH1::AddDirectory(kFALSE);

    TF1 *func = new TF1("func", "cos(2 * TMath::Pi() * (3 * x)) * TMath::Exp(-TMath::Pi() * x * x)", -2, 2);
    func->Draw();

    Int_t n = 25;
    TH1 *hist = 0;
    // funchist = new TH1D("hsin", "hsin", TMath::Sin(), lowval, hival);
    hist = new TH1D("hist", "hist", n + 1, -2, 2);
    Double_t x;

    for (Int_t i = 0; i <= n; i++) {
        x = (Double_t(i) / n) * (4*TMath::Pi());
        hist->SetBinContent(i+1, func->Eval(x));
    }

    return hist;
}

TH1D* hist (Double_t *parm) {
    static Double_t parmsave[6] = {-999.,-999.,-999.,-999.,-999.,-999.};
    static int ncalls=0;
    int n = 8847360;
    int nparms = 6;
    Double_t lowval = 0.0;
    Double_t hival = 50000.0;
    Double_t xtemp;
    int itemp;

    // long double functot[8847360];
    // static Double_t finalfuncval[8847360];
    bool sameparm;

    TH1D *funchist = 0;
    TH1 *ftransform = 0;

    //Calculate and fill histogram for one interaction
    funchist = new TH1D("funchist", "funchist", n+1, lowval, hival);
    for (Int_t i=0; i<=n; i++){
      xtemp = hival*double(i)/double(n);
      funchist->SetBinContent(i+1,
          parm[4]*parm[1]*exp(-parm[1]*xtemp)+(1.-parm[4])*parm[5]*exp(-parm[5]*xtemp));
     }

    return funchist;

}

TH1* fourier (TH1D* hist) {

    TCanvas *fft_canv = new TCanvas("fft_canv", "fft canvas");
    TH1 *hm = 0;
    TVirtualFFT::SetTransform(0);
    hm = hist->FFT(hm, "MAG");
    // hm->SetTitle("Magnitude of the 1st transform");
    hm->Draw();

    return hm;
}
