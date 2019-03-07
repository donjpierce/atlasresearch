#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TMath.h>

void fft_test (Double_t *x) {

    TH1D *funchist = 0; 
    funchist = new TH1D("funchist", "funchist", TMath::Sin(), lowval, hival);

}
