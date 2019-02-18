#include "sumet_func_fft.h"

void sumet_func_compilation() {
    // create array for mu values
    int muvals[12];
    for (int i = 0; i < 12; i++) {
        muvals[i] = i * 5;
    }

    TF1 *fft[60];
    for (Int_t i = 0; i < 60; i++) {
        fft[i] = new TF1("first_fft", sumet_func_fft, 0, 2000, 6);
        fft[i]->SetParNames("Norm","Slope","Mu","Offset","Frac1","Slope2");
        fft[i]->FixParameter(0, 5.0);
        fft[i]->FixParameter(1, 0.065);
        fft[i]->FixParameter(2, muvals[i]);
        fft[i]->FixParameter(3, 50.0 - 6.0 * muvals[i]);
        fft[i]->FixParameter(4, 0.995);
        fft[i]->FixParameter(5, 0.0085);
    }

    return fft;
}
