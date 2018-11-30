#include "TF1.h"

Double_t rayleigh (Double_t varMET, Double_t *parm) {
    /*
        Arguments
        _________

        varMET   :   function varibale   :   varMET
        parm[0]     :   constant            :   slope from alorithm
        parm[1]     :   constant            :   intercept from algorithm
        parm[2]     :   constant            :   SUMET

        Returns
        _______

        rayleigh :  double  :   Rayleigh distribution result for given args
    */

    double rayleigh, sigma;
    // begin Rayleigh part
    sigma = parm[0] * sqrt(parm[2]) + parm[1];
    // if (parm[2] < 16.0) {
        // sigma = parm[0] + 4.0 * parm[1];
    // }
    rayleigh = varMET / (sigma * sigma) * exp(-0.5 * varMET * varMET / (sigma * sigma));

    return rayleigh;
}


Double_t PWGD (Double_t varSUMET, Double_t *parm) {
    /*
        Arguments
        _________

        varSUMET  : function variable :   sumet integration varibale
        parm[0] : constant          :   gamma
        parm[1] : constant          :   mu

        Returns
        _______

        pwgd   : double : Poisson-Weighted Gamma Densities result for given args
    */

    long double sum;
    double buffer, pwgd;
    sum = 1.;
    buffer = 1.;
    for (Int_t n = 2; n < 50 * int(parm[1]); n++)
    {
        buffer *= parm[0] * parm[1] * varSUMET / (n * (n - 1));
        sum += buffer;
    }

    // Poisson-weighted gamma densities function result
    pwgd = double(sum) * parm[1] * parm[0] * exp(-(parm[0] * varSUMET + parm[1]));

    return pwgd;
}

Double_t integration(Double_t *MET, Double_t *parm) {
    /*
        Arguments
        _________
        var[0]  :   variable    : MET
        parm[0] :   Double_t    : number of subintervals for integration
        parm[1] :   Double_t    : lower bound of integration
        parm[2] :   Double_t    : upper bound of integration
        parm[3] :   Double_t    : mu

    */
    Double_t cell15_slope = 0.5;
    Double_t cell15_intercept = -0.6;
    Double_t gamma = 0.042;
    Double_t mu = parm[3];

    int n = parm[0];  // number of subintervals for integration
    double a, b, h;   // limits of integration (a, b) and trapezoidal width h
    a = parm[1];
    b = parm[2];
    h = (b - a) / n;  // get width of the subintervals

    Double_t pwgdParams[2];
    pwgdParams[0] = gamma;
    pwgdParams[1] = mu;
    Double_t rayleighParams[3];
    rayleighParams[0] = cell15_slope;
    rayleighParams[1] = cell15_intercept;

    double SUMET[n+1], R[n+1];  // integration variable = SUMET, result = Rayleigh
    for (int i = 0; i < n; i ++) {
        SUMET[i] = a + i * h;
        rayleighParams[2] = SUMET[i];
        R[i] = PWGD(SUMET[i], pwgdParams) * rayleigh(MET[0], rayleighParams);
        // R[i] /= (1 - exp(-mu)); // corrects for not starting the Poisson sum at 0
    }
    double sum = 0;
    for (int j = 1; j < n; j++) {
        sum += h * R[j];
    }
    double integral;
    integral = h / 2.0 * (R[0] + R[n]) + sum;

    return integral;
}

void fitConvolution() {

    Double_t cell15_slope = 0.5;
    Double_t cell15_intercept = -0.6;

    Int_t n_subint = 2000;
    Double_t upper_bound = 1000.0;

    TF1 *mu1 = new TF1("mu1", integration, 0, 40, 4);
    mu1->SetParNames("number of subintervals", "lower bound", "upper bound", "mu");
    mu1->SetParameters(n_subint, 0.0, upper_bound, 1.0);
    mu1->SetLineColor(3); // green
    mu1->Draw();

    TF1 *mu4 = new TF1("mu4", integration, 0, 40, 4);
    mu4->SetParNames("number of subintervals", "lower bound", "upper bound", "mu");
    mu4->SetParameters(n_subint, 0.0, upper_bound, 4.0);
    mu4->SetLineColor(1); // black
    mu4->Draw("sames");

    TF1 *mu5 = new TF1("mu5", integration, 0, 40, 4);
    mu5->SetParNames("number of subintervals", "lower bound", "upper bound", "mu");
    mu5->SetParameters(n_subint, 0.0, upper_bound, 5.0);
    mu5->SetLineColor(2); // red
    mu5->Draw("sames");

    TF1 *mu7 = new TF1("mu7", integration, 0, 40, 4);
    mu7->SetParNames("number of subintervals", "lower bound", "upper bound", "mu");
    mu7->SetParameters(n_subint, 0.0, upper_bound, 7.0);
    mu7->SetLineColor(4); // blue
    mu7->Draw("sames");

    mu1->SetTitle("Marginal p.d.f. for MET for different values of #mu");
    // mu1->GetXaxis()->SetTitle("MET");
    // mu1->GetYaxis()->SetTitle("Probability");

 }
