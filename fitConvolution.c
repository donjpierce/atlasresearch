#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"


Double_t rayleigh (Double_t varMET, Double_t *parm) {
    /*
        Arguments
        _________

        varMET      :   function varibale   :   varMET
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
    //     sigma = parm[0] + 4.0 * parm[1];
    // }
    rayleigh = varMET / (sigma * sigma) * exp(-0.5 * varMET * varMET / (sigma * sigma));

    return rayleigh;
}


Double_t PWGD (Double_t varSUMET, Double_t *parm) {
    /*
        Arguments
        _________

        varSUMET  : function variable :   sumet integration varibale
        parm[0]   : constant          :   gamma
        parm[1]   : constant          :   mu
        parm[2]   : constant          :   offset

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
    // if (varSUMET > parm[2]) {
    //     pwgd = parm[0] * exp(-parm[0] * (varSUMET - parm[2]));
    // }
    // else {
    //     pwgd = 0;
    // }

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
        parm[4] :   Double_t    : slope
        parm[5] :   Double_t    : intercept

    */

    // gamma for 2011
    // Double_t gamma = 0.042;
    // gamma for 2017
    Double_t gamma = 0.065;
    Double_t mu = parm[3];
    Double_t offset = 92 - 7.4 * parm[3];
    Double_t cell17_slope = parm[4];
    Double_t cell17_intercept = parm[5];

    int n = parm[0];  // number of subintervals for integration
    double a, b, h;   // limits of integration (a, b) and trapezoidal width h
    a = parm[1];
    b = parm[2];
    h = (b - a) / n;  // get width of the subintervals

    Double_t pwgdParams[3];
    pwgdParams[0] = gamma;
    pwgdParams[1] = mu;
    pwgdParams[2] = offset;
    Double_t rayleighParams[3];
    rayleighParams[0] = cell17_slope;
    rayleighParams[1] = cell17_intercept;

    double SUMET[n+1], R[n+1];  // integration variable = SUMET, result = Rayleigh
    for (int i = 0; i < n; i ++) {
        SUMET[i] = a + i * h;
        rayleighParams[2] = SUMET[i];
        R[i] = PWGD(SUMET[i], pwgdParams) * rayleigh(MET[0], rayleighParams);
        R[i] /= (1 - exp(-mu)); // corrects for not starting the Poisson sum at 0
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

    // TF1 *rayleighFit = new TF1("rayleighFit", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
    // rayleighFit->SetParameters(1., 1.);
    // rayleighFit->SetParLimits(0, 0.1, 10000000.);
    // rayleighFit->SetParLimits(1, 0.1, 10000000.);
    // rayleighFit->SetParNames("amplitude", "sigma");
    //
    // // Defining linear fit function
    // TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
    // linfit->SetParameters(0, -80.);
    // linfit->SetParameters(1, -80.);
    // linfit->SetParLimits(0, -80., 80.);
    // linfit->SetParLimits(1, -80., 80.);
    // linfit->SetParNames("slope", "intercept");
    //
    // TFile *jburr17 = TFile::Open("../jburr_data_2017.root");
    // TTree *tree17 = (TTree*)jburr17->Get("METTree");
    //
    // TCanvas *myCanv = new TCanvas("MyCanv", "");
    // TH2F *cell17 = new TH2F("cell17", "", 100, 0., 100., 100, 0., 100.);
    // // tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17", "HLT_noalg_zb_L1ZB.passed>0.1");
    // tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17");
    // cell17->FitSlicesY(rayleighFit, 0, -1, 10, "L");
    // TH1D *cell17_sigma = (TH1D*)gDirectory->Get("cell17_1");
    // cell17_sigma->Fit(linfit, "M", "L", 0., 200.);
    // cell17_sigma->Draw();
    // Double_t cell17_slope, cell17_intercept;
    // cell17_slope = linfit->GetParameter(0);
    // cell17_intercept = linfit->GetParameter(1);
    Double_t cell17_slope = 0.17744; // 2017 cell
    Double_t cell17_intercept = 22.3986; // 2017 cell

    Int_t n_subint = 700;
    Double_t upper_bound = 1000.0;
    Int_t n_curves = 12;
    TCanvas *dists = new TCanvas("dists", "");
    TLegend* legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    TF1 *mu[n_curves];
    char *funcName = new char[5];
    for (Int_t i = 0; i < n_curves; i++) {
        int color = i + 1;
        int muValue = (i + 1) * 5;
        sprintf(funcName, "mu%d", muValue);
        mu[i] = new TF1(funcName, integration, 0, 100, 6);
        mu[i]->SetParNames("number of subintervals", "lower bound", "upper bound", "mu", "slope", "intercept");
        mu[i]->SetParameters(n_subint, 0.0, upper_bound, muValue, cell17_slope, cell17_intercept);
        mu[i]->SetLineColor(color);

        char *legendEntryName = new char[7];
        sprintf(legendEntryName, "#mu = %d", muValue);
        legend->AddEntry(mu[i], legendEntryName);
    }

    int normal = 1;
    if (normal == 1) {
        // draw lrves in normal order
        for (Int_t k = 0; k < n_curves; k++) {
             if (k == 0) { mu[k]->Draw(); }
             else { mu[k]->Draw("sames"); }
        }
    }
    else
    {
        // draw curves in reverse order
        for (Int_t j = 0; j < n_curves; j++) {
            if (j == 0) { mu[n_curves - 1]->Draw(); }
            else { mu[n_curves - j - 1]->Draw("sames"); }
        }
    }

    legend->Draw();

    // ROOT::Math::WrappedTF1 wf1(*mu[1]);
    // ROOT::Math::GaussIntegrator ig;
    // ig.SetFunction(wf1);
    // ig.SetRelTolerance(0.001);
    // cout << ig.Integral(0, 1000.0) << endl;


    // code for just plotting one curve
    // TF1 *mu7 = new TF1("mu7", integration, 0, 100, 6);
    // mu7->SetParNames("number of subintervals", "lower bound", "upper bound", "mu", "slope", "intercept");
    // mu7->SetParameters(n_subint, 0.0, upper_bound, 7.0, cell17_slope, cell17_intercept);
    // mu7->SetLineColor(4); // blue
    // mu7->Draw();
    // TLegend* legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    // legend->AddEntry(mu7, "#mu = 7");
    // legend->Draw();

 }
