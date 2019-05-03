#include <TROOT.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TVirtualFFT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include <iostream>

Double_t frechet (Double_t *x, Double_t *parm) {
    /*
        Arguments
        _________
        x       :  function variable    : MET
        parm[0] :  constant             : norm (function of SUMET)
        parm[1] :  constant             : alpha
        parm[2] :  constant             : s
        parm[3] :  constant             : m

        Returns
        _________
        value   :  Double_t             : value of the Frechet PDF
    */

    double argument = (x[0] - parm[3]) / parm[2];
    double value = parm[0] * (parm[1] / parm[2]) * pow(argument, -1.0 - parm[1]) \
                   * exp(-pow(argument, -parm[1]));
    return value;
}

Double_t rayleigh (Double_t *varMET, Double_t *parm) {
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
    rayleigh = varMET[0] / (sigma * sigma) * exp(-0.5 * varMET[0] * varMET[0] / (sigma * sigma));

    return rayleigh;
}

Double_t PWGD (Double_t *varSUMET, Double_t *parm) {
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
    long double sum, pwgd;

    if (varSUMET[0] < parm[2]) {
        pwgd = 0;
        return pwgd;
    }

    double buffer;
    sum = 1.;
    buffer = 1.;
    for (Int_t n = 2; n < 50 * int(parm[1]); n++)
    {
        buffer *= parm[0] * parm[1] * varSUMET[0] / (n * (n - 1));
        sum += buffer;
    }

    // Poisson-weighted gamma densities function result
    pwgd = double(sum) * parm[1] * parm[0] * exp(-(parm[0] * varSUMET[0] + parm[1]));

    // if (varSUMET > parm[2]) {
    //     pwgd = parm[0] * exp(-parm[0] * (varSUMET - parm[2]));
    // }
    // else {
    //     pwgd = 0;
    // }

    return pwgd;
}

Double_t sumet_func_fft (Double_t *x, Double_t *parm) {
    static Double_t parmsave[6] = {-999.,-999.,-999.,-999.,-999.,-999.};
    static int ncalls=0;
    int n = 8847360;
    int nparms = 6;
    Double_t lowval = 0.0;
    Double_t hival = 50000.0;
    Double_t xtemp;
    int itemp;

    long double functot[8847360];
    static Double_t finalfuncval[8847360];
    bool sameparm;

    TH1D *funchist = 0;
    TH1 *ftransform = 0;

    // ncalls++;
    // if(ncalls%100==0) {
    //   std::cout << "Number of routine calls = " << ncalls << "\n";
    // }

    if(x[0] - parm[3]<0.) {
      return 0.0;
    }

    sameparm=true;
    for(int k=0; k<nparms; k++) {
      if (parm[k] != parmsave[k]) {
        sameparm=false;
        parmsave[k]=parm[k];
      }
    }

    if(!sameparm) {

      //Recompute function
      for (int i=0; i<n; i++) {
        functot[i]=0.000000000000000000000000000;
      }

    //Calculate and fill histogram for one interaction
    funchist = new TH1D("funchist", "funchist", n+1, lowval, hival);
    for (Int_t i=0; i<=n; i++){
      xtemp = hival*double(i)/double(n);
      funchist->SetBinContent(i+1,parm[4]*parm[1]*exp(-parm[1]*xtemp)+(1.-parm[4])*parm[5]*exp(-parm[5]*xtemp));
     }

    //Compute and store fft for one interaction function
    TVirtualFFT::SetTransform(0);
    ftransform = funchist->FFT(ftransform, "MAG");
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();

    //Compute fft for full mu dependent function
    Double_t *re_full = new Double_t[n];
    Double_t *im_full = new Double_t[n];
    fft->GetPointsComplex(re_full,im_full);

    for(int i=0; i<n; i++) {
      Double_t re_full_save = re_full[i];
      Double_t im_full_save = im_full[i];
      re_full[i] = exp(parm[2]*(re_full_save)*hival/double(n))*cos(parm[2]*im_full_save*hival/double(n))-1.0;
      im_full[i] = exp(parm[2]*(re_full_save)*hival/double(n))*sin(parm[2]*im_full_save*hival/double(n));
    }

    //get the back transform of full function
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
    fft_back->SetPointsComplex(re_full,im_full);
    fft_back->Transform();
    TH1 *hb = 0;
    hb = TH1::TransformHisto(fft_back,hb,"Re");
    for (Int_t i=0; i<n; i++){
      finalfuncval[i] = exp(-parm[2])*(hb->GetBinContent(i+1))/hival;
      //following for extra tail not in the function itself
      //finalfuncval[i] = parm[4]*exp(-parm[2])*(hb->GetBinContent(i+1))/hival +(1.0-parm[4])*exp(-parm[5]*hival*double(i)/double(n));
    }

    //Clean up memory
    delete hb;
    hb=0;
    delete ftransform;
    ftransform = 0;
    delete funchist;
    funchist=0;
    delete fft;
    fft=0;
    delete fft_back;
    fft_back=0;
    delete [] re_full;
    delete [] im_full;

  } // finished filling the function


  itemp=int((x[0] - parm[3])*double(n)/hival);
  if(itemp>n-1 || itemp<0){
      return 0.0;
  }
  else {
      //std::cout << " itemp " << itemp << " x = " << x[0] << "  f(x) " << finalfuncval[itemp] <<"\n";
      return parm[0] * finalfuncval[itemp];
  }

}

Double_t integration(Double_t *MET, Double_t *parm) {
    /*
        Arguments
        _________
        MET     :   variable    : MET
        parm[0] :   Double_t    : number of subintervals for integration
        parm[1] :   Double_t    : lower bound of integration
        parm[2] :   Double_t    : upper bound of integration
        parm[3] :   Double_t    : mu
        parm[4] :   Double_t    : slope
        parm[5] :   Double_t    : intercept
        parm[6] :   bool        : TRUE for FFT, FALSE for PWGD
    */

    // gamma for 2011
    // Double_t gamma = 0.042;
    // gamma for 2017
    Double_t gamma = 0.070;
    Double_t mu = parm[3];
    Double_t offset = 25.0 - 5.0 * parm[3];
    Double_t cell17_slope = parm[4];
    Double_t cell17_intercept = parm[5];

    int n = parm[0];  // number of subintervals for integration
    double a, b, h;   // limits of integration (a, b) and trapezoidal width h
    a = parm[1];
    b = parm[2];
    h = (b - a) / n;  // get width of the subintervals

    // parameters for FFT component
    Double_t sumetFFT_params[6];
    sumetFFT_params[0] = 1.0;  // Normalization constant
    sumetFFT_params[1] = gamma; // gamma
    sumetFFT_params[2] = mu;
    sumetFFT_params[3] = offset;
    sumetFFT_params[4] = 0.995; // fraction of first exponential
    sumetFFT_params[5] = 0.0090; // fraction of second exponential
    TF1 *fft = new TF1("fft", sumet_func_fft, 0, 2000, 6);
    fft->SetParameters(sumetFFT_params);

    // parameters for PWGD component
    Double_t pwgdParams[3];
    pwgdParams[0] = gamma;
    pwgdParams[1] = mu;
    pwgdParams[2] = offset;
    TF1 *pwgd = new TF1("pwgd", PWGD, 0, 2000, 3);
    pwgd->SetParameters(pwgdParams);

    // parameters for Rayleigh component
    Double_t rayleighParams[3];
    rayleighParams[0] = cell17_slope;
    rayleighParams[1] = cell17_intercept;
    TF1 *rayleigh_func = new TF1("rayleigh_func", rayleigh, 0, 500, 3);

    Double_t SUMET[n+1], R[n+1];  // SUMET = integration variable, R = result
    for (int i = 0; i < n; i ++) {
        SUMET[i] = a + i * h;

        // Perform convolution of Raylegih and PWDG / FFT
        // set the SUMET-dependent rayleigh parameters
        rayleighParams[2] = SUMET[i];
        rayleigh_func->SetParameters(rayleighParams);
        Double_t r1 = rayleigh_func->Eval(MET[0]);

        Double_t r2;
        if (parm[6] == true) { r2 = fft->Eval(SUMET[i]); }
        else { r2 = pwgd->Eval(SUMET[i]); }

        R[i] = r1 * r2;
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

Double_t linear_combination(Double_t *MET, Double_t *parm) {
  /*
    Arguments
    _________
    MET     :   Double_t    :   MET value
    parm[0] :   Double_t    :   integration 0 (number of subintervals)
    parm[1] :   Double_t    :   integration 1 (lower bound)
    parm[2] :   Double_t    :   integration 2 (upper bound)
    parm[3] :   Double_t    :   integration 3 (mu)
    parm[4] :   Double_t    :   integration 4 (resolution slope)
    parm[5] :   Double_t    :   integration 5 (resolution intercept)
    parm[6] :   bool        :   true for fft, false for PWGD
    parm[7] :   Double_t    :   Frechet 0 (norm)
    parm[8] :   Double_t    :   Frechet 1 (alpha)
    parm[9] :   Double_t    :   Frechet 2 (s)
    parm[10]:   Double_t    :   Frechet 3 (m)
    parm[11]:   Double_t    :   alpha (fitting parameter)

    Returns
    _________
    linsum  :   Double_t    :   result of the linear sum
  */

  Double_t muValue = parm[3];

  // parameters for the integration
  Double_t integrationParams[7] = {parm[0], parm[1], parm[2], parm[3],
                                   parm[4], parm[5], parm[6]};
  TF1 *integration_func = new TF1("integration_funcs", integration, 0, 1000, 7);
  integration_func->SetParameters(integrationParams);

  // parameters for the Frechet distribtuion
  Double_t frechetParams[4] = {parm[7], parm[8], parm[9], parm[10]};
  TF1 *frechet_func = new TF1("frechet_func", frechet, 0, 1000, 4);
  frechet_func->SetParameters(frechetParams);

  Double_t linsum, integration_result, frechet_result;
  integration_result = integration_func->Eval(MET[0]);
  frechet_result = frechet_func->Eval(MET[0]);

  // perform mu-dependent linear sum
  Double_t alpha = parm[11];
  Double_t int_coeff =  1 - (alpha * muValue);
  // Double_t int_coeff = 0;
  Double_t frechet_coeff = alpha * muValue;
  linsum = int_coeff * integration_result + frechet_coeff * frechet_result;
  // linsum = (1 - alpha * parm[3]) * integration_result + alpha * parm[3] * frechet_result;

  return linsum;
}

void fitConvolution() {
  /*
    Description:

    This function constructions the MET distribution from data and fits the
    data to the theoretical MET curve.
  */

  TF1 *func = new TF1("func", linear_combination, 0.0, 300.0, 12);
  /*
  Descrption of func parameters.
  _________
  MET     :   Double_t    :   MET value
  parm[0] :   Double_t    :   integration 0 (number of subintervals)
  parm[1] :   Double_t    :   integration 1 (lower bound)
  parm[2] :   Double_t    :   integration 2 (upper bound)
  parm[3] :   Double_t    :   integration 3 (mu)
  parm[4] :   Double_t    :   integration 4 (resolution slope)
  parm[5] :   Double_t    :   integration 5 (resolution intercept)
  parm[6] :   bool        :   true for fft, false for PWGD
  parm[7] :   Double_t    :   Frechet 0 (norm)
  parm[8] :   Double_t    :   Frechet 1 (alpha)
  parm[9] :   Double_t    :   Frechet 2 (s)
  parm[10]:   Double_t    :   Frechet 3 (m)
  parm[11]:   Double_t    :   alpha (fitting parameter)
  */
  func->FixParameter(0, 1700);
  func->FixParameter(1, 0.0);
  func->FixParameter(2, 2000.0);
  func->FixParameter(4, 0.465);
  func->FixParameter(5, 3.0);
  func->FixParameter(6, true);
  func->FixParameter(7, 1);
  func->FixParameter(8, 18.0);
  func->FixParameter(9, 100.0);
  func->FixParameter(10, -80.0);
  func->SetParameter(11, 0.001); // alpha
  func->SetParLimits(11, 0.0005, 0.01);

  Double_t integrationParams[7] = {};
  Double_t frechetParams[4] = {};

  TFile *file = new TFile("data/mu_analysis.root");

  TLegend *legend = new TLegend(0.37, 0.7, 0.55, 0.88);
  TObjArray* reconstructed_distributions = 0;
  file->GetObject("reconstructed_distributions",reconstructed_distributions);
  TH1F *reconcorrmuxx[7];
  TCanvas *canvMu[7];
  for (int i=0;i<7;i++) {
    reconcorrmuxx[i] = (TH1F*)reconstructed_distributions->At(i)->Clone();
    reconcorrmuxx[i]->SetNormFactor(0.2);
    reconcorrmuxx[i]->SetMarkerColor(6);
    reconcorrmuxx[i]->SetMarkerStyle(21);
    reconcorrmuxx[i]->SetMarkerSize(0.8);
    reconcorrmuxx[i]->SetLineColor(6);
    reconcorrmuxx[i]->GetXaxis()->SetTitle("Cell Algorithm E_{T}^{MISS} [GeV]");
    reconcorrmuxx[i]->GetYaxis()->SetTitle("Fraction of Events / GeV");
    reconcorrmuxx[i]->SetStats(kFALSE);
    reconcorrmuxx[i]->SetTitle("");

    int muValue = (i + 1) * 10 - 5;

    char *canvName = new char[10];
    sprintf(canvName, "canvMu_%i", muValue);

    canvMu[i] = new TCanvas(canvName, canvName);

    func->FixParameter(3, muValue);

    reconcorrmuxx[i]->Draw();
    reconcorrmuxx[i]->Fit(func, "R");
    // func->Draw("sames");
    canvMu[i]->SetLogy();

    char *entryName = new char[10];
    sprintf(entryName, "#mu = %i", muValue);

    char *filename = new char[10];
    sprintf(filename, "mu_%i.png", muValue);

    legend->AddEntry(reconcorrmuxx[i], entryName);

    char *plotTitle = new char[60];
    Double_t alpha = func->GetParameter(11);
    sprintf(plotTitle, "MET Curve Fit for #mu = %i and #alpha = %e", muValue, alpha);
    reconcorrmuxx[i]->SetTitle(plotTitle);

    canvMu[i]->SaveAs(filename);
  }

}

void plotMuCurves() {
    // Defining the Rayleigh function in main()
    TF1 *rayleighFit = new TF1("rayleighFit", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
    rayleighFit->SetParameters(1., 1.);
    rayleighFit->SetParLimits(0, 0.1, 10000000.);
    rayleighFit->SetParLimits(1, 0.1, 10000000.);
    rayleighFit->SetParNames("amplitude", "sigma");

    // Defining linear fit function
    TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
    linfit->SetParameters(-80., -80.);
    linfit->SetParLimits(0, -80., 80.);
    linfit->SetParLimits(1, -80., 80.);
    linfit->SetParNames("slope", "intercept");

    // TFile *jburr17 = TFile::Open("data/jburr_data_2017.root");
    TFile *jburr17 = TFile::Open("data/user.jburr.2017_11_17.data17.ZB.root");
    TTree *tree17 = (TTree*)jburr17->Get("METTree");

    TCanvas *myCanv = new TCanvas("MyCanv", "");
    TH2F *cell17 = new TH2F("cell17", "", 100, 0., 100., 100, 0., 100.);
    tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17", "HLT_noalg_zb_L1ZB.passed>0.1");
    // tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17");
    cell17->FitSlicesY(rayleighFit, 0, -1, 10, "L");
    TH1D *cell17_sigma = (TH1D*)gDirectory->Get("cell17_1");
    cell17_sigma->Fit(linfit, "M", "L", 0., 200.);
    cell17_sigma->Draw();
    Double_t cell17_slope, cell17_intercept;
    cell17_slope = linfit->GetParameter(0);
    cell17_intercept = linfit->GetParameter(1);

    // function parameters
    Int_t n_subint = 1700;
    Double_t upper_bound = 2000.0;
    Double_t lower_bound = 0.0;
    // Double_t slope = cell17_slope;
    Double_t slope = 0.465;
    // Double_t intercept = cell17_intercept;
    Double_t intercept = 3.0;
    Double_t fft = true;

    // set parametes for linear sum
    // NOTE: SKIPPED PARAMS ARE SET AT LOOP LEVEL
    // NOTE: linear_combinationParams[11] is a fitting parameter
    Double_t linear_combinationParams[12];
    linear_combinationParams[0] = n_subint;
    linear_combinationParams[1] = lower_bound;
    linear_combinationParams[2] = upper_bound;
    linear_combinationParams[4] = slope;
    linear_combinationParams[5] = intercept;
    linear_combinationParams[6] = fft;
    linear_combinationParams[8] = 18.0;
    linear_combinationParams[9] = 100.0;
    linear_combinationParams[10] = -80.0;

    Int_t n_curves = 12;
    TCanvas *dists = new TCanvas("dists", "");
    TLegend *legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    TF1 *mu[n_curves];
    char *funcName = new char[5];
    for (Int_t i = 0; i < n_curves; i++) {
        int color = i + 1;
        short int muValue = (i + 1) * 5;
        sprintf(funcName, "mu%i", muValue);

        // set mu parameter for linear sum
        linear_combinationParams[3] = muValue;
        // Frechet NORM set to 1 for almost all purposes
        linear_combinationParams[7] = 1;
        // initialize function which performs convolution
        mu[i] = new TF1("met", linear_combination, 0, 100, 12);
        mu[i]->SetParameters(linear_combinationParams);
        mu[i]->SetParNames("number of subintervals", "lower bound", "upper bound",
                           "mu", "slope", "intercept", "FFT", "frechet_norm",
                           "frechet_alpha", "frechet_s", "frechet_m");
        mu[i]->SetLineColor(color);

        char *legendEntryName = new char[10];
        sprintf(legendEntryName, "#mu = %i", muValue);
        legend->AddEntry(mu[i], legendEntryName);
    }

    int normal = 1;
    if (normal == 1) {
        // draw curves in normal order
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
    mu[0]->SetTitle("Convolution with FFT");
    mu[0]->GetXaxis()->SetTitle("MET [GeV]");
    mu[0]->GetYaxis()->SetTitle("Probability of an Event");
    dists->SetLogy();
    dists->SaveAs("results.png");

}

void plotIndividualSUMET() {
    // This method plots two SUMET plots:
    // (1) the SUMET using one-exponential (PWGD method)
    // (2) the SUMET using two-exponential (FFT method)

    // BEGIN FFT PART
    Int_t n_curves = 12;
    Double_t sumetFFT_params[6];
    sumetFFT_params[0] = 1.0;   // Normalization constant
    sumetFFT_params[1] = 0.070; // gamma
    sumetFFT_params[2] = 0.0;   // mu
    sumetFFT_params[3] = 0.0;   // offset
    sumetFFT_params[4] = 0.995;
    sumetFFT_params[5] = 0.0090;
    TCanvas *fft_dists = new TCanvas("dists", "");
    TLegend* legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    TF1 *fft_mu[n_curves];
    char *funcName = new char[5];
    for (Int_t i = 0; i < n_curves; i++) {
        int color = i + 1;
        short int muValue = (i + 1) * 5;
        double offset = 25.0 - 5.0 * muValue;

        sprintf(funcName, "mu%i", muValue);
        fft_mu[i] = new TF1(funcName, sumet_func_fft, 0, 1600, 6);

        sumetFFT_params[2] = muValue;
        sumetFFT_params[3] = 25.0 - 5.0 * muValue;

        fft_mu[i]->SetParameters(sumetFFT_params);
        fft_mu[i]->SetLineColor(color);

        char *fft_legendEntryName = new char[10];
        sprintf(fft_legendEntryName, "#mu = %i", muValue);
        legend->AddEntry(fft_mu[i], fft_legendEntryName);
    }

    int fft_normal = 1;  // change to 0 if drawings do not fit on canvas
    if (fft_normal == 1) {
        // draw lrves in normal order
        for (Int_t k = 0; k < n_curves; k++) {
             if (k == 0) { fft_mu[k]->Draw(); }
             else { fft_mu[k]->Draw("sames"); }
        }
    }
    else
    {
        // draw curves in reverse order
        for (Int_t j = 0; j < n_curves; j++) {
            if (j == 0) { fft_mu[n_curves - 1]->Draw(); }
            else { fft_mu[n_curves - j - 1]->Draw("sames"); }
        }
    }

    legend->Draw();
    fft_mu[0]->SetTitle("SumEt Distribution using Two-Exponential (FFT method)");
    fft_mu[0]->GetXaxis()->SetTitle("SumEt [GeV]");
    fft_mu[0]->GetYaxis()->SetTitle("Probability of an Event");
    fft_dists->SaveAs("fft_dists.png");


    // BEGIN PWGD PART
    Double_t pwgd_params[3];
    sumetFFT_params[0] = 0.070; // gamma
    sumetFFT_params[1] = 0.0;   // mu
    sumetFFT_params[2] = 0.0;   // offset
    TCanvas *pwgd_dists = new TCanvas("dists", "");
    TF1 *pwgd_mu[n_curves];
    for (Int_t i = 0; i < n_curves; i++) {
        int color = i + 1;
        short int muValue = (i + 1) * 5;
        double offset = 25.0 - 5.0 * muValue;

        sprintf(funcName, "mu%i", muValue);
        pwgd_mu[i] = new TF1(funcName, PWGD, 0, 1600, 3);

        sumetFFT_params[1] = muValue;
        sumetFFT_params[2] = 25.0 - 5.0 * muValue;

        pwgd_mu[i]->SetParameters(sumetFFT_params);
        pwgd_mu[i]->SetLineColor(color);

        char *pwgd_legendEntryName = new char[10];
        sprintf(pwgd_legendEntryName, "#mu = %i", muValue);
        legend->AddEntry(pwgd_mu[i], pwgd_legendEntryName);
    }

    int pwgd_normal = 1;  // change to 0 if drawings do not fit on canvas
    if (pwgd_normal == 1) {
        // draw lrves in normal order
        for (Int_t k = 0; k < n_curves; k++) {
             if (k == 0) { pwgd_mu[k]->Draw(); }
             else { pwgd_mu[k]->Draw("sames"); }
        }
    }
    else
    {
        // draw curves in reverse order
        for (Int_t j = 0; j < n_curves; j++) {
            if (j == 0) { pwgd_mu[n_curves - 1]->Draw(); }
            else { pwgd_mu[n_curves - j - 1]->Draw("sames"); }
        }
    }

    legend->Draw();
    pwgd_mu[0]->SetTitle("SumEt Distribution using Two-Exponential (PWGD method)");
    pwgd_mu[0]->GetXaxis()->SetTitle("SumEt [GeV]");
    pwgd_mu[0]->GetYaxis()->SetTitle("Probability of an Event");
    pwgd_dists->SaveAs("pwgd_dists.png");
}
