#include <TROOT.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TVirtualFFT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

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

    if (varSUMET < parm[2]) {
        pwgd = 0;
        return pwgd;
    }

    long double sum;
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
    //int n=62914560;
    // int n = 8847360;
    int n = 8847360;
    int nparms = 6;
    Double_t lowval = 0.0;
    Double_t hival = 50000.0;
    Double_t xtemp;
    int itemp;

    // long double functot[n];
    // static Double_t finalfuncval[62914560];
    // static Double_t finalfuncval[8847360];

    long double functot[8847360];
    //static Double_t finalfuncval[62914560];
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
    //onetermfunc->SetParameters(parm[0],parm[1],1.,0.);
    //delete funchist;
    //funchist=0;
    funchist = new TH1D("funchist", "funchist", n+1, lowval, hival);
    for (Int_t i=0; i<=n; i++){
      xtemp = hival*double(i)/double(n);
      //gamma*exp(-gamma*x)
      //funchist->SetBinContent(i+1,parm[1]*exp(-parm[1]*xtemp));
      //gamma^2*exp(-gamma*sqrt(x))
      //funchist->SetBinContent(i+1,parm[1]*parm[1]*exp(-parm[1]*sqrt(xtemp)));
      //(gamma-1)*x^-gamma
      //funchist->SetBinContent(i+1,(parm[1]-1)*pow(xtemp,-parm[1]));
      //alpha*gamma*exp(-gamma*x)+(1-alpha)*gamma2*exp(-gamma2*x)
      funchist->SetBinContent(i+1,parm[4]*parm[1]*exp(-parm[1]*xtemp)+(1.-parm[4])*parm[5]*exp(-parm[5]*xtemp));
     }
    //Compute and store fft for one interaction function
    //delete ftransform;
    //ftransform = 0;

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
      //re_full[i] = exp(parm[2]*(re_full_save)*hival/double(n))*cos(parm[2]*im_full_save*hival/double(n));
      re_full[i] = exp(parm[2]*(re_full_save)*hival/double(n))*cos(parm[2]*im_full_save*hival/double(n))-1.0;
      im_full[i] = exp(parm[2]*(re_full_save)*hival/double(n))*sin(parm[2]*im_full_save*hival/double(n));
      //re_full[i] = re_full_save;
      //im_full[i] = im_full_save;
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
    sumetFFT_params[4] = 0.995;
    sumetFFT_params[5] = 0.0090;
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
        rayleighParams[2] = SUMET[i];
        rayleigh_func->SetParameters(rayleighParams);
        Double_t r1 = rayleigh_func->Eval(MET[0]);
        // Double_t r2 = fft->Eval(SUMET[i]);
        Double_t r2 = pwgd->Eval(SUMET[i]);
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

void fitConvolution() {

    // *** JUST FOR TESTING ***
    // TCanvas *canv = new TCanvas("canv", "canv");
    // Double_t params[6] = {5.0, 0.067, 5.0, 20.0, 0.995, 0.0085};
    // TF1 *fft = new TF1("fft", sumet_func_fft, 0, 2000, 6);
    // fft->SetParameters(params);

    TF1 *rayleighFit = new TF1("rayleighFit", "[0]*(1/[1])*(x/[1])*exp(-.5*(x/[1])*(x/[1]))");
    rayleighFit->SetParameters(1., 1.);
    rayleighFit->SetParLimits(0, 0.1, 10000000.);
    rayleighFit->SetParLimits(1, 0.1, 10000000.);
    rayleighFit->SetParNames("amplitude", "sigma");

    // Defining linear fit function
    TF1 *linfit = new TF1("linfit", "[0]*x + [1]");
    linfit->SetParameters(0, -80.);
    linfit->SetParameters(1, -80.);
    linfit->SetParLimits(0, -80., 80.);
    linfit->SetParLimits(1, -80., 80.);
    linfit->SetParNames("slope", "intercept");

    TFile *jburr17 = TFile::Open("data/jburr_data_2017.root");
    TTree *tree17 = (TTree*)jburr17->Get("METTree");

    TCanvas *myCanv = new TCanvas("MyCanv", "");
    TH2F *cell17 = new TH2F("cell17", "", 100, 0., 100., 100, 0., 100.);
    // tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17", "HLT_noalg_zb_L1ZB.passed>0.1");
    tree17->Draw("cell.met:sqrt(cell.sumet)>>cell17");
    cell17->FitSlicesY(rayleighFit, 0, -1, 10, "L");
    TH1D *cell17_sigma = (TH1D*)gDirectory->Get("cell17_1");
    cell17_sigma->Fit(linfit, "M", "L", 0., 200.);
    cell17_sigma->Draw();
    Double_t cell17_slope, cell17_intercept;
    cell17_slope = linfit->GetParameter(0);
    cell17_intercept = linfit->GetParameter(1);

    Int_t n_subint = 1700;
    Double_t upper_bound = 2000.0;

    // TCanvas *dists = new TCanvas("dists", "");
    // TLegend* legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    // TF1 *mu5 = new TF1("mu5", integration, 0, 100, 6);
    // mu5->SetParameters(n_subint, 0.0, upper_bound, 5, cell17_slope, cell17_intercept);
    // legend->AddEntry(mu5, "#mu = 5");
    // mu5->Draw();
    // legend->Draw();

    Int_t n_curves = 12;
    TCanvas *dists = new TCanvas("dists", "");
    TLegend* legend = new TLegend(0.37, 0.7, 0.55, 0.88);
    TF1 *mu[n_curves];
    char *funcName = new char[5];
    for (Int_t i = 0; i < n_curves; i++) {
        int color = i + 1;
        short int muValue = (i + 1) * 5;
        sprintf(funcName, "mu%i", muValue);
        mu[i] = new TF1(funcName, integration, 0, 100, 6);
        mu[i]->SetParNames("number of subintervals", "lower bound", "upper bound", "mu", "slope", "intercept");
        mu[i]->SetParameters(n_subint, 0.0, upper_bound, muValue, cell17_slope, cell17_intercept);
        mu[i]->SetLineColor(color);

        char *legendEntryName = new char[10];
        sprintf(legendEntryName, "#mu = %i", muValue);
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
    mu[0]->SetTitle("Convolution with PWGD");
    mu[0]->GetXaxis()->SetTitle("MET [GeV]");
    mu[0]->GetYaxis()->SetTitle("Probability of an Event");
    // cout << "Integral of mu5:  " << mu[0]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu10:  " << mu[1]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu15:  " << mu[2]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu20:  " << mu[3]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu25:  " << mu[4]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu30:  " << mu[5]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu35:  " << mu[6]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu40:  " << mu[7]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu45:  " << mu[8]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu50:  " << mu[9]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu55:  " << mu[10]->Integral(1, 2000) << "\n";
    // cout << "Integral of mu60:  " << mu[11]->Integral(1, 2000) << "\n";
    // dists->SaveAs("results.png");

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
