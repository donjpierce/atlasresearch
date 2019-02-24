#include <TROOT.h>
#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TF1.h>



Double_t sumet_func_fft ( Double_t *x, Double_t *parm) {
  static Double_t parmsave[6] = {-999.,-999.,-999.,-999.,-999.,-999.};
  static int ncalls=0;
  //int n=62914560;
  int n=8847360;
  int nparms=6;
  Double_t lowval =0.0;
  Double_t hival =50000.0;
  Double_t xtemp;
  int itemp;
  long double functot[n];
  //static Double_t finalfuncval[62914560];
  static Double_t finalfuncval[8847360];
  bool sameparm;

  TH1D *funchist = 0;
  TH1 *ftransform=0;

  ncalls++;
  if(ncalls%10000==0) {
    std::cout << "Number of routine calls = " << ncalls << "\n";
  }

  if(x[0]-parm[3]<0.) {
    return 0.;
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

  itemp=int((x[0]-parm[3])*double(n)/hival);
  if(itemp>n-1 || itemp<0){
    return 0.;
  }
  else {
    //std::cout << " itemp " << itemp << " x = " << x[0] << "  f(x) " << finalfuncval[itemp] <<"\n";
    return parm[0]*finalfuncval[itemp];
  }

}
