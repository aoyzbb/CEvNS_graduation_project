#include "BetaBranch.hh"
#include "TH1D.h"
#include "TGraph.h"
#include <iostream>
#include "TMath.h"
#include "TComplex.h"
#include <numeric>
#include <algorithm>
using namespace std;

double BetaBranch::me = 0.511;
double BetaBranch::alpha = 1/137.035999074;

BetaBranch::BetaBranch()
{
	type = BetaMinus;
	fBranchRatio = 0;
	Qvalue = 0;
	Z=0;
	A=0;
	C=-1;
	verbosity = 0;

}

BetaBranch::BetaBranch(int branchId, double ratio, DecayType newtype, double Qv, int z, int a, int ds, int dp,  vector<double>& gammalist)
{
	//cout<<"Construct beta decay brach!"<<endl;
	verbosity = 0;
	type = newtype;
	fBranchRatio = ratio;
	Qvalue = Qv;
	Z = z;
	A = a;
	dSpin = ds;
	dParity = dp;
	if(type==BetaMinus) C = -1;
	else C = 1;
	fGammaEnergy = gammalist;

	// Here calculate energy spectra
	
	// EC only emits gamma, We don't need to calculate EC spectrum
	if(type == EC)
	{
	//	EnergyHisMax = 0;
	//	fSpectraVector.resize(1);
	//	fSpectraVector[0] = 0;
		return;
	}
	
    Zdaughter = Z-C; 
    R = 0.0029*pow(A, 1./3) + 0.0063*pow(A, -1./3) - 0.0017/A; // radius in me
    // R = 0.42587*alpha * pow(A,1./3.);

    //me = 0.511; // mass of electron
    //alpha = 1/137.035999074; // fine structure constant
    
    // tabulated screening potential
    double x[10] = {  1,    8,    13,    16,   23,    27,    29,    49,    84,    92};
    double y[10] = {1.0, 1.42, 1.484, 1.497, 1.52, 1.544, 1.561, 1.637, 1.838, 1.907};
    gScreeningV = new TGraph(10, x, y);

	double binWidth = 1e-3; // MeV
    double Q0 = Qvalue;
    //TString name = Form("%s_%03i", isotope.c_str(), i);
    //TString title = Form("%s_%.1f_MeV", isotope.c_str(), Q0);
    int nBins = int(floor(Q0/binWidth)) + 1;
    //hBeta[i] = new TH1D(name.Data(), title.Data(), nBins, 0, nBins*binWidth);
	EnergyHisMax = nBins*binWidth;
    hBeta = new TH1D(TString::Format("Z%dA%d_%d", z,a,branchId), "BetaSpectrum", nBins, 0, EnergyHisMax);
    
	//hBetaNoCorrection[i] = new TH1D((name+"_nc").Data(), title.Data(), nBins, 0, Q0);
    bool useForbiddenCorrection = UseForbiddenCorrection(dSpin, dParity);
    for (int j=1; j<=nBins; j++) 
	{
        double Te = (j-0.5)*binWidth;
        double Ee = Te+ me;
        double dNdE =sqrt(Ee*Ee - me*me) * (Q0-Te)*(Q0-Te) * Ee/me/me/me/me;
        //hBetaNoCorrection[i]->SetBinContent(j, dNdE);

        dNdE *= Fermi(Zdaughter, Te, C)*Screening(Zdaughter, Te, C);

        if (useForbiddenCorrection) 
		{
            dNdE *= Forbiddenness(Zdaughter, Te, Q0, dSpin, dParity, C);
        }            
        if (!(dParity==1 && (dSpin==0 || dSpin==1))) 
		{  
			// not first non-unique:not 0-,1-
            dNdE *= FiniteSizeEM(Zdaughter, Te, C)
                    *FiniteSizeWI(Zdaughter, Te, Q0, C)
                    * WeakMagnetism(Te, C);
        }
            
	    //    cout <<dNdE << endl;//<<"--"<< Fermi(Zdaughter, Te, C)
        hBeta->SetBinContent(j, dNdE);
    }
    //hBeta->Scale(1./hBeta->Integral());
    //hBetaNoCorrection[i]->Scale(1./hBetaNoCorrection[i]->Integral());

}

BetaBranch::~BetaBranch()
{
	delete gScreeningV;
	delete hBeta;
}

double BetaBranch::GetRandomEk()
{
	if(type==EC)
	  return 0;

	double result = hBeta->GetRandom();
	//delete hBeta;
	return result;
}

bool BetaBranch::UseForbiddenCorrection(int deltaSpin, int deltaParity)
{
    bool doCorrection = false;
    
    if ((deltaSpin == 0 || deltaSpin == 1) && deltaParity == 0) { // allowed decay
        doCorrection = false;                                     
    } else if ((deltaSpin == 0 || deltaSpin == 1) && deltaParity == 1) { // Non-unique First Forbidden Decay
        doCorrection = false;
    } else if (deltaSpin == 2 && deltaParity == 1) { // Unique First Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 3 && deltaParity == 0) { // Unique Second Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 4 && deltaParity == 1) { // Unique Third Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 5 && deltaParity == 0) { // Unique Fourth Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 2 && deltaParity == 0) { // Non-unique Second Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 3 && deltaParity == 1) { // Non-unique Third Forbidden Decay
        doCorrection = false; //crude approximation
    } else if (deltaSpin == 4 && deltaParity == 0) { // Non-unique Fourth Forbidden Decay
        doCorrection = false; //crude approximation
    } else {  // higher orders
        doCorrection = false; //crude approximation allowed
    }
    
    if (doCorrection) {
        if(verbosity>0)cout << "   Forbidden Shape Correction: (" << deltaSpin << ", " << deltaParity << ")" << endl; 
    }
    else {
        if(verbosity>0)cout << "No Forbidden Shape Correction: (" << deltaSpin << ", " << deltaParity << ")" << endl;
    }
    return doCorrection;
}

double BetaBranch::Fermi(int z, double Te, int charge)
{
   
    double Ee = Te + me;
    double pe = sqrt(Ee*Ee-me*me);
    
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    
    double nu = -charge * alpha * z * Ee / pe;

    
    double result =2 * (1+gamma0) * exp(TMath::Pi()*nu)
                  * pow(2*pe*R/me, 2*(gamma0-1)) 
                  / pow(TMath::Gamma(2*gamma0+1), 2);
    TComplex a(gamma0, nu);

    TComplex b = CGamma(a);
    result *= b.Rho2();

    return result;
}

double BetaBranch::Screening(int z, double Te, int charge)
{
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    
    double W = Te/me + 1;
        
    double p = sqrt(W*W-1);
    double y = - charge * alpha * z * W / p;
 
    double Ztilt = z +charge;//error, parent atom should be "Z+charge".Attention,z is daughter atom here.	
    double V0 = - charge * alpha * alpha * pow(Ztilt,4./3.) * gScreeningV->Eval(Ztilt);
    if (W<V0) return 1;
          
                                                            
    double Wtilt = W - V0 > 1 ? W-V0 : 1.00001;
    double ptilt = sqrt(Wtilt*Wtilt-1);
    double ytilt = - charge * alpha * z * Wtilt / ptilt;
      
    TComplex a(gamma0, ytilt);
    TComplex b = CGamma(a);
    TComplex c(gamma0, y);
    TComplex d = CGamma(c);
      
    double result = Wtilt/W 
        * pow(ptilt/p, 2*gamma0-1) 
        * exp(TMath::Pi() * (ytilt-y))
        * b.Rho2() / d.Rho2();
        

    return result;
    
}



double BetaBranch::Forbiddenness(int z, double Te, double Q0, int deltaSpin, int deltaParity, int charge)
{    
    double tmp = deltaParity; tmp = 0; // currently don't use deltaParity in this fuction
    // shape correction to fobidden decays
    double Ee = Te + me;
    double pe = sqrt(Ee*Ee-me*me);
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    double nu = -charge * alpha * z * Ee / pe;

    double result = 1.0;
    double corr = 0;
    
    /* Get the factorials in the normalization of the forbiddenness correction factor */
    double factorial_1 = 1;
    double dubl_fact_1 = 1;
    for(int k=1; k<=(deltaSpin-1); k++) factorial_1 = factorial_1*k;//factorial_1=(deltaSpin-1)!
    for(int k=1; k<=(2*deltaSpin-1); k++) {dubl_fact_1 = dubl_fact_1*k; k++;} //dubl_fact_1=(2deltaSpin-1)!!
    
    for (int j=1; j<=deltaSpin; j++) {
        /* Get Double Factorials over spin */
        double dubl_fact_num = 1;
        for(int k=1; k<=j; k++) { dubl_fact_num = dubl_fact_num*(2*k-1); }//dubl_fact_num =(2j-1)!!
        double dubl_fact_denom = 1;
        for(int k=1; k<=(2*deltaSpin-2*j+1); k++){ dubl_fact_denom = dubl_fact_denom*(k); k++; }//dubl_fact_denom =(2*deltaSpin-2*j+1)!!

        /* Get Single Factorials over spin */
        double factorial_denom_1 = 1;
        for(int k=1; k<=(j-1); k++){ factorial_denom_1=factorial_denom_1*(k);}//factorial_denom_1=(j-1)!
        double factorial_denom_2=1;
        for(int k=1; k<=(deltaSpin-j); k++) { factorial_denom_2=factorial_denom_2*(k); }//factorial_denom_2=(deltaSpin-j)!

        double gamma1 = sqrt(j*j - alpha*z*alpha*z);

        TComplex a(gamma0, nu);
        TComplex b = CGamma(a);
        TComplex c(gamma1, nu);
        TComplex d = CGamma(c);

        corr += 8 * TMath::Pi()
              * factorial_1/dubl_fact_1
              * pow(R, 2.*((double)(deltaSpin-1)))/ (1.+gamma0)
              * dubl_fact_num / (dubl_fact_denom*factorial_denom_1*factorial_denom_2) 
              * pow((Q0-Te)/me, 2.*(deltaSpin-(2*j-1)/2.)-1.)//2(deltaSpin-j),why write so complex?
              * 0.5 * pow(2.*pe/me, 2.*(gamma1-gamma0)) 
          //    * pow(R, 2.*(gamma1-gamma0)-(2.*j-1.))  //error,no this part in Kai Siegbahn's book 
              * j * (j+gamma1)
              * pow(TMath::Gamma(2*gamma0+1)/TMath::Gamma(2*gamma1+1), 2)
              * d.Rho2()/b.Rho2();
    }
    result *= corr;
    return result;
}


double BetaBranch::FiniteSizeEM(int z, double Te, int charge)
{
    
    double W = Te/me + 1;
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    
    //------------add the last 3 items--------------------
    //----------important for large Z---------------------
   
   //-----for negatron---------------------------
    double a1[7];
    double a2[7];
    double an;
    double data1[7][6]=
      {
	{0.115,-1.8123,8.2498,-11.223,-14.854,32.086},
	{-0.00062,0.007165,0.01841,-0.53736,1.2691,-1.5467},
	{0.02482,-0.5975,4.84199,-15.3374,23.9774,-12.6534},
	{-0.14038,3.64953,-38.8143,172.137,-346.708,288.787},
	{0.008152,-1.15664,49.9663,-273.711,657.629,-603.703},
	{1.2145,-23.9931,149.972,-471.299,662.191,-305.68},
	{-1.5632,33.4192,-255.133,938.53,-1641.28,1095.36}
      };
    
    for(int i=0;i<=6;i++){
      for(int j=0;j<=5;j++){
	a1[i]=data1[i][j]*pow(-1*charge*alpha*z,j+1) ;   
      }
    }
    //------------------------------------------------
    
    //--------for positron------------------------
    double data2[7][6]=
      {
	{-0.0701,-2.572,-27.5971,-128.658,-272.264,-214.925},
	{0.002308,0.066483,0.6407,2.83606,5.6317,4.0011},
	{-0.07936,-2.09284,-18.45462,-80.9375,-160.8384,-124.8927},
	{0.93832,22.02513,197.00221,807.1878,1566.6077,1156.3287},
	{-4.276181,-96.82411,-835.26505,-3355.8441,-6411.3255,-4681.573},
	{8.2135,179.0862,1492.1295,5872.5362,11038.7299,7963.4701},
	{-5.4583,-115.8922,-940.8305,-3633.9181,-6727.6296,-4795.0481}
      };
    for(int i=0;i<=6;i++){
      for(int j=0;j<=5;j++){
	a2[i]=data2[i][j]*pow(-1*charge*alpha*z,j+1) ;   
      }
    } 
    //-------------------------------------------------
    
    double result=0;
    if(charge<0){
      an=0;
      for(int i=1;i<=6;i++){
	an+=a1[i]*pow(W*R,i-1);
      }
      
      result = 1 
	+13. * pow(alpha*z, 2) / 60.
	+ charge * W*R*alpha*z*(41-26*gamma0)/15./(2*gamma0-1) 
	+ charge * alpha*z*R*gamma0*(17-2*gamma0)/30./W/(2*gamma0-1)  //error,neglect the last 3 items
	+ a1[0]*R/W+an+0.41*(R-0.0164)*pow(alpha*z,4.5);//the last 3 items , z is always positive here.       
    }     
    
    
    if(charge>0){
      an=0;
      for(int i=1;i<=6;i++){
	an += a2[i]*pow(W*R,i-1);
      }
      
      result = 1 
	+13. * pow(alpha*z, 2) / 60.
	+ charge * W*R*alpha*z*(41-26*gamma0)/15./(2*gamma0-1) 
	+ charge * alpha*z*R*gamma0*(17-2*gamma0)/30./W/(2*gamma0-1)  //error,neglect the last 3 items
	+ a2[0]*R/W+an+0.22*(R-0.0164)*pow(alpha*z,4.5);//the last 3 items, z is always positive here.       
      
    }     
    
    return result;
}

double BetaBranch::FiniteSizeWI(int z, double Te, double Q0, int charge)
{

    double W = Te/me + 1;
    double W0 = Q0/me+1;//error,should be "W0=Q/me+1",W0 is the maximum energy of electron
    
    double C0, C1, C2;
    C0 = -233/630.*pow(alpha*z, 2) - pow(W0*R, 2)/5. - charge * 2./35.*W0*R*alpha*z;//error,should be "-charge"
    C1 = charge * 21./35.*R*alpha*z + 4./9.*W0*R*R;
    C2 = -4./9.*R*R;
    
    double result = 1 + C0 + C1*W + C2*W*W;
    return result;
}


double BetaBranch::WeakMagnetism(double Te, int charge)
{
    return 1 -charge * 0.0048 * (Te+me);  // use a universal WM slope
}                                           



// Complex Gamma function
TComplex BetaBranch::CGamma(TComplex z)
{
  // complex gamma function
  // modified f2c of CERNLIB cgamma64.F
  // taken from Knat and converted to compile without ROOT libs

  Double_t zr=z.Re();

  if (z.Im()==0.0) {
    if (zr==0.0) {
      TComplex w(0.0,0.0);
      return w;
    } else if (-fabs(zr)==floor(zr)) {
      TComplex w(0.0,0.0);
      return w;
    }
  }

  TComplex u;
  TComplex v;

  if (zr>=1.0) {
    u=1.0;
    v=z;
  } else if (zr>=0.0) {
    u=1.0/z;
    v=1.0+z;
  } else {
    u=1.0;
    v=1.0-z;
  }

  const Double_t c1=2.50662827463100050;
  const Double_t c2[16]={
     41.624436916439068,-51.224241022374774,
     11.338755813488977, -0.747732687772388,
      0.008782877493061, -0.000001899030264,
      0.000000001946335, -0.000000000199345,
      0.000000000008433,  0.000000000001486,
     -0.000000000000806,  0.000000000000293,
     -0.000000000000102,  0.000000000000037,
     -0.000000000000014,  0.000000000000006
  };

  TComplex w(1.0,0.0);
  TComplex s(c2[0],0.0);
  TComplex cpi(3.14159265358979312,0.0);

  for (int i=1;i<16;++i) {
    w*=(v-((Double_t)i))/(v+((Double_t)(i-1)));
    s+=c2[i]*w;
  }

  w=v+4.5;
  TComplex logw = TComplex::Log(w);
  w=c1*s*TComplex::Exp((v-0.5)*TComplex::Log(w)-w);
  if (zr<0.0) w=cpi/(TComplex::Sin(cpi*z)*w);

  return w*u;
}
