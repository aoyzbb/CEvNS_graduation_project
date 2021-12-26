#include "TMath.h"

class SeaLevelMuon
{
public:
	SeaLevelMuon() 
	{
		p1=0.102573;
		p2=-0.068287;
		p3=0.958633;
		p4=0.0407253;
		p5=0.817285;
	}
	~SeaLevelMuon() {}

	Double_t StandardGaisser(Double_t* x, Double_t* /*p*/)
	{
		Double_t E=x[0];
		Double_t costheta=x[1];
		Double_t val=0.14*TMath::Power(E,-2.7)*(1./(1.+1.1*E*costheta/115)+0.054/(1+1.1*E*costheta/850));
		return val;
	}
	
	Double_t logStandardGaisser(Double_t* x, Double_t* /*p*/)
	{
		Double_t E=x[0];
		Double_t Ed = TMath::Exp(E);
		Double_t costheta=x[1];
		Double_t val=0.14*TMath::Power(Ed,-2.7)*(1./(1.+(1.1*Ed*costheta)/(115))+0.054/(1+(1.1*Ed*costheta)/(850)))*Ed;
		return val;
	}
	
	Double_t logDistribution(Double_t* x, Double_t* /*p*/)
	{
		Double_t E=x[0];
		Double_t Ed = TMath::Exp(E);
		Double_t costheta=x[1];
		Double_t costhetastar=fCosthetastar(costheta);
		Double_t E1=fEstar(Ed,costhetastar);
		Double_t val=0.14*TMath::Power(E1,-2.7)*(1./(1.+(1.1*Ed*costhetastar)/(115))+0.054/(1+(1.1*Ed*costhetastar)/(850)))*Ed;
		return val;
	}
	Double_t Distribution(Double_t* x, Double_t* /*p*/)
	{
		Double_t E=x[0];
		Double_t costheta=x[1];
		Double_t costhetastar=fCosthetastar(costheta);
		Double_t E1=fEstar(E,costhetastar);
		Double_t val=0.14*TMath::Power(E1,-2.7)*(1./(1.+(1.1*E*costheta)/115.)+0.054/(1+(1.1*E*costheta)/850.));   
		return val;

	}
	
private:
	Double_t p1,p2,p3,p4,p5;
	Double_t fCosthetastar(Double_t costheta)
	{
		Double_t val=sqrt((costheta*costheta+p1*p1+p2*TMath::Power(costheta,p3)+p4*TMath::Power(costheta,p5))/(1+p1*p1+p2+p4));
		return val;
	}	
	Double_t fEstar(Double_t E, Double_t costhetastar)
	{
		return E*(1+3.64/(E*TMath::Power(costhetastar,1.29)));
	}	

};

