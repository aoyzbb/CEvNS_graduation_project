#include "DecaySpec.hh"
#include <iostream>
#include <time.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

int main()
{
	//gRandom->SetSeed(0);

	TFile* f1 = new TFile("output.root", "recreate");
	TTree* t = new TTree("energy","merged");
	Double_t EkMerged = 0;
	Double_t betaEk;
	Int_t nGammas;
	Double_t GammaE[50];
	t->Branch("EkMerged", &EkMerged);
	t->Branch("betaEk", &betaEk);
	t->Branch("nGammas", &nGammas);
	t->Branch("gammaEk", GammaE);

	DecaySpec* radDecay = new DecaySpec("Bi214", 83, 214);

	for(Int_t j=0; j<10000; j++)
	{
		EkMerged = 0;
		//cout<<"Random particle: "<<endl;
		double ek;
		vector<double> eG;
		DecayType type = radDecay->GetVertex(ek, eG);
		EkMerged += ek;
		betaEk = ek;
		nGammas = eG.size(); 
		//cout<<DecayBranch::DecayTypeString[type]<<": "<<ek<<"MeV"<<endl;
		for(size_t i=0; i<eG.size(); i++)
		{
			//cout<<"Gamma: "<<eG[i]<<"MeV"<<endl;
			EkMerged += eG[i];
			GammaE[i] = eG[i];
		}
		t->Fill();
	}
	t->Write();
	f1->Close();
	delete radDecay;
	return 0;

}
