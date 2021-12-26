#include "DecaySpec.hh"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "AlphaBranch.hh"
#include "BetaBranch.hh"
#include "TRandom3.h"
#include "Randomize.hh"
using namespace std;

const char* DecayBranch::DecayTypeString[] = {"Alpha", "BetaMinus", "BetaPlus", "EC", "Gamma"};

DecaySpec::DecaySpec()
{


}


DecaySpec::DecaySpec(std::string name, int z, int a)
{
	//cout<<"Construct DecaySpec!"<<endl;
	isotope = name;
	Z = z;
	A = a;

	string path = "/mnt/g/books/graduation_project/python_simulation/geant4_simulation/coherent_new_bkg/source/Generator/data";

	// Load Alpha decay data
	string alphaDataFile = path+"/alpha/"+isotope+".txt";
	LoadDataFile(alphaDataFile, Alpha);

	// Load Beta- decay data
	string betaDataFile = path+"/beta/"+isotope+".txt";
	LoadDataFile(betaDataFile, BetaMinus);

	// Load Beta+ decay data
	string beta1DataFile = path+"/betaplus/"+isotope+".txt";
	LoadDataFile(beta1DataFile, BetaPlus);

	// Load EC data
	string beta2DataFile = path+"/EC/"+isotope+".txt";
	LoadDataFile(beta2DataFile, EC);

	/*
	cout<<"BranchList size: "<<DecayBranchList.size()<<endl;
	for(auto&& branch : DecayBranchList)
	{
		cout<<branch->GetType()<<"\t"<<branch->GetBranchRatio()<<endl;

	}*/

}

DecaySpec::~DecaySpec()
{
	for(size_t i=0; i<DecayBranchList.size(); i++)
	  delete DecayBranchList[i];

}

void DecaySpec::LoadDataFile(std::string filename,  DecayType type)
{

	ifstream f(filename.c_str());
	if (f.good()) 
	{
		//cout << "loading " << isotope<<" "<<DecayBranch::DecayTypeString[type]
		//	<< " decay data from " << filename << endl;
		//cout << "Type\tQ value\tB.R.\tdeltaSpin\tdeltaParity\tnGamma\teGammas\n";
		//cout << "=================================================" << endl;
		string line;
		double Ek, br, ds, dp, nG;
		vector<double> gammalist;
		int nlines = 0;
		while (f.good()) 
		{
			getline(f, line);
			// Come across problem in Windows noeol
			//if(f.eof())
			//  break;
			nlines++;

			// Check if the nGamma and eGammas are matched
			string check;
			stringstream inputcheck(line);
			int n = 0;
			int nG1 = 0;
			while(inputcheck>>check)
			{
				n++;
				if(n==5)
				{
					stringstream nG2(check);
					nG2>>nG1;
				}
			}
			if(n<2) continue;
			if( nG1!=0 && n != 5+(nG1) )
			{
				cerr<<"Error: nGammas and eGammas are not matched at "<<filename<<" Line "<<nlines<<endl;
				exit(0);
			}

			// Read alpha decay data table
			stringstream input(line);
			input>>Ek>>br>>ds>>dp>>nG;
			gammalist.clear();
			gammalist.resize(nG);
			for(size_t i=0; i<nG; i++)
			  input>>gammalist[i];

			
			switch(type)
			{
			case Alpha:
				DecayBranchList.push_back(
							new AlphaBranch(br, Ek, gammalist));
				//cout<<"Alpha\t";
				break;
			case BetaMinus:
				{
				DecayBranchList.emplace_back(new BetaBranch
							(DecayBranchList.size(), br, BetaMinus, Ek, Z, A, ds,dp, gammalist));
				//cout<<"Beta-\t";
				break;
				}
			case BetaPlus:
				DecayBranchList.push_back(new BetaBranch
							(DecayBranchList.size()-1, br, BetaPlus, Ek, Z, A, ds, dp, gammalist));
				//cout<<"Beta+\t";
				break;
			case EC:
				DecayBranchList.push_back(
					new BetaBranch(DecayBranchList.size()-1, br, EC, Ek, Z, A, ds, dp, gammalist));
				//cout<<"EC\t";
				break;

			default:
				break;


			}

			//cout<<Ek<<"\t"<<br<<"\t"<<ds<<"\t"<<dp<<"\t"<<nG<<"\t";
			//for(size_t i=0; i<nG; i++)
			//  cout<<gammalist[i];
			//cout<<endl;
		}
		f.close();
		//cout << "=================================================" << endl;
	}



}


DecayType DecaySpec::GetVertex(double& Ek, vector<double>& eGammas)
{
	Ek = 1;
	eGammas.clear();
	size_t nsize = DecayBranchList.size();
	//Double_t rndm = gRandom->Uniform();
	Double_t rndm = G4UniformRand();
	Double_t brt = 0;

	/*size_t iB=138;
	Ek = DecayBranchList[iB]->GetRandomEk();
	for(size_t i=0; i<DecayBranchList[iB]->fGammaEnergy.size(); i++)
		{
			eGammas.push_back(DecayBranchList[iB]->fGammaEnergy[i]);
			//cout<<"Gamma\t"<<DecayBranchList[iB]->fGammaEnergy[i]<<endl;
		}
	return DecayBranchList[iB]->GetType();*/

	for(size_t iB=0;iB<nsize;iB++)
	{
		brt+=DecayBranchList[iB]->GetBranchRatio();
		if(rndm>brt)
			continue;

		Ek = DecayBranchList[iB]->GetRandomEk();
		//cout<<iB<<"\t"<<DecayBranchList[iB]->GetBranchRatio()<<"\t"<<DecayBranch::DecayTypeString[DecayBranchList[iB]->GetType()]<<"\t"<<Ek<<endl;

		for(size_t i=0; i<DecayBranchList[iB]->fGammaEnergy.size(); i++)
		{
			eGammas.push_back(DecayBranchList[iB]->fGammaEnergy[i]);
			//cout<<"Gamma\t"<<DecayBranchList[iB]->fGammaEnergy[i]<<endl;
		}
		
		return DecayBranchList[iB]->GetType();
	}

	return INVALID;
}
