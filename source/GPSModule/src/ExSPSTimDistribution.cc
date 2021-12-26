//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        ExSPSTimDistribution.cc
//
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
//
// 2016, Linyan, Created.
//    To enable time distribution storing and sampling
//
///////////////////////////////////////////////////////////////////////////////
//

#include "ExSPSTimDistribution.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4double ExSPSTimDistribution::fSegmentSize = 1*CLHEP::second;

ExSPSTimDistribution::ExSPSTimDistribution()
  : particle_definition(0), eneRndm(0), Splinetemp(0)
{
	//
	// Initialise all variables
	particle_time = 0.0 ;
	bDependency = false;

	TimeDisType = "Lin";
	weight = 1.;
	MonoTime = 0.0;
	Emin = 0.;
	Emax = 1.e9;
	alpha = 0.;
	biasalpha = 0.;
        prob_norm = 1.0;
	Ezero = 0.;
	SE = 0.;
	grad = 0.;
	cept = 1.;
	Biased = false; // not biased
	TimeSpec = true; // true - time spectra, false - momentum spectra
	DiffSpec = true; // true - differential spec, false integral spec
	IntType = "NULL"; // Interpolation type
	IPDFTimeExist = false;
	IPDFArbExist = false;

	ArbEmin = 0.;
	ArbEmax = 1.e30;

	verbosityLevel = 0;

}

ExSPSTimDistribution::~ExSPSTimDistribution() {
}

void ExSPSTimDistribution::SetTimeDisType(G4String DisType) {
	TimeDisType = DisType;
	if (TimeDisType == "User") {
		UDefTimeH = IPDFTimeH = ZeroPhysVector;
		IPDFTimeExist = false;
	} else if (TimeDisType == "Arb") {
		ArbTimeH = IPDFArbTimeH = ZeroPhysVector;
		IPDFArbExist = false;
	} else if (TimeDisType == "Epn") {
		UDefTimeH = IPDFTimeH = ZeroPhysVector;
		IPDFTimeExist = false;
		EpnTimeH = ZeroPhysVector;
	}
}

void ExSPSTimDistribution::SetEmin(G4double emi) {
	Emin = emi;
}

void ExSPSTimDistribution::SetEmax(G4double ema) {
	Emax = ema;
}

void ExSPSTimDistribution::SetMonoTime(G4double mtime) {
	MonoTime = mtime;
}

void ExSPSTimDistribution::SetBeamSigmaInE(G4double e) {
	SE = e;
}
void ExSPSTimDistribution::SetAlpha(G4double alp) {
	alpha = alp;
}

void ExSPSTimDistribution::SetBiasAlpha(G4double alp) {
	biasalpha = alp;
	Biased = true;
}

void ExSPSTimDistribution::SetEzero(G4double eze) {
	Ezero = eze;
}

void ExSPSTimDistribution::SetGradient(G4double gr) {
	grad = gr;
}

void ExSPSTimDistribution::SetInterCept(G4double c) {
	cept = c;
}

void ExSPSTimDistribution::UserTimeHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In UserTimeHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	UDefTimeH.InsertValues(ehi, val);
	Emax = ehi;
}

void ExSPSTimDistribution::ArbTimeHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In ArbTimeHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	ArbTimeH.InsertValues(ehi, val);
}

void ExSPSTimDistribution::ArbTimeHistoFile(G4String filename) {
	std::ifstream infile(filename, std::ios::in);
	if (!infile)
		G4Exception("ExSPSTimDistribution::ArbTimeHistoFile",
                "Event0301",FatalException,
                "Unable to open the histo ASCII file");
	G4double ehi, val;
	while (infile >> ehi >> val) {
		ArbTimeH.InsertValues(ehi, val);
	}
}

void ExSPSTimDistribution::EpnTimeHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In EpnTimeHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	EpnTimeH.InsertValues(ehi, val);
	Emax = ehi;
	Epnflag = true;
}


void ExSPSTimDistribution::InputTimeSpectra(G4bool value) {
	// Allows user to specifiy spectrum is momentum
	TimeSpec = value; // false if momentum
	if (verbosityLevel > 1)
		G4cout << "TimeSpec has value " << TimeSpec << G4endl;
}

void ExSPSTimDistribution::InputDifferentialSpectra(G4bool value) {
	// Allows user to specify integral or differential spectra
	DiffSpec = value; // true = differential, false = integral
	if (verbosityLevel > 1)
		G4cout << "Diffspec has value " << DiffSpec << G4endl;
}

void ExSPSTimDistribution::ArbInterpolate(G4String IType) {
	if (TimeDisType != "Arb")
		G4cout << "Error: this is for arbitrary distributions" << G4endl;
	IntType = IType;
	ArbEmax = ArbTimeH.GetMaxLowEdgeEnergy();
	ArbEmin = ArbTimeH.GetMinLowEdgeEnergy();

	// Now interpolate points
	if (IntType == "Lin")
		LinearInterpolation();
	if (IntType == "Log")
		LogInterpolation();
	if (IntType == "Exp")
		ExpInterpolation();
	if (IntType == "Spline")
		SplineInterpolation();
}

void ExSPSTimDistribution::LinearInterpolation() {
	// Method to do linear interpolation on the Arb points
	// Calculate equation of each line segment, max 1024.
	// Calculate Area under each segment
	// Create a cumulative array which is then normalised Arb_Cum_Area

	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbTimeH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbTimeH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbTimeH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made time.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (TimeSpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Time**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - time equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to time unit and its value to per time unit
			G4double total_time;
			for (count = 0; count < maxi; count++) {
				total_time = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total time

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_time;
				Arb_x[count] = total_time - mass; // kinetic time
			}
		}
	}
	//
	i = 1;
	Arb_grad[0] = 0.;
	Arb_cept[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	while (i < maxi) {
		// calc gradient and intercept for each segment
		Arb_grad[i] = (Arb_y[i] - Arb_y[i - 1]) / (Arb_x[i] - Arb_x[i - 1]);
		if (verbosityLevel == 2)
			G4cout << Arb_grad[i] << G4endl;
		if (Arb_grad[i] > 0.) {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is positive" << G4endl;
			Arb_cept[i] = Arb_y[i] - (Arb_grad[i] * Arb_x[i]);
		} else if (Arb_grad[i] < 0.) {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is negative" << G4endl;
			Arb_cept[i] = Arb_y[i] + (-Arb_grad[i] * Arb_x[i]);
		} else {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is 0." << G4endl;
			Arb_cept[i] = Arb_y[i];
		}

		Area_seg[i] = ((Arb_grad[i] / 2) * (Arb_x[i] * Arb_x[i] - Arb_x[i - 1]
				* Arb_x[i - 1]) + Arb_cept[i] * (Arb_x[i] - Arb_x[i - 1]));
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		sum = sum + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_x[i] << Arb_y[i] << Area_seg[i] << sum << Arb_grad[i]
					<< G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
		IPDFArbTimeH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbTimeH, needed by Probability()
	ArbTimeH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1) {
		G4cout << "Leaving LinearInterpolation" << G4endl;
		ArbTimeH.DumpValues();
		IPDFArbTimeH.DumpValues();
	}
}

void ExSPSTimDistribution::LogInterpolation() {
	// Interpolation based on Logarithmic equations
	// Generate equations of line segments
	// y = Ax**alpha => log y = alpha*logx + logA
	// Find area under line segments
	// create normalised, cumulative array Arb_Cum_Area
	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbTimeH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbTimeH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbTimeH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made time.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (TimeSpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Time**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - time equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to time unit and its value to per time unit
			G4double total_time;
			for (count = 0; count < maxi; count++) {
				total_time = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total time

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_time;
				Arb_x[count] = total_time - mass; // kinetic time
			}
		}
	}
	//
	i = 1;
	Arb_alpha[0] = 0.;
	Arb_Const[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	if (Arb_x[0] <= 0. || Arb_y[0] <= 0.) {
		G4cout << "You should not use log interpolation with points <= 0."
				<< G4endl;
		G4cout << "These will be changed to 1e-20, which may cause problems"
				<< G4endl;
		if (Arb_x[0] <= 0.)
			Arb_x[0] = 1e-20;
		if (Arb_y[0] <= 0.)
			Arb_y[0] = 1e-20;
	}

	G4double alp;
	while (i < maxi) {
		// Incase points are negative or zero
		if (Arb_x[i] <= 0. || Arb_y[i] <= 0.) {
			G4cout << "You should not use log interpolation with points <= 0."
					<< G4endl;
			G4cout
					<< "These will be changed to 1e-20, which may cause problems"
					<< G4endl;
			if (Arb_x[i] <= 0.)
				Arb_x[i] = 1e-20;
			if (Arb_y[i] <= 0.)
				Arb_y[i] = 1e-20;
		}

		Arb_alpha[i] = (std::log10(Arb_y[i]) - std::log10(Arb_y[i - 1]))
				/ (std::log10(Arb_x[i]) - std::log10(Arb_x[i - 1]));
		Arb_Const[i] = Arb_y[i] / (std::pow(Arb_x[i], Arb_alpha[i]));
		alp = Arb_alpha[i] + 1;
		if (alp == 0.) {
		  Area_seg[i] =	Arb_Const[i] * (std::log(Arb_x[i]) - std::log(Arb_x[i - 1])); 
		} else {
		  Area_seg[i] = (Arb_Const[i] / alp) * (std::pow(Arb_x[i], alp)
				- std::pow(Arb_x[i - 1], alp));
		}
		sum = sum + Area_seg[i];
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_alpha[i] << Arb_Const[i] << Area_seg[i] << G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
		IPDFArbTimeH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbTimeH, needed by Probability()
	ArbTimeH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1)
		G4cout << "Leaving LogInterpolation " << G4endl;
}

void ExSPSTimDistribution::ExpInterpolation() {
	// Interpolation based on Exponential equations
	// Generate equations of line segments
	// y = Ae**-(x/e0) => ln y = -x/e0 + lnA
	// Find area under line segments
	// create normalised, cumulative array Arb_Cum_Area
	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbTimeH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbTimeH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbTimeH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made time.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (TimeSpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Time**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - time equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to time unit and its value to per time unit
			G4double total_time;
			for (count = 0; count < maxi; count++) {
				total_time = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total time

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_time;
				Arb_x[count] = total_time - mass; // kinetic time
			}
		}
	}
	//
	i = 1;
	Arb_ezero[0] = 0.;
	Arb_Const[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	while (i < maxi) {
		G4double test = std::log(Arb_y[i]) - std::log(Arb_y[i - 1]);
		if (test > 0. || test < 0.) {
			Arb_ezero[i] = -(Arb_x[i] - Arb_x[i - 1]) / (std::log(Arb_y[i])
					- std::log(Arb_y[i - 1]));
			Arb_Const[i] = Arb_y[i] / (std::exp(-Arb_x[i] / Arb_ezero[i]));
			Area_seg[i] = -(Arb_Const[i] * Arb_ezero[i]) * (std::exp(-Arb_x[i]
					/ Arb_ezero[i]) - std::exp(-Arb_x[i - 1] / Arb_ezero[i]));
		} else {
			G4cout << "Flat line segment: problem" << G4endl;
			Arb_ezero[i] = 0.;
			Arb_Const[i] = 0.;
			Area_seg[i] = 0.;
		}
		sum = sum + Area_seg[i];
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_ezero[i] << Arb_Const[i] << Area_seg[i] << G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
		IPDFArbTimeH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbTimeH, needed by Probability()
	ArbTimeH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1)
		G4cout << "Leaving ExpInterpolation " << G4endl;
}

void ExSPSTimDistribution::SplineInterpolation() {
	// Interpolation using Splines.
	// Create Normalised arrays, make x 0->1 and y hold
	// the function (Time)
        // 
        // Current method based on the above will not work in all cases. 
        // new method is implemented below.
  
	G4double sum, Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;

	G4int maxi = ArbTimeH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbTimeH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbTimeH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made time.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (TimeSpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
                    G4Exception("ExSPSTimDistribution::SplineInterpolation",
                                "Event0302",FatalException,
			        "Error: particle not defined");
		else {
			// Apply Time**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - time equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to time unit and its value to per time unit
			G4double total_time;
			for (count = 0; count < maxi; count++) {
				total_time = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total time

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_time;
				Arb_x[count] = total_time - mass; // kinetic time
			}
		}
	}

	//
	i = 1;
	Arb_Cum_Area[0] = 0.;
	sum = 0.;
	Splinetemp = new G4DataInterpolation(Arb_x, Arb_y, maxi, 0., 0.);
	G4double ei[101],prob[101];
	while (i < maxi) {
	  // 100 step per segment for the integration of area
	  G4double de = (Arb_x[i] - Arb_x[i - 1])/100.;
	  G4double area = 0.;

	  for (count = 0; count < 101; count++) {
	    ei[count] = Arb_x[i - 1] + de*count ;
	    prob[count] =  Splinetemp->CubicSplineInterpolation(ei[count]);
	    if (prob[count] < 0.) { 
              G4ExceptionDescription ED;
	      ED << "Warning: G4DataInterpolation returns value < 0  " << prob[count] <<" "<<ei[count]<< G4endl;
              G4Exception("ExSPSTimDistribution::SplineInterpolation","Event0303",
              FatalException,ED);
	    }
	    area += prob[count]*de;
	  }
	  Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + area;
	  sum += area; 

	  prob[0] = prob[0]/(area/de);
	  for (count = 1; count < 100; count++)
	    prob[count] = prob[count-1] + prob[count]/(area/de);

	  SplineInt[i] = new G4DataInterpolation(prob, ei, 101, 0., 0.);
	  // note i start from 1!
	  i++;
	}
	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
		IPDFArbTimeH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}
	// now scale the ArbTimeH, needed by Probability()
	ArbTimeH.ScaleVector(1., 1./sum);

	if (verbosityLevel > 0)
	  G4cout << "Leaving SplineInterpolation " << G4endl;
}

void ExSPSTimDistribution::GenerateMonoTime() {
	// Method to generate MonoEnergetic particles.
	particle_time = MonoTime;
//	if (verbosityLevel >= 1)G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::GenerateGaussTimes() {
	// Method to generate Gaussian particles.
	particle_time = G4RandGauss::shoot(MonoTime,SE);
	if (particle_time < 0) particle_time = 0.;
//	if (verbosityLevel >= 1)G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::GenerateLinearTimes(G4bool bArb = false) {
	G4double rndm;
	G4double emaxsq = std::pow(Emax, 2.); //Emax squared
	G4double eminsq = std::pow(Emin, 2.); //Emin squared
	G4double intersq = std::pow(cept, 2.); //cept squared

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	G4double bracket = ((grad / 2.) * (emaxsq - eminsq) + cept * (Emax - Emin));
	bracket = bracket * rndm;
	bracket = bracket + (grad / 2.) * eminsq + cept * Emin;
	// Now have a quad of form m/2 E**2 + cE - bracket = 0
	bracket = -bracket;
	//  G4cout << "BRACKET" << bracket << G4endl;
	if (grad != 0.) {
		G4double sqbrack = (intersq - 4 * (grad / 2.) * (bracket));
		//      G4cout << "SQBRACK" << sqbrack << G4endl;
		sqbrack = std::sqrt(sqbrack);
		G4double root1 = -cept + sqbrack;
		root1 = root1 / (2. * (grad / 2.));

		G4double root2 = -cept - sqbrack;
		root2 = root2 / (2. * (grad / 2.));

		//      G4cout << root1 << " roots " << root2 << G4endl;

		if (root1 > Emin && root1 < Emax)
			particle_time = root1;
		if (root2 > Emin && root2 < Emax)
			particle_time = root2;
	} else if (grad == 0.)
		// have equation of form cE - bracket =0
		particle_time = bracket / cept;

	if (particle_time < 0.)
		particle_time = -particle_time;

//	if (verbosityLevel >= 1)G4cout << "Time is " << Emax<<" "<<Emin<<" "<<cept<<" "<<grad<<" "<<bracket<<" "<<particle_time << G4endl;
}

void ExSPSTimDistribution::GeneratePowTimes(G4bool bArb = false) {
	// Method to generate particle energies distributed as
	// a power-law

	G4double rndm;
	G4double emina, emaxa;

	emina = std::pow(Emin, alpha + 1);
	emaxa = std::pow(Emax, alpha + 1);

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	if (alpha != -1.) {
		particle_time = ((rndm * (emaxa - emina)) + emina);
		particle_time = std::pow(particle_time, (1. / (alpha + 1.)));
	} else {
		particle_time = (std::log(Emin) + rndm * (std::log(Emax) - std::log(
				Emin)));
		particle_time = std::exp(particle_time);
	}
//	if (verbosityLevel >= 1)
//		G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::GenerateBiasPowTimes() {
	// Method to generate particle energies distributed as
	// in biased power-law and calculate its weight

        G4double rndm;
	G4double emina, emaxa, emin, emax;

	G4double normal = 1. ;

	emin = Emin;
	emax = Emax;
	//	if (TimeDisType == "Arb") { 
	//  emin = ArbEmin;
	//  emax = ArbEmax;
	//}

	rndm = eneRndm->GenRandEnergy();

	if (biasalpha != -1.) {
	        emina = std::pow(emin, biasalpha + 1);
	        emaxa = std::pow(emax, biasalpha + 1);
		particle_time = ((rndm * (emaxa - emina)) + emina);
		particle_time = std::pow(particle_time, (1. / (biasalpha + 1.)));
		normal = 1./(1+biasalpha) * (emaxa - emina);
	} else {
		particle_time = (std::log(emin) + rndm * (std::log(emax) - std::log(
				emin)));
		particle_time = std::exp(particle_time);
		normal = std::log(emax) - std::log(emin) ;
	}
	weight = GetProbability(particle_time) / (std::pow(particle_time,biasalpha)/normal);

//	if (verbosityLevel >= 1)
//		G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::GenerateExpTimes(G4bool bArb = false) {
	// Method to generate particle energies distributed according
	// to an exponential curve.
	//G4cout<<"Exp"<<G4endl;
	G4double rndm;

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	//G4cout<<Ezero<<G4endl;
	//particle_time = -Ezero * (std::log(rndm * (std::exp(-Emax / Ezero)
	//		- std::exp(-Emin / Ezero)) + std::exp(-Emin / Ezero)));
	double ExpMax = 60*Emax;
	particle_time = -Ezero * (std::log(rndm * (std::exp(-ExpMax / Ezero)
			- std::exp(-Emin / Ezero)) + std::exp(-Emin / Ezero)));
	if (verbosityLevel >= 1)
		G4cout << "Time is " << particle_time << G4endl;
}


void ExSPSTimDistribution::GenUserHistTimes() {
	// Histograms are DIFFERENTIAL.
	//  G4cout << "In GenUserHistTimes " << G4endl;
	if (IPDFTimeExist == false) {
		G4int ii;
		G4int maxbin = G4int(UDefTimeH.GetVectorLength());
		G4double bins[1024], vals[1024], sum;
		sum = 0.;

		if ((TimeSpec == false) && (particle_definition == NULL))
			G4cout << "Error: particle definition is NULL" << G4endl;

		if (maxbin > 1024) {
			G4cout << "Maxbin > 1024" << G4endl;
			G4cout << "Setting maxbin to 1024, other bins are lost" << G4endl;
		}

		if (DiffSpec == false)
			G4cout << "Histograms are Differential!!! " << G4endl;
		else {
			bins[0] = UDefTimeH.GetLowEdgeEnergy(size_t(0));
			vals[0] = UDefTimeH(size_t(0));
			sum = vals[0];
			for (ii = 1; ii < maxbin; ii++) {
				bins[ii] = UDefTimeH.GetLowEdgeEnergy(size_t(ii));
				vals[ii] = UDefTimeH(size_t(ii)) + vals[ii - 1];
				sum = sum + UDefTimeH(size_t(ii));
			}
		}

		if (TimeSpec == false) {
			G4double mass = particle_definition->GetPDGMass();
			// multiply the function (vals) up by the bin width
			// to make the function counts/s (i.e. get rid of momentum
			// dependence).
			for (ii = 1; ii < maxbin; ii++) {
				vals[ii] = vals[ii] * (bins[ii] - bins[ii - 1]);
			}
			// Put time bins into new histo, plus divide by time bin width
			// to make evals counts/s/time
			for (ii = 0; ii < maxbin; ii++) {
				bins[ii] = std::sqrt((bins[ii] * bins[ii]) + (mass * mass))
						- mass; //kinetic time
			}
			for (ii = 1; ii < maxbin; ii++) {
				vals[ii] = vals[ii] / (bins[ii] - bins[ii - 1]);
			}
			sum = vals[maxbin - 1];
			vals[0] = 0.;
		}
		for (ii = 0; ii < maxbin; ii++) {
			vals[ii] = vals[ii] / sum;
			IPDFTimeH.InsertValues(bins[ii], vals[ii]);
		}

		// Make IPDFTimeExist = true
		IPDFTimeExist = true;
		if (verbosityLevel > 1)
			IPDFTimeH.DumpValues();
	}

	// IPDF has been create so carry on
	G4double rndm = eneRndm->GenRandEnergy();
	particle_time = IPDFTimeH.GetEnergy(rndm);

//	if (verbosityLevel >= 1)
//		G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::GenArbPointTimes() {
	if (verbosityLevel > 0)
		G4cout << "In GenArbPointTimes" << G4endl;
	G4double rndm;
	rndm = eneRndm->GenRandEnergy();
	//      IPDFArbTimeH.DumpValues();
	// Find the Bin
	// have x, y, no of points, and cumulative area distribution
	G4int nabove, nbelow = 0, middle;
	nabove = IPDFArbTimeH.GetVectorLength();
	//      G4cout << nabove << G4endl;
	// Binary search to find bin that rndm is in
	while (nabove - nbelow > 1) {
	  middle = (nabove + nbelow) / 2;
	  if (rndm == IPDFArbTimeH(size_t(middle)))
	    break;
	  if (rndm < IPDFArbTimeH(size_t(middle)))
	    nabove = middle;
	  else
	    nbelow = middle;
	}
	if (IntType == "Lin") {
	  Emax = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow));
	  grad = Arb_grad[nbelow + 1];
	  cept = Arb_cept[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << grad << " " << cept << G4endl;
	  GenerateLinearTimes(true);
	} else if (IntType == "Log") {
	  Emax = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow));
	  alpha = Arb_alpha[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << alpha << G4endl;
	  GeneratePowTimes(true);
	} else if (IntType == "Exp") {
	  Emax = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow));
	  Ezero = Arb_ezero[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << Ezero << G4endl;
	  GenerateExpTimes(true);
	} else if (IntType == "Spline") {
	  Emax = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbTimeH.GetLowEdgeEnergy(size_t(nbelow));
	  particle_time = -1e100;
	  rndm = eneRndm->GenRandEnergy();
	  while (particle_time < Emin || particle_time > Emax) {
	    particle_time = SplineInt[nbelow+1]->CubicSplineInterpolation(rndm);
	    rndm = eneRndm->GenRandEnergy();
	  }
//	  if (verbosityLevel >= 1)
//	    G4cout << "Time is " << particle_time << G4endl;
	} else
		G4cout << "Error: IntType unknown type" << G4endl;
}

void ExSPSTimDistribution::GenEpnHistTimes() {
	//  G4cout << "In GenEpnHistTimes " << Epnflag << G4endl;

	// Firstly convert to time if not already done.
	if (Epnflag == true)
	// epnflag = true means spectrum is epn, false means e.
	{
		// convert to time by multiplying by A number
		ConvertEPNToTime();
		// EpnTimeH will be replace by UDefTimeH.
		//      UDefTimeH.DumpValues();
	}

	//  G4cout << "Creating IPDFTime if not already done so" << G4endl;
	if (IPDFTimeExist == false) {
		// IPDF has not been created, so create it
		G4double bins[1024], vals[1024], sum;
		G4int ii;
		G4int maxbin = G4int(UDefTimeH.GetVectorLength());
		bins[0] = UDefTimeH.GetLowEdgeEnergy(size_t(0));
		vals[0] = UDefTimeH(size_t(0));
		sum = vals[0];
		for (ii = 1; ii < maxbin; ii++) {
			bins[ii] = UDefTimeH.GetLowEdgeEnergy(size_t(ii));
			vals[ii] = UDefTimeH(size_t(ii)) + vals[ii - 1];
			sum = sum + UDefTimeH(size_t(ii));
		}

		for (ii = 0; ii < maxbin; ii++) {
			vals[ii] = vals[ii] / sum;
			IPDFTimeH.InsertValues(bins[ii], vals[ii]);
		}
		// Make IPDFEpnExist = true
		IPDFTimeExist = true;
	}
	//  IPDFTimeH.DumpValues();
	// IPDF has been create so carry on
	G4double rndm = eneRndm->GenRandEnergy();
	particle_time = IPDFTimeH.GetEnergy(rndm);

//	if (verbosityLevel >= 1)
//		G4cout << "Time is " << particle_time << G4endl;
}

void ExSPSTimDistribution::ConvertEPNToTime() {
	// Use this before particle generation to convert  the
	// currently stored histogram from time/nucleon
	// to time.
	//  G4cout << "In ConvertEpntoTime " << G4endl;
	if (particle_definition == NULL)
		G4cout << "Error: particle not defined" << G4endl;
	else {
		// Need to multiply histogram by the number of nucleons.
		// Baryon Number looks to hold the no. of nucleons.
		G4int Bary = particle_definition->GetBaryonNumber();
		//      G4cout << "Baryon No. " << Bary << G4endl;
		// Change values in histogram, Read it out, delete it, re-create it
		G4int count, maxcount;
		maxcount = G4int(EpnTimeH.GetVectorLength());
		//      G4cout << maxcount << G4endl;
		G4double ebins[1024], evals[1024];
		if (maxcount > 1024) {
			G4cout << "Histogram contains more than 1024 bins!" << G4endl;
			G4cout << "Those above 1024 will be ignored" << G4endl;
			maxcount = 1024;
		}
		if (maxcount < 1) {
			G4cout << "Histogram contains less than 1 bin!" << G4endl;
			G4cout << "Redefine the histogram" << G4endl;
			return;
		}
		for (count = 0; count < maxcount; count++) {
			// Read out
			ebins[count] = EpnTimeH.GetLowEdgeEnergy(size_t(count));
			evals[count] = EpnTimeH(size_t(count));
		}

		// Multiply the channels by the nucleon number to give energies
		for (count = 0; count < maxcount; count++) {
			ebins[count] = ebins[count] * Bary;
		}

		// Set Emin and Emax
		Emin = ebins[0];
		if (maxcount > 1)
			Emax = ebins[maxcount - 1];
		else
			Emax = ebins[0];
		// Put time bins into new histogram - UDefTimeH.
		for (count = 0; count < maxcount; count++) {
			UDefTimeH.InsertValues(ebins[count], evals[count]);
		}
		Epnflag = false; //so that you dont repeat this method.
	}
}

//
void ExSPSTimDistribution::ReSetHist(G4String atype) {
	if (atype == "time") {
		UDefTimeH = IPDFTimeH = ZeroPhysVector;
		IPDFTimeExist = false;
		Emin = 0.;
		Emax = 1e30;
	} else if (atype == "arb") {
		ArbTimeH = IPDFArbTimeH = ZeroPhysVector;
		IPDFArbExist = false;
	} else if (atype == "epn") {
		UDefTimeH = IPDFTimeH = ZeroPhysVector;
		IPDFTimeExist = false;
		EpnTimeH = ZeroPhysVector;
	} else {
		G4cout << "Error, histtype not accepted " << G4endl;
	}
}

G4double ExSPSTimDistribution::GenerateOne(G4ParticleDefinition* a) {
	particle_definition = a;
	particle_time = -1.;

	//while ((TimeDisType == "Arb") ? (particle_time < ArbEmin
	//		|| particle_time > ArbEmax) : (particle_time < Emin
	//		|| particle_time > Emax)) 
	//{
		if (Biased) {
			GenerateBiasPowTimes();
		} else {
			if (TimeDisType == "Mono")
				GenerateMonoTime();
			else if (TimeDisType == "Lin")
			{
				GenerateLinearTimes();
			}
			else if (TimeDisType == "Pow")
				GeneratePowTimes();
			else if (TimeDisType == "Exp")
				GenerateExpTimes();
			else if (TimeDisType == "Gauss")
				GenerateGaussTimes();
			else if (TimeDisType == "User")
				GenUserHistTimes();
			else if (TimeDisType == "Arb")
				GenArbPointTimes();
			else if (TimeDisType == "Epn")
				GenEpnHistTimes();
			else
				G4cout << "Error: TimeDisType has unusual value" << G4endl;
		}
	//}
//	if (verbosityLevel >= 1)G4cout << "Time is " << particle_time << G4endl;
	return particle_time;
}

G4double ExSPSTimDistribution::GetProbability(G4double ene) {
	G4double prob = 1.;

	if (TimeDisType == "Lin") {
	  if (prob_norm == 1.) {
	    prob_norm = 0.5*grad*Emax*Emax + cept*Emax - 0.5*grad*Emin*Emin - cept*Emin;
	  }
	  prob = cept + grad * ene;
	  prob /= prob_norm;
	}
	else if (TimeDisType == "Pow") {
	  if (prob_norm == 1.) {
	    if (alpha != -1.) {
	      G4double emina = std::pow(Emin, alpha + 1);
	      G4double emaxa = std::pow(Emax, alpha + 1);
	      prob_norm = 1./(1.+alpha) * (emaxa - emina);
	    } else {
	      prob_norm = std::log(Emax) - std::log(Emin) ;
	    }
	  }
	  prob = std::pow(ene, alpha)/prob_norm;
	}
	else if (TimeDisType == "Exp"){
	  if (prob_norm == 1.) {
	    prob_norm = -Ezero*(std::exp(-Emax/Ezero) - std::exp(Emin/Ezero));
	  }  
	  prob = std::exp(-ene / Ezero);
	  prob /= prob_norm;
	}
	else if (TimeDisType == "Arb") {
	  prob = ArbTimeH.Value(ene);
	  //  prob = ArbEInt->CubicSplineInterpolation(ene);
	  //G4double deltaY;
	  //prob = ArbEInt->PolynomInterpolation(ene, deltaY);
	  if (prob <= 0.) {
	    //G4cout << " Warning:ExSPSTimDistribution::GetProbability: prob<= 0. "<<prob <<" "<<ene << " " <<deltaY<< G4endl;
	    G4cout << " Warning:ExSPSTimDistribution::GetProbability: prob<= 0. "<<prob <<" "<<ene << G4endl;
	    prob = 1e-30;
	  }
	  // already normalised
	}
	else
		G4cout << "Error: TimeDisType not supported" << G4endl;
       
	return prob;
}
