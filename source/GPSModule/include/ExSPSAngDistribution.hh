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
// MODULE:        ExSPSAngDistribution.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
//
// CHANGE HISTORY
// --------------
// 26/10/2004  F Lei
//    Added a "focused" option to allow all primary particles pointing to 
//    a user specified focusing point. 
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//    Based on the G4GeneralParticleSource class in Geant4 v6.0
//
// 2016, Linyan, Modified
// To support ExPos
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// To generate the direction of a primary vertex according to the defined distribution 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// ExSPSAngDistribution ()
//    Constructor: Initializes variables
//
// ~ExSPSAngDistribution ()
//    Destructor: 
//
// void SetAngDistType(G4String)
//    Used to set the type of angular distribution wanted. Arguments
//    are iso, cos, beam  and user for isotropic, cosine-law, beam and user-defined
//    respectively.
//
// void DefineAngRefAxes(G4String, G4ThreeVector)
//    DefineAngRefAxes is used in a similar way as SetPosRot to
//    define vectors, one x' and one in the plane x'y', to create
//    a rotated set of axes for the angular distribution.
//
// void SetMinTheta(G4double)
//    Sets the minimum value for the angle theta.
//
// void SetMinPhi(G4double)
//    Sets the minimum value for phi.  
//
// void SetMaxTheta(G4double)
//    Sets the maximum value for theta.
//
// void SetMaxPhi(G4double)
//    Sets the maximum value for phi.
//
// void UserDefAngTheta(G4ThreeVector)
//    This method allows the user to define a histogram in Theta.
//
// void UserDefAngPhi(G4ThreeVector)
//    This method allows the user to define a histogram in phi.
//
// void GenerateIsotropicFlux()
//    This method generates momentum vectors for particles according
//    to an isotropic distribution.
//
// void GenerateCosineLawFlux()
//    This method generates momentum vectors for particles according
//    to a cosine-law distribution.
//
// void GenerateFocusedFlux()
//    This method generates momentum vectors for particles pointing to
//    an user specified focusing point.
//
// void GenerateUserDefFlux()
//    Controls generation of momentum vectors according to user-defined
//    distributions.
//
// G4double GenerateUserDefTheta()
//    Generates the theta angle according to a user-defined distribution.
//
// G4double GenerateUserDefPhi()
//    Generates phi according to a user-defined distribution.
//
//  void SetBeamSigmaInAngR(G4double);
//    Sets the sigma for 1D beam
//
//  void SetBeamSigmaInAngX(G4double);
//    Sets the first sigma for 2D beam
// 
//  void SetBeamSigmaInAngY(G4double);
//    Sets the second sigma for 2D beam
//
// void SetUserWRTSurface(G4bool)
//    Allows user to have user-defined spectra either with respect to the
//    co-ordinate system (default) or with respect to the surface normal.
//
//  void SetPosDistribution(G4SPSPosDistribution* a) {posDist = a; };
//    Sets the required position generator, required for determining the cosine-law distribution
// 
//  void SetBiasRndm (G4SPSRandomGenerator* a)
//    Sets the biased random number generator
//
//  G4ThreeVector GenerateOne();
//    Generate one random direction
//
//  void ReSetHist(G4String);
//    Re-sets the histogram for user defined distribution
//
//  void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef ExSPSAngDistribution_h
#define ExSPSAngDistribution_h 1

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4DataInterpolation.hh"
#include "G4ParticleMomentum.hh"

#include "ExSPSPosDistribution.hh"
#include "G4SPSRandomGenerator.hh"

class ExSPSAngDistribution 
{
public:
  ExSPSAngDistribution (); 
  ~ExSPSAngDistribution ();
 
  // Angular Distribution Methods
  void SetAngDistType(G4String);
  void DefineAngRefAxes(G4String, G4ThreeVector);
  void SetMinTheta(G4double);
  void SetMinPhi(G4double);
  void SetMaxTheta(G4double);
  void SetMaxPhi(G4double);
  void SetBeamSigmaInAngR(G4double);
  void SetBeamSigmaInAngX(G4double);
  void SetBeamSigmaInAngY(G4double);
  void UserDefAngTheta(G4ThreeVector);
  void UserDefAngPhi(G4ThreeVector);
  void SetFocusPoint(G4ThreeVector);
  inline void SetParticleMomentumDirection
  (G4ParticleMomentum aMomentumDirection)
  { particle_momentum_direction =  aMomentumDirection.unit(); }
  void SetUseUserAngAxis(G4bool);
  void SetUserWRTSurface(G4bool);
  //
  void SetPosDistribution(ExSPSPosDistribution* a) {posDist = a; }
  void SetBiasRndm(G4SPSRandomGenerator* a) {angRndm = a;}
  // method to re-set the histograms
  void ReSetHist(G4String);
  //
  // Set the verbosity level.
  void SetVerbosity(G4int a) {verbosityLevel = a; }
  // some get methods
  G4String GetDistType() { return AngDistType;}
  G4double GetMinTheta() { return MinTheta; }
  G4double GetMaxTheta() { return MaxTheta; }
  G4double GetMinPhi() { return MinPhi; }
  G4double GetMaxPhi() { return MaxPhi; }
  //
  G4ParticleMomentum GenerateOne();
  
private:
  // These methods generate the momentum vectors for the particles.
  void GenerateFocusedFlux();
  void GenerateIsotropicFlux();
  void GenerateCosineLawFlux();
  void GenerateBeamFlux();
  void GeneratePlanarFlux();
  void GenerateUserDefFlux();
  G4double GenerateUserDefTheta();
  G4double GenerateUserDefPhi();

private:

   // Angular distribution variables.
  G4String AngDistType; // String to hold Ang dist type iso, cos, user
  G4ThreeVector AngRef1, AngRef2, AngRef3; // Reference axes for ang dist
  G4double MinTheta, MaxTheta, MinPhi, MaxPhi; // min/max theta/phi
  G4double DR,DX,DY ; // Standard deviations for beam divergence 
  G4double Theta, Phi; // Store these for use with DEBUG
  G4ThreeVector FocusPoint ; // the focusing point in mother coordinates
  G4bool IPDFThetaExist, IPDFPhiExist; // tell whether IPDF histos exist
  G4PhysicsOrderedFreeVector UDefThetaH; // Theta histo data
  G4PhysicsOrderedFreeVector IPDFThetaH; //Cumulative Theta histogram.
  G4PhysicsOrderedFreeVector UDefPhiH; // Phi histo bins
  G4PhysicsOrderedFreeVector IPDFPhiH; // Cumulative phi histogram.
  G4String UserDistType; //String to hold user distributions
  G4bool UserWRTSurface; // G4bool to tell whether user wants distribution wrt
                       // surface normals or co-ordinate system
  G4bool UserAngRef; // Set to true when user defines a new coordinates
  //
  G4ParticleMomentum     particle_momentum_direction;
  //
  ExSPSPosDistribution* posDist;  // need it here for the cosine-law distri 
  G4SPSRandomGenerator* angRndm; // biased random generator

  // Verbosity
  G4int verbosityLevel;
  //
  G4PhysicsOrderedFreeVector ZeroPhysVector ; // for re-set only 
};

#endif
