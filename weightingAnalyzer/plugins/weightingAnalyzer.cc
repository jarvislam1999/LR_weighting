// -*- C++ -*-
//
// Package:    GenStudy/weightingAnalyzer
// Class:      weightingAnalyzer
// 
/**\class weightingAnalyzer weightingAnalyzer.cc GenStudy/weightingAnalyzer/plugins/weightingAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Prakash Thapa
//     Contributor:  Jarvis Lam
//         Created:  Mon, 16 Jul 2018 16:12:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>


//prunedGenParticles: I guess I do not need othere header files because everything is inheritated in GenParticle.h
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <Pythia8/Pythia.h>
#include <cmath>
#include <complex>
#include "TLorentzVector.h"


#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/Math/interface/deltaR.h"
#include <TTree.h>
#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class weightingAnalyzer : public edm::EDAnalyzer {
public:
  explicit weightingAnalyzer(const edm::ParameterSet&);
  ~weightingAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:

  std::auto_ptr<Pythia8::Pythia> fMasterGen;
  //  Pythia8::Pythia *fMasterGen;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  // virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup) override;
  // const reco::Candidate* getBoson( const reco::GenParticleCollection& genParts);

  //  bool isBoson(int pid);
  bool isMuon(int pid);
  //  bool checkBosonStatus(const reco::GenParticleCollection& genParts);
 
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup) override;

  //  double findInvariantMass(const reco::Candidate* particle, const reco::Candidate* antiparticle);

  //Count total number of leptons after event loops
  double muonCountTotal;
  double muonCountStored;
  double antiMuonCountTotal;
  double antiMuonCountStored;

  //Event count
  double eventcount;
  
  double fracLR=1;
  double fracRL=1;
  double sumRL_Plus_LR;
  double sH,tH,uH, sH2, tH2, uH2;

  struct P4Struct {
    float energy,et,eta,phi,pt,mass,theta;
    void fill(const math::XYZTLorentzVector& p4){
      if(p4.Pt()!=0 && p4.Et()!=0){
        energy = p4.E();
        et = p4.Et();
        eta = p4.Eta();
        phi = p4.Phi();
        pt =p4.Pt();
        mass = p4.mag();
        theta = p4.Theta();
      }else clear();
    }
    void clear(){energy=-999;et=-999;eta=-999;phi=-0;pt=-999;mass=-999;theta=-999;}
    static std::string contents(){return "energy/F:et:eta:phi:pt:mass:theta";}
  };

  const reco::Candidate* getLeptonMother(const reco::GenParticle& p, bool second);
  const reco::Candidate* getBosonMother(const reco::GenParticle& p);


  // Histograming variables
  double massInvariant;

  // TTree 
  TTree* tree_;

  TH1F * h_Zmass, *h_Zpt,*h_Zeta,*h_Zphi,*h_Zcharge;
  TH1F *h_muMinusmass,*h_muMinuspt,*h_muMinuseta,*h_muMinusphi,*h_muMinuscharge;
  TH1F *h_muPlusmass,*h_muPluspt,*h_muPluseta,*h_muPlusphi,*h_muPluscharge;
  TH1F *h_dphi,*h_dtheta, *h_dr, *h_thetaMuMinus,*h_thetaMuPlus;
  TH1F *h_massInvar, *h_dimuonPt, *h_dimuonEta, *h_dimuonPhi;
  TH1F *h_cosTheta, *h_tanPhi, *h_csTheta, *h_csPhi;
  TH1F *h_finalInvariantMass;
  TH1F *h_cosThetaPlusInvariantMass, *h_cosThetaMinusInvariantMass;
  TH2F *h_InvariantMassvscosTheta;

  TH2F *h2_pt1_vs_pt2,*h2_eta1_vs_eta2,*h2_phi1_vs_phi2;

  TH1F *h_SHat, *h_T_Hat, *h_U_Hat, *h_THat2,*h_UHat2, *h_fractionLR,*h_fractionRL, *h_sumFrac;
  TH1F *h_alphaRunQED;
  // Delete it after checking it

  TH1F *h_DQuark, *h_UQuark, *h_SQuark, *h_CQuark, *h_BQuark, *h_TQuark;

  P4Struct bosonP4_; // as a sanity check we have the right event...
  P4Struct muMinusP4_;
  P4Struct muPlusP4_;
  int muMinusPID_;
  int muPlusPID_;
  int bosonId_;
  double crossSec, cosTheta, tanPhi, csTheta, csPhi;
  double mCosThetaPlus, mCosThetaMinus;
  double finalInvariantMass;

  // mandelston variables
  double SHat,T_Hat,U_Hat,THat2,UHat2,fractionLR,fractionRL,sumFrac;
  double alphaRunQED;
  //Delete it after checking it
  int DQuark,UQuark,SQuark,CQuark,BQuark,TQuark;


  // Histograms for initial state invariant mass

  TH1F *h_massInvariant;

  // Input tag and Tokens  
  int debug_;
  edm::InputTag prunedGenParticlesTag_;
  int decayParticlePID_;
  edm::InputTag genInfoProduct_;
  edm::EDGetTokenT<GenRunInfoProduct> genInfoProductToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;

  //edm::InputTag genInfoProduct_;
  //edm::EDGetTokenT<GenRunInfoProduct> genInfoProductToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;

  edm::InputTag genEventInfoProduct_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProductToken_;


  // ----------member data ---------------------------
};

//
// constants, enums and typedefs

//

//
// static data member definitions
//

//
// constructors and destructor
//
weightingAnalyzer::weightingAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  debug_=iConfig.getParameter<int>("debug");
  prunedGenParticlesTag_=iConfig.getParameter<edm::InputTag>("prunedGenParticlesTag");

  genInfoProduct_ = iConfig.getParameter<edm::InputTag>("genInfoProduct");
  
  
  genInfoProductToken_ = consumes<GenRunInfoProduct,edm::InRun>(genInfoProduct_);
  prunedGenParticlesToken_= consumes<reco::GenParticleCollection>(prunedGenParticlesTag_);
 
  genInfoProduct_ = iConfig.getParameter<edm::InputTag>("genInfoProduct");
  genInfoProductToken_ = consumes<GenRunInfoProduct,edm::InRun>(genInfoProduct_);
  genEventInfoProduct_ = iConfig.getParameter<edm::InputTag>("genEventInfoProduct");
  genEventInfoProductToken_ = consumes<GenEventInfoProduct,edm::InEvent>(genEventInfoProduct_);

  
}


weightingAnalyzer::~weightingAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
weightingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventcount += 1;
  std::cout << "==========================================================================================================================================================================================================================================================================================================================================================================================" << std::endl;

  std::cout <<"EVENT " <<eventcount<< std::endl;

  using namespace edm;

  edm::Handle<reco::GenParticleCollection> prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticlesHandle);
  const reco::GenParticleCollection& prunedGenParticlesCollection = *prunedGenParticlesHandle;
  // Object declaration using additional infromation
  edm::Handle<GenEventInfoProduct> genEventInfoProductHandle;
  iEvent.getByToken(genEventInfoProductToken_, genEventInfoProductHandle);
  const GenEventInfoProduct& genEventInfoProduct = *genEventInfoProductHandle;


  // const reco::Candidate* firstMother;
  // const reco::Candidate* muonMother;

  bosonId_=0;
  bosonP4_.clear();
  muMinusP4_.clear();
  muPlusP4_.clear();
  muMinusPID_=0;
  muPlusPID_=0;

  //  const reco::Candidate* daughter1;
  //  const reco::Candidate* daughter2;

  //Store boson and initial leptons in these variables
  const reco::Candidate* boson;
  const reco::Candidate* leptons;

  double muonpdgId;
  double muonpt;
  double muoneta;
  double muonphi;
  double muoncharge;

  double antiMuonpdgId;
  double antiMuonpt;
  double antiMuoneta;
  double antiMuonphi;
  double antiMuoncharge;

  double quarkpx,antiquarkpx, leptonPluspx, leptonMinuspx;
  double quarkpy,antiquarkpy, leptonPluspy, leptonMinuspy;
  double quarkpz,antiquarkpz, leptonPluspz, leptonMinuspz;
  double quarkenergy,antiquarkenergy, leptonPlusenergy, leptonMinusenergy;

  //Store final state leptons in these variables
  const reco::Candidate* muMinus;
  const reco::Candidate* muPlus;

  math::XYZTLorentzVectorD dimuon;
  double dimuonPx, dimuonPy, dimuonPz, dimuonPt, pseudorapidity, Phi, mu1Energy, mu2Energy, dimuonQ;
  double thetaCos, thetaCS, phiTan, phiCS;
  double muPlusKPlus, muPlusKMinus, muMinusKPlus, muMinusKMinus, invariantK;

  //  int  motherID;
  // Use mother function

  //  const reco::Candidate* mother1;
  //  const reco::Candidate* mother2;

  //Initialize lambdaSquare and etaLR using default value
  //double qCLambda2 = pow(fMasterGen->settings.parm("ContactInteractions:Lambd\a"),2);
  double qCLambda=16000.0;
  double qCLambda2=qCLambda*qCLambda;
  int qCetaLR =-1;

  if (debug_ >5) {
  std::cout << "######### PRINTING OUT STEVE'S CODE VARIABLES BEFORE THE FOR LOOP:##########" << std::endl;


  std::cout << "Test_lambdaSquare: default value for lambda is 1000 GeV"<< std::endl;
  std::cout <<"Test_lambdaSquare: " <<qCLambda2 << std::endl;
  //int qCetaLR = fMasterGen->settings.mode("ContactInteractions:etaLR");
  std::cout << "Test_etaLR: defalut value for etaLR is 0.00 "<< std::endl;
  std::cout <<"Test_etaLR: " << qCetaLR << std::endl;
  };
  // Running coupling constat-> changes event by event and accessible in miniAOD
  double genEventalphaEM =genEventInfoProduct.alphaQED();
  double alpEM=genEventalphaEM;
  // std::cout << " ########################################" <<std::endl;
  // std::cout << " genEventalphaEM is  alpEM :  " << alpEM <<std::endl;

  //Variables to get quarks and initial leptons in for loop
  const reco::Candidate* quark;
  const reco::Candidate* antiquark;
  const reco::Candidate* leptonPlus;
  const reco::Candidate* leptonMinus;

  if (debug_ > 5){
    std::cout << "#### IMPLEMENTATION OF STEVE'S CODE::STOP HERE AND GO AFTER THE FOR LOOP ##########" <<std::endl;
  };

  //Begin particle loop
  for(const reco::GenParticle& particle : prunedGenParticlesCollection){
     //Continue to next particle if not quark, protons or leptons
    if (abs(particle.pdgId())!= 11 && abs(particle.pdgId()) != 13)
      {
	continue;
      };

    // Check for final state lepton
    if ((particle.pdgId() == 13 || particle.pdgId() == 11) &&  particle.status() == 1){
 
      //Test Mother  Particle tree

      int daughterNumber=particle.numberOfDaughters();
      if (debug_ >5){
	std::cout<< "Number of daughter : " <<  daughterNumber << std::endl;
      };
      int MotherNumber=particle.numberOfMothers();
      if (debug_ >5){
	std::cout<< "Number of Mother : " << MotherNumber  << std::endl;
      };

      // Function call to get quark and boson from final state lepton
      boson = getBosonMother(particle);
      leptons =getLeptonMother(particle, false);// leptons here are quarks
      if (!leptons)
	{
	  std::cout << "Lepton pointer is null!" << std::endl;
	  continue;
	}

      if (leptons != nullptr){
	//Trying to store particle into final state lepton by changing its type from genParticle into Reco::Candidate
	muonCountTotal = muonCountTotal +1;
	muMinus = particle.clone();	
	
	if (muMinus == nullptr){
	  std::cout << "muMinus is null!" << std::endl;
	}
	     
      }

      if (boson != nullptr){
	std::cout<< "Z boson is present!" << std::endl;
	std::cout << "- BOSON INFORMATION: " << std::endl;
	std::cout <<"boson ID: " << boson->pdgId() << std::endl;
	std::cout <<"boson Energy: " << boson->energy()<< std::endl;
	std::cout <<"boson px: " <<boson->px() <<std::endl;
	std::cout <<"boson py: " <<boson->py() << std::endl;
	std::cout <<"boson pz: " << boson->pz() <<std::endl;
	std::cout <<"boson pt: " << boson->pt() <<std::endl;
	std::cout <<"boson eta: " << boson->eta() <<std::endl;
	std::cout <<"boson phi: " << boson->phi() <<std::endl;
	std::cout << std::endl;	
      }
      //Assigning quarks to leptons
      quark = leptons;

      debug_ =1;

      //Printing quark information  
      if (debug_ >0) {
	std::cout << std::endl;
	std::cout << " - QUARK INFORMATION: " << std::endl;
	std::cout <<"quark ID: " << quark->pdgId()<< std::endl;
	std::cout <<"quark Energy: " << quark->energy()<< std::endl;
	std::cout <<"quark px: " <<quark->px() <<std::endl;
	std::cout <<"quark py: " <<quark->py() << std::endl;
	std::cout <<"quark pz: " << quark->pz() <<std::endl;
	std::cout <<"quark pt: " << quark->pt() <<std::endl;
	std::cout <<"quark eta: " << quark->eta() <<std::endl;
	std::cout <<"quark phi: " << quark->phi() <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
      };

      quarkenergy =  quark->energy();
      quarkpx = quark->px();
      quarkpy = quark->py();
      quarkpz =  quark->pz();

      const reco::Candidate* otherLepton = getLeptonMother(particle, true);
      if (!otherLepton)
	{
	  std::cout << "otherLepton pointer is null!" << std::endl;
	  continue;
	}

      antiquark = otherLepton;
      antiquarkenergy =  antiquark->energy();
      antiquarkpx = antiquark->px();
      antiquarkpy = antiquark->py();
      antiquarkpz =  antiquark->pz();

      if (debug_ > 0) {
	//Printing out antiquark information
	std::cout << std::endl;
	std::cout << "- ANTIQUARK INFORMATION: " << std::endl;
	std::cout <<"antiquark ID: " << antiquark->pdgId() << std::endl;
	std::cout <<"antiquark Energy: " << antiquark->energy()<< std::endl;
	std::cout <<"antiquark px: " <<antiquark->px() <<std::endl;
	std::cout <<"antiquark py: " <<antiquark->py() << std::endl;
	std::cout <<"antiquark pz: " << antiquark->pz() <<std::endl;
	std::cout <<"antiquark pt: " << antiquark->pt() <<std::endl;
	std::cout <<"antiquark eta: " << antiquark->eta() <<std::endl;
	std::cout <<"antiquark phi: " << antiquark->phi() <<std::endl;
	std::cout << std::endl;	
	std::cout << std::endl;
      };
      
      //Get initial lepton for contact interaction
      if (boson == nullptr){      
	if (leptons->daughter(0)->pdgId()>0)
	  {
	    leptonMinus = leptons->daughter(0);
	    muonpdgId= (leptons->daughter(0))->pdgId();
	    muonpt= (leptons->daughter(0))->pt();
	    muoneta= (leptons->daughter(0))->eta();
	    muonphi= (leptons->daughter(0))->phi();
	    muoncharge= (leptons->daughter(0))->charge();
	    muonCountStored = muonCountStored +1;

	    leptonMinusenergy = leptonMinus->energy();
	    leptonMinuspx = leptonMinus->px();
	    leptonMinuspy = leptonMinus->py();
	    leptonMinuspz = leptonMinus->pz();
	    
	    if (leptons->daughter(1) !=nullptr){
	      leptonPlus = leptons->daughter(1);	      
	      antiMuonpdgId = (leptons->daughter(1))->pdgId();
	      antiMuonpt = (leptons->daughter(1))->pt();
	      antiMuoneta = (leptons->daughter(1))->eta();
	      antiMuonphi = (leptons->daughter(1))->phi();
	      antiMuoncharge = (leptons->daughter(1))->charge();
	      antiMuonCountStored = antiMuonCountStored +1;
	      
	      leptonPlusenergy = leptonPlus->energy();
	      leptonPluspx = leptonPlus->px();
	      leptonPluspy = leptonPlus->py();
	      leptonPluspz = leptonPlus->pz();
	    }

	  }
	else
	  {
	    leptonPlus = leptons->daughter(0);
	    antiMuonpdgId = (leptons->daughter(0))->pdgId();
	    antiMuonpt = (leptons->daughter(0))->pt();
	    antiMuoneta = (leptons->daughter(0))->eta();
	    antiMuonphi = (leptons->daughter(0))->phi();
	    antiMuoncharge = (leptons->daughter(0))->charge();
	    antiMuonCountStored = antiMuonCountStored +1;
	    
	    leptonPlusenergy = leptonPlus->energy();
	    leptonPluspx = leptonPlus->px();
	    leptonPluspy = leptonPlus->py();
	    leptonPluspz = leptonPlus->pz();
	    
	    if (leptons->daughter(1) !=nullptr){
	      leptonMinus = leptons->daughter(1);
	      muonpdgId= (leptons->daughter(1))->pdgId();
	      muonpt= (leptons->daughter(1))->pt();
	      muoneta= (leptons->daughter(1))->eta();
	      muonphi= (leptons->daughter(1))->phi();
	      muoncharge= (leptons->daughter(1))->charge();
	      muonCountStored = muonCountStored +1;

	      leptonMinusenergy = leptonMinus->energy();
	      leptonMinuspx = leptonMinus->px();
	      leptonMinuspy = leptonMinus->py();
	      leptonMinuspz = leptonMinus->pz();
	    }
	  };
      };

      //Initial lepton for Drell-Yan
      if (boson != nullptr){      
	if (boson->daughter(0)->pdgId()>0)
	  {
	    leptonMinus = boson->daughter(0);
	    muonpdgId= (boson->daughter(0))->pdgId();
	    muonpt= (boson->daughter(0))->pt();
	    muoneta= (boson->daughter(0))->eta();
	    muonphi= (boson->daughter(0))->phi();
	    muoncharge= (boson->daughter(0))->charge();
	    muonCountStored = muonCountStored +1;

	    leptonMinusenergy = leptonMinus->energy();
	    leptonMinuspx = leptonMinus->px();
	    leptonMinuspy = leptonMinus->py();
	    leptonMinuspz = leptonMinus->pz();
	    
	    if (boson->daughter(1) !=nullptr){
	      leptonPlus = boson->daughter(1);	      
	      antiMuonpdgId = (boson->daughter(1))->pdgId();
	      antiMuonpt = (boson->daughter(1))->pt();
	      antiMuoneta = (boson->daughter(1))->eta();
	      antiMuonphi = (boson->daughter(1))->phi();
	      antiMuoncharge = (boson->daughter(1))->charge();
	      antiMuonCountStored = antiMuonCountStored +1;
	      
	      leptonPlusenergy = leptonPlus->energy();
	      leptonPluspx = leptonPlus->px();
	      leptonPluspy = leptonPlus->py();
	      leptonPluspz = leptonPlus->pz();
	    }

	  }
	else
	  {
	    leptonPlus = boson->daughter(0);
	    antiMuonpdgId = (boson->daughter(0))->pdgId();
	    antiMuonpt = (boson->daughter(0))->pt();
	    antiMuoneta = (boson->daughter(0))->eta();
	    antiMuonphi = (boson->daughter(0))->phi();
	    antiMuoncharge = (boson->daughter(0))->charge();
	    antiMuonCountStored = antiMuonCountStored +1;
	    
	    leptonPlusenergy = leptonPlus->energy();
	    leptonPluspx = leptonPlus->px();
	    leptonPluspy = leptonPlus->py();
	    leptonPluspz = leptonPlus->pz();
	    
	    if (boson->daughter(1) !=nullptr){
	      leptonMinus = boson->daughter(1);
	      muonpdgId= (boson->daughter(1))->pdgId();
	      muonpt= (boson->daughter(1))->pt();
	      muoneta= (boson->daughter(1))->eta();
	      muonphi= (boson->daughter(1))->phi();
	      muoncharge= (boson->daughter(1))->charge();
	      muonCountStored = muonCountStored +1;

	      leptonMinusenergy = leptonMinus->energy();
	      leptonMinuspx = leptonMinus->px();
	      leptonMinuspy = leptonMinus->py();
	      leptonMinuspz = leptonMinus->pz();
	    }
	  };
      };

    }//end if muon test

    //Check for antiquark
    if ((particle.pdgId() == -13  || particle.pdgId() == -11) &&  particle.status() == 1){
     
      boson = getBosonMother(particle);
      leptons =getLeptonMother(particle, false);

      if (!leptons)
	{
	  std::cout << "Lepton pointer is null!" << std::endl;
	  continue;
	}

      if (leptons != nullptr){
	antiMuonCountTotal = antiMuonCountTotal + 1;
	muPlus = particle.clone();
	
	if (!muPlus){
	  std::cout << "muPlus is null!" << std::endl;
	}
      }

      quark = leptons;
      const reco::Candidate* otherLepton = getLeptonMother(particle, true);

      if (!otherLepton)
	{
	  std::cout << "Lepton pointer is null!" << std::endl;
	  continue;
	}

      antiquark = otherLepton;
      //Store quark information
      quarkenergy =  quark->energy();
      quarkpx = quark->px();
      quarkpy = quark->py();
      quarkpz =  quark->pz();

      antiquarkenergy =  antiquark->energy();
      antiquarkpx = antiquark->px();
      antiquarkpy = antiquark->py();
      antiquarkpz =  antiquark->pz();

      if (leptons == nullptr){
	std::cout << "leptons is null" << std::endl;
      }

      //Initial leptons for contact interaction
      if (boson == nullptr){
	if (leptons->daughter(0)->pdgId()>0)
	  {
	    leptonMinus = leptons->daughter(0);
	    muonpdgId= (leptons->daughter(0))->pdgId();
	    muonpt= (leptons->daughter(0))->pt();
	    muoneta= (leptons->daughter(0))->eta();
	    muonphi= (leptons->daughter(0))->phi();
	    muoncharge= (leptons->daughter(0))->charge();
	    muonCountStored = muonCountStored +1;
	    
	    leptonMinusenergy = leptonMinus->energy();
	    leptonMinuspx = leptonMinus->px();
	    leptonMinuspy = leptonMinus->py();
	    leptonMinuspz = leptonMinus->pz();
	    
	    if (leptons->daughter(1) !=nullptr){
	      leptonPlus = leptons->daughter(1);
	      antiMuonpdgId = (leptons->daughter(1))->pdgId();
	      antiMuonpt = (leptons->daughter(1))->pt();
	      antiMuoneta = (leptons->daughter(1))->eta();
	      antiMuonphi = (leptons->daughter(1))->phi();
	      antiMuoncharge = (leptons->daughter(1))->charge();
	      antiMuonCountStored = antiMuonCountStored +1;
	      
	      leptonPlusenergy = leptonPlus->energy();
	      leptonPluspx = leptonPlus->px();
	      leptonPluspy = leptonPlus->py();
	    leptonPluspz = leptonPlus->pz();
	    }
	  }
	else
	  {
	    leptonPlus = leptons->daughter(0);
	    antiMuonpdgId = (leptons->daughter(0))->pdgId();
	    antiMuonpt = (leptons->daughter(0))->pt();
	    antiMuoneta = (leptons->daughter(0))->eta();
	    antiMuonphi = (leptons->daughter(0))->phi();
	    antiMuoncharge = (leptons->daughter(0))->charge();
	    antiMuonCountStored = antiMuonCountStored +1;
	    
	    leptonPlusenergy = leptonPlus->energy();
	    leptonPluspx = leptonPlus->px();
	    leptonPluspy = leptonPlus->py();
	    leptonPluspz = leptonPlus->pz();
	    
	    if (leptons->daughter(1) !=nullptr){
	      leptonMinus = leptons->daughter(1);
	      muonpdgId= (leptons->daughter(1))->pdgId();
	      muonpt= (leptons->daughter(1))->pt();
	      muoneta= (leptons->daughter(1))->eta();
	      muonphi= (leptons->daughter(1))->phi();
	      muoncharge= (leptons->daughter(1))->charge();
	      muonCountStored = muonCountStored +1;
	      
	      leptonMinusenergy = leptonMinus->energy();
	      leptonMinuspx = leptonMinus->px();
	      leptonMinuspy = leptonMinus->py();
	      leptonMinuspz = leptonMinus->pz();
	    }
	    
	  }
      };

      //Initial lepton for Drell_Yan
      if (boson != nullptr){      
	if (boson->daughter(0)->pdgId()>0)
	  {
	    leptonMinus = boson->daughter(0);
	    muonpdgId= (boson->daughter(0))->pdgId();
	    muonpt= (boson->daughter(0))->pt();
	    muoneta= (boson->daughter(0))->eta();
	    muonphi= (boson->daughter(0))->phi();
	    muoncharge= (boson->daughter(0))->charge();
	    muonCountStored = muonCountStored +1;

	    leptonMinusenergy = leptonMinus->energy();
	    leptonMinuspx = leptonMinus->px();
	    leptonMinuspy = leptonMinus->py();
	    leptonMinuspz = leptonMinus->pz();
	    
	    if (boson->daughter(1) !=nullptr){
	      leptonPlus = boson->daughter(1);	      
	      antiMuonpdgId = (boson->daughter(1))->pdgId();
	      antiMuonpt = (boson->daughter(1))->pt();
	      antiMuoneta = (boson->daughter(1))->eta();
	      antiMuonphi = (boson->daughter(1))->phi();
	      antiMuoncharge = (boson->daughter(1))->charge();
	      antiMuonCountStored = antiMuonCountStored +1;
	      
	      leptonPlusenergy = leptonPlus->energy();
	      leptonPluspx = leptonPlus->px();
	      leptonPluspy = leptonPlus->py();
	      leptonPluspz = leptonPlus->pz();
	    }
	    if (debug_ > 5){
	      std::cout << "muon pdgID :" <<  muonpdgId << std::endl;
	      std::cout << "muon pt :" <<  muonpt << std::endl;
	      std::cout << "muon eta :" <<  muoneta << std::endl;
	      std::cout << "muon phi :" <<  muonphi << std::endl;
	      std::cout << "muon charge :" <<  muoncharge << std::endl;
	    }; 

	  }
	else
	  {
	    leptonPlus = boson->daughter(0);
	    antiMuonpdgId = (boson->daughter(0))->pdgId();
	    antiMuonpt = (boson->daughter(0))->pt();
	    antiMuoneta = (boson->daughter(0))->eta();
	    antiMuonphi = (boson->daughter(0))->phi();
	    antiMuoncharge = (boson->daughter(0))->charge();
	    antiMuonCountStored = antiMuonCountStored +1;
	    
	    leptonPlusenergy = leptonPlus->energy();
	    leptonPluspx = leptonPlus->px();
	    leptonPluspy = leptonPlus->py();
	    leptonPluspz = leptonPlus->pz();
	    
	    if (boson->daughter(1) !=nullptr){
	      leptonMinus = boson->daughter(1);
	      muonpdgId= (boson->daughter(1))->pdgId();
	      muonpt= (boson->daughter(1))->pt();
	      muoneta= (boson->daughter(1))->eta();
	      muonphi= (boson->daughter(1))->phi();
	      muoncharge= (boson->daughter(1))->charge();
	      muonCountStored = muonCountStored +1;

	      leptonMinusenergy = leptonMinus->energy();
	      leptonMinuspx = leptonMinus->px();
	      leptonMinuspy = leptonMinus->py();
	      leptonMinuspz = leptonMinus->pz();
	    }
	  };
      };

      if (debug_ > 5){
	std::cout << "antimuon pdgID :" <<  antiMuonpdgId << std::endl;
	std::cout << "antimuon pt :" <<  antiMuonpt << std::endl;
	std::cout << "antimuon eta :" <<  antiMuoneta << std::endl;
	std::cout << "antimuon phi :" <<  antiMuonphi << std::endl;
	std::cout << "antimuon charge :" <<  antiMuoncharge << std::endl;
      }; 

    }//end if antimuon

   
    //tree_->Fill();  
  }//end for loop

	std::cout << std::endl;
	//Check for right leptons
	if (boson == nullptr){
	  if ((leptonPlus->status() != 23) || (leptonMinus->status() != 23)){
	    std::cout << "WRONG INITIAL LEPTON!" << std::endl;
	  }
	  if ((leptonPlus->status() == 23) && (leptonMinus->status() == 23)) {
	    std::cout<< "RIGHT INITIAL LEPTON" << std::endl;
	  };
	};

	//Print out initial leptons
  std::cout << " - INITIAL LEPTON INFORMATION:  " << std::endl;
  std::cout << std::endl;
  std::cout << "+ LEPTON:" <<std::endl;
  std::cout <<"leptonMinus ID:"     << leptonMinus->pdgId() <<std::endl;
  std::cout <<"leptonMinus status:"     << leptonMinus->status() <<std::endl;
  std::cout <<"leptonMinus Energy: " << leptonMinus->energy()<< std::endl;
  std::cout <<"leptonMinus px: " <<leptonMinus->px() <<std::endl;
  std::cout <<"leptonMinus py: " <<leptonMinus->py() << std::endl;
  std::cout <<"leptonMinus pz: " << leptonMinus->pz() <<std::endl;
  std::cout <<"leptonMinus pt: " << leptonMinus->pt() <<std::endl;
  std::cout <<"leptonMinus eta: " << leptonMinus->eta() <<std::endl;
  std::cout <<"leptonMinus phi: " << leptonMinus->phi() <<std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "+ ANTILEPTON" << std::endl;
  std::cout <<"leptonPlus ID:" << leptonPlus->pdgId() <<std::endl;
  std::cout <<"leptonMinus status:" << leptonPlus->status() <<std::endl;      
  std::cout <<"leptonPlus Energy: " << leptonPlus->energy()<< std::endl;
  std::cout <<"leptonPlus px: " <<leptonPlus->px() <<std::endl;
  std::cout <<"leptonPlus py: " <<leptonPlus->py() << std::endl;
  std::cout <<"leptonPlus pz: " << leptonPlus->pz() <<std::endl;
  std::cout <<"leptonPlus pt: " << leptonPlus->pt() <<std::endl;
  std::cout <<"leptonPlus eta: " << leptonPlus->eta() <<std::endl;
  std::cout <<"leptonPlus phi: " << leptonPlus->phi() <<std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  
  leptonMinusenergy = leptonMinus->energy();
  leptonMinuspx = leptonMinus->px();
  leptonMinuspy = leptonMinus->py();
  leptonMinuspz = leptonMinus->pz();

  // std::cout << "##### FIGURE OUT HOW TO DEFINE A FOUR VECTOR IN FUTURE ######## " << std::endl;


      //Mandelston varibles goes here

  /*
      std::vector<double> tempShat;
      tempShat.push_back(quark->px()+antiquark->px());
      tempShat.push_back(quark->py()+antiquark->py());
      tempShat.push_back(quark->pz()+antiquark->pz());
      tempShat.push_back(quark->energy()+antiquark->energy());

      sH= tempShat[3]*tempShat[3] -tempShat[0]*tempShat[0]-tempShat[1]*tempShat{1]-tempShat[2]*tempShat[2];
      sH2 =sH*sH;

      // test using final state
      std::vector<double> test_lepton_Shat;
      test_lepton_Shat.push_back(leptonPlus->px()+leptonMinus->px());
      test_lepton_Shat.push_back(leptonPlus->py()+leptonMinus->py());
      test_lepton_Shat.push_back(leptonPlus->pz()+leptonMinus->pz());
      test_lepton_Shat.push_back(leptonPlus->energy()+leptonMinus->energy());

      double test_finalState_sHat=test_lepton_Shat[3]*test_lepton_Shat[3] -test\
_lepton_Shat[0]*test_lepton_Shat[0]-test_lepton_Shat[1]*test_lepton_Shat[1]-tes\
	t_lepton_Shat[2]*tempShat[2];

      std::vector<double> tempThat;
      tempThat.push_back(quark->px()-leptonPlus->px());
      tempThat.push_back(quark->py()-leptonPlus->py());
      tempThat.push_back(quark->pz()-leptonPlus->pz());
      tempThat.push_back(quark->energy()-leptonPlus->energy());
      tH= (tempThat[3]*tempThat[3] -tempThat[0]*tempThat[0]-tempThat[1]*tempTha\
	   t[1]-tempThat[2]*tempThat[2]);
      tH2=tH*tH;
      //double tH2=(tempThat[3]*tempThat[3] -tempThat[0]*tempThat[0]-tempThat[1]*tempThat[1]-tempThat[2]*tempThat[2])*(tempThat[3]*tempThat[3] -tempThat[0]*te					      mpThat[0  -tempThat[1]*tempThat[1]-tempThat[2]*tempThat[2]);
      
      
      std::vector<double> tempUhat;
      tempUhat.push_back(quark->px()-leptonMinus->px());
      tempUhat.push_back(quark->py()-leptonMinus->py());
      tempUhat.push_back(quark->pz()-leptonMinus->pz());
      tempUhat.push_back(quark->energy()-leptonMinus->energy());
      uH= (tempUhat[3]*tempUhat[3] -tempUhat[0]*tempUhat[0]-tempUhat[1]*tempUhat[1]-tempUhat[2]*tempUhat[2]);
      uH2=uH*uH;
  

      if (debug_ > 0) {
	std::cout << "###################   MANDELSTAM VARIABLES  #######################################################" << std::endl;
      };

      if (debug_ > 5){
	std::cout << "sHat I guess:" << sH << std::endl;
	std::cout << "tHat I guess:" << tH << std::endl;
	std::cout << "uHat I guess:" << uH << std::endl;
	std::cout << "sHat square I guess:" << sH2 << std::endl;
	std::cout << "tHat square I guess:" << tH2 << std::endl;
	std::cout << "uHat square I guess:" << uH2 << std::endl;
      };

  */

  /*
      TLorentzVector s_hat_4vector(quark->px() + antiquark->px(), quark->py() +antiquark->py(), quark->pz() + antiquark->pz(), quark->energy() + antiquark->energy());

      double s_hat = s_hat_4vector.E()*s_hat_4vector.E() - s_hat_4vector.Px()*s_hat_4vector.Px() - s_hat_4vector.Py()*s_hat_4vector.Py() - s_hat_4vector.Pz()*s_hat_4vector.Pz();

      TLorentzVector t_hat_4vector(quark->px() - leptonPlus->px(), quark->py() - leptonPlus->py(), quark->pz() - leptonPlus->pz(), quark->energy() - leptonPlus->energy());

      double t_hat = t_hat_4vector.E()*t_hat_4vector.E() - t_hat_4vector.Px()*t_hat_4vector.Px() - t_hat_4vector.Py()*t_hat_4vector.Py() - t_hat_4vector.Pz()*t_hat_4vector.Pz();

      TLorentzVector u_hat_4vector(quark->px() - leptonMinus->px(), quark->py()- leptonMinus->py(), quark->pz() - leptonMinus->pz(), quark->energy() - leptonMinus->energy());

      double u_hat = u_hat_4vector.E()*u_hat_4vector.E() - u_hat_4vector.Px()*u_hat_4vector.Px() - u_hat_4vector.Py()*u_hat_4vector.Py() - u_hat_4vector.Pz()*u_hat_4vector.Pz();
  */
  
  //Get Mandelstam variables
      TLorentzVector s_hat_4vector(quarkpx + antiquarkpx, quarkpy +antiquarkpy, quarkpz + antiquarkpz, quarkenergy + antiquarkenergy);

      double s_hat = s_hat_4vector.E()*s_hat_4vector.E() - s_hat_4vector.Px()*s_hat_4vector.Px() - s_hat_4vector.Py()*s_hat_4vector.Py() - s_hat_4vector.Pz()*s_hat_4vector.Pz();
      // std::cout << "leptonPluspx = " << leptonPluspx << std::endl;
      //std::cout << "leptonPluspy = " << leptonPluspy << std::endl;
      //std::cout << "leptonPluspz = " << leptonPluspz << std::endl;

      TLorentzVector t_hat_4vector(quarkpx - leptonMinuspx, quarkpy - leptonMinuspy, quarkpz - leptonMinuspz, quarkenergy - leptonMinusenergy);

      double t_hat = t_hat_4vector.E()*t_hat_4vector.E() - t_hat_4vector.Px()*t_hat_4vector.Px() - t_hat_4vector.Py()*t_hat_4vector.Py() - t_hat_4vector.Pz()*t_hat_4vector.Pz();
      //std::cout << "leptonMinuspx = " << leptonMinuspx << std::endl;
      //std::cout << "leptonMinuspy = " << leptonMinuspy << std::endl;
      //std::cout << "leptonMinuspz = " << leptonMinuspz << std::endl;

      TLorentzVector u_hat_4vector(quarkpx - leptonPluspx, quarkpy- leptonPluspy, quarkpz - leptonPluspz, quarkenergy - leptonPlusenergy);

      double u_hat = u_hat_4vector.E()*u_hat_4vector.E() - u_hat_4vector.Px()*u_hat_4vector.Px() - u_hat_4vector.Py()*u_hat_4vector.Py() - u_hat_4vector.Pz()*u_hat_4vector.Pz();

      sH = s_hat;
      tH = t_hat;
      uH = u_hat;
      sH2 = sH*sH;
      tH2 = tH*tH;
      uH2 = uH*uH ;
      debug_ = 1;

      if (debug_ > 0){
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "##### MANDELSTAM VARIABLES: " <<std::endl;
	std::cout << "sHat :" << sH << std::endl;
	std::cout << "tHat:" << tH << std::endl;
	std::cout << "uHat:" <<uH << std::endl;
	std::cout << std::endl;
	std::cout << "sHat square :" << sH2 << std::endl;
	std::cout << "tHat square :" << tH2 << std::endl;
	std::cout << "uHat square :" << uH2 << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
      };

      //std::cout << "############### IMPLEMENT STEVE'S CODE HERE #############\      ###" << std::endl;

      //Begin calculating fracLR and fracRL
      std::complex<double> I(0.0, 1.0);
      // Complex amplitudes.
      std::complex<double> meLL(0., 0.);
      std::complex<double> meRR(0., 0.);
      std::complex<double> meLR(0., 0.);
      std::complex<double> meRL(0., 0.);
      std::complex<double> meLR_SM(0., 0.);
      std::complex<double> meRL_SM(0., 0.);

      int quarkId=quark->pdgId();
      int idAbs=quarkId;
      // int idAbs=5;
      int leptonMinusId=leptonMinus->pdgId();
      int idNew=leptonMinusId;

      if (debug_ > 5){
	std::cout << "############### TEST AFTER A FOR LOOP:  ################"		  << std::endl;
	std::cout << "############### Test_etaLR: defalut value for etaLR is 0.00  #############################"<< std::endl;
	std::cout <<"Test_etaLR after for loop : " << qCetaLR << std::endl;
	std::cout << "############################## Test_alpEM #############################################"<< std::endl;
	std::cout << "Test_alpEM  after for loop: " << alpEM << std::endl;
	std::cout << "############### Test_lambdaSquare after for loop : #############################"<< std::endl;
	std::cout <<"Test_lambdaSquare: " <<qCLambda2 << std::endl;
	std::cout << "############### TESST HEPMC TO GET FINE STRUCTURE CONSTANT  ################" << std::endl;

	std::cout << "## ARGUAMENTS (idAbs AND idNew) DO NOT TAKE NEGATIVE VALUE:LETS SEE WHAT HAPPENS !!! ### " << std::endl;
	std::cout << "#fMasterGen->couplingsPtr->vf(-2)&fMasterGen->couplingsPtr->vf(-13)#" <<fMasterGen->couplingsPtr->vf(-2) << " " << fMasterGen->couplingsPtr->vf(-13)<< std::endl;
	std::cout << "#fMasterGen->couplingsPtr->vf(-3)&fMasterGen->couplingsPtr->vf(-11)#" <<fMasterGen->couplingsPtr->vf(-3) << " " << fMasterGen->couplingsPtr->vf(-11)<< std::endl;
	std::cout << "#fMasterGen->couplingsPtr->af(-2)&fMasterGen->couplingsPtr->af(-13)#" <<fMasterGen->couplingsPtr->af(-2) << " " << fMasterGen->couplingsPtr->af(-13) << std::endl;
	std::cout << "#fMasterGen->couplingsPtr->af(-1)&fMasterGen->couplingsPtr->af(-11)#" <<fMasterGen->couplingsPtr->af(-1) << " " << fMasterGen->couplingsPtr->af(-11) << std::endl;
	std::cout << "## fMasterGen->couplingsPtr->vf(13) ### " << fMasterGen->couplingsPtr->vf(13)<< std::endl;
	std::cout << "## fMasterGen->couplingsPtr->af(13) ### " <<fMasterGen->couplingsPtr->af(13) << std::endl;
	std::cout << "## fMasterGen->couplingsPtr->vf(-2) ### " << fMasterGen->couplingsPtr->vf(-2)<< std::endl;
	std::cout << "## fMasterGen->couplingsPtr->af(-2) ### " <<fMasterGen->couplingsPtr->af(-2) << std::endl;
	std::cout << "## fMasterGen->couplingsPtr->vf(2) ### " << fMasterGen->couplingsPtr->vf(2)<< std::endl;
	std::cout << "## fMasterGen->couplingsPtr->af(2) ### " <<fMasterGen->couplingsPtr->af(2) << std::endl;

      };

      //Process name
      // Incoming quarks
      double tmPgvf = 0.25 * fMasterGen->couplingsPtr->vf(idAbs);
      double tmPgaf = 0.25 * fMasterGen->couplingsPtr->af(idAbs);
      //Outgoing fermions
      double tmPgvl = 0.25 * fMasterGen->couplingsPtr->vf(idNew);
      double tmPgal = 0.25 * fMasterGen->couplingsPtr->af(idNew);
      double tmPgLf = tmPgvf + tmPgaf;
      double tmPgLl = tmPgvl + tmPgal;
      double tmPgRf = tmPgvf - tmPgaf;
      double tmPgRl = tmPgvl - tmPgal;

      if (debug_ > 5){
	std::cout <<" ############## process name :" << " tmPgLl: " << tmPgLl <<"tmPgRl: " << tmPgRl <<std::endl;
	std::cout << tmPgLl << " " << tmPgRl << std::endl;
      };

      // Kinematics
      //double qCmNew  = fMasterGen->particleData.m0(idNew);
      //double qCmNew2 = qCmNew * qCmNew;
      double qCmZ    = fMasterGen->particleData.m0(23);
      double qCmZ2   = qCmZ * qCmZ;
      double qCGZ    = fMasterGen->particleData.mWidth(23);
      double qCGZ2   = qCGZ * qCGZ;

      // Necessary variables to ampitude
      // First term
      // double alpEM =1./137;
      double tmPe2QfQl = 4. * M_PI * alpEM * fMasterGen->couplingsPtr->ef(idAbs) * fMasterGen->couplingsPtr->ef(idNew);
      double qCPropGm   = 1./sH;
      if (debug_ > 5) {
	// double alpha= fMasterGen->info.alphaEM();
	std::cout <<"DIRECTLY ACCESSIBLE VARIABLES "<< std::endl;
	std::cout <<" M_P ouputs Value of pi : " <<  M_PI <<std::endl;
	std::cout <<"fMasterGen->couplingsPtr->ef(idAbs) outputs Charge of incoming quark: " << fMasterGen->couplingsPtr->ef(idAbs) <<std::endl;
	std::cout <<" fMasterGen->couplingsPtr->ef(idNew) outputs Charge of outgoing leopton: " << fMasterGen->couplingsPtr->ef(idNew) <<std::endl;
	std::cout <<"alpEM CAN NOT BE ACCESSED DIRECTLY "<<std::endl;
        // std::cout <<" alphaQED() pythia->info.alphaEM();. "<< alpha <<std::endl;
        // std::cout <<"AlphaQED "<< AlphaQED <<std::endl;
      };
      //Second term.Model depended variables are defined using incoming quark and outgoing fermion information
  double tmPe2s2c2 = 4. * M_PI * alpEM;
 double denomPropZ = pow((sH - qCmZ2), 2) + qCmZ2 * qCGZ2;
 double qCrePropZ  = (sH - qCmZ2) / denomPropZ;
 double qCimPropZ  = -qCmZ * qCGZ / denomPropZ;

 //Third term:4. * M_PI * qCetaLR / qCLambda2;
 // Amplitudes, M = gamma + Z + CI.
      meLL = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgLl * (qCrePropZ + I * qCimPropZ);
      meRR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgRl * (qCrePropZ + I * qCimPropZ);
      meLR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ)
        + 4. * M_PI * qCetaLR / qCLambda2;
      meRL = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ)
        + 4. * M_PI * qCetaLR / qCLambda2;

      // According to Steve's idea, add SM matrix elements for sigmaNew.
      // Define standard model matrix elements of RL and LR model

      meLR_SM = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ);

      meRL_SM = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ);


      // Calculate weighting facror
      double sigma0 = 1.0;
      double sigmaOld = sigma0 * uH2 * std::real(meLL*std::conj(meLL));
      sigmaOld += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
      sigmaOld += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
      sigmaOld += sigma0 * tH2 * std::real(meRL*std::conj(meRL));

      double sigmaNewLR = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
      sigmaNewLR += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
      sigmaNewLR += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
      // sigma += sigma0 * tH2 * std::real(meRL*std::conj(meRL));
      sigmaNewLR += sigma0 * tH2 * std::real(meRL_SM *std::conj(meRL_SM));
      double fracLR = sigmaNewLR / sigmaOld;

      double sigmaNewRL = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
      sigmaNewRL += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
      //sigmaNew += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
      sigmaNewRL += sigma0 * tH2 * std::real(meRL*std::conj(meRL));
      sigmaNewRL += sigma0 * tH2 * std::real(meLR_SM*std::conj(meLR_SM));
      double fracRL = sigmaNewRL / sigmaOld;
      // double fracSum += frac;
      //double weight *= frac;
      double sumRL_Plus_LR =fracLR+fracRL;

      //Print out weighting factor
      if (debug_ >0){
	std::cout << std::endl;
	std::cout <<"############ FRACTION LR and RL #############" << std::endl;
	std::cout << "fracLR:  " << fracLR << "  " << "fracRL:  " <<fracRL << std::endl;
	std::cout << std::endl;
	std::cout << "sumRL_Plus_LR:  " << fracLR+fracRL << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout <<"############ WEIGHT  #############" << std::endl;
	// std::cout << weight << std::endl;
      };

      if(debug_ > 5){
	std::cout << "Eta of daughter1 is: " << muMinus->eta() << "\n";
	std::cout << "Eta of daughter2 is: " << muPlus->eta() << "\n";
      }
      //Fill histogram for DY
      if(boson){
        bosonId_=boson->pdgId();
        bosonP4_.fill(boson->p4());

        h_Zmass->Fill(boson->mass(), fracLR);
        h_Zpt->Fill(boson->pt(), fracLR);
        h_Zeta->Fill(boson->eta(),fracLR);
        h_Zphi->Fill(boson->phi(), fracLR);
        h_Zcharge->Fill(boson->charge(), fracLR);
      
      }

      if (muMinus != nullptr){
	std::cout << "muMinus is not empty" << std::endl;
      }

      if (debug_ >0){
	// Print out final state lepton
	std::cout << " - FINAL STATE LEPTON INFORMATION:  " << std::endl;
	std::cout << std::endl;
	std::cout << "+ LEPTON:" <<std::endl;
	std::cout <<"muMinus ID:" << muMinus->pdgId() <<std::endl;
	std::cout <<"muMinus Energy: " << muMinus->energy()<< std::endl;
	std::cout <<"muMinus px: " <<muMinus->px() <<std::endl;
	std::cout <<"muMinus py: " <<muMinus->py() << std::endl;
	std::cout <<"muMinus pz: " << muMinus->pz() <<std::endl;
	std::cout <<"muMinus pt: " << muMinus->pt() <<std::endl;
	std::cout <<"muMinus eta: " << muMinus->eta() <<std::endl;
	std::cout <<"muMinus phi: " << muMinus->phi() <<std::endl;
	std::cout << std::endl;
	std::cout << "+ ANTILEPTON" << std::endl;
	std::cout <<"muPlus ID:" << muPlus->pdgId() <<std::endl;
	std::cout <<"muPlus Energy: " << muPlus->energy()<< std::endl;
	std::cout <<"muPlus px: " <<muPlus->px() <<std::endl;
	std::cout <<"muPlus py: " <<muPlus->py() << std::endl;
	std::cout <<"muPlus pz: " << muPlus->pz() <<std::endl;
	std::cout <<"muPlus pt: " << muPlus->pt() <<std::endl;
	std::cout <<"muPlus eta: " << muPlus->eta() <<std::endl;
	std::cout <<"muPlus phi: " << muPlus->phi() <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
      };
      // Fill histograms
      muMinusP4_.fill(muMinus->p4());
      muMinusPID_=muMinus->pdgId();
      if(debug_ > 3){
	std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = "<< muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " <<muMinus->charge();
	std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << muPlus->phi() << "\tq = " << muPlus->charge();
      }

      h_muMinusmass->Fill(muMinus->mass(),fracLR);
      h_muMinuspt->Fill(muMinus->pt(), fracLR);
      h_muMinuseta->Fill(muMinus->eta(), fracLR);
      h_muMinusphi->Fill(muMinus->phi(), fracLR);
      h_muMinuscharge->Fill(muMinus->charge(), fracLR);
      h_thetaMuMinus->Fill(muMinus->theta(), fracLR);

      muPlusP4_.fill(muPlus->p4());
      muPlusPID_=muPlus->pdgId();
      h_muPlusmass->Fill(muPlus->mass(), fracLR);
      h_muPluspt->Fill(muPlus->pt(), fracLR);
      h_muPluseta->Fill(muPlus->eta(), fracLR);
      h_muPlusphi->Fill(muPlus->phi(), fracLR);
      h_muPluscharge->Fill(muPlus->charge(), fracLR);
      h_thetaMuPlus->Fill(muPlus->theta(), fracLR);

      //Calculating invariant mass for final state leptons
      muPlusKPlus = (1/sqrt(2))*(muPlus->energy() + muPlus->pz());
      muPlusKMinus = (1/sqrt(2))*(muPlus->energy() - muPlus->pz());
      muMinusKPlus = (1/sqrt(2))*(muMinus->energy() + muMinus->pz());
      muMinusKMinus = (1/sqrt(2))*(muMinus->energy() - muMinus->pz());

      invariantK = (muPlusKPlus*muMinusKMinus - muMinusKPlus*muPlusKMinus);
      if (debug_ > 0) {
	std::cout << std::endl;
	std::cout << "\n\nInvariantK is: " << invariantK << std::endl;
	std::cout << std::endl;
      };

      dimuon = muMinus->p4() + muPlus->p4();

      dimuonPt =dimuon.pt();
      dimuonPz = dimuon.pz();
      pseudorapidity = asinh(dimuonPz/dimuonPt);
      dimuonPx = dimuon.px();
      dimuonPy = dimuon.py();
      Phi = acos(dimuonPx/dimuonPt);
      dimuonQ = sqrt(pow(dimuon.energy(),2) - pow(dimuon.pt(),2) - pow(dimuon.pz(),2));
      if (debug_ >0) {
	std::cout << std::endl;
	std::cout << "\n\nDimuon Energy is: " << dimuon.energy() << std::endl << std::endl;
	std::cout << std::endl;
      };

      double denominatorTheta, denominatorPhi1, denominatorPhi2, numeratorPhi1, numeratorPhi2;
      double denominatorPhi, numeratorPhi;
      double deltaX, deltaY;
      double invariantMass;
      denominatorTheta = dimuonQ*sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
      thetaCos = (dimuon.pz()/fabs(dimuon.pz()))*(2/denominatorTheta)*invariantK;
      thetaCS = acos(thetaCos);

      denominatorPhi1 = dimuonQ*dimuon.pt();
      numeratorPhi1 = sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
      deltaX = muPlus->px() - muMinus->px();
      deltaY = muPlus->py() - muMinus->py();
      denominatorPhi2 = ((deltaX*dimuon.px()) + (deltaY*dimuon.py()));
      numeratorPhi2 = ((deltaX*dimuon.py()) + (deltaY*dimuon.px()));
      numeratorPhi = numeratorPhi1*numeratorPhi2;
      denominatorPhi = denominatorPhi1 * denominatorPhi2;

      phiTan = numeratorPhi/denominatorPhi;
      phiCS = atan(phiTan);


      mu1Energy = muPlus->energy();
      mu2Energy = muMinus->energy();

      if (debug_ > 5) {
	std::cout << "\n\nmuon1,2 Energies are: " << mu1Energy << "__" << mu2Energy<< std::endl << std::endl;
	std::cout << "\ndimuon1,2 px_py_pz are: "<< dimuonPx << "_" << dimuonPy << "_" << dimuonPz << std::endl;
      };
      //Exaple TTree weighting
      //      cosTheta=thetaCos*fracLR;
      cosTheta=thetaCos;
      tanPhi=phiTan;
      csTheta=thetaCS;
      csPhi=phiCS;

      //      h_cosTheta->Fill(thetaCos);
      h_cosTheta->Fill(cosTheta);
      h_csTheta->Fill(thetaCS);
      h_tanPhi->Fill(phiTan);
      h_csPhi->Fill(phiCS);

      if (debug_ > 3){
	std::cout << "\n\n\ncos(Theta_CS) = " << thetaCos << "\tThetaCS = " << thetaCS << std::endl;
	std::cout << "\n\n\nTan(phi_CS) = " << phiTan << "\tPhiCS = " << phiCS <<  std::endl;
      };

      invariantMass = sqrt(2 * muMinus->pt() * muPlus->pt() *( cosh(muMinus->eta() - muPlus->eta()) - cos(TVector2::Phi_mpi_pi(muMinus->phi() - muPlus->phi()))));

      finalInvariantMass = invariantMass;

      h_finalInvariantMass->Fill(invariantMass,fracLR);
      h_InvariantMassvscosTheta->Fill(thetaCos,invariantMass);
      //Quark selection: we do not need it, just sanity check
      // Delete it after using it, it is just a garbage

      int dquark,uquark,squark,cquark,bquark,tquark;
      if(idAbs==1){
	dquark=idAbs;
      }
      else if(idAbs==2){
	uquark=idAbs;
      }
      else if(idAbs==3){
	squark=idAbs;
      }
      else if(idAbs==4){
	cquark=idAbs;
      }
      else if(idAbs==5){
	bquark=idAbs;
      }
      else {
	tquark=idAbs;
      }


      // Delete it after checking

      DQuark=dquark;
      UQuark=uquark;
      SQuark=squark;
      CQuark=cquark;
      BQuark=bquark;
      TQuark=tquark;
      h_DQuark->Fill(dquark);
      h_UQuark->Fill(uquark);
      h_SQuark->Fill(squark);
      h_CQuark->Fill(cquark);
      h_BQuark->Fill(bquark);
      h_TQuark->Fill(tquark);


      // Fill the hostograms for Mandelston variables

      SHat=sH;
      THat2=tH2;
      UHat2=uH2;
      fractionLR=fracLR;
      fractionRL=fracRL;

      sumFrac=sumRL_Plus_LR;
      T_Hat=tH;
      U_Hat=uH;

      alphaRunQED=alpEM;

      h_SHat->Fill(sH);
      h_THat2->Fill(tH2);
      h_UHat2->Fill(uH2);
      h_fractionLR->Fill(fracLR);
      h_fractionRL->Fill(fracRL);
      h_sumFrac->Fill(sumFrac);

      h_T_Hat->Fill(tH);
      h_U_Hat->Fill(uH);
      h_alphaRunQED->Fill(alpEM);

      if(thetaCos < 0.0){
	h_cosThetaMinusInvariantMass->Fill(invariantMass);
	mCosThetaMinus = invariantMass;
      }
      else{
	h_cosThetaPlusInvariantMass->Fill(invariantMass);
	mCosThetaPlus = invariantMass;
      };


      if(thetaCos < 0.0){
	h_cosThetaMinusInvariantMass->Fill(invariantMass);
	mCosThetaMinus = invariantMass;
      }
      else{
	h_cosThetaPlusInvariantMass->Fill(invariantMass);
	mCosThetaPlus = invariantMass;
      };

      h_dphi->Fill(TVector2::Phi_mpi_pi(muMinus->phi()- muPlus->phi()));
      h_dtheta->Fill(TVector2::Phi_mpi_pi(muMinus->theta()- muPlus->theta()));
      h_dr->Fill(reco::deltaR(muMinus->p4(),muPlus->p4()));
      h_massInvar->Fill(sqrt(2 * muMinus->pt() * muPlus->pt() *( cosh(muMinus->eta() - muPlus->eta()) - cos(TVector2::Phi_mpi_pi(muMinus->phi() - muPlus->phi())))));
      h_dimuonPt->Fill(dimuonPt);
      h_dimuonEta->Fill(pseudorapidity);
      h_dimuonPhi->Fill(Phi);
      h2_phi1_vs_phi2->Fill(muMinus->phi(),muPlus->phi());
      h2_eta1_vs_eta2->Fill(muMinus->eta(),muPlus->eta());
      h2_pt1_vs_pt2->Fill(muMinus->pt(),muPlus->pt());

      tree_->Fill();

      //  }


      double invariant_Mass;
      
      // Get Invariant Mass of initial leptons
      double product = 2* muonpt* antiMuonpt;
      double diff = cosh(muoneta-antiMuoneta)-cos(muonphi-antiMuonphi);
      double invariantMassinit = product*diff;
      if (invariantMassinit > 0)
	{
	  invariant_Mass =  sqrt(invariantMassinit);
	}
	std::cout << std::endl;
      std::cout << "Invariant mass of initial lepton:" <<  invariant_Mass << std::endl;

	std::cout << std::endl;
      std::cout << "\nInvariant mass of final state leptons: " << invariantMass <<std::endl;

	std::cout << std::endl;

  
 // Fill histograms
 
  massInvariant=invariant_Mass;

  h_massInvariant->Fill(massInvariant);
 
  

 
  std::cout << "\n\n===========================================================================================================" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;  

}// end analyze
      

      
// ------------ method called once each job just before starting event loop  ------------
void weightingAnalyzer::beginJob()
{
  //Acesses pythia parameters in CMSSW environment
  fMasterGen.reset( new Pythia8::Pythia );
  fMasterGen->init();


  edm::Service<TFileService> fs;
  h_Zmass = fs->make<TH1F>("Zmass" , "m", 1000, 0., 600);
  h_Zpt  = fs->make<TH1F>( "Zpt"  , "p_{t}", 500,  0., 2500. );
  h_Zeta = fs->make<TH1F>( "Zeta" , "#eta" , 100, -10., 10.    );
  h_Zphi = fs->make<TH1F>( "Zphi" , "#phi" , 100,  -3.20, 3.20   );
  h_Zcharge = fs->make<TH1F>( "Zcharge" , "Q" ,3,  -1.5, 1.5    );
  h_muMinusmass = fs->make<TH1F>("muMinusmass" , "m", 1000, 0., 500);
  h_muMinuspt  = fs->make<TH1F>( "muMinuspt"  , "p_{t}", 500,  0., 2500. );
  h_muMinuseta = fs->make<TH1F>( "muMinuseta" , "#eta" , 100, -5., 5.    );
  h_muMinusphi = fs->make<TH1F>( "muMinusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muMinuscharge = fs->make<TH1F>( "muMinuscharge" , "Q" ,3,  -1.5, 1.5    );

  h_muPlusmass = fs->make<TH1F>("muPlusmass" , "m", 1000, 0., 500);
  h_muPluspt  = fs->make<TH1F>( "muPluspt"  , "p_{t}", 500,  0., 2500. );
  h_muPluseta = fs->make<TH1F>( "muPluseta" , "#eta" , 100, -5., 5.    );
  h_muPlusphi = fs->make<TH1F>( "muPlusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muPluscharge = fs->make<TH1F>( "muPluscharge" , "Q" ,3,  -1.5, 1.5    );

  h_dphi = fs->make<TH1F>("delta phi", "#delta #phi", 100, -3.15, 3.15 );

  h_dtheta = fs->make<TH1F>("delta theta", "#delta #theta", 100, -3.15, 3.15);
  h_dr = fs->make<TH1F>("delta r", "#delta r", 100, 0, 10);
  h_thetaMuMinus = fs->make<TH1F>("theta muMinus", "#theta", 100, -3.15, 3.15);

  h_thetaMuPlus = fs->make<TH1F>("theta muPlus", "#theta", 100, -3.15, 3.15);
  h_massInvar = fs->make<TH1F>("Invariant mass", "Invariant mass", 350, 0., 3500.);
  h_dimuonPt = fs->make<TH1F>("Dimuon Pt", "Dimuon Pt", 500, 0, 2500);
  h_dimuonEta = fs->make<TH1F>("Dimuon eta", "Dimuon #eta", 100, -5, 5);
  h_dimuonPhi = fs->make<TH1F>("Dimuon Phi", "Dimuon #phi", 100, -3.15, 3.15);
  h_cosTheta = fs->make<TH1F>("cosTheta", "cos #theta", 100, -1.01, 1.01);
  h_tanPhi = fs->make<TH1F>("tanPhi", "tan #phi", 100, -1000.0, 1000.0);
  h_csTheta = fs->make<TH1F>("csTheta", "#theta_{CS}", 100, -3.15, 3.15);
  h_csPhi = fs->make<TH1F>("csPhi", "#phi_{CS}", 100, -3.15, 3.15);
  h_finalInvariantMass = fs->make<TH1F>("finalInvariantMass", "finalInvariantMass", 350, 0., 3500.);
  h_InvariantMassvscosTheta = fs->make<TH2F>("InvariantMassvscosTheta", "InvariantMassvscosTheta", 100, -1.01, 1.01, 50, 800., 1300.);
  h_cosThetaMinusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaMinus", "InvariantMass_cosThetaMinus", 350, 0., 3500.);
  h_cosThetaPlusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaPlus", "InvariantMass_cosThetaPlus", 350, 0., 3500.);
  h2_pt1_vs_pt2   = fs->make<TH2F>( "pt1_vs_pt2"   , "p_{t,1} vs. p_{t,2}"   , 500,  0., 2500., 500,  0., 2500.);
  h2_eta1_vs_eta2 = fs->make<TH2F>( "eta1_vs_eta2" , "#eta_{1} vs. #eta_{2}" , 100, -5., 5.   , 100, -5., 5.   );
  h2_phi1_vs_phi2 = fs->make<TH2F>( "phi1_vs_phi2" , "#phi_{1} vs. #phi_{2}" , 100,  -3.15, 3.15  , 100,  -3.15, 3.15  );

  //Book histogram for Mandelston variables
  h_SHat = fs->make<TH1F>("SHat", "SHat", 100, 0.1, 10000000.01);
  h_THat2 = fs->make<TH1F>("THat2", "THat2", 100, 0.1, 10000000000000.01);
  h_UHat2 = fs->make<TH1F>("UHat2", "UHat2", 100, 0.1, 10000000000000.01);
  h_fractionLR= fs->make<TH1F>("fractionLR", "fractionLR", 150, 0.1, 1.51);
  h_fractionRL= fs->make<TH1F>("fractionRL", "fractionRL", 150, 0.1, 1.51);
  h_sumFrac= fs->make<TH1F>("sumFrac", " sumFrac", 300, 0.1, 3.01);
  h_T_Hat = fs->make<TH1F>("T_Hat", "T_Hat", 100, -100000000.1, -1000.01);
  h_U_Hat = fs->make<TH1F>("U_Hat", "U_Hat", 100, -100000000.1, -1000.01);
  h_alphaRunQED = fs->make<TH1F>("alphaRunQED", "alphaRunQED", 100, 0.0001, 0.00001);

  //Delete it after checking it

  h_DQuark= fs->make<TH1F>("DQuark", "DQuark", 10, 0.1, 10.01);
  h_UQuark= fs->make<TH1F>("UQuark", "UQuark", 10, 0.1, 10.01);
  h_SQuark= fs->make<TH1F>("SQuark", "SQuark", 10, 0.1, 10.01);
  h_CQuark= fs->make<TH1F>("CQuark", "CQuark", 10, 0.1, 10.01);
  h_BQuark= fs->make<TH1F>("BQuark", "BQuark", 10, 0.1, 10.01);
  h_TQuark= fs->make<TH1F>("TQuark", "TQuark", 10, 0.1, 10.01);

  h_massInvariant = fs->make<TH1F>("Invariant Mass", "Invariant Mass", 350, 0.,4000);


  //  h_massInvar = fs->make<TH1F>("Invariant mass", "Invariant mass", 350, 0., 3500.);


  tree_= fs->make<TTree>("pdfTree","PDF Tree");
  // tree_->Branch("evtId",&evtId_,EventId::contents().c_str());
  tree_->Branch("bosonP4",&bosonP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1P4",&muMinusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay2P4",&muPlusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1PID",&muMinusPID_,"decay1PID/I");
  tree_->Branch("decay2PID",&muPlusPID_,"decay2PID/I");
  tree_->Branch("bosonPID",&bosonId_,"bosonPID/I");
  tree_->Branch("crossSec", &crossSec, "crossSec/D");
  tree_->Branch("cosTheta", &cosTheta, "cosTheta/D");
  tree_->Branch("finalInvariantMass", &finalInvariantMass, "finalInvariantMass/D");
  //tree_->Branch("InvariantMassvscosTheta",&InvariantMassthetaCos, "InvariantMassvscosTheta/D");
  tree_->Branch("tanPhi", &tanPhi, "tanPhi/D");
  tree_->Branch("csTheta", &csTheta, "csTheta/D");
  tree_->Branch("csPhi", &csPhi, "csPhi/D");
  tree_->Branch("mCosThetaPlus", &mCosThetaPlus, "mCosThetaPlus/D");
  tree_->Branch("mCosThetaMinus", &mCosThetaMinus, "mCosThetaMinus/D");

  // tree_->Branch("pdfInfo",&pdfInfo_,PDFInfo::contents().c_str());

  // Tree for Mandelston variabales

  tree_->Branch("SHat", &SHat, "SHat/D");
  tree_->Branch("THat2", &THat2, "THat2/D");
  tree_->Branch("UHat2", &UHat2, "UHat2/D");
  tree_->Branch("fractionLR", &fractionLR, "fractionLR/D");
  tree_->Branch("fractionRL", &fractionRL, "fractionRL/D");

  tree_->Branch("sumFrac", &sumFrac, "sumFrac/D");

  tree_->Branch("T_Hat", &T_Hat, "T_Hat/D");
  tree_->Branch("U_Hat", &U_Hat, "U_Hat/D");
  tree_->Branch("alphaRunQED", &alphaRunQED, "alphaRunQED/D");


  // Delete it after checking
  tree_->Branch("DQuark", &DQuark, "DQuark/I");
  tree_->Branch("UQuark", &UQuark, "UQuark/I");
  tree_->Branch("SQuark", &SQuark, "SQuark/I");
  tree_->Branch("CQuark", &CQuark, "CQuark/I");
  tree_->Branch("BQuark", &BQuark, "BQuark/I");
  tree_->Branch("TQuark", &TQuark, "TQuark/I");



};

// ------------ method called once each job just after ending the event loop  ------------
void 
weightingAnalyzer::endJob() 
{

  std::cout << "Total muon:" << muonCountTotal << std::endl;
  std::cout << "Stored muon:" << muonCountStored << std::endl;
  std::cout << "Total antimuon:" << antiMuonCountTotal << std::endl;
  std::cout << "Stored antimuon:" << antiMuonCountStored << std::endl;

}


/*


  const reco::Candidate* weightingAnalyzer::getIncomingQurks(const reco::Candidate* particle,int pid)
  {  
  for(size_t particleNum =0;  particleNum < particle->numberOfMothers(); particleNum++){
  //    if(particle->mother(particleNum)->pdgId()==pid) return particle->mother(particleNum);


  const reco::Candidate* mother= particle->mother(particleNum);
  std::cout<< "mother1 ID: " << mother->pdgId() <<std::endl;


  }
  return nullptr;
  }


*/


/*
bool weightingAnalyzer::isBoson(int pid)
{
  if(pid==23 || abs(pid)==22 || pid==32){
    if(debug_ > 0) std::cout << "\n\nFound Boson\n";
    return true;
  }
  else return false;
}
*/
bool weightingAnalyzer::isMuon(int pid){
  if(abs(pid)==11 || abs(pid) ==13){
    if(debug_ > 0) std::cout << "\n\nFound A Muon!\n";
    return true;
  }
  else return false;
}
/*
bool weightingAnalyzer::checkBosonStatus( const reco::GenParticleCollection& genParts){
  const reco::Candidate* boson = getBoson(genParts);
  if(boson == nullptr){
    if(debug_ > 0) std::cout << "\nBoson is: "  << boson;
    return false;
  }

  else if( boson->status() != 22){
    if(debug_ > 0)  std::cout <<"\nBoson Status is: "<< boson->status();
    return false;
  }

  return true;
}
*/

/*
const reco::Candidate* weightingAnalyzer::getBoson( const reco::GenParticleCollection& genParts)
{
  for(auto &part : genParts){
    if(isBoson(part.pdgId())){
      if(debug_ > 1){
	std::cout << "\npId is: " << part.pdgId();
	std::cout << "\nStatus is: " << part.status();
      }
      return getLastDaughter(&part,part.pdgId());
    }
  }
  return nullptr;
}
*/


const reco::Candidate* weightingAnalyzer::getBosonMother(const reco::GenParticle& p){
  const reco::Candidate* firstMother = nullptr;
  const reco::Candidate* muonMother = nullptr;
  const reco::Candidate* tempMother = nullptr;
  const reco::Candidate* firstMotherdaughter = nullptr;
  //  const reco::Candidate* motherdaughter = nullptr;
  int motherID;
  int firstMotherdaughterID;
  int motherdaughterID;

  
  firstMother = p.mother();
  firstMotherdaughter = firstMother->daughter(0);
   if (!firstMother){
    std::cout << "First mother failed for boson." <<std::endl;
    return nullptr;
  };
   if (!firstMotherdaughter){
    std::cout << "First mother daughter failed for boson." <<std::endl;
    return nullptr;
  };

  double firstMotherID = firstMother->pdgId();
  firstMotherdaughterID = firstMotherdaughter->pdgId();

  // First mother from hard interaction
  if (firstMotherID == 23 && (abs(firstMotherdaughterID) ==  11 || abs(firstMotherdaughterID) == 13)) {
	muonMother= firstMother;
  }//end if 
  
  //More than one mother
  tempMother = p.mother();
  
  motherID =tempMother->pdgId();
  const reco::Candidate* motherdaughter = nullptr;
  // Loops over untill it finds mother from hard interaction
  while(motherID != 23) {
    motherdaughter = tempMother;
    tempMother=tempMother->mother();
    
    if (!tempMother){
      return nullptr;
    }
    if (!motherdaughter){
      return nullptr;
    }

    motherID=tempMother->pdgId();
    motherdaughterID = motherdaughter->pdgId();

    if (motherID == 23 && (abs(motherdaughterID) ==  11 || abs(motherdaughterID) == 13)) {

      //      const reco::GenParticle& lastMuon=p;
      //std::cout<< "lastMuon" << lastMuon.pdgId() << std::endl;
	  muonMother=tempMother;
      }//end if 
  }


  // Checking 
  if(muonMother!=nullptr){

    std::cout<<" bosonMother pointer is  filled " << std::endl;

  }
  else{
    std::cout<<" bosonMother pointer is not filled " << std::endl;
  }// end of checking

  return muonMother;
}


const reco::Candidate* weightingAnalyzer::getLeptonMother(const reco::GenParticle& p, bool second){
  const reco::Candidate*  firstMother=nullptr;
  const reco::Candidate*  muonMother=nullptr;
  const reco::Candidate* tempMother=nullptr;
  int motherID;

  
  firstMother = p.mother();
  if (!firstMother)
    return nullptr;

  // First mother from hard interaction
  if (abs(firstMother->pdgId()) <= 6 && firstMother->status()== 21){

      if (second)
	if (firstMother->pdgId() > 0){
	  muonMother=firstMother->daughter(0)->mother(1);
	}
	else{
	  muonMother = firstMother;
	}
      else
	if (firstMother->pdgId() > 0){
	  muonMother= firstMother;
	}
	else{
	  muonMother=firstMother->daughter(0)->mother(1);
	}
    }


    // Function call for muon

    //    getLastMuon(

    //const reco::GenParticle muon =p;
    // muonMother=p.mother();

    //const reco::Candidate* lastMuon =p;
    //    const reco::GenParticle* lastMuon=p;
    // std::cout<< "lastMuon" << lastMuon->pdgId() << std::endl;
  //end if 
  


  //More than one mother
  tempMother = p.mother();
  motherID =tempMother->pdgId();
 
  // Loops over untill it finds mother from hard interaction
  while(abs(motherID) > 6) {

    tempMother=tempMother->mother();
    if (!tempMother)
      return nullptr;
    motherID=tempMother->pdgId();

    if (abs(motherID) <= 6 && tempMother->status()== 21){

      //      const reco::GenParticle& lastMuon=p;
      //std::cout<< "lastMuon" << lastMuon.pdgId() << std::endl;
      if (second)
	if (motherID > 0){
	  muonMother=tempMother->daughter(0)->mother(1);
	}
	else{
	  muonMother = tempMother;
	}
      else
	if (motherID > 0){
	  muonMother= tempMother;
	}
	else{
	  muonMother=tempMother->daughter(0)->mother(1);
	}	



    }

      //      const reco::GenParticle& lastMuon=p;
      //std::cout<< "lastMuon" << lastMuon.pdgId() << std::endl;  
    


  }//end while


  // Checking 
  if(muonMother!=nullptr){

    std::cout<<" muonMother pointer is  filled " << std::endl;

  }
  else{
    std::cout<<" muonMother pointer is not filled " << std::endl;
  }// end of checking

  return muonMother;
}


/*

double weightingAnalyzer:: findInvariantMass(const reco::Candidate* particle, const reco::Candidate* antiparticle)
{
  //  double product = 2*particle->pt()*antiparticle->pt();
  double product = 2* (*particle).pt()* (*antiparticle).pt();
  double diff = cosh(particle->eta()-antiparticle->eta())-cos(particle->phi()-antiparticle->phi());
  double invariantMass = product*diff;
  if (invariantMass > 0)
    {
      return sqrt(invariantMass);
    }
  else
    {
      return false;
    }
}

*/



  void weightingAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
  {
  edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByToken(genInfoProductToken_, genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();
  std::cout<< "Cross Section is: "  << crossSec << std::endl;  
 
  }




// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void weightingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(weightingAnalyzer);

