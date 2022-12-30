class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
#include <iostream>
using std::cout;
using std::endl;
void ana_D2535(TString InFile="Ds2535_DstarK_fast.root", int nevts=100000, double pbarmom = 9.83)
{
	// *** some variables
	int i=0,j=0, k=0, l=0;
	gStyle->SetOptFit(1011);
	
	// *** the output file for FairRunAna
	TString OutFile="dummy_out.root";  
					
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	
	FairRunAna* fRun = new FairRunAna();
	fRun->SetWriteRunInfoFile(kFALSE);
	fRun->SetInputFile(InFile);
	fRun->SetOutputFile(OutFile); // only dummy; the real output is 
	
	fRun->Init(); 
	
	// *** take constant field; needed for PocaVtx
	RhoCalculationTools::ForceConstantBz(20.0);
	
	// *** create an output file for all histograms
	TFile *out = TFile::Open(InFile+"_ana.root","RECREATE");
	
	// *** create some ntuples
	RhoTuple *ntp1 = new RhoTuple("ntp1", "Ds- analysis");
	RhoTuple *ntp2 = new RhoTuple("ntp2", "Ds2535 analysis");
	RhoTuple *ntp3 = new RhoTuple("ntp3", "D0 analysis");
	RhoTuple *ntp4 = new RhoTuple("ntp4", "Ds*+ analysis");
	RhoTuple *nmc  = new RhoTuple("nmc",  "mctruth info");
	RhoTuple *ntp5 = new RhoTuple("ntp5", "Vertex fitter");
	
	//Histograms for vertex fitters

	
	//PndKinVtxFitter
	TH1F *Ds_vtx_diff_X  = new TH1F("Ds_vtx_diff_X", "Ds-: vertex resolution, X",100,-0.5,0.5);
	TH1F *Ds_vtx_diff_Y = new TH1F("Ds_vtx_diff_Y", "Ds-: vertex resolution, Y",100,-0.5,0.5);
	TH1F *Ds_vtx_diff_Z = new TH1F("Ds_vtx_diff_Z", "Ds-: vertex resolution, Z",100,-0.5,0.5);
	//TH2F *Ds_vtx_2Ddiff = new TH2F("Ds_vtx_diff_Pos","fitted decay vertex poisition",100,-0.5,0.5,100,-0.5,0.5);
	TH1F *Ds_vtx_Pos = new TH1F("Ds_vtx_diff_Pos","Ds_vtx_diff_Pos",100,-0.5,0.5);
	TH1F *Ds_vtx_PosMC = new TH1F("Ds_vtx_diff_PosMC","Ds_vtx_diff_PosMC",100,-4,4);
	TH1F *Ds_chi2_vf  = new TH1F("Ds_chi2_vf", "Ds-: #chi^{2} vertex fit",100,-100,100);
	TH1F *Ds_chi2_prob_vf   = new TH1F("Ds_chi2_prob_vf",  "Ds: prob vertex fit",100,-0.01,1);
	TH1F *Ds_vf_mass   = new TH1F("Ds_vf_mass","Ds- mass (vertex fit)",100,1,3);
	TH2F *hvpos2 = new TH2F("hvpos2","(x,y) projection of fitted decay vertex Ds-",100,-2,2,100,-2,2);
	TH1F *Ds_vtx_z   = new TH1F("Ds_vtx_z",  "Ds: Z",100,-0.5,0.5);
	TH1F *hm_diff2  = new TH1F("hm_diff2","Ds- mass diff to truth after vertex fit",100,-0.5,0.5);
	
	TH1F *Ds_missingmass_ftm  = new TH1F("Ds_missingmass_ftm","Ds1(2535): Ds- missing mass (full truth match)",100,0,5.);
	TH1F *Ds_px_missing_ftm  = new TH1F("Ds_px_missing_ftm","Ds missing px ",100,-1.5,1.5);
	TH1F *Ds_py_missing_ftm  = new TH1F("Ds_py_missing_ftm","Ds missing py ",100,-1.5,1.5);
	TH1F *Ds_pz_missing_ftm  = new TH1F("Ds_pz_missing_ftm","Ds missing pz ",100,-1,9);
	TH1F *Ds_E_missing_ftm  = new TH1F("Ds_E_missing ftm","Ds missing E ",100,0,10);
	TH1F *Ds_mass_diff_ftm  = new TH1F("Ds_mass_diff_ftm","Ds E ",100,-0.2,0.2);
	TH1F *Ds_missingmass_diff_ftm  = new TH1F("Ds_missingmass_diff_ftm","Ds E ",100,-0.2,0.2);

	TH1F *Kaon_missingmass_ftm  = new TH1F("Kaon_missingmass_ftm","D0*mass as missing mass (full truth match)",100,-10XS,10.);
	TH1F *Kaon_missing_px_ftm  = new TH1F("Kaon_missing_px_ftm","K missing px ",100,-1.5,1.5);
	TH1F *Kaon_missing_py_ftm  = new TH1F("Kaon_missing_py_ftm","K missing py ",100,-1.5,1.5);
	TH1F *Kaon_missing_pz_ftm  = new TH1F("Kaon_missing_pz_ftm","K missing pz ",100,-1,9);
	TH1F *Kaon_missing_px_diff_ftm  = new TH1F("Kaon_missing_px_diff_ftm","K missing px res",100,-0.2,0.2);
	TH1F *Kaon_missing_py_diff_ftm  = new TH1F("Kaon_missing_py_diff_ftm","K missing py res",100,-0.2,0.2);
	TH1F *Kaon_missing_pz_diff_ftm  = new TH1F("Kaon_missing_pz_diff_ftm","K missing pz res",100,-2,2);
	TH1F *Kaon_missing_E_diff_ftm  = new TH1F("Kaon_missing_E_diff_ftm","K missing E res",100,-0.2,0.2);
	TH1F *Kaon_missing_E_ftm  = new TH1F("Kaon_missing_E_ftm","K missing E ",100,0,10);
	TH1F *Kaon_missingmass_diff_ftm  = new TH1F("Kaon_missingmass_diff_ftm","K missing mass resolution ",100,-0.2,0.2);
	TH1F *Kaon_missing_charge  = new TH1F("Kaon_missing_charge","charge ",100,-2,2);

	
	TH1F *D0star_chi2_mf  = new TH1F("D0star_chi2_mf", "D0*: #chi^{2} mass fit",100,-1,10);
	TH1F *D0star_mcf  = new TH1F("D0star_mcf","D0* mass (mass constraint fit)",100,0,3);
	TH1F *D0star_mcf_res  = new TH1F("D0star_mcf_res","D0* mass res (mass constraint fit)",100,-0.2,0.2);
	TH1F *D0star_chi2_prob_mf  = new TH1F("D0star_chi2_prob_mf", "D0*: prob mass constraint fit",100,0,1);

	//
	// Now the analysis stuff comes...
	//
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** name of the only PidAlgo TClonesArray in fsim
	TString pidalg = "PidChargedProbability";
	TString pidalg0 ="PidNeutralProbability";

	// *** QA tool for simple dumping of analysis results in RhoRuple
	// *** WIKI: https://panda-wiki.gsi.de/foswiki/bin/view/Computing/PandaRootAnalysisJuly13#PndRhoTupleQA
	PndRhoTupleQA qa(theAnalysis, pbarmom); 
	
	// *** Mass selector for the jpsi cands
	double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi


	RhoMassParticleSelector *D0MassSel=new RhoMassParticleSelector("D0",1.8645,0.5);
	RhoMassParticleSelector *DsMassSel=new RhoMassParticleSelector("D_s-",1.9686,0.5);
	RhoMassParticleSelector *D0starMassSel=new RhoMassParticleSelector("D0*",2.00662,0.2);
	RhoMassParticleSelector *Ds2535MassSel=new RhoMassParticleSelector("D'_s1+",2.53535,0.5);
	RhoMassParticleSelector *D0MassSelPrefit=new RhoMassParticleSelector("D0",1.8645,0.8);
	RhoMassParticleSelector *DsMassSelPrefit=new RhoMassParticleSelector("D_s-",1.9686,0.95);

	RhoMomentumParticleSelector *GammaMomentum=new 	RhoMomentumParticleSelector("gamma",1.0,0.99);	
	// *** the lorentz vector of the initial psi(2S)
	//double m0_p    = TDatabasePDG::Instance()->GetParticle("proton")->Mass();   // Get nominal PDG mass of the proton
	//	TLorentzVector ini(0, 0, pbarmom, sqrt(m0_p*m0_p + pbarmom*pbarmom) + m0_p);
	TLorentzVector ini(0, 0, 9.83,10.813);

	Double_t kppx, kppy, kppz, kpE;
	Double_t kmpx, kmpy, kmpz, kmE;
	Double_t pmpx, pmpy, pmpz, pmE;
	Double_t m12, m23, m13, m2P12, m2P13, m2P23;
	Double_t pxD, pyD, pzD, ED;
	Double_t pxDs, pyDs, pzDs, EDs;
	Double_t Ebeam, ss, myDeltaE,myDeltaEDs,myDeltaEDDbar;
	Double_t pxm,pym,pzm,Em;
	Double_t pxmm,pymm,pzmm,Emm;
	Double_t pxt,pyt,pzt,Et;
	Double_t pKx, pKy, pKz, EK, Kcharge;
	Double_t pKxm, pKym, pKzm, EKm;
        Double_t pKxmt, pKymt, pKzmt, EKmt;


	RhoCandList kplus, kminus, piplus, piminus, D0list, Dslist,Ds2535, D0star, gammas, all, mclist, muplus;

	// ***
	// the event loop
	// ***
		
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
		
		// *** get MC list and store info in ntuple
	 	theAnalysis->FillList(mclist, "McTruth");
		
		nmc->Column("ev", (Int_t) i);
		qa.qaMcList("",   mclist, nmc);
		nmc->DumpData();
	
			
		// *** Setup event shape object
		theAnalysis->FillList(all,  "All", pidalg);
		PndEventShape evsh(all, ini, 0.05, 0.1);//here a cut is set up on momentim;
		                                        //p_gamma>50 MeV and  p_charge >100 MeV

		// *** Select with no PID info ('All'); type and mass are set 
		theAnalysis->FillList(kplus,  "KaonBestPlus",  pidalg);
		theAnalysis->FillList(kminus, "KaonBestMinus", pidalg);
		theAnalysis->FillList(piplus, "PionBestPlus",  pidalg);
		theAnalysis->FillList(piminus,"PionBestMinus", pidalg);
		theAnalysis->FillList(gammas, "Neutral",pidalg0);

		
		// *****
		// ***
		// *****
	
	        
		gammas.SetType(22);
		gammas.Select(GammaMomentum);
		int ngammasmct = theAnalysis->McTruthMatch(gammas);
		
		
		gammas.Select(GammaMomentum);
		// *** combinatorics for D0 to K- pi+
		D0list.Combine(kminus,piplus);
		
		// *** combinatorics for Ds- to K+ K- pi-
		Dslist.Combine(kminus,kplus,piminus);
		// *** combinatorics for Ds- to K+ K- pi-
		Dslist.Combine(kminus,kplus,piminus);
		
		// *** combinatorics for D0* to D0 gamma
		D0star.Combine(D0list,gammas);
		
		// *** combinatorics for Ds1+ to D0* K+
		Ds2535.Combine(D0star,kplus);
		
		
		
		D0list.SetType(421);
		gammas.Select(GammaMomentum);
		
		int nd0mct = theAnalysis->McTruthMatch(D0list);
		
		for (j=0;j<D0list.GetLength();++j) 
		  {
		    // some general info about event (actually written for each candidate!)
		    ntp3->Column("ev",		(Float_t) i);
		    ntp3->Column("cand",	(Float_t) j);
		    ntp3->Column("ncand",       (Float_t) D0list.GetLength());
		    ntp3->Column("nmct",        (Float_t) nd0mct);
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp3);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("D0", D0list[j], ntp3);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp3);
		    
		    
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = D0list[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCD0", lv, ntp3);
		    
		    ntp3->DumpData();
		  }
		
		//ntp2
		Dslist.SetType(-431);
		
		int ndsmct = theAnalysis->McTruthMatch(Dslist);
		
		for (j=0;j<Dslist.GetLength();++j) 
		  {
		    // some general info about event (actually written for each candidate!)
		    ntp2->Column("ev",		(Float_t) i);
		    ntp2->Column("cand",	(Float_t) j);
		    ntp2->Column("ncand",       (Float_t) Dslist.GetLength());
		    ntp2->Column("nmct",        (Float_t) ndsmct);
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp2);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Ds", Dslist[j], ntp2);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp2);
		   
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = Dslist[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDs", lv, ntp2);


	
		//Vertex fitters: PndKinVtxFitter on Ds-

	
		  PndKinVtxFitter vtxfitter(Dslist[j]);	// instantiate a vertex fitter
		  vtxfitter.Fit();
		  
		  double chi2_vtx=vtxfitter.GetChi2();	// access chi2 of fit
		  double prob_vtx=vtxfitter.GetProb();
		  Ds_chi2_vf->Fill(chi2_vtx);
		  Ds_chi2_prob_vf->Fill(prob_vtx);
			  
		  // if (chi2_vtx>0.&&chi2_vtx<1.)				
		       //{
		      RhoCandidate *jfit = Dslist[j]->GetFit();	// access the fitted cand	      
		      TVector3 jVtx=jfit->Pos();		// and the decay vertex position
		      
		      Ds_vf_mass->Fill(jfit->M());    		      
		      hvpos2->Fill(jVtx.X(),jVtx.Y());
		      Ds_vtx_z->Fill(jVtx.Z());
		      
		      RhoCandidate *mcdsm = Dslist[j]->GetMcTruth();
		      if (mcdsm && mcdsm->PdgCode()==-431){
			hm_diff2->Fill( mcdsm->M() - jfit->M());
		      }
		      RhoCandidate *mcvtx = Dslist[j]->GetMcTruth();
		      if (mcdsm && mcdsm->PdgCode()==-431){
			Ds_vtx_diff_X->Fill( mcdsm->Pos().X() - jVtx.X());
			Ds_vtx_diff_Y->Fill( mcdsm->Pos().Y() - jVtx.Y());
			Ds_vtx_diff_Z->Fill( mcdsm->Pos().Z() - jVtx.Z());
			Ds_vtx_diff_PosMC->Fill(  TMath::Sqrt(mcdsm->Pos().X()*mcdsm->Pos().X()+ mcdsm->Pos().Y()*mcdsm->Pos().Y() + mcdsm->Pos().Z()*mcdsm->Pos().Z())) ;
			Ds_vtx_diff_Pos->Fill(TMath::Sqrt(jVtx.X()*jVtx.X()+jVtx.Y()*jVtx.Y()++jVtx.Z()*jVtx.Z()));
		      }	
		      //}

		      //Add here : if (theAnalysis->McTruthMatch(Dslist[j])) if you want teh true values ftm.
		      if (theAnalysis->McTruthMatch(Dslist[j]))
			{ 
			  
			  //missing mass
			  
			  RhoCandidate *kpluscand = Dslist[j]->Daughter(0);
			  RhoCandidate *kminuscand = Dslist[j]->Daughter(1);
			  RhoCandidate *piminuscand = Dslist[j]->Daughter(2);
			  
			  //Ds- bachelor 4mom
			  pxDs = Dslist[j]->Px();
			  pyDs = Dslist[j]->Py();
			  pzDs = Dslist[j]->Pz();
			  EDs  = Dslist[j]->E();
			  
			  
			  //Ds2535+ 4mom as missing mass
			  pxmm = -pxDs;  //ini_x - pDs_x
			  pymm = -pyDs;
			  pzmm = ini.Z()-pzDs;
			  Emm =  ini.E()-EDs;
			  
			  pxt = -lv.X();
			  pyt = -lv.Y();
			  pzt = ini.Z()-lv.Z();
			  Et = ini.E()-lv.E();
			  
		      	  
			  TLorentzVector qD(pxDs, pyDs, pzDs, EDs); //Ds- components    
			  TLorentzVector missingmass(pxmm,pymm,pzmm,Emm); //Ds(2535)+ new components
			  TLorentzVector truevalueDs2535(pxt,pyt,pzt,Et); //Ds(2535)+ true values
			  
			  
			  Ds_missingmass_ftm->Fill(missingmass.M());
			  Ds_px_missing_ftm->Fill(missingmass.Px());
			  Ds_py_missing_ftm->Fill(missingmass.Py());
			  Ds_pz_missing_ftm->Fill(missingmass.Pz());
			  Ds_E_missing_ftm->Fill(missingmass.E());
			  Ds_missingmass_diff_ftm->Fill(truevalueDs2535.M() - missingmass.M());
			  
			  ntp2->DumpData();
			}
		      
		  }
		


		
		//ntp4
		D0star.SetType(423);
		gammas.SetType(22);
		
		D0list.Select(D0MassSelPrefit);
		gammas.Select(GammaMomentum);
		int ndsmct = theAnalysis->McTruthMatch(D0star);
		
		for (j=0;j<D0star.GetLength();++j) 
		  {
		     if (theAnalysis->McTruthMatch(D0star[j])){ 
		    // some general info about event (actually written for each candidate!)
		    ntp4->Column("ev",		(Float_t) i);
		    ntp4->Column("cand",	(Float_t) j);
		    ntp4->Column("ncand",       (Float_t) Dslist.GetLength());
		    ntp4->Column("nmct",        (Float_t) ndsmct);
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp4);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("D0star", D0star[j], ntp4);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp4);
		   
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = D0star[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCD0star", lv, ntp4);




		    
		    //for D*0 we use the mass constraint fit:
		    PndKinFitter mfitter(D0star[j]);        // instantiate the PndKinFitter in Ds+1
		    mfitter.AddMassConstraint(2.00697);
		    mfitter.Fit();	
		    
		    double chi2_m = mfitter.GetChi2();	// get chi2 of fit
		    double prob_m = mfitter.GetProb();
		    D0star_chi2_mf->Fill(chi2_m);
		    
		    RhoCandidate *jfit = D0star[j]->GetFit();// access the fitted cand
		    D0star_chi2_prob_mf->Fill(prob_m);  
		    D0star_mcf->Fill(jfit->M());
		    D0star_mcf_res->Fill(jfit->M() - lv.M());
		     }
		     ntp4->DumpData();
		     
		     
		  }

		
		//ntp1
		Ds2535.SetType(10433);
		gammas.Select(GammaMomentum);
		int nd2535mct = theAnalysis->McTruthMatch(Ds2535);
		
		for (j=0;j<Ds2535.GetLength();++j) 
		  {

		    if (theAnalysis->McTruthMatch(Ds2535[j])){ 
		    // some general info about event (actually written for each candidate!)
		    ntp1->Column("ev",		(Float_t) i);
		    ntp1->Column("cand",	(Float_t) j);
		    ntp1->Column("ncand",       (Float_t) Ds2535.GetLength());
		    ntp1->Column("nmct",        (Float_t) nd2535mct);
		    
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp1);
		   
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Ds2536", Ds2535[j], ntp1);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp1);
		    
		    
		    // store the 4-vector of the truth matched candidate 
		    RhoCandidate *truth = Ds2535[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDs2536", lv, ntp1);


   

		    //missing mass to evaluate p4 of D0* --->D0 K+

		    
		      RhoCandidate *kplusFromD0 = Ds2535[j]->Daughter(1);
		      RhoCandidate *trueKFromD0 = kplusFromD0->GetMcTruth();

		      //cout<<trueKFromD0.M()<<endl;
		      TLorentzVector kv = trueKFromD0->P4();
		      //if (trueKFromD0) kv = trueKFromD0->P4();
		      
		      //Kaon: p4 components 
		      pKx = kplusFromD0->Px();
		      pKy = kplusFromD0->Py();
		      pKz = kplusFromD0->Pz();
		      EK  = kplusFromD0->E();
		      Kcharge = kplusFromD0->Charge();
		      
		      TLorentzVector mykplus(pKx,pKy,pKz,EK);
		      
		      //D0* p4 components, as missing mass of the event Ds(2535)+
		      
		      //reco components
		      pKxm = Ds2535[j]->Px() - pKx;
		      pKym = Ds2535[j]->Py() - pKy;
		      pKzm = Ds2535[j]->Pz() - pzD;
		      EKm =  Ds2535[j]->E() - EK;
		      //true values of D0* as missing mass of the event
		      pKxmt = lv.X() - kv.X();
		      pKymt = lv.Y() - kv.Y();
		      pKzmt = lv.Z() - kv.Z();
		      EKmt =  lv.E() - kv.E();
		      
		      
		      TLorentzVector missingKaon(pKxm,pKym,pKzm,EKm); //new D0* components
		      
		      
		      
		      Kaon_missingmass_ftm->Fill( missingKaon.M() );
		      Kaon_missingmass_diff_ftm->Fill( lv.M() - missingKaon.M() );
		      Kaon_missing_px_ftm->Fill( missingKaon.Px());
		      Kaon_missing_py_ftm->Fill( missingKaon.Py());
		      Kaon_missing_pz_ftm->Fill( missingKaon.Pz());
		      Kaon_missing_E_ftm->Fill( missingKaon.E());
		      Kaon_missing_px_diff_ftm->Fill( pKxmt - pKxm);
		      Kaon_missing_py_diff_ftm->Fill( pKymt - pKym);
		      Kaon_missing_pz_diff_ftm->Fill( pKzmt - pKzm);
		      Kaon_missing_E_diff_ftm->Fill( EKmt - EKm);

		     
		      
		      // cout<<missingD2535.M()<<endl;
		    }
		      ntp1->DumpData();
		  }
			

	}
     
		
	// *** write out all the histos
	out->cd();
	
	hvpos2->Write();
	hm_diff2->Write();
	Ds_vf_mass->Write();
	Ds_chi2_prob_vf->Write();
	Ds_vtx_z->Write();
	Ds_vtx_diff_X->Write();
	Ds_vtx_diff_Y->Write();
	Ds_vtx_diff_Z->Write();
	Ds_vtx_Pos->Write();
	Ds_vtx_PosMC->Write();
	Ds_chi2_vf->Write();
	Ds_missingmass_ftm->Write();
	Ds_missingmass_diff_ftm->Write();
	Ds_px_missing_ftm->Write();
	Ds_py_missing_ftm->Write();
	Ds_pz_missing_ftm->Write();
	Ds_E_missing_ftm->Write();
	Ds_mass_diff_ftm->Write();
	Ds_missingmass_diff_ftm->Write();
	Ds_px_missing_ftm->Write();
        Ds_py_missing_ftm->Write();
        Ds_pz_missing_ftm->Write();
        Ds_E_missing_ftm->Write();
	Ds_mass_diff_ftm->Write();

	D0star_chi2_mf->Write();
	D0star_mcf->Write();
	D0star_mcf_res->Write();
	D0star_chi2_prob_mf->Write();
	
	Kaon_missingmass_ftm->Write();
	Kaon_missingmass_diff_ftm->Write();
	Kaon_missing_px_ftm->Write();
	Kaon_missing_py_ftm->Write();
	Kaon_missing_pz_ftm->Write();
	Kaon_missing_E_ftm->Write();
	Kaon_missing_px_diff_ftm->Write();
	Kaon_missing_py_diff_ftm->Write();
	Kaon_missing_pz_diff_ftm->Write();
	Kaon_missing_E_diff_ftm->Write();


	ntp1->GetInternalTree()->Write();
	ntp2->GetInternalTree()->Write();
	ntp3->GetInternalTree()->Write();
	ntp4->GetInternalTree()->Write();
	nmc->GetInternalTree()->Write();
	ntp5->GetInternalTree()->Write();
	out->Save();
	
}








