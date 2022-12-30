class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
#include <iostream>
using std::cout;
using std::endl;
void ana_D2317(TString InFile="Ds2317_fast.root", int nevts=10000, double pbarmom = 9.83)
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
	RhoTuple *ntp1 = new RhoTuple("ntp1", "pi0 analysis");
	RhoTuple *ntp2 = new RhoTuple("Dsm", "Ds- analysis");
	RhoTuple *ntp3 = new RhoTuple("Dsp", "Ds+ analysis");
	RhoTuple *ntp4 = new RhoTuple("Ds2317", " Ds2317+ analysis");
	RhoTuple *nmc  = new RhoTuple("nmc",  "mctruth info");
	RhoTuple *ntp5 = new RhoTuple("ntp5", "Vertex fitter");
	
	//Histograms for vertex fitters

	
	//PndKinVtxFitter
	TH1F *Dsp_vtx_diff_X  = new TH1F("Dsp_vtx_diff_X", "Ds+: vertex resolution, X",100,-0.5,0.5);
	TH1F *Dsp_vtx_diff_Y = new TH1F("Dsp_vtx_diff_Y", "Ds+: vertex resolution, Y",100,-0.5,0.5);
	TH1F *Dsp_vtx_diff_Z = new TH1F("Dsp_vtx_diff_Z", "Ds+: vertex resolution, Z",100,-0.5,0.5);
	TH1F *Dsp_vtx_Pos = new TH1F("Dsp_vtx_diff_Pos","Ds_vtx_diff_Pos",100,-0.5,0.5);
	TH1F *Dsp_vtx_PosMC = new TH1F("Dsp_vtx_diff_PosMC","Ds_vtx_diff_PosMC",100,-4,4);
	TH1F *Dsp_chi2_vf  = new TH1F("Dsp_chi2_vf", "Ds+: #chi^{2} vertex fit",100,-100,100);
	TH1F *Dsp_chi2_prob_vf   = new TH1F("Dsp_chi2_prob_vf",  "Ds+: prob vertex fit",100,-0.01,1);
	TH1F *Dsp_vf_mass   = new TH1F("Dsp_vf_mass","Ds+ mass (vertex fit)",100,1,3);
	TH2F *hvpos2dsp = new TH2F("hvpos2dsp","(x,y) projection of fitted decay vertex Ds+",100,-2,2,100,-2,2);
	TH1F *Dsp_vtx_z   = new TH1F("Dsp_vtx_z",  "Ds+: Z",100,-0.5,0.5);
	TH1F *hm_diff2dsp  = new TH1F("hm_diff2dsp","Ds+ mass diff to truth after vertex fit",100,-0.5,0.5);


	TH1F *Dsm_vtx_diff_X  = new TH1F("Dsm_vtx_diff_X", "Ds-: vertex resolution, X",100,-0.5,0.5);
	TH1F *Dsm_vtx_diff_Y = new TH1F("Dsm_vtx_diff_Y", "Ds-: vertex resolution, Y",100,-0.5,0.5);
	TH1F *Dsm_vtx_diff_Z = new TH1F("Dsm_vtx_diff_Z", "Ds-: vertex resolution, Z",100,-0.5,0.5);
	TH1F *Dsm_vtx_Pos = new TH1F("Dsm_vtx_diff_Pos","Ds_vtx_diff_Pos",100,-0.5,0.5);
	TH1F *Dsm_vtx_PosMC = new TH1F("Dsm_vtx_diff_PosMC","Ds_vtx_diff_PosMC",100,-4,4);
	TH1F *Dsm_chi2_vf  = new TH1F("Dsm_chi2_vf", "Ds-: #chi^{2} vertex fit",100,-100,100);
	TH1F *Dsm_chi2_prob_vf   = new TH1F("Dsm_chi2_prob_vf",  "Ds: prob vertex fit",100,-0.01,1);
	TH1F *Dsm_vf_mass   = new TH1F("Dsm_vf_mass","Ds- mass (vertex fit)",100,1,3);
	TH2F *hvpos2dsm = new TH2F("hvpos2dsm","(x,y) projection of fitted decay vertex Ds-",100,-2,2,100,-2,2);
	TH1F *Dsm_vtx_z   = new TH1F("Dsm_vtx_z",  "Ds: Z",100,-0.5,0.5);
	TH1F *hm_diff2dsm  = new TH1F("hm_diff2dsm","Ds- mass diff to truth after vertex fit",100,-0.5,0.5);


	TH1F *pi0px  = new TH1F("pi0px ", "pi0, px",100,-0.5,0.5);
	TH1F *pi0py = new TH1F("pi0py", "pi0, py",100,-0.5,0.5);
	TH1F *pi0pz = new TH1F("pi0pz", "pi0, pz",100,-0.5,0.5);
	TH1F *pi0mass = new TH1F("pi0mass","pi0 mass",100,0.,0.5);

	
	TH1F *Ds_missingmass_ftm  = new TH1F("Ds_missingmass_ftm","Ds1(2535): Ds- missing mass (full truth match)",100,2.,2.6);
	TH1F *Ds_px_missing_ftm  = new TH1F("Ds_px_missing_ftm","Ds missing px ",100,-1.5,1.5);
	TH1F *Ds_py_missing_ftm  = new TH1F("Ds_py_missing_ftm","Ds missing py ",100,-1.5,1.5);
	TH1F *Ds_pz_missing_ftm  = new TH1F("Ds_pz_missing_ftm","Ds missing pz ",100,-1,9);
	TH1F *Ds_E_missing_ftm  = new TH1F("Ds_E_missing ftm","Ds missing E ",100,0,10);
	TH1F *Ds_mass_diff_ftm  = new TH1F("Ds_mass_diff_ftm","Ds E ",100,-0.2,0.2);
	TH1F *Ds_missingmass_diff_ftm  = new TH1F("Ds_missingmass_diff_ftm","Ds E ",100,-0.2,0.2);

	TH1F *Kaon_missingmass_ftm  = new TH1F("Kaon_missingmass_ftm","D0*mass as missing mass (full truth match)",100,-0.1,0.5);
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

	RhoMassParticleSelector *pi0MassSel=new RhoMassParticleSelector("pi0",0.1349,0.5);
	RhoMassParticleSelector *DsMassSel=new RhoMassParticleSelector("D_s-",1.9686,0.5);
	RhoMassParticleSelector *Ds2317MassSel=new RhoMassParticleSelector("D_s0*+",2.3178,0.5);//10431
	RhoMassParticleSelector *DsMassSelPrefit=new RhoMassParticleSelector("D_s+",1.9686,0.95);
	RhoMomentumParticleSelector *GammaMomentum=new 	RhoMomentumParticleSelector("gamma",1.0,0.99);	
	// *** the lorentz vector of the initial psi(2S)
	//double m0_p    = TDatabasePDG::Instance()->GetParticle("proton")->Mass();   // Get nominal PDG mass of the proton
	//	TLorentzVector ini(0, 0, pbarmom, sqrt(m0_p*m0_p + pbarmom*pbarmom) + m0_p);
	TLorentzVector ini(0, 0, 8.805,9.793);/////aggiustala

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
	Double_t pDx, pDy, pDz, ED, Dcharge;
	Double_t pPix, pPiy, pPiz, EPi;
        Double_t pPixt, pPiyt, pPizt, EPit;


	RhoCandList kplus, kminus, piplus, piminus,  DslistM, DslistP,D2317,pi0list, gammas, all, mclist, newpi0list;

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

		// *** combinatorics for Ds- to K+ K- pi-
		pi0list.Combine(gammas,gammas);
		// *** combinatorics for Ds- to K+ K- pi-
		DslistM.Combine(kminus,kplus,piminus);
		// *** combinatorics for Ds+ to K+ K- pi+
		DslistP.Combine(kminus,kplus,piplus);
		
		// *** combinatorics for Ds1+ to D0* K+
		D2317.Combine(DslistP,pi0list);

		newpi0list.Combine(DslistP,DslistM);
		
		
		
		DslistM.SetType(-431);
		DslistM.Select(DsMassSel);
		
		int ndsmct = theAnalysis->McTruthMatch(DslistM);
		 


		for (j=0;j<DslistM.GetLength();++j) 
		  {
		    //cout<<"mom2="<<DslistM[j]->PdgCode()<<endl;
		    // some general info about event (actually written for each candidate!)
		    ntp2->Column("ev",		(Float_t) i);
		    ntp2->Column("cand",	(Float_t) j);
		    ntp2->Column("ncand",       (Float_t) DslistM.GetLength());
		    ntp2->Column("nmct",        (Float_t) ndsmct);
		    

		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp2);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Dsm", DslistM[j], ntp2);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp2);
		    
		    
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = DslistM[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDsm", lv, ntp2);
		    
	
		  


		    //Vertex fitters: PndKinVtxFitter on Ds-
		    
		    
		    PndKinVtxFitter vtxfitter(DslistM[j]);	// instantiate a vertex fitter
		    vtxfitter.Fit();
		    
		    double chi2_vtx=vtxfitter.GetChi2();	// access chi2 of fit
		    double prob_vtx=vtxfitter.GetProb();
		    Dsm_chi2_vf->Fill(chi2_vtx);
		    Dsm_chi2_prob_vf->Fill(prob_vtx);


		    
		    // if (chi2_vtx>0.&&chi2_vtx<1.)				
		    //{
		    RhoCandidate *jfit = DslistM[j]->GetFit();	// access the fitted cand	      
		    TVector3 jVtx=jfit->Pos();		// and the decay vertex position
		    
		    Dsm_vf_mass->Fill(jfit->M());    		      
		    hvpos2dsm->Fill(jVtx.X(),jVtx.Y());
		    Dsm_vtx_z->Fill(jVtx.Z());
		    
		    RhoCandidate *mcdsm = DslistM[j]->GetMcTruth();
		    if (mcdsm && mcdsm->PdgCode()==-431){
		      hm_diff2dsm->Fill( mcdsm->M() - jfit->M());
		    }
		    RhoCandidate *mcvtx = DslistM[j]->GetMcTruth();
		    if (mcdsm && mcdsm->PdgCode()==-431){
		      Dsm_vtx_diff_X->Fill( mcdsm->Pos().X() - jVtx.X());
		      Dsm_vtx_diff_Y->Fill( mcdsm->Pos().Y() - jVtx.Y());
		      Dsm_vtx_diff_Z->Fill( mcdsm->Pos().Z() - jVtx.Z());
		      Dsm_vtx_diff_PosMC->Fill(  TMath::Sqrt(mcdsm->Pos().X()*mcdsm->Pos().X()+ 
							     mcdsm->Pos().Y()*mcdsm->Pos().Y() + 
							     mcdsm->Pos().Z()*mcdsm->Pos().Z())) ;
		      Dsm_vtx_diff_Pos->Fill(TMath::Sqrt(jVtx.X()*jVtx.X()+jVtx.Y()*jVtx.Y()++jVtx.Z()*jVtx.Z()));
		    }	
		    

		 
			      //Add here : if (theAnalysis->McTruthMatch(Dslist[j])) if you want the true values ftm.
		    if (theAnalysis->McTruthMatch(DslistM[j]))
		      { 	    
			//missing mass: I will get a better Ds2317+ here!
			
			RhoCandidate *kpluscand = DslistM[j]->Daughter(0);
			RhoCandidate *kminuscand = DslistM[j]->Daughter(1);
			RhoCandidate *piminuscand = DslistM[j]->Daughter(2);
			
			//Ds- bachelor 4mom
			pxDs = DslistM[j]->Px();
			pyDs = DslistM[j]->Py();
			pzDs = DslistM[j]->Pz();
			EDs  = DslistM[j]->E();
			
			//Ds2317+ 4mom as missing mass of Ds-
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
			TLorentzVector truevalueDs2317(pxt,pyt,pzt,Et); //Ds(2535)+ true values
			
			
			Ds_missingmass_ftm->Fill(missingmass.M());
			Ds_px_missing_ftm->Fill(missingmass.Px());
			Ds_py_missing_ftm->Fill(missingmass.Py());
			Ds_pz_missing_ftm->Fill(missingmass.Pz());
			Ds_E_missing_ftm->Fill(missingmass.E());
			Ds_missingmass_diff_ftm->Fill(truevalueDs2317.M() - missingmass.M());   
			
			
			//now I want to get a pi0 as missing mass of Ds+, daughter of Ds2317+.
			//pi0 = ini - Ds- - Ds+
			//not yet ready or this

			//cout<<"Ds2317="<<missingmass.M()<<endl;
			ntp2->DumpData();
		      }
		  }


		//ntp3
		DslistP.SetType(431);
		
		int ndspct = theAnalysis->McTruthMatch(DslistP);
		
		for (j=0;j<DslistP.GetLength();++j) 
		  {
		    // some general info about event (actually written for each candidate!)
		    ntp3->Column("ev2",		(Float_t) i);
		    ntp3->Column("cand2",	(Float_t) j);
		    ntp3->Column("ncand2",       (Float_t) DslistP.GetLength());
		    ntp3->Column("nmct2",        (Float_t) ndspct);
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam2", ini, ntp3);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Dsp", DslistP[j], ntp2);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp2);
		   
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = DslistP[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDsp", lv, ntp2);


	
		//Vertex fitters: PndKinVtxFitter on Ds+

	
		  PndKinVtxFitter vtxfitter(DslistP[j]);	// instantiate a vertex fitter
		  vtxfitter.Fit();
		  
		  double chi2_vtx=vtxfitter.GetChi2();	// access chi2 of fit
		  double prob_vtx=vtxfitter.GetProb();
		  Dsp_chi2_vf->Fill(chi2_vtx);
		  Dsp_chi2_prob_vf->Fill(prob_vtx);
			  
		  // if (chi2_vtx>0.&&chi2_vtx<1.)				
		       //{
		      RhoCandidate *jfit = DslistP[j]->GetFit();	// access the fitted cand	      
		      TVector3 jVtx=jfit->Pos();		// and the decay vertex position
		      
		      Dsp_vf_mass->Fill(jfit->M());    		      
		      hvpos2dsp->Fill(jVtx.X(),jVtx.Y());
		      Dsp_vtx_z->Fill(jVtx.Z());
		      
		      RhoCandidate *mcdsp = DslistP[j]->GetMcTruth();
		      if (mcdsp && mcdsp->PdgCode()==431){
			hm_diff2dsp->Fill( mcdsp->M() - jfit->M());
		      }
		      RhoCandidate *mcvtx = DslistP[j]->GetMcTruth();
		      if (mcdsp && mcdsp->PdgCode()==431){
			Dsp_vtx_diff_X->Fill( mcdsp->Pos().X() - jVtx.X());
			Dsp_vtx_diff_Y->Fill( mcdsp->Pos().Y() - jVtx.Y());
			Dsp_vtx_diff_Z->Fill( mcdsp->Pos().Z() - jVtx.Z());
			Dsp_vtx_diff_PosMC->Fill(  TMath::Sqrt(mcdsp->Pos().X()*mcdsp->Pos().X()+ mcdsp->Pos().Y()*mcdsp->Pos().Y() + mcdsp->Pos().Z()*mcdsp->Pos().Z())) ;
			Dsp_vtx_diff_Pos->Fill(TMath::Sqrt(jVtx.X()*jVtx.X()+jVtx.Y()*jVtx.Y()+jVtx.Z()*jVtx.Z()));
		      }	
		      //}
		      //cout<<"Dsminus mass="<<jfit->M()<<endl;

		      TLorentzVector mypi0(pxmm-jfit->Px(),pymm-jfit->Py(),pzmm-jfit->Pz(),Emm-jfit->E());
		      pi0px ->Fill(mypi0.X());
		      pi0py ->Fill(mypi0.Y());
		      pi0pz ->Fill(mypi0.Z());
		      pi0mass ->Fill(mypi0.M());

		      
			  ntp2->DumpData();

		  }
	

		
		//ntp4
		D2317.SetType(10431);
		gammas.SetType(22);
		pi0list.SetType(111);
		gammas.Select(GammaMomentum);
		D2317.Select(Ds2317MassSel);
		pi0list.Select(pi0MassSel);
	
		int nd2317mct = theAnalysis->McTruthMatch(D2317);

		for (j=0;j<D2317.GetLength();++j) 
		  {

		    if (theAnalysis->McTruthMatch(D2317[j])){ 
		   
		    ntp4->Column("ev",		(Float_t) i);
		    ntp4->Column("cand",	(Float_t) j);
		    ntp4->Column("ncand",       (Float_t) D2317.GetLength());
		    ntp4->Column("nmct",        (Float_t) nd2317mct);
		    
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp4);
		   
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Ds2317", D2317[j], ntp4);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp4);
		    
		    
		    // store the 4-vector of the truth matched candidate 
		    RhoCandidate *truth = D2317[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDs2317", lv, ntp4);


		    //missing mass to evaluate p4 of pi0 

		    
		    RhoCandidate *dsplusFromD2317 = D2317[j]->Daughter(0); //ds+
		    RhoCandidate *pi0FromD2317 = D2317[j]->Daughter(1);  //pi0
		    RhoCandidate *trueDsFromD2317 = dsplusFromD2317->GetMcTruth();
		    RhoCandidate *truepi0FromD2317 = pi0FromD2317->GetMcTruth();
		    
		     
		    TLorentzVector ds = trueDsFromD2317->P4();
		    TLorentzVector pi = truepi0FromD2317->P4();
		   
		      
		    //ds+: p4 components 
		    pDx = dsplusFromD2317->Px();
		    pDy = dsplusFromD2317->Py();
		    pDz = dsplusFromD2317->Pz();
		    ED  = dsplusFromD2317->E();
		    Dcharge = dsplusFromD2317->Charge();
		      
		    TLorentzVector mydsp(pDx,pDy,pDz,ED);
		      
		    
		    
		    //reco components
		    pPix = D2317[j]->Px() - pDx;   //pi0_px
		    pPiy = D2317[j]->Py() - pDy;   //pi0_py
		    pPiz = D2317[j]->Pz() - pDz;   //pi0_pz
		    EPi  =  D2317[j]->E() - ED;     //pi0_E
		    //true values of pi0 as missing mass of the event
		    pPixt = lv.X() - ds.X();
		    pPiyt = lv.Y() - ds.Y();
		    pPizt = lv.Z() - ds.Z();
		    EPit =  lv.E() - ds.E();
		    
		      
		    TLorentzVector missingpi0(pPix,pPiy,pPiz,EPi); //new pi0 reco components
		      
		      
		    
		    Kaon_missingmass_ftm->Fill( missingpi0.M() );
		    Kaon_missingmass_diff_ftm->Fill( lv.M() - missingpi0.M() );
		    Kaon_missing_px_ftm->Fill( missingpi0.Px());
		    Kaon_missing_py_ftm->Fill( missingpi0.Py());
		    Kaon_missing_pz_ftm->Fill( missingpi0.Pz());
		    Kaon_missing_E_ftm->Fill( missingpi0.E());
		    Kaon_missing_px_diff_ftm->Fill( pPixt - pPix);
		    Kaon_missing_py_diff_ftm->Fill( pPiyt - pPiy);
		    Kaon_missing_pz_diff_ftm->Fill( pPizt - pPiz);
		    Kaon_missing_E_diff_ftm->Fill( EPit - EPi);

		    
		    
		    ntp4->DumpData();
		    }
		  }





		int ndspct = theAnalysis->McTruthMatch(DslistP);
		
		for (j=0;j<DslistP.GetLength();++j) 
		  {
		    // some general info about event (actually written for each candidate!)
		    ntp3->Column("ev",		(Float_t) i);
		    ntp3->Column("cand",	(Float_t) j);
		    ntp3->Column("ncand",       (Float_t) DslistP.GetLength());
		    ntp3->Column("nmct",        (Float_t) ndspct);
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp3);
		    
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("Dsp", DslistP[j], ntp3);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp3);
		   
		    // store the 4-vector of the truth matched candidate (or a dummy, if not matched to keep ntuple consistent)
		    RhoCandidate *truth = DslistP[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCDsp", lv, ntp3);


	
		//Vertex fitters: PndKinVtxFitter on Ds+

	
		  PndKinVtxFitter vtxfitter(DslistP[j]);	// instantiate a vertex fitter
		  vtxfitter.Fit();
		  
		  double chi2_vtx=vtxfitter.GetChi2();	// access chi2 of fit
		  double prob_vtx=vtxfitter.GetProb();
		  Dsp_chi2_vf->Fill(chi2_vtx);
		  Dsp_chi2_prob_vf->Fill(prob_vtx);
			  
		  // if (chi2_vtx>0.&&chi2_vtx<1.)				
		       //{
		      RhoCandidate *jfit = DslistP[j]->GetFit();	// access the fitted cand	      
		      TVector3 jVtx=jfit->Pos();		// and the decay vertex position
		      
		      Dsp_vf_mass->Fill(jfit->M());    		      
		      hvpos2dsp->Fill(jVtx.X(),jVtx.Y());
		      Dsp_vtx_z->Fill(jVtx.Z());
		      
		      RhoCandidate *mcdsp = DslistP[j]->GetMcTruth();
		      if (mcdsp && mcdsp->PdgCode()==431){
			hm_diff2dsp->Fill( mcdsp->M() - jfit->M());
		      }
		      RhoCandidate *mcvtx = DslistP[j]->GetMcTruth();
		      if (mcdsp && mcdsp->PdgCode()==431){
			Dsp_vtx_diff_X->Fill( mcdsp->Pos().X() - jVtx.X());
			Dsp_vtx_diff_Y->Fill( mcdsp->Pos().Y() - jVtx.Y());
			Dsp_vtx_diff_Z->Fill( mcdsp->Pos().Z() - jVtx.Z());
			Dsp_vtx_diff_PosMC->Fill(  TMath::Sqrt(mcdsp->Pos().X()*mcdsp->Pos().X()+ mcdsp->Pos().Y()*mcdsp->Pos().Y() + mcdsp->Pos().Z()*mcdsp->Pos().Z())) ;
			Dsp_vtx_diff_Pos->Fill(TMath::Sqrt(jVtx.X()*jVtx.X()+jVtx.Y()*jVtx.Y()+jVtx.Z()*jVtx.Z()));
		      }
		  
		      //}
			  
			  ntp3->DumpData();

		  }
	
	
		
		//ntp1

		int npi0mct = theAnalysis->McTruthMatch(newpi0list);
		
		for (j=0;j<newpi0list.GetLength();++j) 
		  {

		   
		    //   cout<<"mom="<<newpi0list[j]->PdgCode()<<endl;
		    ntp4->Column("ev",		(Float_t) i);
		    ntp4->Column("cand",	(Float_t) j);
		    ntp4->Column("ncand",       (Float_t) newpi0list.GetLength());
		    ntp4->Column("nmct",        (Float_t) npi0mct);
		    
		    
		    // store info about initial 4-vector
		    qa.qaP4("beam", ini, ntp1);
		   
		    // store information about composite candidate tree recursively (see PndTools/AnalysisTools/PndRhoTupleQA)
		    qa.qaComp("pi0", newpi0list[j], ntp1);
		    
		    // store info about event shapes
		    qa.qaEventShapeShort("es",&evsh, ntp1);
		    
		    
		    // store the 4-vector of the truth matched candidate 
		    RhoCandidate *truth = newpi0list[j]->GetMcTruth();		
		    TLorentzVector lv;
		    if (truth) lv = truth->P4();
		    qa.qaP4("MCpi0", lv, ntp4);


		    
		    ntp1->DumpData();
		    
		  }
	}
     
		
	// *** write out all the histos
	out->cd();
	
	Dsp_vtx_diff_X->Write();
	Dsp_vtx_diff_Y->Write();
	Dsp_vtx_diff_Z->Write();
	Dsp_vtx_Pos->Write();
	Dsp_vtx_PosMC->Write();
	Dsp_chi2_vf->Write();
	Dsp_chi2_prob_vf->Write();
	Dsp_vf_mass->Write();
	hvpos2dsp->Write();
	Dsp_vtx_z->Write();
	hm_diff2dsp->Write();

	Dsm_vtx_diff_X->Write();
	Dsm_vtx_diff_Y->Write();
	Dsm_vtx_diff_Z->Write();
	Dsm_vtx_Pos->Write();
	Dsm_vtx_PosMC->Write();
	Dsm_chi2_vf->Write();
	Dsm_chi2_prob_vf->Write();
	Dsm_vf_mass->Write();
	hvpos2dsm->Write();
	Dsm_vtx_z->Write();
	hm_diff2dsm->Write();

	Ds_missingmass_ftm->Write();
	Ds_px_missing_ftm->Write();
	Ds_py_missing_ftm->Write();
	Ds_pz_missing_ftm->Write();
	Ds_E_missing_ftm->Write();
	Ds_mass_diff_ftm->Write();
	Ds_missingmass_diff_ftm->Write();
	Kaon_missingmass_ftm->Write();
	Kaon_missing_px_ftm->Write();
	Kaon_missing_py_ftm->Write();
	Kaon_missing_pz_ftm->Write();
	Kaon_missing_px_diff_ftm->Write();
	Kaon_missing_py_diff_ftm->Write();
	Kaon_missing_pz_diff_ftm->Write();
	Kaon_missing_E_diff_ftm->Write();
	Kaon_missing_charge->Write();
	Kaon_missingmass_diff_ftm->Write();

	pi0px->Write();
	pi0py->Write();
	pi0pz->Write();
	pi0mass->Write();

	ntp1->GetInternalTree()->Write();
	ntp2->GetInternalTree()->Write();
	ntp3->GetInternalTree()->Write();
	ntp4->GetInternalTree()->Write();
	nmc->GetInternalTree()->Write();
	ntp5->GetInternalTree()->Write();
	out->Save();
	
}







