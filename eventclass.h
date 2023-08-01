#include "/work/halla/sbs/adr/GMn_analysis/physics_analysis/ElasticEventsStudy/includes/beam_variables.h"
#include "/work/halla/sbs/adr/GMn_analysis/physics_analysis/ElasticEventsStudy/includes/constants.h"
#include "/work/halla/sbs/adr/GMn_analysis/physics_analysis/ElasticEventsStudy/includes/calc_HCalintersect.h"
#include "/work/halla/sbs/adr/GMn_analysis/physics_analysis/ElasticEventsStudy/includes/HCalConstants.h"
#include <Math/Vector4D.h>
#include "TVector3.h"
#include "TMath.h"

#ifndef EVENTCLASS_H
#define EVENTCLASS_H

class Event
{
private:
	
	// Read-in from the main script.
	TChain* m_C; // TChain that holds the "T" TTree(s) of input replayed ROOT file(s).
	const int m_SBSKineNum; 
	const double m_cut_BBTrVzCut;
	const double m_cut_HCalE;
	const double m_cut_CoinCutLow;
	const double m_cut_CoinCutHigh;
	const double m_cut_RCut;
	
	// Variables needed to copy information from "T" tree.
	double m_bbtrn{0};
	const static int m_MAXNTRACKS {100};
	//double m_bbtrnhits[m_]
	double m_bbtrvz[m_MAXNTRACKS];
	double m_bbtrth[m_MAXNTRACKS];
	double m_bbtrx[m_MAXNTRACKS];
	double m_bbtrp[m_MAXNTRACKS];
	double m_bbtrpx[m_MAXNTRACKS];
	double m_bbtrpy[m_MAXNTRACKS];
	double m_bbtrpz[m_MAXNTRACKS];
	double m_sbshcale{0.};
	double m_sbshcalx{0.};
	double m_sbshcaly{0.};
	double m_bbshe{0.};
	double m_bbpse{0.};
	double m_hcaladctime{0.};
	double m_shadctime{0.};
	double m_adctimediff_hcalsh{0.};
	
	// Kinematic parameters.
	GMnKinInfo m_kininfo;
	const double m_Ebeam;
	const double m_photonecut;
	const double m_HCaldist;
	const double m_HCalangle;

	// Calss to calculate intersection point of the neutron.
	HCalVectors m_hcal;
		
public:

	Event(TChain* C, int num_SBSKine, double cut_BBTrVzCut, double cut_HCalE, double cut_CoinCutLow, double cut_CoinCutHigh, double cut_RCut) 
	: m_C{C}, m_SBSKineNum{num_SBSKine}, m_cut_BBTrVzCut{cut_BBTrVzCut}, m_cut_HCalE{cut_HCalE}, m_cut_CoinCutLow{cut_CoinCutLow}, m_cut_CoinCutHigh{cut_CoinCutHigh}, m_cut_RCut{cut_RCut}, m_kininfo{num_SBSKine}, m_Ebeam{m_kininfo.return_BeamEnergy()}, m_photonecut{m_kininfo.return_BeamEnergy()+m_kininfo.return_BeamEnergy()*0.03}, m_HCaldist{m_kininfo.return_HCalDis()}, m_HCalangle{m_kininfo.return_HCalTheta()}
	{
		m_C->SetBranchStatus("*", 0);
		m_C->SetBranchStatus("bb.tr.n",1);
		m_C->SetBranchStatus("bb.tr.vz",1);
		m_C->SetBranchStatus("bb.tr.th", 1);
		m_C->SetBranchStatus("bb.tr.x", 1);
		m_C->SetBranchStatus("bb.tr.p",1);
		m_C->SetBranchStatus("bb.tr.px", 1);
		m_C->SetBranchStatus("bb.tr.py", 1);
		m_C->SetBranchStatus("bb.tr.pz", 1);
		m_C->SetBranchStatus("sbs.hcal.x", 1);
		m_C->SetBranchStatus("sbs.hcal.y", 1);
		m_C->SetBranchStatus("sbs.hcal.e",1);
		m_C->SetBranchStatus("bb.sh.e",1);
		m_C->SetBranchStatus("bb.ps.e", 1);
		m_C->SetBranchStatus("sbs.hcal.atimeblk",1);
		m_C->SetBranchStatus("bb.sh.atimeblk",1);

		m_C->SetBranchAddress("bb.tr.n", &m_bbtrn);
		m_C->SetBranchAddress("bb.tr.vz", m_bbtrvz);
		m_C->SetBranchAddress("bb.tr.th", m_bbtrth);
		m_C->SetBranchAddress("bb.tr.x", m_bbtrx);
		m_C->SetBranchAddress("bb.tr.p", m_bbtrp);
		m_C->SetBranchAddress("bb.tr.px", m_bbtrpx); 
		m_C->SetBranchAddress("bb.tr.py", m_bbtrpy); 
		m_C->SetBranchAddress("bb.tr.pz", m_bbtrpz);
		m_C->SetBranchAddress("sbs.hcal.x", &m_sbshcalx);
		m_C->SetBranchAddress("sbs.hcal.y", &m_sbshcaly);
		m_C->SetBranchAddress("sbs.hcal.e", &m_sbshcale);
		m_C->SetBranchAddress("bb.sh.e", &m_bbshe);
		m_C->SetBranchAddress("bb.ps.e", &m_bbpse);
		m_C->SetBranchAddress("sbs.hcal.atimeblk", &m_hcaladctime);
		m_C->SetBranchAddress("bb.sh.atimeblk", &m_shadctime);

		m_hcal.make_HCal_vectors_NoVerticalOffset(m_HCaldist, m_HCalangle);
	}

	int getEntry(int n) //Copies the enries of the "T" to the above defined member variables.
	{
		return m_C->GetEntry(n);
	}

private:	

	bool pass_HCalSHADCtime_coincut()
	{
		m_adctimediff_hcalsh = m_hcaladctime - m_shadctime;

	    if ( m_adctimediff_hcalsh < m_cut_CoinCutLow || m_adctimediff_hcalsh > m_cut_CoinCutHigh ) return false;

		return true;
	}
	
	bool pass_EndPointCut() 
	{
		double E_gamma = m_Ebeam;

		// Mandelstam s variable for the photon(with beam energy) + proton(at rest = target) pair.
		double s = std::pow(Constants::p_mass,2) + 2*E_gamma*Constants::p_mass;

		// Establish the CM frame.
		double v_CMframe = E_gamma/(E_gamma + Constants::p_mass); // Velocity of the CM frame w.r.t the lab frame (in the direction of the beam i.e. +Z).
		double beta = v_CMframe;
		double gamm_fac = 1/std::sqrt(1 - std::pow(beta,2));
		double E_tot_CM = std::sqrt(s);

		// photon + p --> pi+ + n
		// Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
		double p_piplus_CM_1pi = std::sqrt( std::pow((s + std::pow(Constants::n_mass,2) - std::pow(Constants::p_mass,2)),2)/(4*s) - std:: pow(Constants::n_mass,2) );
		// Transform the pi+ momentum into lab frame.
		double E_piplus_CM_1pi = std::sqrt(std::pow(p_piplus_CM_1pi,2) + std::pow(Constants::piplus_mass,2));
		double p_piplus_lab_1pi = gamm_fac*(p_piplus_CM_1pi + v_CMframe*E_piplus_CM_1pi); // Maximum possible pi+ momentum in the lab frame.
		double E_piplus_lab_1pi = gamm_fac*(E_piplus_CM_1pi + v_CMframe*p_piplus_CM_1pi);

		// photon + p --> pi+ + pi0 + n 
		// Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
		double m_T = Constants::n_mass + Constants::pizero_mas; // Treat the pi0 and n as a composite system.
		double p_piplus_CM_2pi = std::sqrt( std::pow((s + std::pow(m_T,2) - std::pow(Constants::piplus_mass,2)),2)/(4*s) - std::pow(m_T,2) ); // Maximum pi+ momentum possible for the given beam energy in the CM frame.
		// Transform the pi+ momentum into lab frame.	
		double E_piplus_CM_2pi = std::sqrt(std::pow(p_piplus_CM_2pi,2) + std::pow(Constants::piplus_mass,2));
		double p_piplus_lab_2pi = gamm_fac*(p_piplus_CM_2pi + v_CMframe*E_piplus_CM_2pi); // Maximum possible pi+ momentum in the lab frame.
		double E_piplus_lab_2pi = gamm_fac*(E_piplus_CM_2pi + v_CMframe*p_piplus_CM_2pi);

		// Now we want to extend the calculation to the case where we have non-zero polar angle for theta/theta'
		//m_bbpoltgt = acos(m_bbtrpz[0]/m_bbtrp[0]); //Polar scattering angle of pi+ in radians.
		double tan_theta = tan(m_bbpoltgt); 

		double cos_thetaprime_1pi_pos = ( -v_CMframe*E_piplus_CM_1pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_1pi,2) - std::pow(v_CMframe*E_piplus_CM_1pi,2)) + std::pow(p_piplus_CM_1pi,2))) / ( p_piplus_CM_1pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
		double p_piplus_lab_1pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_1pi,2)*std::pow(cos_thetaprime_1pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_1pi*v_CMframe*E_piplus_CM_1pi*cos_thetaprime_1pi_pos + std::pow(p_piplus_CM_1pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_1pi,2));
		double p_piplus_lab_1pi_theta_max = p_piplus_lab_1pi_theta_pos;  

		double cos_thetaprime_2pi_pos = ( -v_CMframe*E_piplus_CM_2pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_2pi,2) - std::pow(v_CMframe*E_piplus_CM_2pi,2)) + std::pow(p_piplus_CM_2pi,2))) / ( p_piplus_CM_2pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
		double p_piplus_lab_2pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_2pi,2)*std::pow(cos_thetaprime_2pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_2pi*v_CMframe*E_piplus_CM_2pi*cos_thetaprime_2pi_pos + std::pow(p_piplus_CM_2pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_2pi,2));
		double p_piplus_lab_2pi_theta_max = p_piplus_lab_2pi_theta_pos;
		double piplus_momentum_limit = p_piplus_lab_2pi_theta_max*(1 + 1.5/100.0) - 0.2; //Substract 200 MeV/c from the momentum threshold give some room. 

		if ( m_bbtrp[0] < piplus_momentum_limit ) return false; // End-point cut.

		if ( m_bbtrp[0] > p_piplus_lab_1pi_theta_max+0.2 ) return false; // Max possible momentum for single pion production.

		return true; 
	}

public:

	bool passCuts()
	{	
		
		//if ( m_bbtrn > 1 ) return false; // Number of tracks found in BigBite.

		if ( abs(m_bbtrvz[0]) > m_cut_BBTrVzCut ) return false; // BigBite track vertex Z cut.

		//if ( m_bbpse > 0.25 ) return false; // PS energy anti cut: To reject electrons. 

		//// *** Cuts that uses HCal data - DO NOT USE/COMMENT OUT for HCal neutron detection efficiency analysis. *** ////
		// if ( m_sbshcale < m_cut_HCalE ) return false; // HCal energy cut. Only usef for preliminary analysis. NOT TO BE USED IN NEUTRON DET. EFF. CALC.

		if ( !pass_HCalSHADCtime_coincut() ) return false; //BBCal and HCal ADC coincidence time cut. Only usef for preliminary analysis. NOT TO BE USED IN NEUTRON DET. EFF. CALC.
		//pass_HCalSHADCtime_coincut();
		// //// **** ////

		//if ( m_photone > m_photonecut ) return false; // Discard events with reconstructed photon energy larger than beam energy. Possibly the events from different processes?

		if ( !pass_EndPointCut() ) return false; // ***End-point cut***

		if ( !pass_WouldHitHCalCut() ) return false; // Predicted hit position on the neutron must be within the HCal active area (with safety margins considered) for the analysis to progress.
 
		return true;
	}

	double return_BBTrn()
	{
		return m_bbtrn;
	}

	double return_BBTrVz() 
	{
		return m_bbtrvz[0];
	}
	
	double return_BBTrth()
	{
		return m_bbtrth[0];
	}

	double return_BBTrx() 
	{
		return m_bbtrx[0];
	}
	
	double return_BBTrP()
	{
		return m_bbtrp[0];
	}
	
	double return_BBTrPx()
	{
		return m_bbtrpx[0];
	}

	double return_BBtrPy()
	{
		return m_bbtrpy[0];
	}
	
	double return_BBSHe()
	{
		return m_bbshe;
	}
	double return_BBPSe()
	{
		return m_bbpse;
	}

	double return_SBSHCalx()
	{
		return m_sbshcalx;
	}

	double return_SBSHCaly()
	{
		return m_sbshcaly;
	}

	double return_SBSHCale()
	{
		return m_sbshcale;
	}

	double return_ADCTimeDiffHCalSH()
	{
		return m_adctimediff_hcalsh;
	}

private: 	

	// Variables to hold the secondary calculation results.
	double m_bbshpse{0.}; // bb.sh.e + bb.ps.e = Total energy deposited in the sh and ps.
	double m_eoverp{0.};
	double m_bbphitgt{0.}; // Azimuthal angle of pi+ detected by the BigBite spectrometer.
	double m_bbpoltgt{0.}; // Polar angle of pi+ detected by the BigBite spectrometer.

	double m_bbtrpipluse{0.}; // Calculated Pi+'s total relativistic energy from bb.tr.p[0]
	double m_photone{0.}; // Reconstructed photon energy of the reaction: gamma + p --> pi+ + n
	double m_sbsneutronpx{0.}; // Calculated x component of the neutron momentum vector.
	double m_sbsneutronpy{0.}; // Calculated y component of the neutron momentum vector.
	double m_sbsneutronpz{0.}; // Calculated z component of the neutron momentum vector.
	double m_sbsneutrone{0.}; // Calculated energy (total relativistic energy) of the neutron.
	double m_sbsneutronp{0.} ; // Calculated magnitude of the neutron momentum vector.
	double m_sbsneutronhcalx{0.}; // Calculated x hit position of the neutron in HCal local detector coordinates.
	double m_sbsneutronhcaly{0.}; // Calculated y hit position of the neutron in HCal local detector coordinates.
	double m_sbsneutronhcaldx{0.}; // X coordinate difference between the measured and prediceted neutron position values. 
	double m_sbsneutronhcaldy{0.}; // Y coordinate difference between the measured and prediceted neutron position values.	
	
public: 
 
	double return_BBSHPSe()
	{
		m_bbshpse = m_bbshe + m_bbpse;
		return m_bbshpse;
	}

	double return_EoverP() //Make sure to call the above "return_BBSHPSe()" function before calling this function.
	{
		m_eoverp = m_bbshpse / m_bbtrp[0];
		return m_eoverp;
	}

	void calc_BBTrackAngles() //Function to calculate bb track polar and azimuthal scattering angles.
	{
		m_bbphitgt = atan2(m_bbtrpy[0], m_bbtrpx[0]); //Azimuthal scattering angle of pi+ in radians.
		m_bbpoltgt = acos(m_bbtrpz[0]/m_bbtrp[0]); //Polar scattering angle of pi+ in radians.
	}

	double return_BBphitgt()
	{
		double bbphitgt_deg = m_bbphitgt*TMath::RadToDeg();
		return bbphitgt_deg;
	}	 

	double return_BBpoltgt()
	{
		double bbpoltgt_deg = m_bbpoltgt*TMath::RadToDeg();
		return bbpoltgt_deg;
	}

	void calc_PhotonE()	
	{
		double piplus_p = m_bbtrp[0];
		double piplus_E = std::sqrt(std::pow(piplus_p,2) + std::pow(Constants::piplus_mass,2)); // Use the relativistic energy-momentum relation to obtain piplus's energy.
		m_bbtrpipluse = piplus_E; 
		// Formula for recontructed photon energy for the reaction, gamma + p --> pi+ + n, when the initial proton is at rest and the final pi+'s momentum and energy is known.
		m_photone = ( 2*Constants::p_mass*piplus_E + pow(Constants::n_mass,2) - pow(Constants::piplus_mass,2) - pow(Constants::p_mass,2) ) / ( 2*( Constants::p_mass + piplus_p*cos(m_bbpoltgt) - piplus_E ) );
	}	

	double return_PhotonE()
	{
		return m_photone;
	}

	void calc_NeutronKin() // Calculates the Energy/Momentum four-vector of the neutron generated from pion photoproduction: photon + p --> pi+ + n
	{
		ROOT::Math::PxPyPzEVector pPhoton; // 4-momentum of the incident photon.
		ROOT::Math::PxPyPzEVector pProton; // 4-momentum of the target proton.
		ROOT::Math::PxPyPzEVector pPiplus; // 4-momentum of the pi+ as measured by BigBite.
		ROOT::Math::PxPyPzEVector pNeutron; // 4-momentum of the neutron, calculated using the above 3, 4-momentum vectors.

		double gamma_e = m_photone;
		double piplus_px = m_bbtrpx[0];
		double piplus_py = m_bbtrpy[0];
		double piplus_pz = m_bbtrpz[0];
		double piplus_e = m_bbtrpipluse;

		pPhoton.SetPxPyPzE( 0, 0, gamma_e, gamma_e ); // Assumption: High energy photons are moving only in the +Z direction.
		pProton.SetPxPyPzE( 0, 0, 0, Constants::p_mass ); // Target proton is at rest and has only its rest mass energy.
		pPiplus.SetPxPyPzE( piplus_px, piplus_py, piplus_pz, piplus_e ); // 3 - Componets of the momentum as measured from BigBite. Energy calculated from the energy-momentum relationship.

		pNeutron = pPhoton + pProton - pPiplus; // Conversation from energy-momentum four-vector.

		m_sbsneutronpx = pNeutron.Px();
		m_sbsneutronpy = pNeutron.Py();
		m_sbsneutronpz = pNeutron.Pz();
		m_sbsneutrone = pNeutron.E();
	}

	double return_NeutronPx()
	{
		return m_sbsneutronpx;
	}  
	
	double return_NeutronPy()
	{
		return m_sbsneutronpy;
	}   

	double return_NeutronPz()
	{
		return m_sbsneutronpz;
	}

	double return_NeutronE()
	{
		return m_sbsneutrone;
	}

	void calc_NeutronHCalIntersect() // Calculates the point of intersection of the neutron with HCal in HCal local coordinates. To be compared with sbs.hcal.dx and sbs.hcal.dy variables.
	{
		ROOT::Math::PxPyPzEVector pNeutron;
		pNeutron.SetPxPyPzE( m_sbsneutronpx, m_sbsneutronpy, m_sbsneutronpz, m_sbsneutrone );

		m_hcal.calc_expected_NeutronxyonHCal( pNeutron, m_bbtrvz );
		m_sbsneutronhcalx = m_hcal.return_xexpected();
		m_sbsneutronhcaly = m_hcal.return_yexpected();
	}

	double return_HCalNeutronPredXpos()
	{
		return m_sbsneutronhcalx;   
	}
	
	double return_HCalNeutronPredYpos()
	{
		return m_sbsneutronhcaly;
	}

	void calc_NeutronHCaldxdy()
	{
		m_sbsneutronhcaldx = m_sbshcalx - m_sbsneutronhcalx;
		m_sbsneutronhcaldy = m_sbshcaly - m_sbsneutronhcaly;
	}

	double return_HCalNeutrondx()
	{
		return m_sbsneutronhcaldx;
	}

	double return_HCalNeutrondy()
	{
		return m_sbsneutronhcaldy;
	}

private:

	// HCal boundaries with safety margins in HCal coordinate system.
	// Safety margin = Exclude the two outer columns and two outer blocks.
	constexpr static double m_hcal_safetymargin_factor = 1;   
	const double m_hcal_active_xlow_safe =  HCalConst::hcal_topXpos + m_hcal_safetymargin_factor*HCalConst::hcalblk_h;
	const double m_hcal_active_xhigh_safe = HCalConst::hcal_botXpos - m_hcal_safetymargin_factor*HCalConst::hcalblk_h;
	const double m_hcal_active_ylow_safe = HCalConst::hcal_rightYpos + m_hcal_safetymargin_factor*HCalConst::hcalblk_w;
	const double m_hcal_active_yhigh_safe = HCalConst::hcal_leftYpos - m_hcal_safetymargin_factor*HCalConst::hcalblk_w;

public:

	bool pass_WouldHitHCalCut()
	{
		if ( m_sbsneutronhcalx < m_hcal_active_xlow_safe || m_sbsneutronhcalx > m_hcal_active_xhigh_safe || m_sbsneutronhcaly < m_hcal_active_ylow_safe || m_sbsneutronhcaly > m_hcal_active_yhigh_safe )
		{
			return false;
		}

		return true;
	}

// Variables and functions to calculate HCal neutron detection efficiency.
private:

	long m_nhcal_shouldhit{0};
	long m_nhcal_didhit{0};
	bool m_passrcut{false};

public:

	void calc_HCalNDE_Rmthd()
	{
		m_nhcal_shouldhit++; // Increment the should hit count by 1.

		double delta_x = m_sbsneutronhcalx - m_sbshcalx;
		double delta_y = m_sbsneutronhcaly - m_sbshcaly;
		double r = std::sqrt( std::pow(delta_x,2) + std::pow(delta_y,2) );

		if ( r <= m_cut_RCut )
		{
			m_nhcal_didhit++; // Increment the did hit count by 1 if the HCal hit is within the search region.
			m_passrcut = true;
		}
		else
		{
			m_passrcut = false;
		}
	}

	bool return_DidPassRcut()
	{
		return m_passrcut;
	}

	void print_HCalNDE_Rmthd()
	{	
		std::cout << "Should hit: " << m_nhcal_shouldhit << '\n';
		std::cout << "Did hit: " << m_nhcal_didhit << '\n';
		std::cout << "*** HCal NDE from the R method: " << std::fixed << std::setprecision(1) << ((double)m_nhcal_didhit/(double)m_nhcal_shouldhit)*100 << "% ***\n";
	}

};

#endif