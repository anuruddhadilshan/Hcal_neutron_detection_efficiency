#include <iostream>
#include <iomanip>

// Standalone code to calculate pion production thresholds.

constexpr double piplus_mass = 0.13957039; //GeV
constexpr double pizero_mas = 0.1349768;   //GeV
constexpr double pimin_mass = 0.13957039;  //GeV
constexpr double n_mass = 0.93956542052;   //GeV
constexpr double p_mass = 0.93827208816;   //GeV

void calc_pionprod_thres_withpolangle(double E_beam = 4.0268, double theta_deg = 49.0)
{
	double E_gamma = E_beam;
	// Mandelstam s variable for the photon(with beam energy) + proton(at rest = target) pair.
	double s = std::pow(p_mass,2) + 2*E_gamma*p_mass;

	// Establish the CM frame.
	std::cout << '\n' << "*** Finding the Center of Momentum frame ***" << '\n';
	double v_CMframe = E_gamma/(E_gamma + p_mass); // Velocity of the CM frame w.r.t the lab frame (in the direction of the beam i.e. +Z).
	double beta = v_CMframe;
	double gamm_fac = 1/std::sqrt(1 - std::pow(beta,2));
	double E_tot_CM = std::sqrt(s);
	std::cout << "Velocity of the CM frame w.r.t lab frame (in the direction of the beam i.e. +Z): " << v_CMframe << "c \n";
	std::cout << "Total Energy in the CM frame (sqrt(mandelstam s)): " << E_tot_CM << " GeV\n";
	////


	// photon + p --> pi+ + n
	// Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
	std::cout << '\n' << "*** Reaction: photon + p --> pi+ + n ***" << '\n';
	double p_piplus_CM_1pi = std::sqrt( std::pow((s + std::pow(n_mass,2) - std::pow(piplus_mass,2)),2)/(4*s) - std:: pow(n_mass,2) );
	// Transform the pi+ momentum into lab frame.
	double E_piplus_CM_1pi = std::sqrt(std::pow(p_piplus_CM_1pi,2) + std::pow(piplus_mass,2));
	double p_piplus_lab_1pi = gamm_fac*(p_piplus_CM_1pi + v_CMframe*E_piplus_CM_1pi); // Maximum possible pi+ momentum in the lab frame.
	//double E_piplus_lab_1pi = std::sqrt(std::pow(p_piplus_lab_1pi,2) + std::pow(piplus_mass,2)); 
	double E_piplus_lab_1pi = gamm_fac*(E_piplus_CM_1pi + v_CMframe*p_piplus_CM_1pi);

	std::cout << "Max possible pi+ momentum in the CM frame: " << p_piplus_CM_1pi << " GeV/c" << '\n';
	std::cout << "Max possible pi+ energy in the CM frame: " << E_piplus_CM_1pi << " GeV" << '\n';
	std::cout << "Maximum possible pi+ momentum in the lab frame: " << p_piplus_lab_1pi << " GeV/c" << '\n';
	std::cout << "Maximum possible pi+ energy in the lab frame: " << E_piplus_lab_1pi << " GeV" << '\n';
	////
	

	// photon + p --> pi+ + pi0 + n
	// Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
	std::cout << '\n' << "*** Reaction: photon + p --> pi+ + pi0 + n ***" << '\n';
	double m_T = n_mass + pizero_mas; // Treat the pi0 and n as a composite system.
	double p_piplus_CM_2pi = std::sqrt( std::pow((s + std::pow(m_T,2) - std::pow(piplus_mass,2)),2)/(4*s) - std::pow(m_T,2) ); // Maximum pi+ momentum possible for the given beam energy in the CM frame.
	// Transform the pi+ momentum into lab frame.	
	double E_piplus_CM_2pi = std::sqrt(std::pow(p_piplus_CM_2pi,2) + std::pow(piplus_mass,2));
	double p_piplus_lab_2pi = gamm_fac*(p_piplus_CM_2pi + v_CMframe*E_piplus_CM_2pi); // Maximum possible pi+ momentum in the lab frame.
	//double E_piplus_lab_2pi = std::sqrt(std::pow(p_piplus_lab_2pi,2) + std::pow(piplus_mass,2));
	double E_piplus_lab_2pi = gamm_fac*(E_piplus_CM_2pi + v_CMframe*p_piplus_CM_2pi);
	
	std::cout << "Max possible pi+ momentum in the CM frame: " << p_piplus_CM_2pi << " GeV/c" << '\n';
	std::cout << "Max possible pi+ energy in the CM frame: " << E_piplus_CM_2pi << " GeV" << '\n';
	std::cout << "Maximum possible pi+ momentum in the lab frame: " << p_piplus_lab_2pi << " GeV/c" << '\n';
	std::cout << "Maximum possible pi+ energy in the lab frame: " << E_piplus_lab_2pi << " GeV" << '\n';
	////	
	
	std::cout << "###############################################################################################################" << '\n' << '\n';

	// Now we want to extend the calculation to the case where we have non-zero polar angle for theta/theta'
	std::cout << '\n' << "*** Extending the calcualtion to the case where the polar scattering angle of pi+ is non-zero ***" << '\n';
	double tan_theta = tan(TMath::DegToRad()*theta_deg); 

	// photon + p --> pi+ + n
	std::cout << '\n' << "*** Reaction: photon + p --> pi+ + n ***" << '\n';
	double cos_thetaprime_1pi_pos = ( -v_CMframe*E_piplus_CM_1pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_1pi,2) - std::pow(v_CMframe*E_piplus_CM_1pi,2)) + std::pow(p_piplus_CM_1pi,2))) / ( p_piplus_CM_1pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
	double cos_thetaprime_1pi_neg = ( -v_CMframe*E_piplus_CM_2pi*std::pow(gamm_fac*tan_theta,2) - std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_1pi,2) - std::pow(v_CMframe*E_piplus_CM_1pi,2)) + std::pow(p_piplus_CM_1pi,2))) / ( p_piplus_CM_1pi*(std::pow(gamm_fac*tan_theta,2) + 1) );

	std::cout << "cos_thetaprime_1pi_pos: " << cos_thetaprime_1pi_pos << '\n';
	std::cout << "Pos sol thetaprime: " << acos(cos_thetaprime_1pi_pos)*TMath::RadToDeg() << " degrees\n";
	std::cout << "cos_thetaprime_1pi_neg: " << cos_thetaprime_1pi_neg << '\n';
	std::cout << "Neg sol thetaprime: " << acos(cos_thetaprime_1pi_neg)*TMath::RadToDeg() << " degrees\n" << '\n';

	// Calculate the two values for pi+ momentum using the above two solutions (pos/neg) for thetaprime.
	double p_piplus_lab_1pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_1pi,2)*std::pow(cos_thetaprime_1pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_1pi*v_CMframe*E_piplus_CM_1pi*cos_thetaprime_1pi_pos + std::pow(p_piplus_CM_1pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_1pi,2));
	double p_piplus_lab_1pi_theta_neg = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_1pi,2)*std::pow(cos_thetaprime_1pi_neg,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_1pi*v_CMframe*E_piplus_CM_1pi*cos_thetaprime_1pi_neg + std::pow(p_piplus_CM_1pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_1pi,2));

	//double E_piplus_lab_1pi_theta_pos = std::sqrt( std::pow(p_piplus_lab_1pi_theta_pos,2) + std::pow(piplus_mass,2) );
	//double E_piplus_lab_1pi_theta_neg = std::sqrt( std::pow(p_piplus_lab_1pi_theta_neg,2) + std::pow(piplus_mass,2) );
	double E_piplus_lab_1pi_theta_pos = gamm_fac*( E_piplus_CM_1pi + v_CMframe*p_piplus_CM_1pi*cos_thetaprime_1pi_pos );
	double E_piplus_lab_1pi_theta_neg = gamm_fac*( E_piplus_CM_1pi + v_CMframe*p_piplus_CM_1pi*cos_thetaprime_1pi_neg );

	std::cout << "Max possible pi+ momentum in the lab frame for a scatterin angle of " << theta_deg << " degrees using the pos sol: " << p_piplus_lab_1pi_theta_pos <<" GeV/c\n";
	std::cout << "Max possible pi+ energy in the lab frame for a scatterin angle of " << theta_deg << " degrees using the pos sol: " << E_piplus_lab_1pi_theta_pos <<" GeV\n";
	std::cout << "Max possible pi+ momentum in the lab frame for a scatterin angle of " << theta_deg << " degrees using the neg sol: " << p_piplus_lab_1pi_theta_neg <<" GeV/c\n";
	std::cout << "Max possible pi+ energy in the lab frame for a scatterin angle of " << theta_deg << " degrees using the neg sol: " << E_piplus_lab_1pi_theta_neg <<" GeV\n";

	// Lab frame calculation //
	// double p_piplus_labframecalc_1pi_max = ( -(s+std::pow(piplus_mass,2)-std::pow(n_mass,2))*E_gamma*cos(theta_deg) + (E_gamma + p_mass)*std::sqrt( 4*std::pow(piplus_mass*E_gamma*cos(theta_deg),2) + (s+std::pow(piplus_mass,2)-std::pow(n_mass,2)) -4*std::pow(piplus_mass,2)*std::pow(E_gamma+p_mass,2) ) ) / ( 4*(std::pow(E_gamma*cos(theta_deg),2) - std::pow(E_gamma+p_mass,2)) );
	// std::cout << "Max possible pi+ momentum in the lab frame from pure lab frame only calcualtion: " << p_piplus_labframecalc_1pi_max << " GeV/c\n";
	////

	// photon + p --> pi+ + pi0 + n
	std::cout << '\n' << "*** Reaction: photon + p --> pi+ + pi0 + n ***" << '\n';
	double cos_thetaprime_2pi_pos = ( -v_CMframe*E_piplus_CM_2pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_2pi,2) - std::pow(v_CMframe*E_piplus_CM_2pi,2)) + std::pow(p_piplus_CM_2pi,2))) / ( p_piplus_CM_2pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
	double cos_thetaprime_2pi_neg = ( -v_CMframe*E_piplus_CM_2pi*std::pow(gamm_fac*tan_theta,2) - std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_2pi,2) - std::pow(v_CMframe*E_piplus_CM_2pi,2)) + std::pow(p_piplus_CM_2pi,2))) / ( p_piplus_CM_2pi*(std::pow(gamm_fac*tan_theta,2) + 1) );

	std::cout << "cos_thetaprime_2pi_pos: " << cos_thetaprime_2pi_pos << '\n';
	std::cout << "Pos sol thetaprime: " << acos(cos_thetaprime_2pi_pos)*TMath::RadToDeg() << " degrees\n";
	std::cout << "cos_thetaprime_2pi_neg: " << cos_thetaprime_2pi_neg << '\n';
	std::cout << "Neg sol thetaprime: " << acos(cos_thetaprime_2pi_neg)*TMath::RadToDeg() << " degrees\n" << '\n';

	// Calculate the two values for pi+ momentum using the above two solutions (pos/neg) for thetaprime.
	double p_piplus_lab_2pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_2pi,2)*std::pow(cos_thetaprime_2pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_2pi*v_CMframe*E_piplus_CM_2pi*cos_thetaprime_2pi_pos + std::pow(p_piplus_CM_2pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_2pi,2));
	double p_piplus_lab_2pi_theta_neg = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_2pi,2)*std::pow(cos_thetaprime_2pi_neg,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_2pi*v_CMframe*E_piplus_CM_2pi*cos_thetaprime_2pi_neg + std::pow(p_piplus_CM_2pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_2pi,2));

	//double E_piplus_lab_2pi_theta_pos = std::sqrt( std::pow(p_piplus_lab_2pi_theta_pos,2) + std::pow(piplus_mass,2) );
	//double E_piplus_lab_2pi_theta_neg = std::sqrt( std::pow(p_piplus_lab_2pi_theta_neg,2) + std::pow(piplus_mass,2) );
	double E_piplus_lab_2pi_theta_pos = gamm_fac*( E_piplus_CM_2pi + v_CMframe*p_piplus_CM_2pi*cos_thetaprime_2pi_pos );
	double E_piplus_lab_2pi_theta_neg = gamm_fac*( E_piplus_CM_2pi + v_CMframe*p_piplus_CM_2pi*cos_thetaprime_2pi_neg );

	std::cout << "Max possible pi+ momentum in the lab frame for a scatterin angle of " << theta_deg << " degrees using the pos sol: " << p_piplus_lab_2pi_theta_pos <<" GeV/c\n";
	std::cout << "Max possible pi+ energy in the lab frame for a scatterin angle of " << theta_deg << " degrees using the pos sol: " << E_piplus_lab_2pi_theta_pos <<" GeV\n";
	std::cout << "Max possible pi+ momentum in the lab frame for a scatterin angle of " << theta_deg << " degrees using the neg sol: " << p_piplus_lab_2pi_theta_neg <<" GeV/c\n";
	std::cout << "Max possible pi+ energy in the lab frame for a scatterin angle of " << theta_deg << " degrees using the neg sol: " << E_piplus_lab_2pi_theta_neg <<" GeV\n";


	// Printing out the final results 
	std::cout << '\n' << "*** Final Threshold Values ***" << '\n';
	double E_piplus_lab_1pi_theta_max = E_piplus_lab_1pi_theta_pos;
	double E_piplus_lab_2pi_theta_max = E_piplus_lab_2pi_theta_pos;
	std::cout << "Max pi+ energy (gamma, pi): " << std::fixed << std::setprecision(3)  << E_piplus_lab_1pi_theta_max << " GeV \n";
	std::cout << "Max pi+ energy (gamma, 2pi): " << std::fixed << std::setprecision(3) << E_piplus_lab_2pi_theta_max << " GeV \n";
	double piplus_energy_limit = E_piplus_lab_2pi_theta_max*(1 + 1.5/100.0);
	std::cout << "pi+ energy limit (gamma, pi): " << std::fixed << std::setprecision(3) << piplus_energy_limit << " GeV \n";

	// Calculating the minimum photon energy required
	double piplus_momentum_limit = std::sqrt( std::pow(piplus_energy_limit,2) - std::pow(piplus_mass,2) );
	double piplus_energy_limit_CM = gamm_fac*( piplus_energy_limit - v_CMframe*piplus_momentum_limit*cos(theta_deg*TMath::DegToRad()) );
	double piplus_momentum_limit_CM = std::sqrt( std::pow(piplus_energy_limit_CM,2) - std::pow(piplus_mass,2) );
	double E_gamma_min = ( std::pow( std::sqrt( std::pow(n_mass,2) + std::pow(piplus_momentum_limit_CM,2) ) + std::sqrt( std::pow(piplus_mass,2) + std::pow(piplus_momentum_limit_CM,2) ) ,2) - std::pow(p_mass,2) ) / (2*p_mass);
	std::cout << "Min photon energy giving pi+ from (gamma, pi) above pi+ energy limit: " << std::fixed << std::setprecision(3) << E_gamma_min << " GeV\n";
}