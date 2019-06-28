// Solid.cpp
//
// 
//This code computes the solid angle subtended by a 
// rapidly rotating star.
//
//
// Based on code written by Coire Cadeau.
//
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

// Full spectral parameters and 2 energy bands included.

#include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include "OblDeflectionTOA.h"
#include "Chi.h"
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "PolyOblModelBase.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"




int main(int argc, char** argv) try {

  int nmu(30); //set smaller nmu when testing
 
  double incl, theta(90.0), mass, rspot(0.0), omega, req, bbrat, ts;
  double distance, Temperature, E0(1.0);
  double north, south;
  double Solid_Angle;
  double omega_cgs;
  unsigned int modelchoice(1), quadmodel(0), fdmodel(0);
  int run;
  unsigned int numbins(59);
  bool incl_is_set(false), theta_is_set(false), bbrat_is_set(false), numbins_is_set(false),
    mass_is_set(false), rspot_is_set(false), omega_is_set(false),
    model_is_set(true), infile_is_set(false), ignore_time_delays(false);
  char out_file[256] = "oblflux.txt"; 

  int power=0;

  char data_file[80];
  char line[80];
  char vert[180];
  char min_file[180];
  //char min_store[180];
  char out_dir[80];
  // New variables added by SMM
 
  double aniso(0.586), Gamma(2.0);
  double shift_t(0.0);    /* Time off-set from data */
  double ftol = 1e-3;
  double b_mid;
  std::ifstream data;
  std::ofstream out;
  std::ofstream angles;
  std::ofstream dOmega;
  std::ofstream boost;
  std::ofstream cosbeta;
  
  //Define an array of energies (keV) for plotting a spectrum
  const int E_len = 100; //length of energy array
  double E_n = 0.05; //starting energy value (keV)
  const double E_bin = 0.05; //bin width
  double E_spec[E_len] = {0.0}; //initialize array of energies for spectrum file (keV)
 
  // Populate E_spec array
  int n = 0;
  for(n = 0; n < E_len; n++)
   {
   E_spec[n] = E_n;
   //std::cout << "Energies (keV):" << E_spec[n] << std::endl;
   E_n = E_n + E_bin;
   }

  
  // Create LightCurve data structure
  class LightCurve curve;
  double mu, cosg;
 
  curve.para_read_in = false;


  // Read in the command line options
  for(int i(1); i < argc; i++) {
    if(argv[i][0]=='-') {
      switch(argv[i][1]) {

      case 'o':
	// Name of output file
	sscanf(argv[i+1], "%s", out_file);
	break;

      case 's':
	// This option controls the use of an external file
	// to define the shape of the spot.
	infile_is_set = true;
	break;

      case 'q':
	// Oblateness model (default is 1)
	sscanf(argv[i+1], "%u", &modelchoice);
	model_is_set = true;
	break;

      case 'Q':
	// Quadrupole model (default is 0)
	sscanf(argv[i+1], "%u", &quadmodel);
	break;

      case 'w':
	// Frame dragging (default is 0)
	sscanf(argv[i+1],"%u",&fdmodel);
	break;

      case 'n':
	//number of bins
	sscanf(argv[i+1], "%u", &numbins);
	numbins_is_set = true;
	break;

      case 'p':
	//power that eta is raised to
	sscanf(argv[i+1], "%d", &power);
	break;

      case 'i':
	// Inclination angle of the observer (degrees)
	sscanf(argv[i+1], "%lf", &incl);
	incl_is_set = true;
	break;

      case 'e':
	// Emission angle of the spot (degrees)
	sscanf(argv[i+1], "%lf", &theta);
	theta_is_set = true;
	break;

      case 'm':
	// Mass of the star (solar mass units)
	sscanf(argv[i+1], "%lf", &mass);
	mass_is_set = true;
	break;

      case 'r':
	// Radius of the star at the equator(km)
	sscanf(argv[i+1], "%lf", &req);
	rspot_is_set = true;
	break;

      case 'f': 
	// Spin frequency (Hz)
	sscanf(argv[i+1], "%lf", &omega);
	omega_cgs = omega;
	omega_is_set = true;
	break;

      case 't': //toggle ignore_time_delays (only affects output)
	ignore_time_delays = true;
	break;

      case 'T':  // Temperature of the spot, in the star's frame, in keV
	sscanf(argv[i+1], "%lf", &Temperature);
	break;

      case 'D':  // Distance to NS
	sscanf(argv[i+1], "%lf", &distance);
	break;

      case 'h':
      
	std::cout << "solid help: " << std::endl
		  << "-o filename: (optional, default is oblflux.txt)" << std::endl
		  << "-s filename: (optional)" <<std::endl
		  << "-n numbins: number of datapoints" << std::endl
		  << "-i inclination of observer in degrees, between 0 and 180." << std::endl
		  << "-e location of emission region in degrees, between 0 and 180." << std::endl
		  << "-m mass of star in Msun." << std::endl
		  << "-r radius of star (at the spot) in km." << std::endl
		  << "-f rotation frequency of star in Hz." << std::endl
		  << "-t ignore time delays in output (see source)." << std::endl
		  << "-p azimuthal location (phi_0) of spot [0.0]." << std::endl
		  << "-q #, where #=:" << std::endl
		  << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		  << "      2 for CFL quark star poly model" << std::endl
		  << "      3 for spherical model" << std::endl
		  << "-Q [0] Quadrupole model 0 or 1" << std::endl
		  << std::endl;
	return 0;
      } // end switch	
    } // end if
  } // end for

  // Check that necessary numbers are set.
  if( !( incl_is_set 
	 && numbins_is_set
	 && mass_is_set
	 && rspot_is_set
	 && omega_is_set
	 ) ) {
    std::cout << "Not all required parameters were specified. Exiting."
	      << std::endl;
    return -1;
  }

  // Print out information about the model
  std::cout << std::endl << "Pulse: Stellar Model Parameters" << std::endl << std::endl;
  std::cout << "Mass = " << mass << " Msun" << std::endl;
  std::cout << "Req = " << req << " km" <<  std::endl;
  std::cout << "2GM/Rc^2 = " << 2.0*Units::G*mass*Units::MSUN/(req*1e+5*Units::C*Units::C)  <<  std::endl;
  std::cout << "spin = " << omega << " Hz" << std::endl;
  std::cout << "v/c = " << req * 1e5 * omega * 2.0*Units::PI/Units::C << std::endl;
  std::cout << "inclination = " << incl << " degrees" << std::endl;

  double M_over_R = Units::G*mass*Units::MSUN/(req*1e+5*Units::C*Units::C);
    
  std::cout << std::endl;
 
  // Units conversions.
  incl *= Units::PI / 180.0;
  theta *= Units::PI / 180.0;
  mu = cos(theta);
  mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
  req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
  omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
  distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );

  double baromega; // dimensionless omega
  baromega = omega * sqrt( pow(req,3)/mass);
  std::cout << "baromega = " << baromega <<std::endl;

  north = incl;
  south = incl; //was 180 - incl for previous version of code

  std::cout << "mass/radius = " << mass/req << std::endl;

  // curve is a structure holding parameters describing the star and emission properties
  curve.para.mass = mass;
  curve.para.omega = omega;
  //  curve.para.theta = theta;
  //curve.para.incl = incl;
  curve.para.aniso = aniso;
  curve.para.Gamma = Gamma;
  curve.para.bbrat = bbrat;
  curve.para.ts = ts;

 
  // Set up model describing the oblate shape of the star
  OblModelBase* model;
  if (modelchoice == 1 ) {
    // Default model 
    model = new PolyOblModelNHQS( rspot, req,
				  PolyOblModelBase::zetaparam(mass,req),
				  PolyOblModelBase::epsparam(omega, mass, req)
				  );
  }
  else if(modelchoice == 2 ) {
    // Alternative model for quark stars (not very different)
    model = new PolyOblModelCFLQS( rspot, req,
				   PolyOblModelBase::zetaparam(mass,req),
				   PolyOblModelBase::epsparam(omega, mass, req)
				   );
  }
  else if(modelchoice == 3 ) {
    // Use a spherical star
    model = new SphericalOblModel( rspot );
  }
  else {
    std::cerr << "Invalid modelchoice parameter. Exiting." << std::endl;
    return -1;
  }

  // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
  // used to compute deflection angles and times of arrivals 
  OblDeflectionTOA* defltoa = new OblDeflectionTOA( model, mass );

 

  curve.junk.shift_t = shift_t;
  curve.junk.infile_is_set = infile_is_set;
  curve.numbins = numbins;

  //  std::cout << "i" << "\t" << "a1" << "\t" << "phase" << "\t" << "Chi-Square" << std::endl;

  double dOmega_mu(0.0), dFlux_mu[E_len] = {0.0};
  double Omega_s(0.0), Flux[E_len] = {0.0};
  double Omega_sQ(0.0), Flux_Q(0.0);
  double dBoloFlux_mu(0.0), BoloFlux(0.0);
  double dArea_mu(0.0), Area(0.0);
  double BB[E_len] = {0.0};

  double red; // Gravitational redshift red=1+z= (1-2M/R)^{-1/2}
  double quad; // Quadrupole Moment
  double riso; // Isotropic radius
  double framedragging; 
  double inertia;
  double speed;

 

  if (quadmodel == 1){
    // quad = - 0.1 * pow(baromega,2) * pow( req/mass, 2);
    quad = -0.2;
  }
  else 
    quad = 0.0;
  std::cout << "mass/req = " << mass/req 
    << "Quadrupole = " << quad  << std::endl;

  double a_kerr;

  if (fdmodel==1)
    a_kerr=0.24;
  else
    a_kerr=0.0;


  if (fdmodel==1)
    inertia = sqrt(mass/req) * (1.14 - 2.53*mass/req + 5.6*pow(mass/req,2));
  else
    inertia = 0;

  double dmu,dphi,dtheta, P2;
  // nmu = 10; //set equal to spacing in python code
  // numbins = 20;
  dmu = 1.0/(nmu*1.0);
  dphi = 2 * Units::PI / (numbins*1.0);
  dtheta = (Units::PI - 1e-6)/(1.0*nmu);
  
  char angle_file[80];
  char dOmega_file[80];
  char boost_file[80];
  char cosbeta_file[80];
  
  std::cout<<"incl = " << incl*180/Units::PI<<"spin =" << omega_cgs<<"M/R = "<< 100.0*M_over_R<<std::endl;
    if (modelchoice == 3){
      if (omega_cgs < 1) {
	if(incl*180/Units::PI != 0){
	  std::cout<<"spherical, spin <1, nonzero incl"<<std::endl;
	  sprintf(angle_file, "angles_sph_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(dOmega_file, "dOmega_sph_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(boost_file, "boost_sph_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(cosbeta_file, "cosbeta_sph_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	}
	else {
	  std::cout<<"spherical, spin <1, zero incl"<<std::endl;
	  sprintf(angle_file, "angles_sph_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(dOmega_file, "dOmega_sph_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(boost_file, "boost_sph_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	  sprintf(cosbeta_file, "cosbeta_sph_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	}
      
      }
     else{
      if(incl*180/Units::PI != 0){
	std::cout<<"spherical, spin > 1 , nonzero incl"<<std::endl;
	sprintf(angle_file, "angles_sph_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_sph_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_sph_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_sph_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
      else{
	std::cout<<"spherical, spin > 1, zero incl"<<std::endl;
	sprintf(angle_file, "angles_sph_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_sph_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_sph_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_sph_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
     }
    }
  else{
    if (omega_cgs < 1) {
      if (incl*180/Units::PI != 0){
        std::cout<<"oblate, spin < 1, nonzero incl"<<std::endl;
	sprintf(angle_file, "angles_obl_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_obl_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_obl_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_obl_spin%1.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
      else{
	std::cout<<"oblate, spin < 1, zero incl"<<std::endl;
	sprintf(angle_file, "angles_obl_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_obl_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_obl_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_obl_spin%1.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
     }
    else{
      if (incl*180/Units::PI != 0){
        std::cout<<"oblate, spin > 1, nonzero incl"<<std::endl;
	sprintf(angle_file, "angles_obl_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_obl_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_obl_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_obl_spin%3.0f_MR%2.0f_incl%2.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
      else{
	std::cout<<"oblate, spin > 1, zero incl"<<std::endl;
	sprintf(angle_file, "angles_obl_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(dOmega_file, "dOmega_obl_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(boost_file, "boost_obl_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
	sprintf(cosbeta_file, "cosbeta_obl_spin%3.0f_MR%2.0f_incl%1.0f.txt", omega_cgs, 100.0*M_over_R, incl*180/Units::PI);
      }
    }
  }
 
   //open file for angles
      angles.open(angle_file);
      angles<<"theta "<< "phi "<<"dOmega "<<"dtheta "<<std::endl;
  //open file for dOmega
      dOmega.open(dOmega_file);
      dOmega<<"theta = rows; phi = columns"<<std::endl;
  //open file for doppler boost factor
      boost.open(boost_file);
      boost<<"theta = rows; phi = columns"<<std::endl;
  //open file for cosbeta
      cosbeta.open(cosbeta_file);
      cosbeta<<"theta = rows; phi = columns"<<std::endl;

    std::cout<<" dtheta = "<<dtheta<<std::endl;
  // Loop through the hemispheres
    for (int k(0);k<=1; k++){
      
       // Loop through the star's latitudes
      for (int j(0); j<nmu/2; j++){
      // for (int j(3); j<4; j++){
	 
	 if(k==0){
	   //if in northern hemisphere
	   //mu = (nmu - 1 - j)*dmu;
	   theta = j*dtheta + 1e-06 + 0.5*dtheta;
	   
	 }
	 else{
	   //if in southern hemisphere
	   //mu = j*dmu;
	   theta = Units::PI/2 + ((j + 1)*dtheta + 1e-06) - 0.5*dtheta; 
	 }
	 // theta = acos(mu);
	 // theta = j*dtheta;
	 mu = cos(theta);
     
	 curve.para.theta = theta;
	 //     P2 = 0.5 * ( 3.0*mu*mu - 1.0);

	 // Values we need in some of the formulas.

	 if (modelchoice==1){
	     rspot = calcrspot( omega, mass, theta, req);
	   // rspot = R_at_costheta(mu);
	   cosg = cosgamma(mu,req,rspot,mass,omega);
	 }
	 if (modelchoice==3){
	   rspot = req;
	   cosg = 1.0;
	 }

	 //     curve.para.dS = pow(rspot/distance,2)/cosg * dmu * dphi;  

	 //curve.para.dS = dmu * dphi/cosg;
	 //removed dtheta dphi May 27th
	 curve.para.dS = (sin(theta)*(dtheta*dphi))/cosg;
	 curve.para.cosgamma = cosg;
	 curve.para.radius = rspot;

	 riso = rspot * pow( 1 + 0.5*mass/rspot, -2);
	 //     red = 1.0/sqrt(1 - 2*mass/rspot) * exp( quad * P2 * pow( mass/riso,3)) ;

	 double x;
	 x = rspot/mass;

	 double F1;
	 F1 = -5.0/8.0*(x-1.0)/(x*(x-2)) * (2.0+6.0*x - 3.0*x*x)
	   - 15.0/16.0 * x*(x-2) * log(x/(x-2.0));

	 double Sigma;
	 Sigma = pow(x,2) + pow(a_kerr*mu,2);

	 double Nsquare;
	 Nsquare = (1.0 - 2.0*x/Sigma) * ( 1.0 + 2.0*quad*F1*P2);

	 red = 1.0/sqrt(Nsquare);

	 // std::cout
	 //  << "j = " << j
	 //  << " mu = " << mu 
	 //  << " dS = " << curve.para.dS
	 //  <<  std::endl;

	 framedragging = 2.0*inertia * pow(req/mass,2) * pow(mass/riso,3) * (1.0 - 3.0*mass/riso);
	 speed = omega*rspot*sin(theta)* red * (1.0 - framedragging);

	 //     std::cout << "mu = " << mu << " P2(mu) = " << P2 << " red^{-4} = " << pow(red*sqrt(1 - 2/x),-4) << std::endl;
	 //std::cout << "mass/radius = " << mass/req << " Quad = " << quad << std::endl;
	 //std::cout << "inertia = " << inertia << " fd = " << framedragging << " speed = " << speed << std::endl;


    
       
	 curve.defl.b_max =  defltoa->bmax_outgoing(rspot);
	 curve.defl.psi_max = defltoa->psi_max_outgoing(curve.defl.b_max,rspot,&curve.problem);

	 // Compute b vs psi lookup table good for the specified M/R and mu
	 b_mid = curve.defl.b_max*0.9;
	 curve.defl.b_psi[0] = 0.0;
	 curve.defl.psi_b[0] = 0.0;
	 for (unsigned int i(1); i < NN+1; i++){ // compute table of b vs psi points 
	   curve.defl.b_psi[i] = b_mid * i/(NN*1.0);
	   curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i],rspot,curve.defl.b_max, curve.defl.psi_max, &curve.problem);
	 }
	 // For arcane reasons, the table is not evenly spaced.
	 for (unsigned int i(NN+1); i < 3*NN; i++){ // compute table of b vs psi points 
	   curve.defl.b_psi[i] = b_mid + (curve.defl.b_max-b_mid)/2.0 * (i-NN)/(NN*1.0);
	   curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i],rspot,curve.defl.b_max, curve.defl.psi_max, &curve.problem);
	 }
	 curve.defl.b_psi[3*NN] = curve.defl.b_max;
	 curve.defl.psi_b[3*NN] = curve.defl.psi_max;
	 // Finished computing lookup table
     
	 //std::cout << "Finished computing lookup table" << std::endl;

      
 
	 if (k==1) curve.para.incl = south;
	 else curve.para.incl = north;

	 curve =  ComputeAngles( &curve, defltoa);

	 dOmega_mu = 0.0;
	 dBoloFlux_mu = 0.0;
	 dArea_mu = 0.0;

	 double eta;
	 double redshift;

	 redshift = 1.0/sqrt(1 - 2/x);
       
	 //Loop through azimuthal angles
	 // std::cout<<"dphi = "<<dphi<<std::endl;
	 for (int i(0);i<numbins;i++){
	   // std::cout<<"i = "<<i<<" phi = "<<curve.t[i]*2*Units::PI<<std::endl;
	   /* std::cout<<"theta ="<<theta<<" mu ="<<mu<<" phi ="<<curve.t[i]*2*Units::PI - Units::PI
		    << " cospsi = " << sin(theta)*cos(curve.t[i]*2*Units::PI)
		    <<" cosbeta ="<<curve.cosbeta[i]
		    << " dOmega = " << curve.dOmega_s[i]
		    <<std::endl; */
	   angles<<theta<< "  "<<curve.t[i]*2*Units::PI<<" "<<curve.dOmega_s[i]<<" "<<dtheta<<std::endl;
	   dOmega<<" "<<curve.dOmega_s[i];
	   boost<<" "<<curve.eta[i];
	   cosbeta<<" "<<curve.cosbeta[i];

	   /*  if (j==1)
	     std::cout<<" t = "<<curve.t[i] 
	     <<" dOmega[i] = "<<curve.dOmega_s[i]<<std::endl; */

	   
	   dArea_mu += curve.para.dS * pow( Units::nounits_to_cgs( rspot*1.0e-5, Units::LENGTH ),2) ;

	   if (curve.dOmega_s[i] != 0.0){

	     eta = sqrt( 1 - speed*speed)/(1.0 - speed * curve.cosxi[i]);
	   
	     dOmega_mu += curve.dOmega_s[i];

	     dBoloFlux_mu += curve.dOmega_s[i]*pow(eta,4)/Units::PI;
	   
	     Solid_Angle += curve.dOmega_s[i];
	   
	     for (int n(0); n<E_len; n++)
	       {
		 dFlux_mu[n] += curve.dOmega_s[i] * pow(eta,3) * BlackBody(Temperature,E_spec[n]*redshift/curve.eta[i]);
	       }
	   
	   }
       
	 } //end of azimuthal angle loop
	 std::cout
	   <<"j = "<<j
	   <<"theta = "<< theta
	   <<" costheta = "<< mu
	   <<" R = "<<Units::nounits_to_cgs( rspot*1.0e-5, Units::LENGTH )
	   <<" dS = "<< curve.para.dS
	   <<" dArea_mu = "<< dArea_mu
	   <<" dOmega = "<< dOmega_mu*pow(rspot/distance,2)
	   <<std::endl;
	 dOmega<<std::endl;
	 boost<<std::endl;
	 cosbeta<<std::endl;

	 Area += dArea_mu  ;

	 Omega_s += dOmega_mu * pow(rspot/distance,2);
	 for (int n(0); n<E_len; n++)
	   {
	     BB[n] = BlackBody(Temperature,E_spec[n]*redshift);
	     //Flux[n] = dFlux_mu[n] * pow( sqrt(1 - 2/x), 3) * (1.0 / ( E_spec[n] * Units::H_PLANCK )) * pow(rspot/distance,2)*E_bin;
	     Flux[n] = dFlux_mu[n]* pow(rspot/distance,2) * pow( sqrt(1 - 2/x), 3);
	     //std::cout << "BB = " << BB[n] << std::endl;
	   }
	 BoloFlux += dBoloFlux_mu * pow(red,-4) * pow(rspot/distance,2);

       } //end of latitude loop
  } //end of hemispheres loop


  std::cout << "Solid Angle = " << Omega_s << std::endl;

  // std::cout<<"Flux = " << BoloFlux<< std::endl;
  
  std::cout << "distance = " << distance <<std::endl;
  
  std::cout<< "rspot = " <<rspot <<std::endl;
  
  std::cout << "Area = " << Area << std::endl;

  std::cout << "pi R^2/D^2 (1+z)^2 = " <<   Units::PI * pow(req,2)/(1.0-2.0*mass/req) << std::endl;

  std::cout << "R^2/D^2 (1+z)^{-2} = " << pow(req,2)*(1.0-2.0*mass/req) << std::endl;
  
  
 // Print out information about the model
  /*
  std::cout << std::endl << "Pulse: Stellar Model Parameters" << std::endl << std::endl;
  std::cout << "Mass = " << Units::nounits_to_cgs( mass/Units::MSUN, Units::MASS ) << " Msun" << std::endl;
  std::cout << "Req = " << Units::nounits_to_cgs( req*1.0e-5, Units::LENGTH ) << " km" <<  std::endl;
  std::cout << "2GM/Rc^2 = " << 2.0*Units::G*(Units::nounits_to_cgs( mass, Units::MASS ))/((Units::nounits_to_cgs( req, Units::LENGTH ))*Units::C*Units::C)  <<  std::endl;
  std::cout << "spin = " << omega_cgs << " Hz" << std::endl;
  std::cout << "v/c = " << rspot * 1e6 * omega_cgs * 2.0*Units::PI/Units::C << std::endl;
  std::cout << "inclination = " << incl << " degrees" << std::endl; 
  
  std::cout << std::endl; */

  out.open(out_file, std::ios::app);
  
  //write energy/fluxes to a file
  std::ofstream myfile;
  myfile.open ("spectra.txt");
  myfile << "Energy (keV)"<<"   Flux (#photons/cm^2)" <<"   BB (erg/cm^2)"<<std::endl;
  for (int n(0); n<E_len; n++)
    {
      //myfile<<E_spec[n]<<"       "<<Flux[n]<< "       "<<E_spec[n]*Flux[n]<<"        "<<BB[n]<<std::endl;
      myfile<<E_spec[n]<<"       "<<Flux[n]<< "       "<<BB[n]<<std::endl;
    }
    myfile.close(); 

    //angles.close();
  //  out = fopen(out_file, "a");
	
 // Print out information about the model
  out << "#M = " << Units::nounits_to_cgs( mass/Units::MSUN, Units::MASS ) << " Msun" 
      << " #Req = " << Units::nounits_to_cgs( req*1.0e-5, Units::LENGTH ) << " km" << std::endl;
    
	  
  out << "#spin      "
      << "Flux(1keV) "  
      << "BolFlux/Fs "
      << "v_{eq}/c   "
      << "Area       "
      << "A/4piR^2   "
      << "SolidAng   "
      << "SolidAng/Ss" << std::endl;

  out  << omega_cgs
       << "          " <<  Flux 
       << "          " <<  BoloFlux/(pow(req/distance,2)*(1.0-2.0*mass/req)) 
       << "          " <<     req * 1e6 * omega_cgs * 2.0*Units::PI/Units::C * red
       << "          " << Area
       << "          " << Area/(4.0 * Units::PI * pow( Units::nounits_to_cgs( rspot*1.0e-5, Units::LENGTH ),2))
       << "          " << Omega_s
       << "          " << Omega_s/(Units::PI * pow(req/distance,2)/(1.0-2.0*mass/req))
      << std::endl;



  delete defltoa;
  delete model;
  return 0;

 } catch(std::exception& e) {
  std::cerr << "Top-level exception caught in application oblflux: " << std::endl
	    << e.what() << std::endl;
  return -1;
 }
