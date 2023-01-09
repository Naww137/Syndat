/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: mappings.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// header for the mappings class

#ifndef DEF_MAPPINGS
#define DEF_MAPPINGS

#include "dd_vector.hpp"  // for Ddvec::dd_entry

namespace Maps
{
// ----------- class particleInfo -----------------
//! mass information needed for discrete 2-body reactions
class particleInfo
{
public:
  double mProj;     //   Projectile's mass 
  double mTarg;     //   Target's mass 
  double mProd;     //   Product's mass 
  double mRes;      //   Residual's mass

  particleInfo( ): mProj(-1.0), mTarg(-1.0), mProd(-1.0), mRes(-1.0) {}
  ~particleInfo( ) {}

  //! Check that we have masses
  void check_data( ) const;

  //! Make a copy
  //! \param to_copy partidcle masses to copy
  void copy( const Maps::particleInfo &to_copy );

  //! Returns true if any particle is a photon
  bool photon_in_out( );
};

// ----------- class map_cm_lab -----------------
//! Class that maps energy-angle pairs between Newtonian reference frames
class map_cm_lab
{
protected:

public:
  // cm energy of ejectum
  // Eout_cm = alpha*E_in + beta_eject*Q
  double alpha;
  double beta_targ;
  double beta_proj;
  double beta_eject;
  double Q_value;     // net energy of the reaction

  // translational energy from initial collision
  // E_transl = gamma*E_in
  double gamma;

  //! Default constructor
  inline map_cm_lab( ): Q_value(0.0) {}

  //! Default destructor
  inline ~map_cm_lab() {}

  //! Sets alpha, beta, and gamma
  //! \param mProj mass of the projectile
  //! \param mTarg mass of the target
  //! \param mEject mass of the ejected particle
  //! \param mRes mass of the residual
  void setup_params( double mProj, double mTarg, double mEject, double mRes );

  //! Sets beta and gamma mass ratios for compound reactions
  //! \param mProj mass of the projectile
  //! \param mTarg mass of the target
  //! \param mEject mass of the ejected particle
  void setup_ratios( double mProj, double mTarg, double mEject );

  //! Gets the translational energy from initial collision
  //! \param E_in energy of incident particle
  inline double get_Etrans( double E_in ) const
  { return gamma*E_in; }

  //! Gets the laboratory energy of outgoing particle
  //! \param E_in energy of incident particle
  //! \param E_cm center-of-mass energy of ejected particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  double get_Elab( double E_in, double E_cm, double mu_cm );

  //! Gets the laboratory energy and cosine of outgoing particle
  //! \param E_in energy of incident particle
  //! \param E_cm center-of-mass energy of ejected particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param Eout_lab computed lab-frame energy of ejected particle
  //! \param mu_lab computed lab-frame direction cosine of ejected particle
  void get_E_mu_lab( double E_in, double E_cm, double mu_cm, double *Eout_lab,
    double *mu_lab );

  //! center-of-mass kinetic energy of incident particle plus target
  //! \param E_in energy of incident particle
  inline double incident_cm_KE( double E_in )
  { return beta_proj*E_in;}

  //! Gets the center-of-mass cosine for given incident energy, E_out_lab, and E_cm 
  //! \param E_in energy of incident particle
  //! \param E_out_lab computed lab-frame energy of ejected particle
  //! \param E_cm center-of-mass energy of ejected particle
  double get_mu_cm( double E_in, double E_out_lab, double E_cm );

  //! The total emission energy
  //! \param E_cm center-of-mass energy of ejected particle
  inline double total_Eout( double E_cm )
  { return E_cm/beta_eject; }
};

// ----------- class two_body_map -----------------
//! Class for cm-lab mappings for discrete 2-body reactions, Newtonian
class two_body_map : public Maps::map_cm_lab
{
private:

public:
  double threshold;
  double minimal_Ein;  // incident enrgy for minimum outgoing lab-frame energy for this cosine
  double minimal_Eout;  // minimum outgoing lab-frame energy for this cosine

  //! center-of-mass energy of outgoing particle when the residual is a photon
  double Ecm_photon_residual;

  //! Default constructor
  two_body_map( ) {}

  //! Default destructor
  ~two_body_map() {}

  //! Sets up the map from center-of-mass to laboratory coordinates
  //! \param particle_info the particle masses
  //! \param Q the energy of the reaction
  void set_map( Maps::particleInfo &particle_info, double Q );

  //! Calculates the threshold energy
  void get_threshold( );

  //! Gets the center-of-mass energy of outgoing particle
  //! Eout_cm = alpha*E_in + beta_eject*Q
  //! \param E_in energy of incident particle
  double get_Ecm( double E_in );

  //! Gets the value of mu_cm given the incident energy and the laboratory energy of outgoing particle.
  //! \param E_in energy of incident particle
  //! \param E_lab lab-frame energy of ejected particle
  double get_mu_cm( double E_in, double E_lab );

  //! Gets the laboratory energy and cosine of outgoing particle for discrete 2-body reactions
  //! \param E_in energy of incident particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param Eout_lab computed lab-frame energy of ejected particle
  //! \param mu_lab computed lab-frame direction cosine of ejected particle
  void two_body_get_E_mu_lab( double E_in, double mu_cm, double *Eout_lab, double *mu_lab );

  //! Returns the incident energy for zero outgoing lab energy when this exists
  double zero_Eout( );

  //! For endothermic reactions, reReturns the incident energy for minimal outgoing lab energy when this minimum is positive
  double min_Eout( );
};

// ----------- class Newton_map_param -----------------
//! Class for 2-body map as parameter for the Newtonian_F functions
class Newton_map_param
{
public:
  Maps::two_body_map *masses;  // the mass ratios
  double mu_cm;       // the centger-of-mass direction cosine

  Newton_map_param( ) {}
  ~Newton_map_param( ) {}

  //! For mu < 0 find the incident energy which minimizes Eout.
  Ddvec::dd_entry find_bottom( double mu );

  //! Find the incident energy which minimizes Eout for mu = -1
  Ddvec::dd_entry find_lowest_bottom(  );

  //! Finds the incident energy for given E_lab and mu_cm.
  //! Returns the number of solutions.
  //! \param E_lab lab-frame energy of ejected particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_lab, double mu_cm, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 );
};

// ----------- class two_step_map -----------------
//! Class for cm-lab mappings for discrete 2-step reactions, Newtonian
class two_step_map : public Maps::two_body_map
{
private:

public:
  //! data for the second step of the reaction
  Maps::two_body_map step_2_map;
  
   //! the energy of the product in the frame of the second step
  double E_cm2_prod;

  //! Default constructor
  two_step_map( ) {}

  //! Default destructor
  ~two_step_map() {}

  //! Sets up the map from center-of-mass to laboratory coordinates
  //! \param step1_particles the particle masses for step 1
  //! \param Q the energy of the first step
  void set_map( Maps::particleInfo &step1_particles, double Q,
		Maps::particleInfo &step2_particles );

  //! Gets the laboratory energy of the outgoing particle for 2-step reactions
  //! \param E_in, energy of incident particle
  //! \param mucm1, center-of-mass direction cosine of ejected particle
  //! \param mucm2, the direction cosine for the second step
  double two_step_get_E_lab( double E_in, double mucm1, double mucm2 );

  //! Returns the direction cosine of outgoing particle for 2-step reactions
  //! \param E_in energy of incident particle
  //! \param mucm1 center-of-mass direction cosine for first step
  //! \param mucm2, the direction cosine for the second step
  //! \param Etrans2, contribution from the ejected particle from step 1
  //! \param Eout_2, lab-frame energy of ejected particle from step 2
  //! \param w, the longitude for the second step, $-\pi/2 \le w \le \pi/2$.
  double two_step_get_mu_lab( double E_in, double mucm1, double mucm2,
			      double Etrans2, double Eout_2, double w );

  //! Returns the translational energy for step 2 from the outgoing energy for step 1
  //! \param Eout_lab1, the outgoing energy for step 1
  double two_body_get_Etrans2( double Eout_lab1 ) const
  { return step_2_map.gamma*Eout_lab1; }

  //! Returns the direction cosine for the first step for given incident energy and desired outgoing energy
  //! \param E_in, the incident energy
  //! \param mucm2, direction cosine for step 2 (+/- 1)
  //! \param Eout_lab, the desired outgoing energy, lab frame
  double get_mucm1( double E_in, int mucm2, double Eout_lab );

  //! Returns the direction cosine for the second step for given translational energy and desired outgoing energy
  //! \param Etrans2, the energy due to motion of the center of mass for step 2
  //! \param Eout_lab, the desired outgoing energy, lab frame
  double get_mucm2( double Etrans2, double Eout_lab );

};

// ----------- class two_step_map_param -----------------
//! Class for cm-lab mappings for discrete 2-step reactions as a function parameter, Newtonian
class two_step_map_param
{
private:

public:
  //! data for the mapping
  Maps::two_step_map *map;
  
  //! the incident energy
  double E_in;

  //! the direction cosines for the 2 steps
  double mucm_1;
  double mucm_2;

  //! Default constructor
  two_step_map_param( ) {}

  //! Default destructor
  ~two_step_map_param() {}

  //! Finds the incident energy for given E_lab.
  //! \param E_lab lab-frame energy of ejected particle
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_lab, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 );
  
  //! Find the incident energy which minimizes Eout for mucm_1 = mucm_2 = -1
  Ddvec::dd_entry find_lowest_bottom(  );
};

// ----------- class Ecm_intersect -----------------
//! Class for data used in the intersection of an E_cm = const. surface with an E_out cylinder
class Ecm_intersect : public Maps::map_cm_lab
{
private:

public:
  double E_in;     // incident energy
  double E_cm;     // center-of-mass energy of outgoing particle
  double G;        // value of hit_G
  double dG;       // value of d_hit_G( );
  double V_trans;  // velocity of translation of the center of mass
  double V_cm;     // center-of-mass velocity of outgoing particle

  //! Default constructor
  Ecm_intersect( ) {}

  //! Default destructor
  ~Ecm_intersect( ) {}

  //! Gets the translational energy from initial collision
  inline double get_Etrans( ) const
  { return Maps::map_cm_lab::get_Etrans( E_in ); }

  //! Sets the energies
  //! \param Ein energy of the incident particle
  //! \param Ecm_out energy of the outgoing particle in the center-of-mass frame
  void set_energies( double Ein, double Ecm_out );

  //! This Ecm curve hits the E_lab=const. curve if and only if hit_G(E_lab) >= 0
  //! \param E_lab the outgoing energy bin boundary that we may intersect
  void set_hit_G( double E_lab );

  //! Sets up the parameters for an intermediate incident energy
  //! \param Ein energy of the incident particle
  //! \param low_data intgersection information at a lower incident energy
  //! \param high_data intgersection information at a higher incident energy
  void interpolate_Ein( double Ein, double E_lab, const Maps::Ecm_intersect &low_data,
			const Maps::Ecm_intersect &high_data );
};

// ----------- class phase_space_map -----------------
//! Class for data used in the intersection of an E_cm = const. surface with an E_out cylinder
class phase_space_map : public Maps::Ecm_intersect
{
private:
  double out_mass_ratio;
  double exponent;  // the phase-space exoponent, (3*n/2) - 4
  double normalize;  // the phase-space normalization constant

public:
  int num_particles;

  //! Default constructor
  phase_space_map( ) {}

  //! Default destructor
  ~phase_space_map( ) {}

  //! sets the private members
  //! \param numParticles the number of emitted particles in the phase space model
  //! \param mEject mass of the ejected particle
  //! \param totalMass the total mass of emitted particles in the phase space model
  //! \param Q net energy of the reaction
  void set_data( int numParticles, double mEject, double totalMass, double Q );

  //! Returns the reaction threshold for the phase-space model
  double get_threshold( );

  //! Returns the maximum center-of-mass energy of the outgoing particle for the phase-space model
  //! \param E_in energy of incident particle
  double get_Ecm_max( double E_in );

  //! Returns the center-of-mass energy probability density
  //! \param E_cm center-of-mass energy of ejected particle
  //! \param Ecm_max maximal center-of-mass energy of ejected particle
  double get_prob( double E_cm, double Ecm_max );
};
} // end of namespace Maps

// ************* functions *************
namespace Newtonian_F
{
  // ------------------ T_out_lab -----------------------
  //! Returns the kinetic energy of the emitted particle in the lab frame
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the Newtonian boost
  double T_out_lab( double T_in_lab, void *params );

  // ------------------ two_step_T_out_lab -----------------------
  //! For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the Newtonian boost
  double two_step_T_out_lab( double T_in_lab, void *params );
  
}

#endif
