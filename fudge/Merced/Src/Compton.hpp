/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Compton.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for Compton scattering

#ifndef COMPTON_CLASS
#define COMPTON_CLASS

#include <iostream>
#include "data_parser.hpp"
#include "transfer.hpp"
#include "param_base.hpp"
#include "x_vector.hpp"

namespace Comp
{
  
class Compton_Ein_param;  // foreward declaration

//! Class for intersections of a data curve with an integration box
// ---------------- class Eout_curve ------------------
class Eout_curve
{
private:

public:
  //! Energy of outgoing gamma
  double E_out;
  //! Energy of incident gamma giving this E_out in backscatter
  double Ein_back;
  //! Do we get backscatter?
  bool backscatter;

  Eout_curve( ) {}
  ~Eout_curve( ) {}

  //! Calculates the incident energy for backscatter
  //! \param Eout energy of outgoing gamma
  void setup( double Eout );

  //! Finds (Ein, mu) where this Eout curve hits the curve x = constant.
  //! \param x ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2} 
  Ddvec::dd_entry hit_x( double x ) const;

  //! Finds mu where this Eout curve hits the line Ein = constant.
  //! \param E_in energy of incident gamma
  double hit_Ein( double E_in ) const;

};

//! Identify the top and bottom of the region of integration
// ---------------- Boundary_ID -----------------------------
enum Boundary_ID{ USE_MU1,     // the top is mu = 1
                  USE_EOUT,    // the boundary is an E_out curve
                  USE_X,       // the boundary is an x curve
                  USE_MU_MINUS1 }; // the bottom is mu = -1

//! Class for a corner of the top or bottom of a quadrature box
// ---------------- class Boundary_corner ------------------
class Boundary_corner
{
public:
  double E_in;
  Boundary_ID boundary_ID;

  Boundary_corner( );
  Boundary_corner( const Boundary_corner& boundary_corner );
  ~Boundary_corner( ){}
};

//! Class for all corners of the top or bottom of a quadrature box
// ---------------- class corner_list ------------------
class corner_list : public std::list< Boundary_corner >
{
private:

public:
  corner_list( )  {}
  ~corner_list( ){}

  //! Append a corner to the list
  //! \param x ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2} 
  //! \param Eoutcurve intersections with the quadrature box
  //! \param boundary_ID identify this boundary of the quadrature box
  void new_corner( double x, const Eout_curve &Eoutcurve, Boundary_ID boundary_ID );

  //! Append a corner to the list
  //! \param E energy of incident gamma
  //! \param boundary_ID identify this boundary of the quadrature box
  void new_corner( double E, Boundary_ID boundary_ID );

  //! Copy the last entry from copy_from but change the boundary ID
  //! Used to copy the first and last bottom corners to the top corners
  //! \param copy_from current list of bottom corners
  //! \param boundary_ID identify this boundary of the quadrature box
  void copy_last( const corner_list &copy_from, Boundary_ID boundary_ID );
};

//! Class for the geometry of the quadrature region
// ---------------- class Compton_region ------------------
class Compton_region
{
private:

public:
  //! incident energy range
  double Ein_left;
  double Ein_right;

  //! the range of x values
  double x0;
  double x1;

  //! Curve for the bottom of the outgoing energy bin
  Eout_curve lower_Eout;
  //! Curve for the top of the outgoing energy bin
  Eout_curve upper_Eout;

  //! The bottom of the quadrature region
  corner_list bottom_corners;
  //! The top of the quadrature region
  corner_list top_corners;
  corner_list::iterator left_corner;  // intersection of the lower E_out with the lower x
  corner_list::iterator right_corner;  // intersection of the upper E_out with the upper x

  //! The bin counts
  int Ein_count;
  int Eout_count;

  Compton_region( )  {}
  ~Compton_region( ) {}

  //! constructs upper_Eout_hits and lower_Eout_hits
  //! \param Eout_ptr the lower bound of the outgoing energy bin
  void setup_Eout( std::vector< double >::const_iterator Eout_ptr );

  //! Identifies the top and bottom of the quadrature region
  //! \param x_0 lower value of ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2} 
  //! \param x_1 higher value of ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2} 
  void setup_region( double x_0, double x_1 );
};

//! Class for parameters for the 1-d quadrature of the Compton scattering factor over mu
//! class for (x_0, SF_0) and (x_1, SF_1) scattering factors
// ---------------- class Compton_mu_param ------------------
class Compton_mu_param : public Ddvec::dd_pair, public Pbase::param_base
{
public:
  double E_in;

  Compton_mu_param( ) {}
  inline ~Compton_mu_param( ) {}
};

//! Class for Compton scattering data
//--------------- class Compton ----------------
class Compton : public Xvec::x_vector
{
private:
  //! Adds the result of one integration
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param data for quadrature over incident energy   
  void update_T( Trf::T_matrix &transfer, Compton_Ein_param *Ein_param );

  //! Calculates the cross section
  //! \param mu_quad_rule the method for integration over direction cosine
  //! \param sigma the computed cross section data
  void get_xsec( Qmeth::Quadrature_Rule mu_quad_rule, dd_vector& sigma );

public:
  //! scattering factor in the file
  dd_vector file_data;

  Compton( ) {}

  ~Compton( ) {}

  // Calculates the cross section and transfer matrix
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param xsec the cross section data
  void get_T( Trf::T_matrix& transfer, dd_vector& xsec );

  // Climbs up the incident energy bins
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param data for quadrature over incident energy   
  void Ein_ladder( Trf::T_matrix& transfer, Compton_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param data for quadrature over incident energy   
  void set_Ein_range( Compton_Ein_param *Ein_param );

  //! Integrates over one x-E box; loop over the (E_in, mu) Compton region
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param data for quadrature over incident energy   
  void one_Ebox( Trf::T_matrix& transfer, Compton_Ein_param *Ein_param );

};

//! Class for parameters for the 2-d quadrature of the Compton scattering factor
// ---------------- class Compton_Ein_param ------------------
class Compton_Ein_param : public Pbase::param_base
{
public:
  Compton_region quad_box;

  // number of 2-d quadratures
  long int quad_count;
  // number of calls to Compton_F::Ein_F
  long int Ein_F_count;
  long int mu_F_count;  // number of calls to Compton_F::mu_F

  corner_list::const_iterator bottom_ptr; // where we are along the bottom
  corner_list::const_iterator next_bottom_ptr; // where we are along the bottom
  corner_list::const_iterator top_ptr; // where we are along the top
  corner_list::const_iterator next_top_ptr; // where we are along the top

  Xvec::x_vector::const_iterator x_ptr;  // where we are in the x-scattering function data
  Xvec::x_vector::const_iterator next_x;

  Eout_curve *Eout_ptr;  // where we are in the outgoing energies
  Eout_curve *next_Eout;

  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for cosine integral

  Compton_Ein_param( ): quad_count( 0 ), Ein_F_count( 0 ), mu_F_count( 0 ) {}
  ~Compton_Ein_param( ) {}
};

} // end of namespace Comp

// ----------------- functions to integrate --------------------
namespace Compton_F
{
  //! Function for computing the cross section
  //! Returns false if there are problems
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool sigma_F( double mu, Qparam::QuadParamBase *FF_param,
		Coef::coef_vector *value );

  //! Function for computing the cross section with singularity sqrt( flip_mu )
  //! Returns false if there are problems
  //! \param flip_mu = 1 - mu
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool flip_sigma_F( double flip_mu, Qparam::QuadParamBase *FF_param,
		     Coef::coef_vector *value );

  //! Function for the 1-d quadrature over mu
  //! Returns false if there are problems
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param SF_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu, Qparam::QuadParamBase *SF_param,
	     Coef::coef_vector *value );

  //! Function for the 1-d quadrature over mu with singularity sqrt( flip_mu )
  //! Returns false if there are problems
  //! \param flip_mu = 1 - mu
  //! \param SF_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool flip_mu_F( double flip_mu, Qparam::QuadParamBase *SF_param,
		  Coef::coef_vector *value );

  //! \param Ecm_quad_param the function parameters
  //! Returns false if there are problems
  //! Function for the 2-d quadrature over E_in
  //! \param E_in the energy of the incident particle
  //! \param compton_mu_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *compton_param,
	      Coef::coef_vector *value );
}

#endif
