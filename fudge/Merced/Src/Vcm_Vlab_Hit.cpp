/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: Vcm_Vlab_Hit.cpp 1 2010-10-05 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle the geometry of the cm-to-lab map

#include <cmath>

#include "Vcm_Vlab_Hit.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class Vhit::Vcm_quadBox_Hit *****************
// ---------------- Vhit::Vcm_quadBox_Hit::print ------------------
// Prints a hit
void Vhit::Vcm_quadBox_Hit::print( ) const
{
  std::cout << " V_cm: " << V_cm;
  switch( hit_corner )
  {
  case Vhit::V_INSIDE:
    std::cout << " Vhit::V_INSIDE";
    break;
  case Vhit::BOTTOM_FORWARD:
    std::cout << " Vhit::BOTTOM_FORWARD";
    break;
  case Vhit::BOTTOM_BACKWARD:
    std::cout << " Vhit::BOTTOM_BACKWARD";
    break;
  case Vhit::TOP_FORWARD:
    std::cout << " Vhit::TOP_FORWARD";
    break;
  case Vhit::TOP_BACKWARD:
    std::cout << " Vhit::TOP_BACKWARD";
    break;
  case Vhit::V_BELOW:
    std::cout << " Vhit::V_BELOW";
    break;
  case Vhit::V_ABOVE:
    std::cout << " Vhit::V_ABOVE";
    break;
  }
  std::cout << std::endl;
}
// ---------------- Vhit::Vcm_quadBox_Hit::set_region ------------------
// Identifies the relative locations
void Vhit::Vcm_quadBox_Hit::set_region( double V_cm, double V_trans, double V_lab_min,
  double V_lab_max )
{
  if( ( V_cm <= V_lab_min - V_trans ) ||
      ( V_cm <= V_trans - V_lab_max ) )
  {
    hit_corner = Vhit::V_BELOW;
  }
  else if( V_cm < V_trans + V_lab_max )
  {
    hit_corner = Vhit::V_INSIDE;
  }
  else
  {
    hit_corner = Vhit::V_ABOVE;
  }
}

// ************* Vcm_quadBox_Hit_F::lessthan_F *****************
// Function to order Vhit::Vcm_quadBox_Hit by its V_cm member
bool Vcm_quadBox_Hit_F::lessthan_F( Vhit::Vcm_quadBox_Hit first, Vhit::Vcm_quadBox_Hit second )
{
  return ( first.V_cm < second.V_cm );
}

// ************* class Vhit::Vcm_Vlab_hit_list *****************
// ----------- Vhit::Vcm_Vlab_hit_list::set_d_hit_G --------------
// Set the derivative of hit_G with respect to alpha
// where Ein = (1 - alpha)*Ein_0 + alpha*Ein_1.
void Vhit::Vcm_Vlab_hit_list::set_d_hit_G( double E_lab )
{
  double delta_Etrans = G1_data.get_Etrans( ) - G0_data.get_Etrans( );
  double delta_Ecm = G1_data.E_cm - G0_data.E_cm;
  if( ( G0_data.E_cm == 0.0 ) && ( G1_data.E_cm == 0.0 ) )
  {
    G0_data.dG = 2*E_lab*delta_Etrans -
      2*G0_data.get_Etrans( ) * delta_Etrans;
    G1_data.dG = 2*E_lab*delta_Etrans -
      2*G1_data.get_Etrans( ) * delta_Etrans;
  }
  else
  {
    G0_data.dG = 2*E_lab*( delta_Etrans + delta_Ecm ) -
      2*( G0_data.get_Etrans( ) - G0_data.E_cm )*( delta_Etrans - delta_Ecm );
    G1_data.dG =  2*E_lab*( delta_Etrans + delta_Ecm ) -
      2*( G1_data.get_Etrans( ) - G1_data.E_cm )*( delta_Etrans - delta_Ecm );
  }
}
// ----------- Vhit::Vcm_Vlab_hit_list::set_G_data --------------
// Calculates hit_G and d_hit_G for given E_out_lab
void Vhit::Vcm_Vlab_hit_list::set_G_data( double E_out_lab )
{
  G0_data.set_hit_G( E_out_lab );
  G1_data.set_hit_G( E_out_lab );
  set_d_hit_G( E_out_lab );
}
// ----------- Vhit::Vcm_Vlab_hit_list::do_check_roots ---------------
// Checks the G value for the computed incident energies giving intersections
void Vhit::Vcm_Vlab_hit_list::do_check_roots( int num_roots, double ein_1, double ein_2 )
{
  Maps::Ecm_intersect G_mid_data;
  G_mid_data.interpolate_Ein( ein_1, E_lab, G0_data, G1_data );
  //  std::cout << "ein_1: " << ein_1 << "  G: " << G_mid_data.G << std::endl;
  if( num_roots > 1 )
  {
    G_mid_data.interpolate_Ein( ein_2, E_lab, G0_data, G1_data );
    //    std::cout << "ein_2: " << ein_2 << "  G: " << G_mid_data.G << std::endl;
  }
}
// ----------- Vhit::Vcm_Vlab_hit_list::solve_hit_G0 --------------
// Solves hit_G = 0 using Taylor series about Ein_0; returns the number of real roots
int Vhit::Vcm_Vlab_hit_list::solve_hit_G0( double *ein_1, double *ein_2 )
{
  // This routine uses the fact that hit_G is quadratic in alpha
  // to find the zeros of G(alpha) = G0_data.G + G0_data.dG*alpha + A*alpha^2
  // with A such that G(1) = G1_data.G
  double A = G1_data.G - ( G0_data.G + G0_data.dG );
  double alpha1;
  double alpha2;
  int num_roots;
  bool check_roots = true;

  if( ( G0_data.E_cm == 0.0) && ( G1_data.E_cm == 0.0) )
  {
    num_roots = 1;
    alpha1 = -G0_data.dG/( 2*A );
    alpha2 = 0.0;
  }
  else
  {
    num_roots = math_F::quadratic( A, G0_data.dG, G0_data.G, &alpha1, &alpha2 );
  }

  if( num_roots == 1 )
  {
    if( ( alpha1 < 0.0 ) || ( alpha1 > 1.0 ) )
    {
      num_roots = 0;
    }
    else
    {
      *ein_1 = ( 1 - alpha1 )*G0_data.E_in + alpha1* G1_data.E_in;
    }
  }
  else if( num_roots == 2 )
  {
    if( alpha1 > 1.0 )
    {
      num_roots = 0;  // because alpha2 > alpha1
    }
    else if( alpha1 >= 0.0 )
    {
      *ein_1 = ( 1 - alpha1 )*G0_data.E_in + alpha1* G1_data.E_in;
      if( alpha2 > 1.0 )
      {
        num_roots = 1;  // because alpha2 > alpha1
      }
      else
      {
        *ein_2 = ( 1 - alpha2 )*G0_data.E_in + alpha2* G1_data.E_in;
      }
    }
    else if ( ( alpha2 >= 0.0 ) && ( alpha2 <= 1.0 ) )
    {
      num_roots = 1;
      *ein_1 = ( 1 - alpha2 )*G0_data.E_in + alpha2* G1_data.E_in;
    }
    else
    {
      num_roots = 0;
    }
  }
  if( check_roots && ( num_roots > 0 ) )
  {
    do_check_roots( num_roots, *ein_1, *ein_2 );
  }
  return num_roots;
}
// ----------- Vhit::Vcm_Vlab_hit_list::solve_hit_G1 --------------
// Solves hit_G = 0 using Taylor series about G1_data.E_in; returns the number of real roots
int Vhit::Vcm_Vlab_hit_list::solve_hit_G1( double *ein_1, double *ein_2 )
{
  // use beta = 1 - alpha
  // Find the zeros of G(beta) = G1_data.G - G1_data.dG*beta + B*beta^2 with 
  // B such that B(1) = G0_data.G
  double B = G0_data.G - ( G1_data.G - G1_data.dG );
  double beta1;
  double beta2;
  int num_roots;
  bool check_roots = true;

  if( ( G0_data.E_cm == 0.0) && ( G1_data.E_cm == 0.0) )
  {
    num_roots = 1;
    beta1 = G1_data.dG/( 2*B );
    beta2 = 0.0;
  }
  else
  {
    num_roots = math_F::quadratic( B, -G1_data.dG, G1_data.G, &beta1, &beta2 );
  }

  if( num_roots == 1 )
  {
    if( ( beta1 < 0.0 ) || ( beta1 > 1.0 ) )
    {
      num_roots = 0;
    }
    else
    {
      *ein_1 = ( 1 - beta1 )*G1_data.E_in + beta1* G0_data.E_in;
    }
  }
  else if( num_roots == 2 )
  {
    if( beta1 > 1.0 )
    {
      num_roots = 0;  // because beta2 > beta1
    }
    else if( beta1 >= 0.0 )
    {
      *ein_1 = ( 1 - beta1 )*G1_data.E_in + beta1* G0_data.E_in;
      if( beta2 > 1.0 )
      {
        num_roots = 1;  // because beta2 > beta1
      }
      else
      {
        *ein_2 = *ein_1;
        *ein_1 = ( 1 - beta2 )*G1_data.E_in + beta2* G0_data.E_in;
      }
    }
    else if ( ( beta2 >= 0.0 ) && ( beta2 <= 1.0 ) )
    {
      num_roots = 1;
      *ein_1 = ( 1 - beta2 )*G1_data.E_in + beta2* G0_data.E_in;
    }
    else
    {
      num_roots = 0;
    }
  }
  if( check_roots && ( num_roots > 0 ) )
  {
    do_check_roots( num_roots, *ein_1, *ein_2 );
  }
  return num_roots;
}
// ----------- Vhit::Vcm_Vlab_hit_list::find_bottom_hits --------------
// Finds the intersections with the bottom of a box
void Vhit::Vcm_Vlab_hit_list::find_bottom_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // don't bother with E_out_lab = 0
  if( E_out_lab <= 0.0 )
  {
    return;
  }
  E_lab = E_out_lab;
  double V_lab = sqrt( E_out_lab );
  set_G_data( E_out_lab );
  // for new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  double ein_1;
  double ein_2;
  int num_roots;

  // check for intersections
  // We intersect the bottom of the energy bin at the lower incident energy iff G0 > 0.
  // We intersect the bottom of the energy bin at the higher incident energy iff G1 > 0.
  if( std::abs( G0_data.G ) < std::abs( G1_data.G ) )
  {
    // use the Taylor series about alpha = 0
    num_roots = solve_hit_G0( &ein_1, &ein_2 );
  }
  else
  {
    num_roots = solve_hit_G1( &ein_1, &ein_2 );
  }
  if( num_roots == 2 )
  {
    if( G0_data.V_trans + G0_data.V_cm < V_lab )
    {
      // We're below the bottom of the bin at the lower incident energy
      // We hit the bottom of the bin for mu = 1 at incident energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::BOTTOM_IN;
      Ein_hits->push_back( Ein_eta_hit );
      // we hit the bottom of the bin for mu = -1 at incident energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else if( G1_data.V_trans + G1_data.V_cm < V_lab )
    {
      // We're below the bottom of the bin at the higher incident energy
      // We hit the bottom of the bin for mu = -1 at incident energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
      // we hit the bottom of the bin for mu = 1 at incident energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::BOTTOM_OUT;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else if( G0_data.G > 0.0 )
    {
      // We intersect the bottom of the bin at the lower incident energy
      // transition to above the bottom at lab mu = -1 and energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
      // transition to hitting the bottom at lab mu = -1 and energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else
    {
      // We are above the bottom of the bin at the lower incident energy
      // transition to hitting the bottom at lab mu = -1 and energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
      // transition to above the bottom at lab mu = -1 and energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
    }
  }
  else if( num_roots == 1 )
  {
    Ein_eta_hit.E_in = ein_1;
    if( G0_data.V_trans + G0_data.V_cm < V_lab )
    {
      // We are below the bottom of the energy bin at the lower incident energy
      Ein_eta_hit.hit_edge = Box::BOTTOM_IN;
    }
    else if( G1_data.V_trans + G1_data.V_cm < V_lab )
    {
      // We are below the bottom of the energy bin at the higher incident energy
      Ein_eta_hit.hit_edge = Box::BOTTOM_OUT;
    }
    else if( G0_data.G > 0.0 )
    {
      // transition from hitting the bottom of the bin at lab mu < 0 to being above
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_IN;
    }
    else
    {
      // transition from being above the bottom of the bin at lab mu < 0 to hitting it
      Ein_eta_hit.hit_edge = Box::BOTTOM_CORNER_OUT;
    }
    Ein_hits->push_back( Ein_eta_hit );
  }
}
// ----------- Vhit::Vcm_Vlab_hit_list::find_top_hits --------------
// Finds the intersections with the top of a box
void Vhit::Vcm_Vlab_hit_list::find_top_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  E_lab = E_out_lab;
  double V_lab = sqrt( E_out_lab );
  set_G_data( E_out_lab );
  // for new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  double ein_1;
  double ein_2;
  int num_roots;

  // check for intersections
  if( std::abs( G0_data.G ) < std::abs( G1_data.G ) )
  {
    // use the Taylor series about alpha = 0
    num_roots = solve_hit_G0( &ein_1, &ein_2 );
  }
  else
  {
    num_roots = solve_hit_G1( &ein_1, &ein_2 );
  }
  if( num_roots == 2 )
  {
    if( std::abs( G0_data.V_trans - G0_data.V_cm ) > V_lab )
    {
      // We're above the top of the bin at the lower incident energy
      // e.g., forward emission near the threshold
      // We hit the top of the bin for lab mu = 1 at incident energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::TOP_IN;
      Ein_hits->push_back( Ein_eta_hit );
      // we hit the top of the bin for lab mu = -1 at incident energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else if( std::abs( G1_data.V_trans - G1_data.V_cm ) > V_lab )
    {
      // We're above the top of the bin at the higher incident energy
      // transition from below the top to hitting it at mu = 1 and incident energy ein_2
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
      // We exit the top of the bin for mu = -1 at incident energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::TOP_OUT;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else if( G0_data.G > 0.0 )
    {
      // We intersect the top of the bin at the lower incident energy
      // transition to below the top at lab mu = 1 and energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
      // transition to hitting the top at lab mu = 1 and energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
    }
    else
    {
      // We are below the top of the bin at the lower incident energy
      // transition to hitting the top at lab mu = 1 and energy ein_1
      Ein_eta_hit.E_in = ein_1;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_IN;
      Ein_hits->push_back( Ein_eta_hit );
      // transition to below the top at lab mu = 1 and energy ein_2
      Ein_eta_hit.E_in = ein_2;
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_OUT;
      Ein_hits->push_back( Ein_eta_hit );
    }
  }
  else if( num_roots == 1 )
  {
    Ein_eta_hit.E_in = ein_1;
    if( std::abs( G0_data.V_trans - G0_data.V_cm ) > V_lab )
    {
      // We are above the top of the energy bin at the lower incident energy
      Ein_eta_hit.hit_edge = Box::TOP_IN;
    }
    else if( std::abs( G1_data.V_trans - G1_data.V_cm ) > V_lab )
    {
      // We are above the top of the energy bin at the higher incident energy
      Ein_eta_hit.hit_edge = Box::TOP_OUT;
    }
    else if( G0_data.G > 0.0 )
    {
      // transition from hitting the top of the bin at lab mu = 1 to being below
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_IN;
    }
    else
    {
      // transition from being below the top of the bin at lab mu = 1 to hitting it
      Ein_eta_hit.hit_edge = Box::TOP_CORNER_OUT;
    }
    Ein_hits->push_back( Ein_eta_hit );
  }
}
// ----------- Vhit::Vcm_Vlab_hit_list::test_left --------------
// Where do we hit the left-hand side of the box?
void Vhit::Vcm_Vlab_hit_list::test_left( double E_in, std::vector< double >::const_iterator Eout_ptr )
{
  // for new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  Ein_eta_hit.E_in = E_in;

  Vhit::Vcm_Vlab_hit_list::iterator first_hit = begin( );
  if( first_hit->hit_edge == Box::TOP_IN )
  {
    Maps::Ecm_intersect G_data;  // for checking the values at incident energy E_in
    G_data.interpolate_Ein( E_in, *Eout_ptr, G0_data, G1_data );
    std::vector< double >::const_iterator next_Eout = Eout_ptr;
    ++next_Eout;
    double V_lab = sqrt( *next_Eout );
    if( G_data.V_cm - G_data.V_trans > V_lab )
    {
      Ein_eta_hit.hit_edge = Box::ABOVE;
    }
    else if( G_data.V_trans - G_data.V_cm > V_lab )
    {
      Ein_eta_hit.hit_edge = Box::ABOVE_FORWARD;
    }
  }
  else if( first_hit->hit_edge == Box::BOTTOM_IN )
  {
    Ein_eta_hit.hit_edge = Box::BELOW;
  }
  else
  {
    Ein_eta_hit.hit_edge = Box::INSIDE;
  }
  // prepend this entry
  push_front( Ein_eta_hit );
}
// ----------- Vhit::Vcm_Vlab_hit_list::test_right --------------
// Where do we hit the right-hand side of the box?
void Vhit::Vcm_Vlab_hit_list::test_right( double E_in, std::vector< double >::const_iterator Eout_ptr )
{
  // for new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  Ein_eta_hit.E_in = E_in;
  // we already have at least one entry (at the left edge)
  Vhit::Vcm_Vlab_hit_list::iterator last_hit = end( );
  --last_hit;

  if( ( last_hit->hit_edge == Box::TOP_OUT ) ||
      ( last_hit->hit_edge == Box::ABOVE ) )
  {
    Maps::Ecm_intersect G_data;  // for checking the values at incident energy E_in
    G_data.interpolate_Ein( E_in, *Eout_ptr, G0_data, G1_data );
    std::vector< double >::const_iterator next_Eout = Eout_ptr;
    ++next_Eout;
    double V_lab = sqrt( *next_Eout );
    if( G_data.V_cm - G_data.V_trans > V_lab )
    {
      Ein_eta_hit.hit_edge = Box::ABOVE;
    }
    else
    {
      Ein_eta_hit.hit_edge = Box::ABOVE_FORWARD;
    }
  }
  else if( ( last_hit->hit_edge == Box::BOTTOM_OUT ) ||
      ( last_hit->hit_edge == Box::BELOW ) )
  {
    Ein_eta_hit.hit_edge = Box::BELOW;
  }
  else
  {
    Ein_eta_hit.hit_edge = Box::INSIDE;
  }

  push_back( Ein_eta_hit );
}
// ----------- Vhit::Vcm_Vlab_hit_list::test_sides --------------
// No top or bottom intersections
void Vhit::Vcm_Vlab_hit_list::test_sides( std::vector< double >::const_iterator Eout_ptr,
  double E_in_left, double E_in_right )
{
  double V_lab = sqrt( *Eout_ptr );
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  Maps::Ecm_intersect G_data;  // for checking intermediate values at the average incident energy
  G_data.interpolate_Ein( 0.5*(E_in_left + E_in_right), *Eout_ptr, G0_data, G1_data );
  if( G_data.V_trans + G_data.V_cm < V_lab )
  {
    Ein_eta_hit.hit_edge = Box::BELOW;
  }
  else
  {
    std::vector< double >::const_iterator next_Eout = Eout_ptr;
    ++next_Eout;
    V_lab = sqrt( *next_Eout );
    G_data.interpolate_Ein( 0.5*(E_in_left + E_in_right), *next_Eout, G0_data, G1_data );
    if( G_data.V_cm - G_data.V_trans > V_lab )
    {
      Ein_eta_hit.hit_edge = Box::ABOVE;
    }
    else if( G_data.V_trans - G_data.V_cm > V_lab )
    {
      Ein_eta_hit.hit_edge = Box::ABOVE_FORWARD;
    }
    else
    {
      Ein_eta_hit.hit_edge = Box::INSIDE;
    }
  }
  Ein_eta_hit.E_in = E_in_left;
  push_back( Ein_eta_hit );
  Ein_eta_hit.E_in = E_in_right;
  push_back( Ein_eta_hit );
}
// ----------- Vhit::Vcm_Vlab_hit_list::consistency --------------
// Checks for inconsistencies caused by peculiarities of real arithmetic.
// Returns true if there are no inconsistencies.
bool Vhit::Vcm_Vlab_hit_list::consistency( )
{
  bool Above = false;
  bool Below = false;
  bool Inside = false;
  bool TopCornerAbove = false;
  bool BottomCornerBelow = false;
  bool is_OK = true;

  // loop through the intersections
  for( Vhit::Vcm_Vlab_hit_list::const_iterator this_link = begin( );
       this_link != end( ); ++this_link )
  {
    if( this_link->hit_edge == Box::INSIDE )
    {
      if( Above || Below )
      {
        is_OK = false;
      }
      Inside = true;
    }
    else if( ( this_link->hit_edge == Box::ABOVE ) || ( this_link->hit_edge == Box::ABOVE_FORWARD ) )
    {
      Above = true;
      TopCornerAbove = true;
      Inside = false;
    }
    else if( this_link->hit_edge == Box::TOP_IN )
    {
      if( Above )
      {
        Above = false;
        TopCornerAbove = true;
        Inside = true;
      }
      else
      {
        is_OK = false;
      }
    }
    else if( this_link->hit_edge == Box::BELOW )
    {
      Below = true;
      BottomCornerBelow = true;
      Inside = false;
    }
    else if( this_link->hit_edge == Box::BOTTOM_IN )
    {
      if( Below )
      {
        Below = false;
        BottomCornerBelow = true;
        Inside = true;
      }
      else
      {
        is_OK = false;
      }
    }
    else if( this_link->hit_edge == Box::TOP_OUT )
    {
      if( Inside )
      {
        Above = true;
        Inside = false;
      }
      else
      {
        is_OK = false;
      }
    }
    else if( this_link->hit_edge == Box::BOTTOM_OUT )
    {
      if( Inside )
      {
        Below = true;
        Inside = false;
      }
      else
      {
        is_OK = false;
      }
    }
    else if( this_link->hit_edge == Box::TOP_CORNER_OUT )
    {
      if( TopCornerAbove || !Inside )
      {
        is_OK = false;
      }
      else
      {
        TopCornerAbove = true;
      }
    }
    else if( this_link->hit_edge == Box::BOTTOM_CORNER_OUT )
    {
      if( BottomCornerBelow || !Inside )
      {
        is_OK = false;
      }
      else
      {
        BottomCornerBelow = true;
      }
    }
    else if( this_link->hit_edge == Box::TOP_CORNER_IN )
    {
      if( !TopCornerAbove && !Inside )
      {
        is_OK = false;
      }
      else
      {
        TopCornerAbove = false;
      }
    }
    else if( this_link->hit_edge == Box::BOTTOM_CORNER_IN )
    {
      if( !BottomCornerBelow && !Inside )
      {
        is_OK = false;
      }
      else
      {
        BottomCornerBelow = false;
      }
    }
  }
  if( !is_OK )
  {
    Msg::Warning( "Vhit::Vcm_Vlab_hit_list::consistency", "Check coding" );
  }
  return is_OK;
}
