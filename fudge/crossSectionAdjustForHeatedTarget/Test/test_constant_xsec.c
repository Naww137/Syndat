/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*
*    This routine test the results for a constant cross-section to theory.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../Src/crossSectionAdjustForHeatedTarget.h"

static long nPoints = 200;
static double Temperature = 1e-3, EMin = 1e-10, EMax = 20., massRatio = 8.;

static double Theory( double E );
static void PrintMsg( char *s, int Usage );
/*
***************************************
*/
int main( int argc, char **argv ) {

    int err, checker = 0;
    long i;
    double f, *E_cs, E, *p, t, rErr, absrErr, rErrMax = 0.;
    char *e;
    crossSectionAdjustForHeatedTarget_info info;
    FILE *fOut;

    info.mode = 0;
    info.verbose = 0;

    if( argc > 6 ) PrintMsg( "too many parameters", 1 );
    if( argc > 1 ) {
        if( strcmp( argv[1], "-check" ) == 0 ) {
            checker = 1; }
        else {
            if( argv[1][0] == '?' ) PrintMsg( NULL, 1 );
            f = strtod( argv[1], &e );
            if( ( e == argv[1] ) || ( *e != 0 ) ) PrintMsg( "cannot convert temperature parameter to double", 1 );
            Temperature = f;
            if( argc > 2 ) {
                f = strtod( argv[2], &e );
                if( ( e == argv[2] ) || ( *e != 0 ) ) PrintMsg( "cannot convert EMin parameter to double", 1 );
                EMin = f;
                if( argc > 3 ) {
                    f = strtod( argv[3], &e );
                    if( ( e == argv[3] ) || ( *e != 0 ) ) PrintMsg( "cannot convert EMax parameter to double", 1 );
                    EMax = f;
                    if( argc > 4 ) {
                        f = strtod( argv[4], &e );
                        if( ( e == argv[4] ) || ( *e != 0 ) ) PrintMsg( "cannot convert massRatio parameter to double", 1 );
                        massRatio = f;
                        if( argc > 5 ) {
                            i = strtol( argv[5], &e, 10 );
                            if( ( e == argv[5] ) || ( *e != 0 ) ) PrintMsg( "cannot convert nPoints parameter to long", 1 );
                            nPoints = i;
                        }
                    }
                }
            }
        }
    }
    if( ( EMin >= EMax ) || ( EMin <= 0. ) ) PrintMsg( "EMax must be greater than EMin and EMin must be greater than 0.", 0 );
    if( nPoints < 2 ) PrintMsg( "nPoints must be greater than 2", 0 );
    if( !checker ) printf( "# T = %e  EMin = %e  EMax = %e  massRatio = %e  nPoints = %ld\n", Temperature, EMin, EMax, massRatio, nPoints );
    f = pow( EMax / EMin, 1. / ( nPoints - 1 ) );
    E_cs = malloc( 2 * nPoints * sizeof( double ) );
    if( E_cs == NULL ) PrintMsg( "cannot allocate memory for E_cs", 0 );

    E = EMin;
    p = E_cs;
    for( i = 0; i < nPoints; i++, E *= f ) {
        *(p++) = E;
        *(p++) = 1.;
    }
    p[-2] = EMax;

    if( !checker ) {
        fOut = fopen( "test_constant_xsec.data", "w" );
        for( i = 0, p = E_cs; i < nPoints; i++, p += 2 ) fprintf( fOut, "%13.5e %13.5e\n", *p, *(p+1) );
        fclose( fOut );
    }

    err = crossSectionAdjustForHeatedTarget( crossSectionAdjustForHeatedTarget_limit_constant, crossSectionAdjustForHeatedTarget_limit_constant,
        &info, 2e-10, massRatio, Temperature, 1e-4, nPoints, E_cs, &p );
    if( err < 0 ) {
        fprintf( stderr, "\ncrossSectionAdjustForHeatedTarget returned err = %d, bytes = %d\n", err, info.bytes );
        PrintMsg( "bad crossSectionAdjustForHeatedTarget return value", 0 );
    }
    if( !checker ) printf( "# WarningStats = %d evaluations = %d, pairs = %d\n", info.WarningStats, info.evaluations, err );
    for( i = 0; i < err; i++, p += 2 ) {
        t = Theory( p[0] );
        rErr = p[1] / t - 1.;
        absrErr = fabs( rErr );
        if( absrErr > rErrMax ) rErrMax = absrErr;
        if( !checker ) printf( "%16.8e %16.8e %16.8e %12.4e %s\n", p[0], p[1], t, rErr, ( absrErr > 1e-6 ? " ****" : "" ) );
    }
    if( rErrMax > 1e-8 ) exit( EXIT_FAILURE );
    exit( EXIT_SUCCESS );
}
/*
***************************************
*/
static double Theory( double E ) {

    double x, x2 = massRatio * E / Temperature, sqrtpi = 1.7724538509055160273, cs;

    x = sqrt( x2 );
    cs =  erf( x ) * ( 1. + 0.5 / x2 ) + exp( -x2 ) / ( sqrtpi * x );
    return( cs );
}
/*
***************************************
*/
static void PrintMsg( char *s, int Usage ) {

    if( s != NULL ) {
        fprintf( stderr, "\nError - %s\n\n", s );
    }

    if( Usage ) {
        printf( "\nUSAGE:\n" );
        printf( "    test_constant_xsec Temperature, EMin, EMax, massRatio, nPoints\n" );
        printf( "\n" );
        printf( "PARAMETERS\n" );
        printf( "    Temperature    Target's temperature (default = %e).\n", Temperature );
        printf( "    EMin           Lowest cross section energy point (default = %e).\n", EMin );
        printf( "    EMax           Highest cross section energy point (default = %e).\n", EMax );
        printf( "    massRatio      Target to incident particle mass ratio (default = %e).\n", massRatio );
        printf( "    nPoints        Number of pairs of (E, xsec) points (default = %ld).\n", nPoints );
        printf( "\n" );
    }
    exit( EXIT_FAILURE );
}
