/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cbdfls.h>

#include "crossSectionAdjustForHeatedTarget.h"

#define EOD "                                                                       1\n"
#define bufferSize 256

static int line = 0, HeatedDataSets = 0;

void PrintMsg( char *s1, char *s2 );
/*
***************************************************
*/
int main( int argc, char **argv ) {

    int i, n_pairs, n_allocated = 0, err, targetZA, yiZA;
    char buffer[bufferSize], newheader1[bufferSize], header1[bufferSize], header2[bufferSize];
    char Str[bufferSize], outputFileName[1024], *c, *p1, *p2;
    double *E_cs = NULL, E, cs, T_MeV, targetT_MeV, *E_cs_out, massTarget, massyi, massRatio, f_interpolation = 0.01, EMin = 1e-11;
    FILE *fIn, *fOut;
    cbdfls_file *bdfls;
    cbdflsErrors Error;
    crossSectionAdjustForHeatedTarget_limit lowerlimit = crossSectionAdjustForHeatedTarget_limit_constant;
    crossSectionAdjustForHeatedTarget_info info;

    info.mode = 0;
    info.verbose = 2;

    bdfls = cbdflsOpen( "", &Error );
    if( bdfls == NULL ) PrintMsg( "cannot open bdfls file", NULL );

    if( ( argc < 3 ) || ( argc > 6 ) ) PrintMsg( "need temperature (MeV), endl input file path [, f_interpolation, EMin, thin]", NULL );
    if( sscanf( argv[1], "%le", &T_MeV ) != 1 ) PrintMsg( "cannot convert temperature to double", argv[1] );

    if( ( c = strrchr( argv[2], '/' ) ) == NULL ) c = argv[2];
    if( *c == '/' ) c++;
    if( T_MeV < 1e-3 ) {
        sprintf( outputFileName, "%s_T%.3feV", c, T_MeV * 1e6 ); }
    else if( T_MeV < 1. ) {
        sprintf( outputFileName, "%s_T%.3fkeV", c, T_MeV * 1e3 ); }
    else {
        sprintf( outputFileName, "%s_T%.3fMeV", c, T_MeV );
    }
    printf( "%s\n", outputFileName );
    if( ( fOut = fopen( outputFileName, "w" ) ) == NULL ) PrintMsg( "cannot open output file", outputFileName );
    if( ( fIn = fopen( argv[2], "r" ) ) == NULL ) PrintMsg( "cannot open endl file", argv[2] );

    if( argc > 3 ) {
        if( sscanf( argv[3], "%le", &f_interpolation ) != 1 ) PrintMsg( "cannot convert f_interpolation to double", argv[3] );
        if( argc > 4 ) {
            if( sscanf( argv[4], "%le", &EMin ) != 1 ) PrintMsg( "cannot convert EMin to double", argv[4] );
            if( argc > 5 ) info.mode |= crossSectionAdjustForHeatedTarget_mode_do_not_thin;
        }
    }

    while( 1 ) {
        line++;
        if( fgets( header1, bufferSize, fIn ) == NULL ) break;
        HeatedDataSets++;
        if( strlen( header1 ) < 71 ) PrintMsg( "invalid header1", header1 );
        strcpy( Str, header1 );
        Str[6] = 0;
        if( sscanf( Str, "%d", &targetZA ) != 1 ) PrintMsg( "error converting target ZA to integer", header1 );
        Str[6] = header1[6];
        Str[9] = 0;
        if( sscanf( &(Str[7]), "%d", &yiZA ) != 1 ) PrintMsg( "error converting target ZA to integer", header1 );
        Str[9] = header1[9];
        Str[70] = Str[69];
        Str[69] = Str[68];
        if( Str[69] == ' ' ) Str[69] = '0';
        Str[68] = Str[67];
        Str[67] = 'e';
        if( sscanf( &(Str[59]), "%le", &targetT_MeV ) != 1 ) PrintMsg( "error converting target temperature to double", header1 );
        if( targetT_MeV > T_MeV ) {
            sprintf( Str, "target temperature = %e  temperature to heat to = %e", targetT_MeV, T_MeV );
            PrintMsg( "target temperature greater than temperature to heat to", Str );
        }
        strcpy( newheader1, header1 );
        sprintf( Str, "%12.5e\n", T_MeV );
        for( p1 = &(Str[8]), p2 = &(Str[9]); *p2 != 0; p1++, p2++ ) *p1 = *p2;
        *p1 = 0;
        strcpy( &(newheader1[59]), Str );

        massTarget = cbdflsGetMass( bdfls, targetZA );
        if( massTarget < 0. ) {
            sprintf( Str, "targetZA = %d", targetZA );
            PrintMsg( "target mass not in bdfls", Str );
        }
        massyi = cbdflsGetMass( bdfls, yiZA );
        if( massyi < 0. ) {
            sprintf( Str, "incident particle ZA = %d", targetZA );
            PrintMsg( "incident particle mass not in bdfls", Str );
        }
        massRatio = massTarget / massyi;
        printf( "massRatio = %.12e\n", massRatio );

        line++;
        if( fgets( header2, bufferSize, fIn ) == NULL ) PrintMsg( "cannot read second header line", NULL );
        n_pairs = 0;
        while( ( line++, fgets( buffer, bufferSize, fIn ) ) != NULL ) {
            if( strcmp( buffer, EOD ) == 0 ) break;
            if( n_pairs == n_allocated ) {
                n_allocated += 10000;
                E_cs = (double *) realloc( E_cs, 2 * n_allocated * sizeof( double ) );
                if( E_cs == NULL ) {
                    sprintf( Str, "cannot allocate memory; n_allocated = %d", n_allocated );
                    PrintMsg( Str, buffer );
                }
            }
            if( sscanf( buffer, "%le %le", &E, &cs ) != 2 ) PrintMsg( "converting E and cs to doubles", buffer );
            E_cs[2 * n_pairs] = E;
            E_cs[2 * n_pairs + 1] = cs;
            n_pairs++;
        }
        printf( "n_pairs = %d\n", n_pairs );

        if( E_cs[0] > 1.1e-10 ) lowerlimit = crossSectionAdjustForHeatedTarget_limit_threshold;
        err = crossSectionAdjustForHeatedTarget( lowerlimit, crossSectionAdjustForHeatedTarget_limit_constant, &info, EMin, 
            massRatio, T_MeV, f_interpolation, n_pairs, E_cs, &E_cs_out );
        if( err < 0 ) {
            sprintf( Str, "\nError - crossSectionAdjustForHeatedTarget returned %d\n\n", err );
            PrintMsg( Str, NULL );
        }

        printf( "n_pairs(out) = %d\n", err );
        fprintf( fOut, "%s%s", newheader1, header2 );
        for( i = 0; i < err; i++ ) fprintf( fOut, " %13.7E %12.5E\n", E_cs_out[2 * i], E_cs_out[2 * i + 1] );
        fprintf( fOut, EOD );
    }
    if( HeatedDataSets == 0 ) PrintMsg( "cannot read first header line", NULL );
    cbdflsRelease( bdfls );
    fclose( fIn );

    exit( EXIT_SUCCESS );
}
/*
***************************************************
*/
void PrintMsg( char *s1, char *s2 ) {

    fprintf( stderr, "\nError: %s\n", s1 );
    fprintf( stderr, "  line = %d  HeatedDataSets = %d\n", line, HeatedDataSets );
    if( s2 != NULL ) fprintf( stderr, s2 );
    fprintf( stderr, "\n" );
    exit( EXIT_FAILURE );
}
