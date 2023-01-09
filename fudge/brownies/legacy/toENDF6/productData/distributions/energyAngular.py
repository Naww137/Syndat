# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import multiD_XYs as multiD_XYsModule

from fudge.productData.distributions import energyAngular as energyAngularModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    subform = self.energyAngularSubform
    frame = self.productFrame
    if( hasattr( subform, 'toENDF6' ) ) :
        LAW, MF6 = subform.toENDF6( flags, targetInfo )
        gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
    else :
        print( 'WARNING: energyAngular subform "%s" has no toENDF6 method' % subform.moniker )

energyAngularModule.Form.toENDF6 = toENDF6

#
# XYs3d
#
def toENDF6( self, flags, targetInfo ) :

    EInInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    EpInterpolation0 = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self[0].interpolation, self[0].interpolationQualifier )
    if( EpInterpolation0 == 1 ) :       # flat interpolation
        LEP = 1
    elif( EpInterpolation0 == 2 ) :     # lin-lin interpolation
        LEP = 2
    else :
        raise 'hell - fix me'
    ENDFDataList = endfFormatsModule.endfInterpolationList( [ len( self ), EInInterpolation ] )
    energyInFactor = PQUModule.PQU( 1, self.axes[3].unit ).getValueAs( 'eV' )
    energyPFactor = PQUModule.PQU( 1, self.axes[2].unit ).getValueAs( 'eV' )
    if( not( isinstance( self[0], multiD_XYsModule.XYs2d ) ) ) : raise 'hell - fix me'
    if( isinstance( self[0][0], energyAngularModule.Legendre ) ) :
        for energyIn in self :
            EpInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( energyIn.interpolation, energyIn.interpolationQualifier )
            if( EpInterpolation != EpInterpolation0 ) : raise 'hell - fix me'
            NA, data = 0, []
            for energy_p in energyIn :
                if( not( isinstance( energy_p, energyAngularModule.Legendre ) ) ) : raise 'hell - fix me'
                NA = max( len( energy_p ), NA )
                coefficients = [ coefficient / energyPFactor for coefficient in energy_p ]
                data += [ energy_p.outerDomainValue * energyPFactor ] + coefficients
            ENDFDataList.append( endfFormatsModule.endfContLine( 0, energyIn.outerDomainValue * energyInFactor, 0, NA - 1, len( data ), len( data ) / ( NA + 1 ) ) )
            ENDFDataList += endfFormatsModule.endfDataList( data )
        LAW = 1
        LANG = 1
    elif( isinstance( self[0][0], energyAngularModule.XYs1d ) ) :
        muInterpolation0 = gndsToENDF6Module.gndsToENDFInterpolationFlag(self[0][0].interpolation)
        for energyIn in self :
            EpInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( energyIn.interpolation, energyIn.interpolationQualifier )
            if( EpInterpolation != EpInterpolation0 ) : raise 'hell - fix me'
            data = []
            NAp = 2 * max([len(energy_p) for energy_p in energyIn])
            data = []
            for energy_p in energyIn :
                if( not( isinstance( energy_p, energyAngularModule.XYs1d ) ) ) : raise 'hell - fix me'
                muInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(energy_p.interpolation)
                if muInterpolation != muInterpolation0:
                    raise Exception("Inconsistent P(mu | E, E') interpolation not supported by ENDF-6")
                energy_p = energy_p.convertAxisToUnit( 0, '1/eV' )
                NA = 2 * len( energy_p )
                while NA != NAp:
                    delta = 0.0
                    for index, muP in enumerate(energy_p):      # Find the widest gap and put a point in the middle of it.
                        mu2 = muP[0]
                        if index != 0:
                            if mu2 - mu1 > delta:
                                delta = mu2 - mu1
                                mu = 0.5 * (mu1 + mu2)
                        mu1 = mu2
                    energy_p.setValue(mu, energy_p.evaluate(mu))
                    NA += 2
                f0 = energy_p.integrate()
                data += [ energy_p.outerDomainValue * energyPFactor, f0 ]
                if( f0 == 0 ) :                 # Special case for when f0 == 0 in the ENDF file.
                    for mu, f1 in energy_p :
                        data.append( mu )
                        data.append( 0.5 )
                else :
                    for mu, f1 in energy_p.normalize( ).copyDataToXYs( ) :
                        data.append( mu )
                        data.append( f1 )
            ENDFDataList.append( endfFormatsModule.endfContLine( 0, energyIn.outerDomainValue * energyInFactor, 0, NA, len( data ), len( energyIn ) ) )
            ENDFDataList += endfFormatsModule.endfDataList( data )
        LANG = 10 + muInterpolation0
        LAW = 1
    ENDFDataList.insert( 0, endfFormatsModule.endfContLine( 0, 0, LANG, LEP, 1, len( self ) ) )
    return( LAW, ENDFDataList )

energyAngularModule.XYs3d.toENDF6 = toENDF6
