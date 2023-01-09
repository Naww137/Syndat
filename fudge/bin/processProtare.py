#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import shutil
import subprocess

import argparse

from PoPs import IDs as IDsPoPsModule
from PoPs.families import nuclide as nuclidePoPsModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from xData import enums as xDataEnumsModule

from pqu import PQU as PQUModule

from LUPY import times as timesModule
from LUPY import locateFudgeBin as locateBinFolderModule

from LUPY import argumentsForScripts as argumentsForScriptsModule

import fudge

from fudge import enums as enumsModule
from fudge import physicalQuantity as physicalQuantityModule
from fudge import styles as stylesModule

from fudge.processing import flux as fluxModule
from fudge.processing import transportables as transportablesModule
from fudge.processing import group as groupModule

timer = timesModule.Times( )

dateTimeStr = None

legendreMaxDefault = 9
reconAccuracyDefault = 1e-3
CoulombPlusNuclearMuCutOffDefault = 0.94                                            # This is the mu cutoff used for ENDL data.
temperatureUnitDefault = 'MeV/k'
energyUnitDefault = 'MeV'
outputDefault = 'SRC'
tagDefault = 'proc'

RCprefixDefault = 'recon'
muCutoffPrefixDefault = 'muCutoff'
AEPprefixDefault = 'apd'
HTprefixDefault = 'heated'
MGprefixDefault = 'MultiGroup'
MCprefixDefault = 'MonteCarlo'
UPprefixDefault = 'UpScatter'

description = """This script processes all data in a GNDS file as requested by input arguments. In addition to entering arguments on the comannd
line, arguments can also be read from an input file by specifying an input file on the command line. The character '@' must prefix
the input file's name. For example, to include the arguments in a file named pp.input execute this script as

    processProtare.py -mc @pp.input -vvv -up eval.xml proc.xml

Multiple input files and arguments can exists on the same command line. The lines in the input file must contain only valid argumments 
and/or a comment. On each line, all character at and after a "#" character are treated as a comment and are ignored.  The following is 
an example of a valid input file:

    -t 2.586e-8 -mg -g
    -t 1e-7
    -t 1e-6
    n=LLNL_gid_7 --groupFile    # This is a comment.
    groups.xml                  # Another comment.
    --fluxFile ./fluxes.xml --energyUnit eV
"""

parserPreview = argparse.ArgumentParser(fromfile_prefix_chars='@', add_help=False)
parserPreview.add_argument('-mg', '--MultiGroup', action='store_true')
argsPreview, dummy = parserPreview.parse_known_args()
multigroupPresent = argsPreview.MultiGroup

class ProcessProtareArgumentParser( argparse.ArgumentParser ) :
    """Sub-classes the argparse.ArgumentParser class to replace the convert_arg_line_to_args function."""

    def convert_arg_line_to_args( self, line ) :
        """Returns the arguments on line. The line is read from a file and passed to this method."""

        return( line.split( '#' )[0].split( ) )

parser = ProcessProtareArgumentParser(description=description, fromfile_prefix_chars='@', formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments( parser )

parser.add_argument( 'outputFile', nargs = '?', default = None,                                         help = 'Output file.' )

parser.add_argument( '--exit', action = 'store', type = int, default = 9999,                            help = 'Exits after phase n where is the value and is an integer.' )
parser.add_argument( '--tag', type = str, default = tagDefault,                                         help = 'Tag to indicate data has been processed.' )
parser.add_argument( '-o', '--output', default = outputDefault,                                         help = 'Directory to write output file. If "SRC", writes to source directory. Default is "%s".' % outputDefault )
parser.add_argument( '--writeConvertedUnits', action = 'store_true', default = False,                   help = 'Write data to file after units converted.' )

parser.add_argument( '--energyUnit', type = str, default = energyUnitDefault,                           help = 'Energy unit to convert to. Default is "%s."' % energyUnitDefault )

# parser.add_argument( '--reconstruct', type = str, choices = [ 'all','crossSection','angular' ], default = 'crossSection', 
#                                                                                                       help = 'What kind of reconstructing should be done ("all", "crossSection", "angular").' )
parser.add_argument( '--CoulombPlusNuclearMuCutOff', type = float, default = CoulombPlusNuclearMuCutOffDefault,
                                                                                                        help = 'For Coulomb + nuclear elastic scattering mu is limited to [ -1, muCutOff ] to make the cross section finite. \nFor identical particles mu is limited to [ -muCutOff, muCutOff ]. muCutOff must be in the range ( -1 to 1 ). Default is "%s".' % CoulombPlusNuclearMuCutOffDefault )
parser.add_argument( '--doNotAddNuclearPlusInterference', action = 'store_true',                        help = 'If not present and projectile is a charged particles, an elastic reaction without Rutherford scattering is added to the incompleteReactions node.' )

parser.add_argument( '-t', '--temperatures', type = float, action = 'append', default = None,           help = 'Temperatures for heating. Use one -t option for each temperature.' )
parser.add_argument( '--temperatureUnit', type = str, default = temperatureUnitDefault,                 help = 'Temperature unit to convert to. Default is "%s".' % temperatureUnitDefault )

parser.add_argument( '--legendreMax', type = int, default = legendreMaxDefault,                         help = 'Maximum Legendre order for Sn prcessed data. Default is "%s".' % legendreMaxDefault )
parser.add_argument( '-mg', '--MultiGroup', action = 'store_true',                                      help = 'Flag to turn on multi-group processing.' )
parser.add_argument( '-up', '--UpScatter', action = 'store_true',                                       help = 'Flag to turn on multi-group upscatter processing.' )

parser.add_argument( '-mc', '--MonteCarlo', action = 'store_true',                                      help = 'Flag to turn on Monte Carlo processing.' )

parser.add_argument( '--fluxFile',  type=str, required=multigroupPresent,                               help = 'File containing a "fluxes" node with a suite of 3d flux nodes (i.e., a flux f given as f(T,E,mu) where T is temperature, \nE is projectile energy and mu is cos of angle). If specified, bdfls file is ignored.' )
parser.add_argument( '--fluxID',    type=str, required=multigroupPresent,                               help = 'Flux ID from flux or bdfls file' )
parser.add_argument( '--groupFile', type=str, required=multigroupPresent,                               help = 'File containing a "groups" node with a suite of group nodes. If specified, bdfls file is ignored.' )
parser.add_argument( '-g', '--gid', type=str, required=multigroupPresent, action='append',              help = 'Specifies a multi-group ID for a specific particle as <particle>=multiGroupID (e.g., --gid n=LLNL_gid_7).' )

parser.add_argument( '--formatVersion', default = GNDS_formatVersionModule.default, choices = GNDS_formatVersionModule.allowed,
                                                                                                        help = 'Specifies the GNDS format for the outputted file. Default = "%s".' % GNDS_formatVersionModule.default )

parser.add_argument( '--cullProcessedData', action = 'store_true',                                      help = 'If set, all existing processed data are removed before data are processed.' )

parser.add_argument( '--threads',type = int, default = 1,                                               help = 'Number of threads to use for temperatures loop.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,                                  help = 'Enable verbose output.' )
parser.add_argument( '--reconAccuracy', type = float, default = reconAccuracyDefault,                   help = 'Accuracy for reconstructing resonances. Default is "%.1e".' % reconAccuracyDefault )
parser.add_argument( '--restart', action = 'store_true',                                                help = 'Continue previous incomplete processProtare run. If enabled, code checks for Merced output files from previous \nruns and if found reads them instead of rerunning Merced. Only impacts "-mg" option processing.' )

parser.add_argument( '--preProcessLabel', type = str, default = None,                                   help = 'Label for style to process. If None, proc will pick the latest from the evaluated, and Realization styles.' )

parser.add_argument( '--printArgsFile', type = str, default = None,                                     help = 'File to which input arguments are written.')

parser.add_argument( '--prefixRecon', type = str, default = RCprefixDefault,                            help = 'Prefix for Resonance Reconstruction styles. Default is "%s".' % RCprefixDefault )
parser.add_argument( '--prefixMuCutoff', type = str, default = muCutoffPrefixDefault,                   help = 'Prefix for Coulomb + nuclear elastic scattering mu cutoff style. Default is "%s".' % muCutoffPrefixDefault )
parser.add_argument( '--prefixAEP', type = str, default = AEPprefixDefault,                             help = 'Prefix for Average Energy to Product styles. Default is "%s".' % AEPprefixDefault )
parser.add_argument( '--prefixHeated', type = str, default = HTprefixDefault,                           help = 'Prefix for heated styles. Default is "%s".' % HTprefixDefault )
parser.add_argument( '--prefixMultiGroup', type = str, default = MGprefixDefault,                       help = 'Prefix for MultiGroup styles. Default is "%s".' % MGprefixDefault )
parser.add_argument( '--prefixMonteCarlo', type = str, default = MCprefixDefault,                       help = 'Prefix for MonteCarlo styles. Default is "%s".' % MCprefixDefault )
parser.add_argument( '--prefixUpScatter', type = str, default = UPprefixDefault,                        help = 'Prefix for UpScatter styles. Default is "%s".' % UPprefixDefault )

args = parser.parse_args( )

# TODO Find a better way to write the command line arguments to file ... the current approach is specific for processProtare.py 
if args.printArgsFile is not None:
    parserActions = dict([(x.dest, x) for x in parser._actions])
    currentArguments = []
    excludeArguments = ['mapOrProtareFileName', 'outputFile', 'printArgsFile']
    for arg in vars(args):
        if arg in excludeArguments:
            continue

        if len(parserActions[arg].option_strings) == 0:
            currentArguments.append(getattr(args, arg))
        else:
            argValue = getattr(args, arg)
            if isinstance(parserActions[arg].const, bool):
                if parserActions[arg].default != argValue:
                    currentArguments.append(parserActions[arg].option_strings[0])
            elif isinstance(argValue, list):
                if arg == 'gid':
                    # for repeated particle definitions ... us the last one
                    gidIndices = {}
                    i = 0
                    for pid_gid in args.gid :
                        pid, gid = pid_gid.split( '=' )
                        gidIndices[pid] = i
                        i += 1

                    for argIndex in sorted(gidIndices.values()):
                        currentArguments.append('%s %s' % (parserActions[arg].option_strings[0], argValue[argIndex]))

                else:
                    for singleArg in argValue:
                        currentArguments.append('%s %s' % (parserActions[arg].option_strings[0], singleArg))
            elif arg == 'verbose':
                if argValue > 0:
                    currentArguments.append(' '.join([parserActions[arg].option_strings[0]] * argValue))
            else:
                if argValue is not None:
                    currentArguments.append('%s %s' % (parserActions[arg].option_strings[0], argValue))

    with open(args.printArgsFile, 'w') as fileObject:
        fileObject.write('\n'.join(currentArguments))

unitMap = { 'MeV' : args.energyUnit, 'eV' : args.energyUnit }

logFile = open( 'logFile.%s' % args.tag, 'w' )

# Write properties for files used in execution
def listFileProperties(_fileName):
    shellCommand = ['ls', '-l', _fileName]
    result = subprocess.run(shellCommand, stdout=subprocess.PIPE)
    return result.stdout.decode('utf-8')

logFile.write('File listing for various files used in execution:\n')
logFile.write('GNDS file: %s\n' % listFileProperties(args.mapOrProtareFileName))
logFile.write('processProtare file: %s\n' % listFileProperties(os.path.realpath(__file__)))
logFile.write('Merced executable: %s\n\n' % listFileProperties( locateBinFolderModule.locateMerced( exceptionIfNotFound=True ) ))

MonteCarlo_cdf = args.prefixMonteCarlo + '_cdf'

if( args.MultiGroup ) :

    fluxes = fluxModule.Fluxes.readXML_file( args.fluxFile )
    flux = fluxes[args.fluxID]
    flux.convertUnits( unitMap )

    gids = {}
    for pid_gid in args.gid :
        pid, gid = pid_gid.split( '=' )
        gids[pid] = gid

    transportables = []
    groups = groupModule.Groups.readXML_file( args.groupFile )
    for pid in gids :
        multiGroup = groups[gids[pid]]
        multiGroup.convertUnits( unitMap )
        transportables.append( transportablesModule.Transportable( pid, enumsModule.Conserve.number, multiGroup ) )
    if( len( transportables ) == 0 ) : raise UserWarning( 'No --gid specified, must have at least 1.' )

reactionSuite = singleProtareArguments.protare( args )
isThermalNeutronScatteringLaw = reactionSuite.isThermalNeutronScatteringLaw( )

if( ( reactionSuite.projectile != IDsPoPsModule.neutron ) or ( reactionSuite.interaction == enumsModule.Interaction.LLNL_TNSL ) ) :
    args.UpScatter = False
    if( args.temperatures is not None ) : raise ValueError( 'Can only heat neutron as projectile cross sections.' )
if( args.temperatures is None ) :
    args.temperatures = [ reactionSuite.styles[0].temperature.getValueAs( args.temperatureUnit ) ]
    if( isThermalNeutronScatteringLaw ) :
        thermalNeutronScatteringLawTemperatures = reactionSuite.thermalNeutronScatteringLawTemperatures( )
        key = 'incoherent-inelastic'
        if( key not in thermalNeutronScatteringLawTemperatures ) : key = list( thermalNeutronScatteringLawTemperatures.keys( ) )[0]
        unit, temperatures = reactionSuite.thermalNeutronScatteringLawTemperatures( )[key]
        args.temperatures = [ PQUModule.PQU( temperature, unit ).getValueAs( args.temperatureUnit ) for temperature in temperatures ]
args.temperatures.sort( )

processedStyles = ( stylesModule.Heated, stylesModule.AverageProductData, stylesModule.MultiGroup, stylesModule.HeatedMultiGroup, 
        stylesModule.MonteCarlo_cdf, stylesModule.SnElasticUpScatter, stylesModule.GriddedCrossSection, stylesModule.URR_probabilityTables )

if( args.cullProcessedData ) :
    stylesToRemove = []
    for style in reactionSuite.styles :
        if( isinstance( style, processedStyles ) ) :
            stylesToRemove.append( style.label )
    reactionSuite.removeStyles( stylesToRemove )

for style in reactionSuite.styles :                 # Fail on detection of existing processed data.
    if( isinstance( style, processedStyles ) ) : 
        raise Exception( 'File already contains processed data. Please use a clean file!' )

reactionSuite.convertUnits( unitMap )
if( args.writeConvertedUnits ) :                    # Convert units if necessary.
    convertedUnitsFileName = args.outputFile
    if( convertedUnitsFileName is None ) : convertedUnitsFileName = reactionSuite.sourcePath
    convertedUnitsFileName = convertedUnitsFileName + '.convertedUnits'
    if( outputDefault != args.output ) : convertedUnitsFileName = os.path.join( args.output, os.path.basename( convertedUnitsFileName ) )
    reactionSuite.saveToFile( convertedUnitsFileName, xs_pdf_cdf1d_singleLine = True, formatVersion = args.formatVersion )

if( args.exit < 2 ) : sys.exit( 0 )

try :
    if( args.preProcessLabel is None ) :
        preProcessingChains = reactionSuite.styles.preProcessingChains( ends = True )
        if( len( preProcessingChains ) != 1 ) : raise Exception( 'Must have only 1 preProcessing chain, have %s' % len( preProcessingChains ) )
        args.preProcessLabel = preProcessingChains[0][0].label
    preLoopStyle = reactionSuite.styles[args.preProcessLabel]
except :
    raise Exception( 'Cannot find the requested processing label "%s"!' % ( args.preProcessLabel ) )

logFile.write( 'Pre resonance reconstruction style label "%s".\n' % preLoopStyle.label )
if( reactionSuite.supportsResonanceReconstruction( ) ) :            # Reconstruct resonances. FIXME: add logic for choice of reconstruction type (reconstructed angular vs. crossSection).
    if( args.verbose > 0 ) : print( 'Reconstructing resonances' )
    reconStyle = reactionSuite.styles.getStyleOfClass( stylesModule.CrossSectionReconstructed )
    if( reconStyle is None ) : 
        reconStyle = stylesModule.CrossSectionReconstructed( args.prefixRecon, preLoopStyle.label, date = dateTimeStr )
        reactionSuite.reconstructResonances(reconStyle, args.reconAccuracy)                     # FIXME need logfile arguments in this process call
        preLoopStyle = reconStyle

additionalReactions = []
if( isinstance( reactionSuite.PoPs[reactionSuite.projectile], nuclidePoPsModule.Particle ) ) :  # Check if Coulomb + nuclear mu cutoff stuff is needed.
    CoulombPlusNuclearMuCutoffs = reactionSuite.CoulombPlusNuclearMuCutoffs( )
    if( CoulombPlusNuclearMuCutoffs is not None ) :
        if( args.CoulombPlusNuclearMuCutOff not in CoulombPlusNuclearMuCutoffs ) :
            if( args.verbose > 0 ) : print( 'Processing CoulombPlusNuclearMuCutoff data' )
            muCutoffStyle = stylesModule.CoulombPlusNuclearElasticMuCutoff( args.prefixMuCutoff, preLoopStyle.label, args.CoulombPlusNuclearMuCutOff, date = dateTimeStr )
            reactionSuite.styles.add( muCutoffStyle )
            nuclearPlusCoulombInterference = reactionSuite.processCoulombPlusNuclearMuCutoff( muCutoffStyle, excludeRutherfordScattering = not( args.doNotAddNuclearPlusInterference ) )
            if( nuclearPlusCoulombInterference is not None ) : additionalReactions.append( nuclearPlusCoulombInterference.reaction )
            preLoopStyle = muCutoffStyle

# Calculate Average Product Data  ( Energy and Momenta )
logFile.write( 'Pre average product data style label "%s".\n' % preLoopStyle.label )
if( not( isThermalNeutronScatteringLaw ) ) :   # TNSL average product data are temperature dependent and must be done in the temperature loop.
    AEPStyle = reactionSuite.styles.getStyleOfClass( stylesModule.AverageProductData )
    if( AEPStyle is None ) :
        if( args.verbose > 0 ) : print( 'Processing average product data' )
        AEPStyle = stylesModule.AverageProductData( args.prefixAEP, preLoopStyle.label, date = dateTimeStr )
        reactionSuite.styles.add( AEPStyle )
        reactionSuite.calculateAverageProductData( AEPStyle, indent = '  ', verbosity = args.verbose - 2,
                additionalReactions = additionalReactions )  # FIXME need logfile arguments in this process call
    preLoopStyle = AEPStyle

    if( args.MonteCarlo ) :
        if( args.verbose > 0 ) : print( 'Processing Monte Carlo' )
        MonteCarloStyle = stylesModule.MonteCarlo_cdf( '%s' % MonteCarlo_cdf, AEPStyle.label, date = dateTimeStr )
        reactionSuite.styles.add( MonteCarloStyle )
        reactionSuite.processMC_cdf( MonteCarloStyle, indent = '  ', verbosity = args.verbose - 2, additionalReactions = additionalReactions )
        preLoopStyle = MonteCarloStyle

logFile.write( 'Initial work (before temperatures loop):\n%s\n\n' % timer )

logFile.write( 'Pre heating loop style label "%s".\n' % preLoopStyle.label )
transportables_href = None
for temperatureIndex, temperatureValue in enumerate( args.temperatures ) : # Heat the cross sections and AEPs with a style for each temperature.

    if( args.verbose > 0 ) : print( 'Processing for temperature %s %s (%s of %s)' % ( temperatureValue, args.temperatureUnit, temperatureIndex + 1, len( args.temperatures ) ) )

    suffix = '_%03d' % temperatureIndex

    temperature = physicalQuantityModule.Temperature( temperatureValue, args.temperatureUnit )
    logFile.write( 'Heating to %s\n' % temperature )

    timerHeating = timesModule.Times( )
    heatStyle = stylesModule.Heated( args.prefixHeated + suffix, preLoopStyle.label, temperature, date = dateTimeStr )
    reactionSuite.styles.add( heatStyle )
    if( isThermalNeutronScatteringLaw ) :
        if( args.verbose > 0 ) : print( '  Processing TNSL' )
        reactionSuite.processThermalNeutronScatteringLaw( heatStyle, indent = '  ', verbosity = args.verbose - 2 )
    else :
        if( args.verbose > 0 ) : print( '  Heating cross sections' )
        reactionSuite.heatCrossSections( heatStyle, setThresholdToZero = True, heatBelowThreshold = False, verbose = args.verbose - 2 )  # FIXME need logfile arguments in this process call
    logFile.write( '    %s\n' % timerHeating )

    loopStyle = heatStyle

    if( args.MonteCarlo ) :
        MonteCarloStyle = loopStyle
        if( isThermalNeutronScatteringLaw ) :
            if( args.verbose > 1 ) : print( '  Processing Monte Carlo' )
            MonteCarloStyle = stylesModule.MonteCarlo_cdf( '%s' % ( MonteCarlo_cdf + suffix ), loopStyle.label, date = dateTimeStr )
            reactionSuite.styles.add( MonteCarloStyle )
            reactionSuite.processMC_cdf( MonteCarloStyle, indent = '    ', verbosity = args.verbose - 2 )

        if( args.verbose > 1 ) : print( '  Processing gridded cross sections' )
        timerMC = timesModule.Times( )
        griddedCrossSectionStyle = stylesModule.GriddedCrossSection( args.prefixMonteCarlo + suffix, MonteCarloStyle.label, date = dateTimeStr )
        reactionSuite.styles.add( griddedCrossSectionStyle )
        reactionSuite.processGriddedCrossSections( griddedCrossSectionStyle, verbosity = args.verbose - 2, indent = '    ',
                additionalReactions = additionalReactions )

        logFile.write( '  Processing Monte Carlo\n    %s\n' % timerMC )

    if( args.MultiGroup ) :
        if( args.verbose > 1 ) : print( '  Processing multi-group' )
        logFile.write( '  Processing multi-group\n' )
        fluxAtTemperature = fluxModule.Flux(args.fluxID, flux.evaluate(temperatureValue, extrapolation=xDataEnumsModule.Extrapolation.flat))
        fluxAtTemperature.data.outerDomainValue = None
        fluxAtTemperature.data.axes.axes.pop( -1 )
        heatedMultiGroupStyle = stylesModule.HeatedMultiGroup( args.prefixMultiGroup + suffix, loopStyle.label, fluxAtTemperature, date = dateTimeStr )
        if( transportables_href is None ) :
            for transportable in transportables : heatedMultiGroupStyle.transportables.add( transportable )
        else :
            heatedMultiGroupStyle.transportables.set_href( transportables_href )
        reactionSuite.styles.add( heatedMultiGroupStyle )
        if( transportables_href is None ) : transportables_href = heatedMultiGroupStyle.transportables.toXLink( )
        workDir = os.path.join( 'Merced.work', reactionSuite.projectile, reactionSuite.target, 'T_%s_%s' % 
                ( temperatureValue, args.temperatureUnit.replace( ' ', '' ).replace( os.sep, '_' ) ) )
        if( os.path.exists( workDir ) and not args.restart ) : shutil.rmtree( workDir )
        reactionSuite.processMultiGroup( heatedMultiGroupStyle, args.legendreMax, verbosity = args.verbose - 2, logFile = logFile, indent = '    ', 
                additionalReactions = additionalReactions, workDir = workDir, restart = args.restart )
        if( args.UpScatter and not( isThermalNeutronScatteringLaw or ( reactionSuite.target in reactionSuite.PoPs.unorthodoxes ) ) ) :
            if( args.verbose > 1 ) : print( '  Processing multi group upscattering' )
            timerUp = timesModule.Times( )
            heatUpscatterStyle = stylesModule.SnElasticUpScatter( args.prefixUpScatter + suffix, heatedMultiGroupStyle.label, date = dateTimeStr ) 
            reactionSuite.styles.add( heatUpscatterStyle )
            if( temperatureValue > 0.0 ) :
                reactionSuite.processSnElasticUpScatter( heatUpscatterStyle, args.legendreMax, verbosity = args.verbose - 2 )  # FIXME need logfile arguments in this process call
            logFile.write( '  Upscatter correction\n    %s\n' % timerUp )

processedFileName = args.outputFile
if( processedFileName is None ) :
    filetype = reactionSuite.sourcePath.split( '.' )[-1]
    processedFileName = reactionSuite.sourcePath[:-len(filetype)] + '%s.%s' % ( args.tag, filetype )
if( outputDefault != args.output ) : processedFileName = os.path.join( args.output, os.path.basename( processedFileName ) )
reactionSuite.saveToFile( processedFileName, xs_pdf_cdf1d_singleLine = True, formatVersion = args.formatVersion )

logFile.write( '\nTotal elapsed time:\n  %s\n' % timer )
logFile.close( )
