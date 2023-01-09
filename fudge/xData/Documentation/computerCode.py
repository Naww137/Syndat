# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node computerCodes and its child nodes classes.
"""

import datetime

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule
from .. import date as dateModule

class Version( textModule.Text ) :

    moniker = 'version'

    def __init__( self, date=None ) :

        textModule.Text.__init__( self, markup = textModule.Markup.none )

        self.__date = dateModule.Date(date)

    @property
    def date( self ) :

        return( self.__date )

    @date.setter
    def date( self, date ) :
        """Set version date."""

        self.__date = dateModule.raiseIfNotDate( date )

    def XML_extraAttributes( self, **kwargs ) :

        if not kwargs['addExtraAttributes'] : return ''
        return ' date="%s"' % self.date

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        Parses a version node.
        """

        textModule.Text.parseNode(self, node, xPath, linkData, **kwargs)
        xPath.append( node.tag )

        self.date = datetime.datetime.strptime( node.get( 'date' ), '%Y-%m-%d' )

        xPath.pop( )

class CodeRepo(ancestryModule.AncestryIO_bare):

    moniker = 'codeRepo'

    def __init__(self, _revisionSystem, _href, _revisionID):

        ancestryModule.AncestryIO_bare.__init__(self)

        self.__revisionSystem = textModule.raiseIfNotString( _revisionSystem, 'revisionSystem' )
        self.__href = textModule.raiseIfNotString( _href, 'href' )
        self.__revisionID = textModule.raiseIfNotString( _revisionID, 'revisionID' )

    @property
    def revisionSystem( self ) :

        return( self.__revisionSystem )

    @property
    def href( self ) :

        return( self.__href )

    @property
    def revisionID( self ) :

        return( self.__revisionID )

    def toXML_strList( self, indent = '', **kwargs ) :

        if( len( self.__revisionSystem + self.__href + self.__revisionID ) == 0 ) : return( [] )

        return( [ '%s<%s revisionSystem="%s" href="%s" revisionID="%s"/>' % ( indent, self.moniker, self.__revisionSystem, self.__href, self.__revisionID ) ] )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append( node.tag )

        self.__revisionSystem = textModule.raiseIfNotString( node.get( 'revisionSystem' ), 'revisionSystem' )
        self.__href = textModule.raiseIfNotString( node.get( 'href' ), 'href' )
        self.__revisionID = textModule.raiseIfNotString( node.get( 'revisionID' ), 'revisionID' )

        xPath.pop( )

class ExecutionArguments( textModule.Text ) :

    moniker = 'executionArguments'

    def __init__( self ) :

        textModule.Text.__init__( self, markup = textModule.Markup.none )

class Note( textModule.Text ) :

    moniker = 'note'

class InputDeck( textModule.Text ) :

    moniker = 'inputDeck'
    keyName = 'label'

    def __init__( self, label, filename, text = None ) :

        textModule.Text.__init__( self, text, markup = textModule.Markup.none, label = label )

        self.__filename = textModule.raiseIfNotString( filename, 'filename' )

    @property
    def filename( self ) :

        return( self.__filename )

    def XML_extraAttributes( self, **kwargs ) :

        if( self.filename == '' ) : return ''

        return ' filename="%s"' % self.filename

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        _label = node.get( 'label' )
        _filename = node.get('filename', '')

        inputDeck = cls(_label, _filename)
        inputDeck.parseNode(node, xPath, linkData, **kwargs)

        return inputDeck

class InputDecks( suiteModule.Suite ) :

    moniker = 'inputDecks'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ InputDeck ] )

class OutputDeck( textModule.Text ) :

    moniker = 'outputDeck'
    keyName = 'label'

    def __init__( self, filename, text = None ) :

        textModule.Text.__init__( self, text, markup = textModule.Markup.none )

        self.__filename = textModule.raiseIfNotString( filename, 'filename' )

    @property
    def filename( self ) :

        return( self.__filename )

    def XML_extraAttributes( self, **kwargs ) :

        if( self.filename == '' ) : return ''

        return ' filename="%s"' % self.filename

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        _filename = node.get('filename', '')

        return cls( _filename)

class OutputDecks( suiteModule.Suite ) :

    moniker = 'outputDecks'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ OutputDeck ] )

class ComputerCode(ancestryModule.AncestryIO):
    """A class representing a GNDS documentation/computerCodes/computerCode node."""

# Also need buildParameters as a text node.

    moniker = 'computerCode'
    keyName = 'label'
    ancestryMembers = ( 'codeRepo', 'executionArguments', 'note', 'inputDecks', 'outputDecks' )

    def __init__( self, label, name, version ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__label = textModule.raiseIfNotString( label, 'label' )
        self.__name = textModule.raiseIfNotString( name, 'name' )
        self.__version = textModule.raiseIfNotString( version, 'version' )

        self.codeRepo = CodeRepo( '', '', '' )

        self.__executionArguments = ExecutionArguments( )
        self.__executionArguments.setAncestor( self )

        self.__note = Note( )
        self.__note.setAncestor( self )

        self.__inputDecks = InputDecks( )
        self.__inputDecks.setAncestor( self )

        self.__outputDecks = OutputDecks( )
        self.__outputDecks.setAncestor( self )

    @property
    def label( self ) :
        """Returns the label instance."""

        return( self.__label )

    @property
    def name( self ) :
        """Returns the name instance."""

        return( self.__name )

    @property
    def version( self ) :
        """Returns the versioninstance."""

        return( self.__version )

    @property
    def codeRepo( self ) :

        return( self.__codeRepo )

    @codeRepo.setter
    def codeRepo( self, _codeRepo ) :

        if( not( isinstance( _codeRepo, CodeRepo ) ) ) : raise TypeError( 'Invalid codeRepo instance.' )
        self.__codeRepo = _codeRepo
        self.__codeRepo.setAncestor( self )

    @property
    def executionArguments( self ) :

        return( self.__executionArguments )

    @property
    def note( self ) :

        return( self.__note )

    @property
    def inputDecks( self ) :

        return( self.__inputDecks )

    @property
    def outputDecks( self ) :

        return( self.__outputDecks )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s label="%s" name="%s" version="%s">' % ( indent, self.moniker, self.__label, self.__name, self.__version ) ]
        XMLList += self.__codeRepo.toXML_strList( indent2, **kwargs )

        XMLList += self.__executionArguments.toXML_strList( indent2, **kwargs )
        XMLList += self.__note.toXML_strList( indent2, **kwargs )
        XMLList += self.__inputDecks.toXML_strList( indent2, **kwargs )
        XMLList += self.__outputDecks.toXML_strList( indent2, **kwargs )

        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        label = node.get( 'label' )
        name = node.get( 'name' )
        version = node.get( 'version', '' )

        computerCode = cls(label, name, version)
        computerCode.parseAncestryMembers(node, xPath, linkData, **kwargs)

        return computerCode

class ComputerCodes( suiteModule.Suite ) :

    moniker = 'computerCodes'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ ComputerCode ] )
