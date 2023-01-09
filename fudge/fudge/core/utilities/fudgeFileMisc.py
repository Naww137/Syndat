# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class FudgeTempFile which simplies the functions in
the tempfile module.
"""

import os
import tempfile
from sys import version_info


class FudgeTempFile :
    """This class creates a temporary file using tempfile.mkstemp and supports common functionallity
    for the file (e.g., write). Currently, reading from the temporary file is not supported."""

    def __init__( self, prefix = None, suffix = "", dir = None, deleteOnClose = False ) :
        """Contructor for the FudgeTempFile class. Calls tempfile.mkstemp to open a temporary file.
        If deleteOnClose is 'True', file will be deleted when it is close."""

        if( dir is None ) :
            dir = tempfile.gettempdir( )
        else :
            if( not os.path.exists( dir ) ) : os.makedirs( dir )
        if( prefix is None ) : prefix = tempfile.gettempprefix( )

        if version_info.major == 3:
            self.fd, self.name = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir, text=True)
        else:
            self.fd, self.name = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir)
        self.deleted = False
        self.deleteOnClose = deleteOnClose

    def __del__( self ) :
        """If class instance is deleted, the file is properly closed."""

        if( self.fd is not None ) : self.close( )

    def close( self, raiseIfClosed = True ) :
        """Closes the file if still opened. If raiseIfClosed is 'True' and file is already
        closed, a raise is executed."""

        if( self.fd is not None ) :
            os.close( self.fd )
            self.fd = None
            if( self.deleteOnClose ) : self.delete( )
        elif( raiseIfClosed ) :
            raise Exception( 'Error from FudgeTempFile.close: file already closed' )

    def delete( self ) :
        """Deletes file if it still exist.  If required, this method calls close first. If file has already
        been deleted, a raise is executed."""

        if( self.deleted ) :
            raise Exception( 'Error from FudgeTempFile.delete: file already deleted' )
        else :
            if( self.fd is not None ) : self.close( )
            os.remove( self.name )
            self.deleted = True

    def getName( self ) :
        """Returns self's name which is the full path name of the temporary file."""

        return( self.name )

    def getFileDescriptor( self ) :
        """Returns self's file descriptor."""

        return( self.fd )

    def isOpened( self ) :
        """Returns 'True' if the file is still opened."""

        return( self.fd is not None )

    def isDeleted( self ) :
        """Returns 'True' if the file has been deleted."""

        return( self.deleted )

    def write( self, str ) :
        """Write str to file. If file is closed or not all characters were written then a raise is executed."""

        if( self.fd is not None ) :
            if version_info.major == 3:
                n = os.write(self.fd, str.encode())
            else:
                n = os.write( self.fd, str )
            if( n != len( str ) ) : raise Exception( 'Error from FudgeTempFile.write: only %d of %d characters written' % ( n, len( str ) ) )
        else :
            raise Exception( 'Error from FudgeTempFile.write: attempted to write to a closed file.' )
