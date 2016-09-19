package Motif_db;
################################################################################
# $Id: Motif_db.pm,v 1.1.2.7 2012/12/20 14:16:36 gribskov Exp $
#
# Class Motif_db
#
# Manages a database of XIOS graphs represented as minimum DFS codes.  Each 
# motif is indexed by a serialized, packed version of the minimum DFS.  As long 
# as the number of stems is seven or less, each row of the DFS can be backed in
# a single byte.
#
# Encoding
#   bit 0 - 1   edge type, i=0, j=1, o=2, x=3
#   bit 2 - 4   to vertex number (0-7)
#   bit 5 - 8   from vertex number (0-7)
#
# TODO need to make functions callable as either oo or serial
#
# usage
#   create a new database from serialized DFS code
#   my $motifdb = new Motif_db( {serial=>'2_to_7_stems_topologies.txt.with_label' } );
#
# Revision log at end of file
################################################################################
use strict;
use Data::Dumper;
use Storable;

use vars qw( @ISA @EXPORT @EXPORT_OK $VERSION );
@ISA       = qw( Exporter );
@EXPORT    = qw( );
@EXPORT_OK = qw( );

my $REVISION = '$Revision: 1.1.2.7 $';
( $VERSION ) = $REVISION =~ /Revision: ([\d.]*) /;

#-------------------------------------------------------------------------------
# new
#
# Motif_db constructor.  Class attributes:
#   db  => ref to hash.  keys are packed dfs, values are motif name as 
#          nstem_motifnumber
#    n  => number of motifs available
#
# Usage
#-------------------------------------------------------------------------------
sub new{
	my ( $class, $param ) = @_;
	
	my $motif = {};
	bless $motif, $class;
	
	$motif->{db} = {};
	$motif->{n}  = 0;
	
	if ( defined $$param{serial} ) {
	    my $nmotif = $motif->loadSerialDfs( $$param{serial} );
	    $motif->n( $nmotif );
	}
	
	if ( defined $$param{db} ) {
	    my $nmotif = $motif->retrieveDb( $$param{db} );
	    $motif->n( $nmotif );
	}
    
	return $motif;    
}

# end of new

#------------------------------------------------------------------------------
# loadSerialDfs
#
# Load motif data in serialized text form:
#2_stem_motif_1: 0 1 2,
#2_stem_motif_2: 0 1 0,
#3_stem_motif_1: 0 1 0, 0 2 0,
#3_stem_motif_2: 0 1 0, 1 2 2, 2 0 2,
#
# USAGE
#    $n_motif_loaded = $motifdb->loadSerialDfs( $filename );
#------------------------------------------------------------------------------
sub loadSerialDfs{
	my ( $motif, $filename ) = @_;
    open( my $dfs, "<", $filename ) or 
        die "Motif_db::loadSerialDfs - Unable to open DFS file ($filename)\n\n";
        
    my %db;
    my $n_motif = 0;
    while ( my $line = <$dfs> ) {
        next if ( $line =~ /^#|^\s*$/ );    # skip comments and blank lines
        chomp $line;
        
        my ( $name, $dfs ) = split /\s*:\s*/, $line, 2;
        $name =~ s/_stem_motif//i;
        
        my $index;
        foreach my $dfsrow ( split /\s*,\s*/, $dfs ) {
            my $bytecode = encodeDfsRow( $dfsrow );
            $index .= sprintf "%02x", $bytecode;
        }
        $db{$index} = $name;
        $n_motif++;
        
    }
    $motif->db( \%db );
	close $dfs;
	
	return $n_motif;
}

# End of loadSerialDfs

#------------------------------------------------------------------------------
# lookupMotif
#
# Return the name of the motif based on the index
#
# USAGE
#    $motif_name = $motif->lookupMotif( $index );
#------------------------------------------------------------------------------
sub lookupMotif{
	my ( $motif, $index ) = @_;
	
	return $motif->{db}->{$index};
}

# End of lookupMotif

#------------------------------------------------------------------------------
# storeDb
#
# Store the current database using Storable.pm.  The resulting file can be read
# with $db->retrieveDb( $filename )
#
# USAGE
#    $n_stored = $motif->storeDb( $filename );
#------------------------------------------------------------------------------
sub storeDb{
    my ( $motif, $filename ) = @_;
    
    print STDERR "Motif_db::storeDb - storing to $filename\n\n";
    store( $motif, $filename );
    
    return $motif->n;
}

# End of storeDb

#------------------------------------------------------------------------------
# retrieveDb  
#
# Read in a database stored using Storable.pm.  Normally this file will be
# written with the function $motif->storeDb( $filename).  Note that this will 
# completely overwrite anything currently in the database.
#
# USAGE
#    $n_graphs = $motif->retrieveDb( $filename );
#------------------------------------------------------------------------------
sub retrieveDb{
    my ( $motif, $filename ) = @_;
    
    #print "Motif_db::retrieveDB - filename:$filename:\n";
    #$motif = retrieve( $filename );
    $_[0] = retrieve( $filename ); # need to change the global $motif, not the local copy
    #print Dumper( $_[0] );
    
    return $_[0]->n;
}

# End of retrieveDb

#------------------------------------------------------------------------------
# encodeDfsRow
#
# Encode a DFS row in eight bits.  Can only be used with graphs that have seven 
# or fewer vertices.
#
# Encoding
#   bit 0 - 1   edge type, i=0, j=1, o=2, x=3
#   bit 2 - 4   to vertex number (0-7)
#   bit 5 - 8   from vertex number (0-7)
#
# USAGE
#    $bytecode = encodeDfsRow( $row_text );
#------------------------------------------------------------------------------
sub encodeDfsRow{
    my ( $dfs ) = @_;
    my $byte = 0;
    
    my ( $v1, $v2, $e ) = split " ", $dfs;
    $byte += $v1 << 5;
    $byte += $v2 << 2;
    $byte += $e;
    return $byte;
}

# End of encodeDfsRow

#------------------------------------------------------------------------------
# encodeDfsRowArray
#
# Encode and look up a DFS code represented as an array or DFS rows, where each row is an 
# array of (v1,v2,e).  Each row is encoded  in eight bits, and the entire DFS 
# encoded as a hexadecimanl string.  Can only be used with graphs that have seven 
# or fewer vertices.
#
# Encoding
#   bit 0 - 1   edge type, i=0, j=1, o=2, x=3
#   bit 2 - 4   to vertex number (0-7)
#   bit 5 - 8   from vertex number (0-7)
#
# USAGE
#    $hex_string = $motif->encodeDfsRowArray( @dfs );
#------------------------------------------------------------------------------
sub encodeDfsRowArray{
    my ( $motif, @dfs ) = @_;
    my $hexstring = "";
    
    foreach my $dfsrow ( @dfs ) {
        my $byte;
        my ( $v1, $v2, $e ) = @{$dfsrow};
        $byte += $v1 << 5;
        $byte += $v2 << 2;
        $byte += $e;
        
        $hexstring .= sprintf "%02x", $byte;
    }
    
    return $hexstring;
}

# End of encodeDfsRowArray

#------------------------------------------------------------------------------
# decodeDfs
#
# Convert a hexadecimal string back to the original DFS code.
#
# USAGE
#    $dfs = decodeDfs( $hex );
#------------------------------------------------------------------------------
sub decodeDfs{
    if ( ref $_[0] eq 'Motif_db' ) {
        shift @_;  
    }
    my ( $hex ) = @_;
    my $dfs = "";
    
    my @row = split /(\w{2})/, $hex;
    while ( @row ) {
        my $value = shift @row;
        next unless $value;
        
        my $dec = hex $value;
        my $v1 = ($dec & 224) >> 5;
        my $v2 = ($dec &  28) >> 2;
        my $e  = ($dec &   3) ;
        
        $dfs .= sprintf "%2d  %2d  %2d\n", $v1, $v2, $e;
    }
    
    return $dfs ;
}

# End of decodeDfs

#-------------------------------------------------------------------------------
# AUTOLOAD
#
# Returns the contents of the class data structure indexed by 'hashkey'
# This function will not create new hash keys.  This keeps typos from creating
# crazy new attributes for objects.  Anything that is not a hash key or 
# a function in this package fails.
#
# Usage
#   $info = $db->{hashkey};
#-------------------------------------------------------------------------------
sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;
  
    my $unknown = $AUTOLOAD;                # The name of the unknown function
    $unknown =~ s/.*:://;                   # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;    # Skip all-caps methods like DESTROY
  
    if ( exists $self->{ $unknown } ) {     # do not create new fields
        $self->{ $unknown } = shift if @_;  # Set new value if one was supplied
        return( $self->{ $unknown } );      # Return current value
    } else {
        warn "Unknown function ($unknown)\n";
    }
        
    return;        
}

# end of AUTOLOAD

################################################################################
# $Log: Motif_db.pm,v $
# Revision 1.1.2.7  2012/12/20 14:16:36  gribskov
# made decodeDfs callable as either oo or serial.
#
# Revision 1.1.2.6  2012/12/20 13:56:00  gribskov
# works with fingerprint_random.pl 1.1.2.4 2012/12/20
#
# Revision 1.1.2.5  2012/12/19 19:54:23  gribskov
# includes motif lookup.  basic strategy works but the hex indices don't match.
#
# Revision 1.1.2.4  2012/12/19 14:46:56  gribskov
# Tested store and retrieve from binary database using Storable.pm.
#
# Revision 1.1.2.3  2012/12/19 14:21:07  gribskov
# Added storeDB and retrieveDB (modified functions from NH_database.pm).
#
# Revision 1.1.2.2  2012/12/19 14:16:07  gribskov
# Successful load of db from serialized file.
#
# Revision 1.1.2.1  2012/12/19 13:47:36  gribskov
# Initial version, builds dictionary of motifs from serialized minimum DFS codes.
#
################################################################################    
1;
