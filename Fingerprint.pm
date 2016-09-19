package Fingerprint;
################################################################################
# $Id: Fingerprint.pm,v 1.5.4.3 2015/10/12 17:38:27 huang147 Exp $
#
# Class Fingerprint
#
# For reading, writing, and manipulating RNA fingerprints.  See end of file for
# sample and compact format.
#
# Revision log at end of file
################################################################################
use strict;
use Data::Dumper;
use XML::Simple;
use XML::LibXML;

use Exporter;
use vars qw( @ISA @EXPORT @EXPORT_OK $VERSION );
@ISA       = qw( Exporter );
@EXPORT    = qw( );
@EXPORT_OK = qw( );

my $REVISION = '$Revision: 1.5.4.3 $';
( $VERSION ) = $REVISION =~ /Revision: ([\d.]*) /;

my $INDENT         =  4;   # number of spaces per indentation level
my $MOTIF_PER_LINE = 20;   # number of motifs per line in a motif_block

#-------------------------------------------------------------------------------
# new
#
# Fingerprint constructor.
#
# Usage
#-------------------------------------------------------------------------------
sub new{
	my ( $class ) = @_;
	
	my $self = {};
	bless $self, $class;
	
	# query info
	$self->{query_id}       = "";	
    $self->{query_vertex}   = undef;
    $self->{query_edge}     = undef;
    
    # fingerprint info
    $self->{type}          = undef;      # exhaustive | random
    $self->{iteration}     = undef;      # only for random fingerprint
    $self->{total_mapping} = undef;      # total number of mappings found
    $self->{program}       = undef;      # program used to calculate fingerprint
    $self->{time_elapsed}  = undef;      # calculation time (sec)
    
    # motif info
	$self->{database_id}  = "";
	$self->{motif_list}   = [];
    
	return $self;    
}

# end of new

#------------------------------------------------------------------------------
# addMotif
#
# Add a motif to the motif list.  The motif is passed as a hash, the 'id' key
# must be present, no other checking is done.  This means that arbitrary new
# fields can be created in the XML.
#
# USAGE
#    my $nmotif = $fingerprint->addMotif( 
#                           {   id             => $motifdb->lookupMotif( $k ),
#                               count          => $count{$k},
#                               first_observed => $first{$k},
#                               encoded_dfs    => $motifdb->decodeDfs( $k )
#                           }
#                 );
#------------------------------------------------------------------------------
sub addMotif{
	my ( $self, $motif ) = @_;
	
	my $motiflist = $self->motifList;
	my $nmotif = @{$motiflist};
	
	unless ( defined $$motif{id} ) {
	    print STDERR "Fingerprint::addMotif - no ID defined for motif\n";
	    print Dumper($motif), "\n";
	    return $nmotif;
	}
	
	push @{$motiflist}, $motif;
	$nmotif++;
	
	return $nmotif;
}

# End of addMotif

#------------------------------------------------------------------------------
# asCompact
#
# Return a string with the compact format, or if a filename is provided, write
# to a file.
#
# USAGE
#   $str = $fingerprint->asCompact( );
#   $str = $fingerprint->asCompact( $filename );
#------------------------------------------------------------------------------
sub asCompact{
	my ( $self, $filename ) = @_;
	my $str;
	
    my $n_motif          = @{$self->motifList};
    my $query_id         = $self->queryId;
    my $query_vertex     = $self->queryVertex;
    my $query_edge       = $self->queryEdge;
    my $type             = $self->type;
    my $iteration        = $self->iteration;
    my $time             = $self->timeElapsed;
    my $total_mapping    = $self->total_mapping;
    my $program          = $self->program;
    my $database_id      = $self->databaseId;

	$str .= stringIf( 'query_id',      $query_id      );
	$str .= stringIf( 'query_vertex',  $query_vertex  );
    $str .= stringIf( 'query_edge',    $query_edge    );
    $str .= stringIf( 'type',          $type          );
    $str .= stringIf( 'iteration',     $iteration     );
    $str .= stringIf( 'time_elapsed',  $time          );
    $str .= stringIf( 'total_mapping', $total_mapping );
    $str .= stringIf( 'program',       $program       );
	$str .= stringIf( 'database_id',   $database_id   );
	$str .= stringIf( 'motif_n',       $n_motif       );
	
	$str .= qq{motif_list:\n};
    my $count = 1;
    foreach my $motif ( @{$self->motifList} ) {
        $str .= qq/$$motif{id}/;
        if ( defined $$motif{count} ) {
            $str .= qq/,c:$$motif{count}/;
        }
        if ( defined $$motif{first_observed} ) {
            $str .= qq/,f:$$motif{first_observed}/;
        }
        if ( defined $$motif{encoded_dfs} ) {
            $str .= qq/,d:$$motif{encoded_dfs}/;
        }
        if ( defined $$motif{mapping} ) {
            $str .= qq/,m:$$motif{mapping}/;
        }
        
        if ( $count % $MOTIF_PER_LINE ) {
            $str .= qq{ };
        } else {
            $str .= qq{\n};
        }
        $count++;
    }  
    if ( ($count-1) % $MOTIF_PER_LINE ) { $str .= qq{\n} };
	$str .= qq{\n};
	$str .= qq{\n};
	$str .= qq{\n};
		
    if ( $filename ) {
        open ( my $file, ">", $filename ) or 
           die "Fingerprint::asCompact unable to open $filename for output\n\n";
        print $file $str;
        close $file;
        
    } 
	return $str;
}

# End of asCompact

#------------------------------------------------------------------------------
# asXML
#
# Return a string with the XML representation of the, or if $filename is 
# provided, write to a file.
#
# USAGE
#   $xml = $fingerprint->asXML;
#   $xml = $fingerprint->asXML( $file );
#------------------------------------------------------------------------------
sub asXML{
	my ( $self, $filename ) = @_;
	my $xml;
	my @istack;
	
	my $n_motif       = @{$self->motifList};
	my $query_id      = $self->queryId;
	my $query_vertex  = $self->queryVertex;
    my $query_edge    = $self->queryEdge;
    my $type          = $self->type;
    my $iteration     = $self->iteration;
    my $total_mapping = $self->totalMapping;
    my $elapsed       = $self->timeElapsed;
    my $program       = $self->program;
	my $database_id   = $self->databaseId;
	
	my $indent = "";
	$xml  = qq{<?xml version="1.0"?>\n};
	$xml .= qq{$indent<XIOS_fingerprint>\n};	
	push @istack, $indent;
    $indent .= " " x $INDENT;
	
	# query info
	$xml .= qq{$indent<query>\n};	
    push @istack, $indent;
    $indent .= " " x $INDENT;	

	$xml .= qq{$indent<query_id>$query_id</query_id>\n};
	$xml .= qq{$indent<query_vertex>$query_vertex</query_vertex>\n};
	$xml .= qq{$indent<query_edge>$query_edge</query_edge>\n};

    $indent = pop @istack;
    $xml .= qq{$indent</query>\n\n};
    
	
	# fingerprint info
	if ( $type ) {
        $xml .= qq{$indent<fingerprint>\n};
        push @istack, $indent;
        $indent .= " " x $INDENT;
        
        $xml .= qq{$indent<type>$type</type>\n};
        if ( $type =~ /rand/i && $iteration ) {
            $xml .= qq{$indent<iteration>$iteration</iteration>\n};
        }
        
        if ( $total_mapping ) {
            $xml .= qq{$indent<total_mapping>$total_mapping</total_mapping>\n};
        }

        if ( $program ) {
            $xml .= qq{$indent<program>$program</program>\n};
        }
                
        if ( $program ) {
            $xml .= qq{$indent<time_elapsed>$elapsed</time_elapsed>\n};
        }
	    $indent = pop @istack;
	    $xml .= qq{$indent</fingerprint>\n};	    
	}
	
    # database info
    $xml .= qq{$indent<database>\n};
    
    push @istack, $indent;
    $indent .= " " x $INDENT;   	
	$xml .= qq{$indent<database_id>$database_id</database_id>\n};
	
	$indent = pop @istack;
	$xml .= qq{$indent</database>\n\n};
	
	# motifs
	$xml .= qq{$indent<motif_list>\n};
	
    push @istack, $indent;
    $indent .= " " x $INDENT;       
	$xml .= qq{$indent<motif_n>$n_motif</motif_n>\n};
#	$xml .= qq{$indent<motif_block>\n};
#
#    push @istack, $indent;
#    $indent .= " " x $INDENT;  
#    my $count = 1;
#    foreach my $motif_id ( @{$self->motifList} ) {
#        if ( $count % $MOTIF_PER_LINE == 1 ) {
#            $xml .= qq{$indent};
#        }
#        $xml .= qq{$motif_id};
#        if ( $count % $MOTIF_PER_LINE ) {
#            $xml .= qq{ };
#        } else {
#            $xml .= qq{\n};
#        }
#        $count++;
#    }  
#    if ( ($count-1) % $MOTIF_PER_LINE ) { $xml .= qq{\n} };
#	
#	$indent = pop @istack;
#	$xml .= qq{$indent<\/motif_block>\n};

    foreach my $motif ( sort { $b->{count}<=>$a->{count} }  @{$self->motifList} ) {
        $xml .= qq{$indent<motif>\n};
    
        push @istack, $indent;
        $indent .= " " x $INDENT;  
        if ( $$motif{id} ) { 
            $xml .= qq{$indent<id>$$motif{id}</id>\n}; }
        if ( $$motif{count} ) { 
            $xml .= qq{$indent<count>$$motif{count}</count>\n}; }
        if ( $$motif{first_observed} ) { 
            $xml .= qq{$indent<first_observed>$$motif{first_observed}</first_observed>\n}; }
        if ( $$motif{encoded_dfs} ) { 
            $xml .= qq{$indent<encoded_dfs>$$motif{encoded_dfs}</encoded_dfs>\n}; }
        if ( $$motif{mapping} ) { 
            $xml .= qq{$indent<mapping>$$motif{mapping}</mapping>\n}; }
        $indent = pop @istack;
        
        $xml .= qq{$indent<\/motif>\n};        
    }
	
	$indent = pop @istack;
	$xml .= qq{$indent</motif_list>\n\n};
	
	$indent = pop @istack;
	$xml .= qq{$indent</XIOS_fingerprint>\n};	
	
	if ( $filename ) {
	    open ( my $file, ">", $filename ) or 
	       die "Fingerprint::asXML unable to open $filename for output\n\n";
        print $file $xml;
        close $file;
	    
	} 
	
	return $xml;
}

# End of asXML

#------------------------------------------------------------------------------
# asXMLSimple
#
# Return a string with the XML representation of the, or if $filename is 
# provided, write to a file.
#
# USAGE
#   $xml = $fingerprint->asXMLSimple;
#   $xml = $fingerprint->asXMLSimple( $file );
#------------------------------------------------------------------------------
sub asXMLSimple{
    my ( $self, $filename ) = @_;
    my $xml;
    
    $xml->{query}->{query_id}        = $self->queryId;
    $xml->{query}->{query_vertex}    = $self->queryVertex;
    $xml->{query}->{query_edge}      = $self->queryEdge;
    $xml->{fingerprint}->{type}      = $self->type;
    $xml->{fingerprint}->{iteration} = $self->iteration;
    $xml->{fingerprint}->{program}   = $self->program;
    $xml->{database}->{database_id}  = $self->databaseId;
    $xml->{database}->{motif_n}      = @{$self->motifList};
    
#    my $count = 1;
#    my $block = "";
#    foreach my $motif_id ( @{$self->motifList} ) {
#        $block .= qq{$motif_id};
#        if ( $count % $MOTIF_PER_LINE ) {
#            $block .= qq{ };
#        } else {
#            $block .= qq{\n};
#        }
#        $count++;
#    }  
#    if ( ($count-1) % $MOTIF_PER_LINE ) { $block .= qq{\n} };
    
    #$xml->{database}->{motif_block} = $block;
    $xml->{motif_list}->{motif} = $self->motifList;
    my $simple = new XML::Simple;
    my $xmlout = $simple->XMLout( $xml, AttrIndent => 1, 
                                        XMLDecl => '<?xml version="1.0"?>',                                        
                                        NoAttr => 1, 
                                        RootName => 'XIOS_fingerprint'   );
           
    if ( $filename ) {
        open ( my $file, ">", $filename ) or 
           die "Fingerprint::asXMLSimple unable to open $filename for output\n\n";
        print $file $xmlout;
        close $file;
        
    } 
    
    return $xmlout;
}

# End of asXMLSimple

#------------------------------------------------------------------------------
# motifFromCompact
#
# The original compact format was simply a list of Motif IDs separated by white
# space.  The extended format supports two formats. the first is simply an 
# ordered comma-separated list
#
#   motif_id[,count[,first_observed[,endocded_dfs]]] ...
#
# with multiple blocks possible per line. The second is a tag value format
#
#   motif_id[,c:count][,f:first_observed][,d:encoded_dfs]
#
# USAGE
#    $ref_motif_list = motifFromCompact( $line );
#------------------------------------------------------------------------------
sub motifFromCompact{
	my ( $line ) = @_;
	my @motif_list;
	
	chomp $line;
	my @token = split " ", $line;
	
	foreach my $token ( @token ) {
	    my %motif;
	    if ( $token =~ /,(f|c|d):/i ) {
	        # tag:value format
	        my @field = split ",", $token;
	        $motif{id} = $field[0];
	        foreach my $field ( @field ) {
                my ( $tag, $value ) = split ":", $field;                
                if ( $tag =~ /c/i ) {
                    # motif count
                    $motif{count} = $value;
	               
                } elsif ( $tag =~ /d/i ) {
                    # encoded dfs
                    $motif{encoded_dfs} = $value;
	               
                } elsif ( $tag =~ /f/i ) {
                    # iteration first observed
                    $motif{first_observed} = $value;
	               
                } elsif ( $tag =~ /m/i ) {
                    # iteration first observed
                    $motif{mapping} = $value;
                   
                } else {
                    print STDERR "Fingerprint::motifFromCompact - Unknown motif tag\n";
                    print STDERR "     line:$line\n";
                }
	        }
	           	           
        } else {
	        # simple comma separated list
	        my @field = split ",", $token;
            $motif{id} = $field[0];
            foreach my $field ( @field ) {
                if ( defined $field[1] ) { $motif{count}          = $field[1]; }
                if ( defined $field[2] ) { $motif{first_observed} = $field[2]; }
                if ( defined $field[3] ) { $motif{encoded_dfs}    = $field[3]; } 
                if ( defined $field[4] ) { $motif{mapping}        = $field[4]; }               
            }
	    }
	    push @motif_list, \%motif;	    
	}
	
	return \@motif_list;
}

# End of motifFromCompact

#------------------------------------------------------------------------------
# motifFromXML
#
# In the original format, motifs were presented in the tag <motif_block>, now
# they are represented as 
#   
#   <motif_list>
#       <motif>
#           <id>
#           <count>
#           <first_observed>
#           <encoded_dfs>
#       </motif>
#   </motif_list>
#
# USAGE
#    $ref_motif_list = motifFromXML( $line );
#------------------------------------------------------------------------------
sub motifFromXML{
	my ( $motiflistblock ) = @_;
	my @motif_list;
	
	foreach my $k ( keys %$motiflistblock ) {
	    
	}
	
	return \@motif_list;
}

# End of motifFromXML

#------------------------------------------------------------------------------
# motifIndex
#
# Return a hash that gives the index number of each motif with the motif ID as
# the key
#
# USAGE
#    %motif_index = $fingerprint->motifIndex;
#------------------------------------------------------------------------------
sub motifIndex{
	my ( $self ) = @_;
	my %index;
	
	my $motiflist = $self->motifList;
	
	my $count = 0;
	foreach my $motif ( @{$motiflist} ) {
	    $index{$$motif{id}} = $count;
	    $count++;
	}
	
	return %index;
}

# End of motifIndex

#------------------------------------------------------------------------------
# motifPerLine
#
# Sets the number of motifs per line for compact/block format.  Modifies global
# variable $MOTIF_PER_LINE
#
# USAGE
#    $motifs_per_line = motifPerLine( n );
#------------------------------------------------------------------------------
sub motifPerLine{
	my ( $self, $n) = @_;
	
	$MOTIF_PER_LINE = $n;
	
	return $n;
}

# End of motifPerLine



#------------------------------------------------------------------------------
# mappingMin
#
# Return the minimum mapping for motifs in the motif list.  
#
# USAGE
#    $min = $fingerprint->mappingMin;
#------------------------------------------------------------------------------
sub mappingMin{
	my ( $self ) = @_;
	my $min= 1000000000000;
	
	foreach my $motif ( @{ $self->motifList } ) {
	    my $mapping = $$motif{mapping} || 0;
	    if ( $mapping < $min ) { $min = $mapping; }
	}
	
	return $min;
}

# End of mappingMin



#------------------------------------------------------------------------------
# motifMin
#
# Return the minimum count for motifs in the motif list.  If no counts are 
# present, value is zero.  If no motifs are present, value is undef.
#
# USAGE
#    $min = $fingerprint->motifMin;
#------------------------------------------------------------------------------
sub motifMin{
	my ( $self ) = @_;
	my $min= 1000000000000;
	
	foreach my $motif ( @{ $self->motifList } ) {
	    my $count = $$motif{count} || 0;
	    if ( $count < $min ) { $min = $count; }
	}
	
	return $min;
}

# End of motifMin

#------------------------------------------------------------------------------
# motifN
#
# Return the number of motifs in motif list.
#
# USAGE
#    $nmotif = $fingerprint->motifN;
#------------------------------------------------------------------------------
sub motifN{
	my ( $self ) = @_;
	
	my $n = @{$self->motifList};
	
	return $n;
}

# End of motifN

#------------------------------------------------------------------------------
# readCompact
#
# Read in compact format fingerprint ifrom a file.  See asCompact for details.
#
# USAGE
#    $n_motif = $fingerprint->readCompact( $filename );
#------------------------------------------------------------------------------
sub readCompact{
	my ( $self, $filename ) = @_;
	
    open ( my $file, "<", $filename ) or 
       die "Fingerprint::readCompact unable to open $filename for input\n\n";
       
    my @motiflist;
    while ( my $line = <$file> ) {
        chomp $line;
        
        if ( $line =~ /^query_id:/ ) {
            my ( $tag, $qid ) = split ":", $line, 2;
            $self->queryId( $qid );
            
        } elsif ( $line =~ /^query_vertex:/ ) {
            my ( $tag, $q_vnum ) = split ":", $line, 2;
            $self->queryVertex( $q_vnum );
            
        } elsif ( $line =~ /^query_edge:/ ) {
            my ( $tag, $q_enum ) = split ":", $line, 2;
            $self->queryEdge( $q_enum );
            
        } elsif ( $line =~ /^type:/ ) {
            my ( $tag, $type ) = split ":", $line, 2;
            $self->type( $type );
            
        } elsif ( $line =~ /^program:/ ) {
            my ( $tag, $program) = split ":", $line, 2;
            $self->iteration( $program );
            
        } elsif ( $line =~ /^iteration:/ ) {
            my ( $tag, $iteration ) = split ":", $line, 2;
            $self->iteration( $iteration );
            
        } elsif ( $line =~ /^database_id:/ ) {
            my ( $tag, $dbid ) = split ":", $line, 2;
            $self->databaseId( $dbid );
            
        } elsif ( $line =~ /^motif_n:/ ) {
            # do nothing, the number of motifs is easily obtained from the 
            # motif list
                        
        } elsif ( $line =~ /^motif_list:/ ) {
            # do nothing, all following values are the motif IDs
            
        } else {
            # any line without a tag is part of the motif list
            my $motifs = motifFromCompact( $line );
            push @motiflist, @{$motifs};
            
        }
    }
    close $file;
    $self->motifList( \@motiflist );
    my $nmotif = @motiflist;	
	
	return $nmotif;
}

# End of readCompact

#------------------------------------------------------------------------------
# readXML_old
#
# Read in XML format fingerprint from a file.  See asXML for details.
# Old version, reads using XML::Simple
#
# USAGE
#    $n_motif = $fingerprint->readXML_old( $filename );
#------------------------------------------------------------------------------
sub readXML_old{
	my ( $self, $filename) = @_;
    open ( my $file, "<", $filename ) or 
       die "Fingerprint::readXML_old unable to open $filename for input\n\n";
       
    my $xml = new XML::Simple;
    my $d   = $xml->XMLin( $filename, SuppressEmpty => 1, XMLDecl => 1 );
    #print Dumper( $d );
    my @motiflist;
    
    $self->queryId    ( $d->{query}->{query_id} );
    $self->queryVertex( $d->{query}->{query_vertex} );
    $self->queryEdge  ( $d->{query}->{query_edge} );
    $self->type       ( $d->{fingerprint}->{type} );
    $self->iteration  ( $d->{fingerprint}->{iteration} );
    $self->program    ( $d->{fingerprint}->{program} );
    $self->databaseId ( $d->{database}->{database_id} );
    my $motif_str  =    $d->{motif_list};
    push @motiflist, @{motifFromXML( $d->{motif_list} )};
    $self->motifList( \@motiflist );
		
	my $n_motif = @motiflist;
	return $n_motif;
}

# End of readXML_old

#------------------------------------------------------------------------------
# readXML
#
# Read in XML format fingerprint from a file.  See asXML for details.
# New version uses XML::LibXML instead of XML::Simple
#
# USAGE
#    $n_motif = $fingerprint->readXML( $filename );
#------------------------------------------------------------------------------
sub readXML{
    my ( $self, $filename) = @_;
    unless ( -r $filename ) { 
       die "Fingerprint::readXML unable to open $filename for input\n\n";
    }
       
    my $parser = XML::LibXML->new();
    my $xml = $parser->parse_file( $filename );
    my @motiflist;
    
    my %query = ( queryId     => '//query/query_id',
                  queryVertex => '//query/query_vertex',
                  queryEdge   => '//query/query_edge',
                  type        => '//fingerprint/type',
                  iteration   => '//fingerprint/iteration',
                  program     => '//fingerprint/program',
                  timeElapsed => '//fingerprint/time_elapsed',
                  databaseId  => '//database/database_id'
                );
    my $node;
    foreach my $q ( keys %query ) {
        $node = $xml->findnodes( $query{$q} )->string_value;
        $self->$q( $node );
    }
    
    foreach my $block ( $xml->findnodes('//motif') ) {
        my %motif;
        $motif{id}             = $block->findnodes( 'id' )->string_value;
        $motif{count}          = $block->findnodes( 'count' ) ->string_value;
        $motif{first_observed} = $block->findnodes( 'first_observed' ) ->string_value;
        unless ( $motif{first_observed} ) { delete $motif{first_observed}; }
        $motif{encoded_dfs}    = $block->findnodes( 'encoded_dfs' ) ->string_value;
        unless ( $motif{encoded_dfs} ) { delete $motif{encoded_dfs}; }
        $motif{mapping}        = $block->findnodes( 'mapping' ) ->string_value;
        unless ( $motif{mapping} ) { delete $motif{mapping}; }
       push @motiflist, \%motif;
    }
    $self->motifList( \@motiflist );
        
    my $n_motif = @motiflist;
    return $n_motif;
}

# End of readXML


#------------------------------------------------------------------------------
# stringIf
#
# Returns a string with the tag:value pair separated by a colon if the variable 
# is defined, otherwise an empty string.
#
# USAGE
#    $str .= stringIf( 'query_id', $query_id );
#------------------------------------------------------------------------------
sub stringIf{
	my ( $tag, $value ) = @_;

	my $str = "";
	if ( defined $value ) {
	    $str .= qq{$tag:$value\n};
	}
	
	return $str;
}

# End of stringIf

#------------------------------------------------------------------------------
# sum
#
# Add a second fingerprint to the object.  
# summed values:
#   elapsed time
#   iteration
#   motif count
#
# mappings cannot be checked so they are set to undef
#
# USAGE
#    n_motif = $fingerprint->sum( $fprint2 );
#------------------------------------------------------------------------------
sub sum{
	my ( $self, $fprt ) = @_;
	
	unless ( $self->queryId eq $fprt->queryId ) {
	    my $id1 = $self->queryId;
        my $id2 = $fprt->queryId;
	    print STDERR "Fingerprint::sum - query IDs do not match ";
	    print STDERR "($id1,$id2)\n"
	}
	
	$self->iteration( $self->iteration + $fprt->iteration );
	$self->timeElapsed( $self->timeElapsed + $fprt->timeElapsed );
	$self->totalMapping( undef );
	    
    my $oldmotiflist  = $self->motifList;
    my %oldmotifindex = $self->motifIndex;
    my $newmotiflist  = $fprt->motifList;
    
    foreach my $newmotif ( @{$newmotiflist} ) {
        my $newid = $$newmotif{id};
        my $oi = $oldmotifindex{$newid};
        if ( defined $oldmotifindex{$newid} ) {
            # new motif exists in old - update
            my $oldmotif = $$oldmotiflist[ $oi ];
            $$oldmotif{count} += $$newmotif{count};
            #delete $$oldmotif{mapping};
                        
        } else {
            # new motif does not exist in old - create new entry
            my $oldmotif = $$oldmotiflist[ $oi ];
            my %this_motif = %{$newmotif};
            $this_motif{first_observed} += $$oldmotif{iteration};
            push @{$oldmotiflist}, \%this_motif;
            #delete $this_motif{mapping};
        }
    }
    
    my $nmotif = @{$oldmotiflist};
	return $nmotif;
}

# End of sum

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
    
    $unknown =~ s/([A-Z])/_\l$1/g;          # convert camel case to _lowercase
  
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
# $Log: Fingerprint.pm,v $
# Revision 1.5.4.3  2015/10/12 17:38:27  huang147
# added mappingMin subroutine.
#
# Revision 1.5.4.2  2014/02/06 20:44:59  huang147
# Keep the mapping information.
#
# Revision 1.5.4.1  2013/06/27 03:20:30  huang147
# Added a sort of motifs according to counts in the subroutine asXML.
#
# Revision 1.5.6.8  2013/03/14 17:58:26  gribskov
# Debugged.
#
# Revision 1.5.6.7  2013/03/14 12:45:05  gribskov
# Added function motifMin to get count of least frequent motif.
# Added function motifN to get count of motifs in motif list.
#
# Revision 1.5.6.6  2013/03/14 12:22:37  gribskov
# Added time_elapsed to compact output.
# Added function motifIndex to allow lookup of motif by name.
# Added function sum to add a second fingerprint to the object.
#
# Revision 1.5.6.5  2013/03/13 18:06:41  gribskov
# Added total mapping and motif mappings to data structure.
# Added time_elapsed.
# Added function motifPerLine to control motifs in compact format.
#
# Revision 1.5.6.4  2013/03/13 12:55:05  gribskov
# Finished update to LibXML and extension of format.
# Included example of format at end of file.
#
# Revision 1.5.6.3  2013/03/12 22:24:49  gribskov
# Could not make new format work with XML::Simple.  converting
# to XML::LibXML with Xpath queries (probably should use SAX for this)
#
# Revision 1.5.6.2  2013/03/11 13:37:21  gribskov
# Added program.
# Moved program, type, iteration to fingerprint info.
#
# Revision 1.5.6.1  2013/03/11 12:55:15  gribskov
# Added query edges, fingerprint type, and number of iterations to object.
#
# Revision 1.5  2012/03/27 21:42:01  gribskov
# Added another XML writing function, asXMLSimple, that uses XML::Simple.
#
# Revision 1.4  2012/03/27 20:49:52  gribskov
# Added readXML to read the XML format.
#
# Revision 1.3  2012/03/27 20:31:36  gribskov
# Added functions to read and write compact format: readCompact, asCompact.
#
# Revision 1.2  2012/03/27 20:04:53  gribskov
# Add asXML function to print out XML format.
#
# Revision 1.1  2012/03/27 17:51:35  gribskov
# Initial version
#
################################################################################    
1;
__END__
<?xml version="1.0"?>
<XIOS_fingerprint>   
    <!-- hand constructed -->
    <query>
        <query_id>trna_2ZZM.xios</query_id>
        <query_vertex>6</query_vertex>
        <query_edge>8</query_edge>
    </query>   

    <fingerprint>
        <type>random</type>
        <iteration>200</iteration>
        <time_elapsed>3.85</time_elapsed>
        <program>fingerprint_shred.pl 1.1.2.4</program>
    </fingerprint>

    <database>
        <database_id>db.storable</database_id>
    </database>
    
        
    <motif_list>
        <motif_n>4</motif_n>
        <motif>
            <id>4_15</id>
            <count>40</count>
            <first_observed>10</first_observed>
            <encoded_dfs>aa1e</encoded_dfs>
        </motif>
        <motif>
            <id>4_17</id>
            <count>47</count>
        </motif>
        <motif>
            <id>4_20</id>
            <count>89</count>
        </motif>
        <motif>
            <id>4_27</id>
            <count>24</count>
        </motif>        
    </motif_list>
    
</XIOS_fingerprint>

compact format

query_id:trna_2ZZM.xios
query_vertex:6
query_edge:8
type:random
iteration:200
program:fingerprint_shred.pl 1.1.2.4
database_id:db.storable
motif_n:4
motif_list:
4_15,c:40,f:10,d:aa1e 4_17,c:47 4_20,c:89 4_27,c:24 
    	
