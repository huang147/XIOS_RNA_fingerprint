package DfsGraph;
#-----------------------------------------------------------------------------
# $Id: DfsGraph.pm,v 1.26.4.1 2015/10/13 17:52:55 huang147 Exp $
# 
# This package holds multiple graphs (structures) in a single object.
# Contains the information needed for DFS graph matching.  Graphs are stored 
# in a single data structure of edges
#
# in addition to subroutines, the following functions are available
#   vertexIndex( $n )   return pointer to first vertex in graph $n
#   edgeIndex( $n )     return pointer to first edge in graph $n
#
# note that all hash keys can also be used as accessor functions
#
# $Author: huang147 $     $Date: 2015/10/13 17:52:55 $
#
# Revision log at end of file.
#-----------------------------------------------------------------------------
use strict;
use Stem;
use Topology;
use warnings;
# export global variable definitions

use vars qw{ $VERSION @ISA @EXPORT };
use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( $MAX_EDGE_TYPE $NO $YES $INCLUDED $JNCLUDED 
                  $OVERLAP $FAILURE $SUCCESS );

#construct version from CVS tag

my $REVISION = '$Revision: 1.26.4.1 $';
( $VERSION ) = $REVISION =~ /Revision: ([\d.]*) /;

# Global Variables
 
our $MAX_EDGE_TYPE = 3;      # index of the highest numbered edge type
our $NO            = 0;      # symbol for false
our $YES           = 1;      # symbol for true

our $INCLUDED      = 0;      # numeric code for included edges I
our $J			   = 1;      # numeric code for included edges J
our $OVERLAP       = 2;      # numeric code for overlap (pseudoknot) edges
our $EXCLUDED      = 3;      # numeric code for exclusive edges

# for converting between = edge type letter to numeric 
our %EDGENUM       = ( i => $INCLUDED, j => $J, o => $OVERLAP, x => $EXCLUDED );
our @EDGETYPE      = ( 'i', 'j', 'o', 'x' );

our $FAILURE       = 0;      # edge extension had insufficient support
our $SUCCESS       = 1;      # edge extension has sufficient support

#-----------------------------------------------------------------------------
# new:
#
# constructor, currently set up to read RNA structure in "pair" format, e.g.
#   1 3 2 5 4 6
#
# All parameters are passed in a single hash.  hash keys are
#   type = pair (default)
#   file = name of file if reading from file
#   string = if entering graph as a string
#
# USAGE
#   graph = new DfsGraph( {parameters} );
#   graph = new DfsGraph->new( {type => "pair", file => "random.pair"} );
#   or
#   my $dfsgraph = new DfsGraph(  );
#   $dfsgraph->addPairCoordinate($st); #Add pair files
#   $dfsgraph->addDiffgraph( $file, $SKIPX  ); #add diffgraph files
#
#-----------------------------------------------------------------------------
sub new{
    my ( $class, $param ) = @_;
    # print "\ncoming before here\n";
    # print STDERR "DfsGraph:new param:",$param,"\n";
    # foreach my $k ( keys %$param ) {
        # print sTDERR "DfsGraph:new parameter$k =>  $$param{$k}\n";
    # }
    my $self = {};
    bless  $self, $class;
    $self->_init;
    
    # Defaults

    my $type = "pair";

    if ( defined $param->{'type'} ) {
        $type = $param->{'type'};
    }

    # add types here for other formats

    if ( $type =~ /pair/i ) {

        my @pairs;

        if ( defined $param->{'file'} ) {
            open( PAIR, "<$param->{'file'}" ) || 
                  die "DfsGraph:new - Cannot open pair file $param->{'file'}\n\n";
	        @pairs = <PAIR>;
	        close PAIR;
        }

        if ( defined $param->{'string'} ) {
            my @graph = split /[,;.|]+/, $param->{'string'};
	        #print "\ncoming inside defind string\n";
	        foreach my $line (@graph) {
		        #print "$line\n";
		        push @pairs, $line;
	        }
	        # push @pairs, @graph;
        }
	
        if ( @pairs ) { $self->addPair( @pairs ); }

    } elsif ( $type =~ /diffgraph/i ) {
        # load diffgraph format 
        if ( defined $param->{'file'} ) {
            $self->addDiffgraph( $param->{'file'} );
        }
    	if ( defined $param->{'description'} ) {
    	    my $descriptions = $self->description;
    	    push @$descriptions, $param->{'description'};
    	}

    } elsif ( $type =~ /xios/i ) {
        # add xios format structure
        $self->addXios( $param->{file} );

    } elsif ( $type =~ /edgedata/i) {
    	if ( defined $param->{'edgelist'} ) {
    	    my $each_struct = $param->{'edgelist'};
    	    my ($edge_arr, $total_v) = @$each_struct;
    	    # my $text='';
    	    my $n_structure = $self->addStructure;
    	    
    	    # check to see if total vertex is greater than zero
    	    my $offset = $self->addVertices( $total_v ); 
            # print "Total vertex: $offset\n";
    	    foreach my $e (@$edge_arr) {
    		my $this_edge;
    		my ($v1, $v2, $edge_type) = @$e;
    		$this_edge = $self->addEdge($v1, $v2, $edge_type);
    		# print "THIS Edge: $this_edge\n";
    		$self->addEdgeToVE( $v1, $v2, $this_edge );
    	    }
    	}
    	
    	if ( defined $param->{'vertexloc'} ) {
        	# if vertexloc is passed, then use it
    	    $self->{'vertex_loc'} = $param->{'vertexloc'};  
    	}
    	
    	if ( defined $param->{'description'} ) {
    	    my $descriptions = $self->description;
    	    push @$descriptions, $param->{'description'};
    	}
	
    }

    return $self;
}

# end of new

#-----------------------------------------------------------------------------
# _init:
#
# create the attributes for the DfsGraph object.  Try to keep all initializ-
# ation here.
#
# The internal structure of the object is a list of vertices and edges for all  
# loaded graphs.  The first edge in each distinct graph is listed in the edge_index 
# hash.  For each edge, the two vertices and edge type are stored in the edge 
# hash as [v1,v2,numeric_edge_type].
#
# vertices are also stored as a single numbered series, the first vertex in each
# graph is stored in the vertex_index hash.  For each vertex, a list of incident 
# edges is stored in the vertex array.  The value of vertex is a reference to
# an array of references to edges.
#
# USAGE:
#   $dfs->_init;
#-----------------------------------------------------------------------------
sub _init{
    my ( $self ) = @_;

    $self->{'revision'} = $REVISION;

    # whole graph.  The description can be used to store any kind of information
    # about the graph, for instance, filename containing the graph
    $self->{'n_structure'} = 0;                 # total number of graphs
    $self->{'description'} = [];                # description string
    $self->{'sorted'} = 0;                      # VE list is unsorted

    # vertices
    $self->{'n_vertex'} = 0;                    # total number of vertices
    $self->{'vertex_index'}[0] = 0;             # first vertex for each graph
    $self->{'vertex'} = ();                     # edge list for each vertex
    
    # vertex_stem is a hash which maps a vertex to its corresponding stem object
    # this hash is indexed by the vertex number and contains a reference to the 
    # stem object corresponding to a vertex.  Since it is not sparse, it is 
    # crazy to use a hash (mrg)
    
    $self->{'vertex_loc'} = {};
   
    # edges
    $self->{'n_edge'} = 0;                      # total number of edges
    $self->{'edge_index'}[0] = 0;               # first edge in each graph
    $self->{'edge'} = ();                       # vertices and edge type

    # DFS - may not be used
    $self->{'dfs'} = [];
    $self->{'dfsLen'} = 0;
    $self->{'dfs_max_vertex'} = undef;
    $self->{'dfs_path_pos'} = 0;
    $self->{'main_path_len'} = 0;

    return;
}

# end of _init

#------------------------------------------------------------------------------
# edgeTypeToNumber
#
# Look up the letter for the edge type and return the corresponding number.  
# Conversion is based on the global definitions at the top of the package.
#
# USAGE
#    $edge_num = $G->edgeTypeToNumber( $edge_letter );
#------------------------------------------------------------------------------
sub edgeTypeToNumber{
	my ( $G, $edge_letter ) = @_;
	
	return $EDGENUM{$edge_letter};
}

# End of edgeTypeToNumber

#-----------------------------------------------------------------------------
# addDiffgraph:
#
# add a Diffgraph structure from a file. X edges can be skipped if the second 
# argument is TRUE
#
# TODO should return number of vertices added
#
# USAGE
#   $dfs->addDiffgraph( $filename);
#   $dfs->addDiffgraph( $filename, $skip_edges); #(if $skip_edges = 1 || TRUE, x edges are skipped)
#-----------------------------------------------------------------------------
sub addDiffgraph{
    my ( $self, $filename, $skip_xedge) = @_;
 
    # if skip_edge is true, then x edges are skipped,  
    # else if skip_edges if false or undefined, they are included. 
    my @vertex;
    open( DIFF, "<$filename" ) || 
	die "DfsGraph:addDiffgraph - Cannot open diffgraph file $filename\n\n";
    
    my @first_part;
    while ( my $line = <DIFF> ) {
        next if ( $line =~ /^\s*$/ );
        
        $line =~s/^\s+//;
    	chomp $line;
    	$line =~s/\s+/ /;
    	
    	if ( $line =~ /\[/ ) {
    	    push @first_part, $line; 
    	} else {
    	    
    	    my ( $v, $edgelist ) = split ":", $line;
    	    #print "   vertex[$v] => $edgelist\n";
    	    $vertex[$v]->{edges} = $edgelist;
    	}
    }
    
    my $n_vertices = @vertex;
    my $n_structure = $self->addStructure;
    $self->addVertices( $n_vertices );
    # print "   $n_vertices vertices read\n";
    
    foreach my $v1 ( 0 .. $#vertex ) {
        # the string $vertex[$v1]->{edges} looks like
        # 3x  4i  5x  6x  7i  8i  9i 10i 11i 12i 13x 14x
        my @edges = split " ", $vertex[$v1]->{edges};
        # print "   vertex $v1\n";
        while ( @edges ) {
    	    my $v2   = shift @edges;
    	    my $type = shift @edges;
    	    #print "      edge $v2, $type\n";
            if ( $type =~ /i/ ) { $type = $INCLUDED; }
            if ( $type =~ /j/ ) { 
                $type = $INCLUDED; 
                die "DfsGraph:AddDiffgraph - j edges types should not be present in diffgraph files\n";
            }
            if ( $type =~ /o/ ) { $type = $OVERLAP; }
            if ( $type =~ /x/ ) { 
    		    #if $skip_xedge is true, then skip x edges
    		    next if ( $skip_xedge ); 
    		    
                $type = $EXCLUDED; 
    	    }
            my $thisedge = $self->addEdge( $v1, $v2, $type );
            $self->addEdgeToVE( $v1, $v2, $thisedge );
        }
    }
    
    my $offset = $self->vertexIndex( $self->n_structure - 1 );
    my $vertex_loc_hash = $self->{'vertex_loc'}; 
    foreach my $line (@first_part) {
        # last if ( $line =~ /^\s*$/ );
        # next unless ( $line =~ /\[/ );
	    # chomp $line;
	    my ( $n, $center, $b1, $l1, $l2, $r1, $r2, $b2, $lvienna, $rvienna ) = split " ", $line;	
       	my $stem_obj = new Stem( {  'left1'       => $l1,
				                    'left2'       => $l2,
				                    'right1'      => $r1,
			                        'right2'      => $r2,
				                    'vienna_left' => $lvienna,
				                    'vienna_right'=> $rvienna
				                 }
				);
        my $newvertex = $offset + $n;
        $vertex_loc_hash->{$newvertex} = $stem_obj;
    }
    close DIFF;
}

# end of addDiffgraph

#------------------------------------------------------------------------------
# addXios
#
# Add a structure in XIOS format using Topology.pm.  j edges, if present, are
# are added as i edges with the vertices reversed. Returns the number of vertices
# added.
#
# USAGE
#    $v_added = $G->addXios( $filename );
#------------------------------------------------------------------------------
sub addXios{
	my ( $dfs, $filename ) = @_;
	
	my $rna = new Topology;
	$rna->XIOSRead( $filename );
	
    my $n_structure = $dfs->addStructure;
    
    # add the edges and vertices from the adjacency matrix.  the number of
    # elements in @adj is the number of vertices
    
    my @adj = @{ $rna->adjacency };
    my $nv = @adj;
    $dfs->addVertices( $nv );
    
    # add the edges
    foreach my $v1 ( 0 .. $nv-1 ) {
        foreach my $v2 ( $v1+1 .. $nv-1 ) {
            
            my $edge = $adj[$v1]->[$v2];
            next unless ( $edge =~ /[ijox]/i );     # skip anything but i, j, o, and x
            
            my $edgenum = $dfs->edgeTypeToNumber( $edge );

            
            if ( $edgenum == $J ) {
                # add j edges as i edges with v1 and v2 reversed
                $edgenum = $INCLUDED;
                my $thisedge = $dfs->addEdge( $v2, $v1, $edgenum );
                $dfs->addEdgeToVE( $v2, $v1, $thisedge );
                
            } else {
                # add other edges in with v1 first
                my $thisedge = $dfs->addEdge( $v1, $v2, $edgenum );
                $dfs->addEdgeToVE( $v1, $v2, $thisedge );
            }                       
        }
    }
    
    # TODO store stems and fill in vertex_loc 
      
 	return $nv;
}

# End of addXios

#-----------------------------------------------------------------------------
# addXiosText:
#
# add a Xios graph text. X edges can be skipped if the second argument is TRUE.
#
# (MRG)
# Its not clear to me what format this is supposed to add.  Many of the comments
# refer to DiffGraph format.  It's not clear if these have just been copied from 
# addDiffgraph when this function was created. Since Diffgraph format is now 
# superceded by the better specified XIOS format.
#
# USE OF THIS FUNCTION IS DEPRECATED
#
# USAGE
#   $dfs->addDiffgraph( $filename);
#   $dfs->addDiffgraph( $filename, $skip_edges); #(if $skip_edges = 1 || TRUE, x edges are skipped)

#-----------------------------------------------------------------------------
sub addXiosText {
    my ( $self, $filename, $skip_xedge) = @_;
    my @vertex;
    open( DIFF, "<$filename" ) || 
	die "DfsGraph:addDiffgraph - Cannot open diffgraph file $filename\n\n";
    
    my @first_part;
    while ( my $line = <DIFF> ) {
        next if ( $line =~ /^\s*$/ );
    	$line =~s /^\s+//;
    	chomp $line;
    	$line =~s /\s+/ /;
    	
    	if ( $line =~ /\[/ ) {
    	    push @first_part, $line; 
    	} else {
    	    my ( $v, $edgelist ) = split ":", $line;
    	    #print "   vertex[$v] => $edgelist\n";
    	    $vertex[$v]->{edges} = $edgelist;
    	}
    }
    
    my $n_vertices = @vertex;
    my $n_structure = $self->addStructure;
    $self->addVertices( $n_vertices );
    #print "   $n_vertices vertices read\n";
    
    foreach my $v1 ( 0 .. $#vertex ) {
        my @edges = split " ", $vertex[$v1]->{edges};
        # print "vertex $v1: |$vertex[$v1]->{edges}|||@edges|\n";
        while ( @edges ) {
            # my $type = pop @edges;
            # my $v2 = pop @edges;
    	    my $v2 = shift @edges;
    	    my $type = chop $v2; #shift @edges;
            # my $type = shift @edges;
    	    # print "      edge $v2, $type\n";
            if ( $type eq 'i' || $type eq 'i,') { $type = 0; }
            if ( $type eq 'j' || $type eq 'j,') { 
                $type = 1; die "j edges types should not be present in diffgraph files\n";}
            if ( $type eq 'o' || $type eq 'o,' ) { $type = 2; }
            if ( $type eq 'x' || $type eq 'x,') { 
    		     # if $skip_edge is true then skip x edges
    		     next if ( $skip_xedge ); 
    		     
    		     $type = 3; 
            }
            my $thisedge = $self->addEdge( $v1, $v2, $type );
            $self->addEdgeToVE( $v1, $v2, $thisedge );
        }
    }
    
    my $offset = $self->vertexIndex( $self->n_structure - 1 );
    # print "the offset fo the graph is: $offset\n";
    my $vertex_loc_hash = $self->{'vertex_loc'}; 
    foreach my $line (@first_part) {
        # last if ( $line =~ /^\s*$/ );
        # next unless ( $line =~ /\[/ );
        # print "$line\n";
    	# chomp $line;
    	my ( $n, $center, $b1, $l1, $l2, $r1, $r2, $b2, $lvienna, $rvienna ) = split " ", $line;	
           	my $stem_obj = new Stem(    {'left1'       => $l1,
    				                     'left2'       => $l2,
    				                     'right1'      => $r1,
    				                     'right2'      => $r2,
    				                     'vienna_left' => $lvienna,
    				                     'vienna_right'=> $rvienna
    				                    }
    				                );
    	my $newvertex = $offset + $n;
    	$vertex_loc_hash->{$newvertex} = $stem_obj;
    }
    close DIFF;
    return;
}

# end of addDiffgraph

#-----------------------------------------------------------------------------
# addPair:
#
# add edges in pair format to the edge list.  Each element of the @pairs
# array should be a complete structure in a string.  
#
# USAGE
#   $dfs->addPair( @pairs );
#-----------------------------------------------------------------------------
    sub addPair{

    my ( $self, @structure ) = @_;

    my $vertex_type = 0;       # currently all vertices are the same
    # my $transaction = 0;

    #print STDERR "DfsGraph:addPair adding\n   structure:@structure\n";

    foreach my $line ( @structure ) {

        my @vertex;
        my @edge;
        my $v = 0;

        # clean input string
        $line =~ s/[^ \d]*//g;   # remove everything but spaces and digits

        my $n_structure = $self->addStructure( $line );

        # commented out print statements would write structure for gspan or gaston
        # print "# $line";
        # print "t # $transaction\n";

        # calculate vertices

        my @halfstem = split " ", $line;
        my $have_stems = @halfstem;
        while (  $have_stems ) {
	    $vertex[$v]{begin} = shift @halfstem;
	    $vertex[$v]{end} = shift @halfstem;
	    
            # my @sorted_halfstem = sort { $a <=> $b } @halfstem;
            # $vertex[$v]{begin} = shift @sorted_halfstem;
            # $vertex[$v]{end} = shift @sorted_halfstem;
            # print "v $v $vertex_type\n";
            $v++;
            # print " $v  begin:end ",$vertex[$v]{begin}.":".$vertex[$v]{end}."\n";
            $have_stems = @halfstem;
        }

        my $offset = 0;
        if ( $v > 0 ) {     # some vertices were found
            $offset = $self->addVertices( $v );

            #print STDERR "Structure $n_structure $v vertices\n";

        } else {
            next;       # skip to next structure if this one was empty
        }

        # calculate edges

        my $e = 0;
        foreach my $v1 ( 0 .. $v - 1 ) {        
            # stop at $v to make sure you don't get vertices from previous structure
            
            my $this_edge;
            foreach my $v2 ( $v1+1 .. $v - 1 ) {

                # skip serial edges
                if ( $vertex[$v1]{end} < $vertex[$v2]{begin} ) { next; }  
                
                my $edge_type = "";
		
                if ( $vertex[$v1]{begin} < $vertex[$v2]{begin} &&
		            $vertex[$v1]{end} > $vertex[$v2]{end} ) {
                    $edge_type = $INCLUDED;
                    # print "e $v1 $v2 0\n";    # edge type 0 is i
                } else {
                    $edge_type = $OVERLAP;
                    # print "e $v1 $v2 1\n";     # if its not i it must be o
                }
                $this_edge = $self->addEdge( $v1, $v2, $edge_type );
		        $self->addEdgeToVE( $v1, $v2, $this_edge );
		        $e++;
	        }
        }
        #print STDERR "     $e edges\n\n";

    }   # end of loop over structures in @pairs

    return;
}

# end of addPair

#-----------------------------------------------------------------------------
# addPairCoord:
#
# add edges in pair format to the edge list.  Each element of the @pairs
# array should be a complete structure in a string. This is similar to addPair except   
# the fact that it adds fake coordinates to each vertex
#
# USAGE
#   $dfs->addPair( @pairs );
#-----------------------------------------------------------------------------
    sub addPairCoordinate{

    my ( $self, @structure ) = @_;

    my $vertex_type = 0;       # currently all vertices are the same
    # my $transaction = 0;
    # print STDERR "DfsGraph:addPair adding\n @structure\n";
    
    foreach my $line ( @structure ) {
        my @vertex;
        my @edge;
        my $v = 0;

        $line =~ s/[^ \d]*//g;   # remove everything but spaces and digits
        # my $n_structure = $self->addStructure( $line );
        my $n_structure = $self->addStructure;
        # print "strucutre: $n_structure, $line";  
        # calculate vertices
        my @halfstem = split " ", $line;

        my $have_stems = @halfstem;
        while (  $have_stems ) {
    	    $vertex[$v]{begin} = shift @halfstem;
    	    $vertex[$v]{end} = shift @halfstem;
            $v++;
            # print " $v  begin:end ",$vertex[$v]{begin}.":".$vertex[$v]{end}."\n";
            $have_stems = @halfstem;
        }
    	my $offset = 0;
        if ( $v > 0 ) {     # some vertices were found, add them
            $offset = $self->addVertices( $v );
	        # print STDERR "Structure $n_structure $v vertices\n";
	    
        } else {
            next;       # skip to next structure if this one was empty
        }
		
        # calculate edges
        my $e = 0;
        foreach my $v1 ( 0 .. $v - 1 ) {
            # stop at $v to make sure you don't get vertices from previous structure
                    
    	    my $this_edge;
            foreach my $v2 ( $v1+1 .. $v - 1 ) {
    		if ( $vertex[$v1]{end} < $vertex[$v2]{begin} ) { next; }  # skip serial edges
                my $edge_type = "";
          		if ( $vertex[$v1]{begin} < $vertex[$v2]{begin} &&
        		     $vertex[$v1]{end} > $vertex[$v2]{end} ) {
                    $edge_type = $INCLUDED;
                    # print "e $v1 $v2 0\n";    # edge type 0 is i
                } else {
                    $edge_type = $OVERLAP;
                    # print "e $v1 $v2 1\n";     # if its not i it must be o
                }
                $this_edge = $self->addEdge( $v1, $v2, $edge_type );
                $self->addEdgeToVE( $v1, $v2, $this_edge );
                $e++;
	       }
        }
    	my $factor = 5;
    	my $new_offset = $self->vertexIndex( $self->n_structure - 1 );
    	my $vertex_loc_hash = $self->{'vertex_loc'};
    	# setting fake coordinates in vertex-loc
    	for my $i ( 0 .. $v - 1) {
    	    # print "index: $i| VAL: $vertex[$i]\n"; 
    	    # print "KEY: begin, Value: $vertex[$i]{$begin}\n";
    	    my $left = (($vertex[$i]{'begin'}) * $factor) + int(rand(4));
    	    my $right = (($vertex[$i]{'end'}) * $factor) + int(rand(4));
    	    
    	    my $stem_hash = {  'left1'       => $left,
    			               'left2'       => $left+1,
    			               'right1'      => $right,
    			               'right2'      => $right+1,
    			               'vienna_left' => '((',
    			               'vienna_right'=> '))'
    			            };
    	    my $stem_obj = new Stem($stem_hash);
    	    my $newvertex = $i + $new_offset;
    	    # print "Adding vertex: $newvertex, and stem obj: $stem_obj|\n";
    	    $vertex_loc_hash->{$newvertex} = $stem_obj;
    	}
    
        # print STDERR "     $e edges\n\n";
    }   # end of loop over structures in @pairs
        return;
}
# end of addPairCoordinate

#-----------------------------------------------------------------------------
# addEdgeToVE
# 
# adds an edge to these vertices in the vertex-edge list 
#
# Usage
#   $self->addEdgeToVE( $v1, $v2, $this_edge );
#-----------------------------------------------------------------------------
sub addEdgeToVE{
    my( $self, $v1, $v2, $edge ) = @_;

    # print STDERR "DfsGraph:AddEdgeToVE v1:$v1  v2:$v2  edge:$edge\n";

    my $offset = $self->vertexIndex( $self->n_structure - 1 );
    $v1 += $offset;
    $v2 += $offset;
    # print STDERR "   offset:$offset v1:$v1   v2:$v2\n";

    # add edge to list of edges for each vertex
    push @{ $self->vertex( $v1 ) }, $edge;
    push @{ $self->vertex( $v2 ) }, $edge;

    # print STDERR "DfsGraph:AddEdgeToVE finished\n";
    return;
}

# end of addEdgeToVE

#-----------------------------------------------------------------------------
# sortVE
#
# sort the VE list for all structures by edge type.  Small edgetypes come 
# first.  This must be done before manipulating the graphs.  Functions
# that need a sorted list should check the $graph->sorted flag.
#
# The VE list, stored in $G->{vertex} has a list of the incident edges for each 
# vertex.  This functions rearranges the list in lexicographic order.  Note that
# because i edges are directed, v1 will not always be the vertex that the list
# corresponds to.
#
# Usage
#   $graph->sortVE;
#-----------------------------------------------------------------------------
sub sortVE{
    my ( $G ) = @_;

    my $min_type = $MAX_EDGE_TYPE;

    my $vertex = 0;
    foreach my $ve ( @{$G->vertex} ) {
        # $v is a reference to an array to refs to edges

        sub lexSort{
            # lexicographic sorting function, assumes i edges are numbered 0
            # regardless of direction
            # compile error warns that $vertex will not stay shared.  This is OK
            # it doesn't need to stay shared.
            if ( $a->[2] == 0 && $b->[2] == 0 ) {                
                # when both edges are type 0, direction must be checked
                our $vertex;
                if ( $a->[0] == $vertex && $b->[0] == $vertex ) {
                    # both edges are have vertex as v2, a == b
                    return 0;
                } elsif ( $a->[0] == $vertex ) {
                    # only a has vertex as v2, a < b
                    return -1;
                } elsif ( $b->[0] == $vertex ) {
                    # only b has vertex as v2, a > b
                    return 1;
                } else {
                    # neither edge has vertex as v2, a == b
                    return 0;
                }
            } else {
                # when edge types are different, just sort by edge type
                return $a->[2] <=> $b->[2];
            }
        }

        @{$ve} = sort lexSort @{$ve};
        $vertex++;
    }
	
    $G->sorted( 1 );
    return;
}

# end of sortVE

#-----------------------------------------------------------------------------
# addStructure:
#
# increment number of structures and set index pointers to the begining of
# the vertex and edge lists for the next structure.  If $pair_string is 
# provided it is saved in the structure. 
#
# It also creates the next structure and update its edgeIndex. 
# It updates the vertex Index for the current structure  
#
# USAGE
#   $n_structure = $dfs->addStructure( $pair_string );
#-----------------------------------------------------------------------------
sub addStructure{
    my ( $self, $pair_string ) = @_;

    # print STDERR "DfsGraph:addStructure   description:$pair_string\n";
    my $descriptions = $self->description;
    push @{$descriptions}, ($pair_string || "");
    $self->description( $descriptions );
    
    # the new structure is number $nstruc
    my $nstruc = $self->n_structure;    
    # print STDERR "   nstruc:$nstruc\n";

    my $n_edge = $self->n_edge;
    # print STDERR "   DfsGraph:addStructure   n edge:$n_edge\n";
    my $edge_index = $self->edgeIndex;
    push @{$edge_index}, $n_edge;
    $self->edgeIndex( $edge_index );
    # print STDERR "   struct:$nstruc  edge_index saved:$edge_index\n";
    # print STDERR "edge index:@{$edge_index}\n";

    # Update the vertex Index for the current structure
    my $n_vertex = $self->n_vertex;
    # print STDERR "   DfsGraph:addStructure   n vertex:$n_vertex\n";
    my $v_index = $self->vertexIndex;
    push @{$v_index}, $n_vertex;
    $self->vertexIndex( $v_index );

    # Creates the next structure in advance, and updates the edgeIndex for the 
    # next structure
    $nstruc = $self->n_structure( $self->n_structure+1 );
    # print STDERR "   nstruc:$nstruc\n";
    
    # initializing the edgeIndex for the next structure
    $self->{'edge_index'}->[$nstruc] = $self->{'edge_index'}->[$nstruc-1];
    # print STDERR "DfsGraph::addStructure finished\n";

    return $self->n_structure;
}

# end of addStructure

#-----------------------------------------------------------------------------
# addVertices:
#
# Add n vertices to the vertex list.  Vertices are currently unlabelled
# so the vertex list doesn't need to be constructed, all we need is the 
# offset needed to get to  a unique range of vertex numbers.
#
# $offset = $dfs->addVertices( $n );
#-----------------------------------------------------------------------------
sub addVertices{
    my ( $self, $new ) = @_;

    # print STDERR "DfsGraph::addVertices begin\n";
    my $offset =  $self->n_vertex;
    my $offset_new = $offset + $new;
    $self->n_vertex( $offset_new );
    # print STDERR "DfsGraph:addVertices offset:$offset   new:$new\n";
    $self->{'vertex_index'}[ $self->n_structure ] += $new;
    for my $v ( $offset .. $offset_new-1 ) {
        # initialize new ve list to anon array
        $self->vertex( $v, [] );           
    }

    # print STDERR "DfsGraph::addVertices end\n";
    return $offset_new;
}

# end of addVertices

#-----------------------------------------------------------------------------
# addEdge
# 
# add an edge of type $edge_type between $vertex1 and $vertex2 to the edge 
# list.  Note that the vertex numbers should be in terms of the  the local
# structure not the overall list.  An offset is automatically added to 
# the vertex number.  Vertices should be consecutive beginning at zero.
#
# USAGE
#   $n_edge_added = $dfs->addEdge( $v1, $v2, $edge_type )
#-----------------------------------------------------------------------------
sub addEdge{
    my ( $self, $v1, $v2, $type ) = @_;

    # print STDERR "DfsGraph:addEdge  v1:$v1  v2:$v2  type:$type\n";
    # get the offset for vertex numbers.  
    # offset is the index of the first edge in this structure

    my $offset = $self->vertex_index( $self->n_structure-1 );
    my $this_edge = [ $v1+$offset, $v2+$offset, $type ];
    push @{$self->{'edge'}}, $this_edge;
    my $n_edge = $self->n_edge( $self->n_edge+1 );
    
    # After an edge is added updating the starting edge index (offset) for the 
    # next structure
    $self->{'edge_index'}->[$self->n_structure] += 1;
    
    # print STDERR "DfsGraph:addEdge   new edge: $this_edge\n";

    return $this_edge;
}

# end of addEdge

#-----------------------------------------------------------------------------
# edgeList
# return a list of indices of edges for all structures, or a specified set. 
# 
#
# Usage:
#   $edge_list = $dfs->edgeList;         # retrieve all
#   $e1 =      = $dfs->edgeList( 2 );    # retrieve edges for structure 2
#   $e24       = $dfs->edgeList( 2, 4 )  # retrieve edges for structures 2 and 4
#
#-----------------------------------------------------------------------------
sub edgeList{
    my ( $self, @struct_id ) = @_;
    my $elist;
    my $n_structure = $self->n_structure;
    
    if ( @struct_id ) {
        # get vertices for selected structures
	
        foreach my $i ( @struct_id ) {

    	    # skip structures that don't exist
    	    next unless ( $i >= 0 and $i < $n_structure );
	    
            my $begin = $self->edgeIndex( $i );
            my $end;
            if ( my $next = $self->edge_index( $i+1 )) {
                $end = $next - 1;
            }  else {
                $end = $self->n_edge - 1;
            }

	        # print "Getting edgelist foreach structures: $i | $begin .. $end\n";
	    
            foreach my $i ( $begin..$end ) {
                my $e = $self->edge( $i );
                push @{$elist}, $e;
		        # print STDERR "DfsGraph:edgeList   retrieving edge $e\n";
            }
        }
    } else {
        # all vertices
        # print STDERR "DfsGraph:edgeList   retrieving all edges\n";
        $elist = $self->edge;
    }
    
    return $elist;
}

# end of edgeList

#-----------------------------------------------------------------------------
# edgeMinList
# return a ponter to a list of indices of the lexically minimal edges for all 
# structures, or of a specified set. 
# 
#
# Usage:
#   $edge_list = $dfs->edgeMinList;         # retrieve all
#   $e1 =      = $dfs->edgeMinList( 2 );    # retrieve edges for structure 2
#   $e24       = $dfs->edgeMinList( 2, 4 )  # retrieve edges for structures 2 and 4
#
#-----------------------------------------------------------------------------
sub edgeMinList{
    my ( $self, @struct_id ) = @_;

    my ( $all, $elist );
    my $n_structure = $self->n_structure;

    my $min_type = $MAX_EDGE_TYPE;
    if ( @struct_id ) {
        # print STDERR "DfsGraph:edgeMinList   edges for @struct_id\n";
        $all = $self->edgeList( @struct_id );
        # print STDERR "DfsGraph:edgeMinList   edgelist: ",@{$all},"\n";
    } else {
        $all = $self->edgeList;
    }

    foreach my $e ( @{$all} ) {
        my ( $v1, $v2, $type ) = @{$e};
        # print STDERR "DfsGraph:edgeMinList   e:$e   v1:$v1   v2:$v2   t:$type\n";
        if ( $type < $min_type ) {
            # print STDERR "DfsGraph:edgeMinList   clearing edge list\n";
            undef @{$elist};
            $min_type = $type;
        }
        if ( $type == $min_type ) { push @{$elist}, $e; }
    }

    return $elist;
}

# end of edgeMinList

#-----------------------------------------------------------------------------
# dump
#
# print the contents of all graph structures
#
# USAGE
#   $dfs->dump( $some_string_to_print );
#
#-----------------------------------------------------------------------------
sub dump{
    my ( $self, $text ) = @_;

    if ( defined $text ) {
        print "$text\n";
    }

    my $n_structures = $self->{'n_structure'};
    print "structures: $n_structures\n";

    foreach my $i ( 0 .. $n_structures - 1 ) {
        print "  structure $i:\n";
        my $description = $self->description;
        print "    description: $description->[$i]\n";
	
        my $vertex_begin = $self->vertexIndex( $i );
        my $vertex_end   = $self->vertexIndex( $i+1 ) - 1;
        print "    vertex: $vertex_begin - $vertex_end\n";
        foreach my $v ( $vertex_begin .. $vertex_end ) {
            my @thisv = @{ $self->vertex( $v ) };
            print "        vertex $v\n";
            # the vertex-edge list (indexed by connected vertex)
            foreach my $ve ( @thisv ) {
                my @ee = @{$ve};
                print "            $ve => @ee\n";;
            }
        }

        my $edge_begin = $self->edgeIndex( $i );
        my $edge_end   = $self->edge_index( $i+1 ) - 1;
        # Another way of getting edge end
        # my $edge_end;
        # if ( my $next = $self->edge_index( $i+1 )) {
        #     $edge_end = $next - 1;
        # }  else {
        #     $edge_end = $self->n_edge - 1;
        # }
        print "    Edge list (offset:$edge_begin)\n";
        foreach my $j ( $edge_begin .. $edge_end ) {
            my @edge = @{$self->edge( $j )};
            print "    $j   v1:$edge[0]   v2:$edge[1]   type:$edge[2]\n";
        }
        print "\n";
    }

    return;
}

# end of dump

#---------------------------------------------------------
# statsEachStructure
#
# Prints statistics of each structure
# TODO better if a string is returned IMO
#
# usage
#   $graph->statsEachStructure;
#----------------------------------------------------
sub statsEachStructure {
    my ($dfsgraph) = @_;
    
    my $n_structures = $dfsgraph->{'n_structure'};
    print "Structures: $n_structures\n";
    my ($included, $overlap, $exclusive) = (0, 2, 3);
    my ($isize, $osize, $xsize);
    
    foreach my $struct ( 0 .. $n_structures - 1 ) { 
    	my $vertex_begin = $dfsgraph->vertexIndex( $struct );
    	my $vertex_end   = $dfsgraph->vertexIndex( $struct + 1 ) - 1;
    	# my $offset  = $dfsgraph->vertexIndex( $struct ); #GET  THE offset
    	# my $vertex_loc_hash = $dfsgraph->{'vertex_loc'}; #vertex loc hash ref of dfsgraph
    	($isize, $osize, $xsize) = (0,0,0);
    	my $edge_hash = {};
    	foreach my $num ($vertex_begin .. $vertex_end) {
    	    my $ve = $dfsgraph->vertex( $num );
    	    foreach my $edge (@$ve) {
    		if (defined($edge_hash->{$edge})) {
    		    next;
    		} else {
    		    $edge_hash->{$edge}++;
    		    my ( $v1, $v2, $type ) = @{$edge};
    		    if ($type==$included) {
    			$isize++;
    		    } elsif ($type==$overlap) {
    			$osize++;
    		    } elsif ($type==$exclusive) {
    			$xsize++;
    		    } 
    		}
    	    }
    	}
    	my $total_vertex = $vertex_end - $vertex_begin + 1;
    	my $total_edges = keys %$edge_hash;
    	print "      Structure $struct: Vertex = $total_vertex; Edge = $total_edges; Included = $isize; Overlap = $osize; Exclusive = $xsize;\n";
    }
    return;
}

# end of statsEachStructure

#---------------------------------------------------
# dfsDiffGraphPrint
#
# given a dfsgraph object and structure number, this subroutine returns a scalar 
# with the structure in diffgraph format
#
# usage
#   my $text = $dfsgraph->dfsDiffGraphPrint($struct);
#-----------------------------------------------------------------------
sub dfsDiffGraphPrint {
    my ($self, $struct) = @_;
    
    my $vertex_begin = $self->vertexIndex( $struct );
    my $vertex_end   = $self->vertexIndex( $struct + 1 ) - 1;
    
    # get the offset
    my $offset  = $self->vertexIndex( $struct );
    my $scalar;
    my $vertex_loc_hash = $self->{'vertex_loc'};

    foreach my $num ($vertex_begin .. $vertex_end) {
    	my $v = $num-$offset;
    	my $stem = $vertex_loc_hash->{$num};
    	my ($l1, $l2, $r1, $r2, $lvienna, $rvienna) = $stem->getStemAsArray;
    	my ($b1, $b2) = ('[', ']');
    	my $center = ($l1+$l2+$r1+$r2)/4 ;
        # print  "$v $center $b1 $l1 $l2 $r1 $r2 $b2 $lvienna $rvienna\n";
    	if (!(defined $lvienna)) { $lvienna = ''; }
    	if (!(defined $rvienna)) { $rvienna = ''; }
    	$scalar .= sprintf( "%4d %5.1f [ %4d %4d %4d %4d ] %s %s\n", 
    	                    $v, $center, $l1, $l2, $r1, $r2, $lvienna, $rvienna);
        # $scalar .= $text;
        # $scalar .= "$v $center $b1 $l1 $l2 $r1 $r2 $b2 $lvienna $rvienna\n";
        # $stem->dump;
    }
    
    $scalar .= "\n";

    foreach my $num ($vertex_begin .. $vertex_end) {
        # $num is the vertex label in Dfsgraph object
    	my $v = $num-$offset; 
    	my $ve = $self->vertex( $num );
    	$scalar .= "$v:";
    	foreach my $edge (@$ve) {
    	    my ( $nv1, $nv2, $type ) = @{$edge};
    	    my ( $v1, $v2 ) =  ($nv1-$offset, $nv2-$offset); 
            # print "$v1, $v2, $type\n";
    	    next if ($v != $v1);
    	    next if ($type == 1);
    	    next if ($v2 < $v1);
            # next if (
    	    $scalar .= " ";
    	    if ($type == 0) {
                $scalar .= "$v2 i"; 
    	    } elsif ($type == 2) {
                $scalar .= "$v2 o"; 
    	    } elsif ($type == 3) {
                $scalar .= "$v2 x"; 
    	    }
    	}
    	$scalar .= "\n";
        # print "\n\n";
    } 
    return $scalar;
}
# end of dfsDiffGraphPrint

#-----------------------------------------------------------
# reduceEdge
#
# This function reduces the edges of a structure based on the embed level number.
# The structure is already a part of the dfsgraph data strucutre.
# 
# Returns the structure with reduced edges in scalar. 
# my $diffgraph = reduceEdge($dfsgraph, $struct, $embedlevel); #embedlevel = 3; 
# my $diffgraph = $dfsgraph->reduceEdge($struct, $embedlevel);
#
# usage
#
# sample code for doing the embed on each structure
# foreach my $i ( 0 .. $n_structures - 1 ) {
#     my $Ndiffgraph = reduceEdge($dfsgraph, $i, 3);
#     print $Ndiffgraph;
# }    
#-----------------------------------------------------------
sub reduceEdge {
    my ($dfsgraph, $struct, $embedlevel) = @_;
    
    my ($included, $overlap, $exclusive) = (0, 2, 3);
    my ($isize, $osize, $xsize) = (0,0,0);
    my $vertex_begin = $dfsgraph->vertexIndex( $struct );
    my $vertex_end   = $dfsgraph->vertexIndex( $struct + 1 ) - 1;
    my $offset  = $dfsgraph->vertexIndex( $struct ); 
    my $vertex_loc_hash = $dfsgraph->{'vertex_loc'}; # vertex loc hash ref of dfsgraph

    my $edge_hash = {};
    my ($ihash, $jhash, $ohash, $xhash) = ({}, {}, {}, {});
    my $new_vertex_loc = {};
   
    ######## Initializing #############
    my $jvert_hash = {};
    foreach my $vnum ($vertex_begin .. $vertex_end) {
        $jvert_hash->{$vnum} ={};
    }
    ###############

    foreach my $num ($vertex_begin .. $vertex_end) {
    	my $ve = $dfsgraph->vertex( $num );
    	my $v = $num-$offset; # will not need to use it
    	$new_vertex_loc->{$v} = $vertex_loc_hash->{$num};
    	foreach my $edge (@$ve) {
    	    if (defined($edge_hash->{$edge})) {
                next;
    	    } else {
                $edge_hash->{$edge}++;
                my ( $v1, $v2, $type ) = @{$edge};
                if ($type==$included) {
        		    $isize++;
        		    $ihash->{$v1}++;
        		    $jhash->{$v2}++;
        		    ########
        		    my $jvert_hash_ref = $jvert_hash->{$v2};
        		    #incrementing count
        		    $jvert_hash_ref->{'count'}++;
        		    ############################ get refercence for neighbor list
        		    my $extracount = 0;
        		    my $nbor_hash;
        		    if (exists $jvert_hash_ref->{'nbor'}) {
        			    $nbor_hash = $jvert_hash_ref->{'nbor'};
                    } else {
            			$nbor_hash = {};
            			$jvert_hash_ref->{'nbor'} = $nbor_hash;
                    }
    		    
                    # go through all edges in v1, to check it has x or o 
                    # relationship with the existing neighbors of i $v2
        		    my $ve = $dfsgraph->vertex($v1);
        		    # print "v1: $v1, v2:$v2\n";
        		    foreach my $edge (@$ve) {
        			my ( $ver1, $ver2, $vertype ) = @{$edge};
            			if ( $v1 == $ver2 ) {
            			   # my $tmp = $ver1;
            			    $ver2 = $ver1;
            			    $ver1 = $v1;
            			    if ($vertype == $included) {
                                $vertype = 1;
            			    }
            			}
                        # print "looping: $ver1, $ver2, $vertype\n";
                        # if ver2 is present in nbor hash... 
                        if (exists $nbor_hash->{$ver2}) {
        			        # print "      neighbor present: $ver2 ||| $nbor_hash->{$ver2}\n";
                            if ( ($vertype == $exclusive) || ($vertype == $overlap)) {
                                $extracount=1;
                                # print "---------------  EXIT vertex: $ver2\n";
                                last;
                            }
                        }
        		    }
        		    $nbor_hash->{$v1}++; #$v1 to nbor hash, list vertices that makes an i edge to $v2
                    # print "   neighbor hash populated: $nbor_hash\n";
        		    ############## taking care of extra count
        		    if (exists $jvert_hash_ref->{'extracount'}) {
                        $jvert_hash_ref->{'extracount'} = $jvert_hash_ref->{'extracount'} + $extracount;
        		    } else {
                        $jvert_hash_ref->{'extracount'} = $extracount;
        		    }
        		    #############
        		} elsif ($type==$overlap) {
        		    #print "0";
        		    $osize++;
        		    $ohash->{$v1}++;
        		    $ohash->{$v2}++;
        		} elsif ($type==$exclusive) {
        		    $xsize++;
        		    $xhash->{$v1}++;
        		    $xhash->{$v2}++;
        		}
    	    }
    	}           # end of loop over edges
    }               # end of loop over vertices
    
#--------------------------- code for printing the jhash structure 
#   print "\n\nFinal hash struct\n";
#   foreach my $vert (keys %$jvert_hash) {
#	    print  "FOR vertex: $vert\n";
#       foreach my $nkey (keys %{$jvert_hash->{$vert}}) {
#           print "   KEY: $nkey; val: $jvert_hash->{$vert}->{$nkey}\n";
#           if ($nkey eq 'nbor') {
#               print "       ";
#               my $nbor_hash = $jvert_hash->{$vert}->{$nkey};
#               foreach my $kk (keys %{$jvert_hash->{$vert}->{$nkey}}) {
#                   print "v: $kk,"; 
#               }
#               print "\n";
#           } 
#       }
#	
#   }
#   print "\n\n";
#----------------------------------------


    # sort the vertex based on number of highest number of i-edges
    $edge_hash = {};
    my $store_edges = [];
    #my $new_vertex_loc = {}; 
    #   foreach my $vnum (sort { $ihash->{$b}<=>$ihash->{$a} }  keys %$ihash) {
    #	my $val = $ihash->{$vnum};
    #	my $jval = $jhash->{$vnum};
    #	print "Vertex: $vnum; ihash Value= $val | jhash val: $jval\n";
    #    }
    
    foreach my $vnum (sort { $ihash->{$b}<=>$ihash->{$a} }  keys %$ihash) {
        # print "looping over vertex: $vnum\n";
    	my $ve = $dfsgraph->vertex( $vnum );
    	my $v = $vnum-$offset;
    	# $new_vertex_loc->{$v} = $vertex_loc_hash->{$vnum};
    	
    	foreach my $edge (@$ve) {	
    	    if (defined($edge_hash->{$edge})) {
    		    next;
    	    } else {
        		$edge_hash->{$edge}++;
        		my ( $v1, $v2, $type ) = @{$edge};
        		if ( ($type==$exclusive) ||  ($type==$overlap) ) {
        		    my ($nv1, $nv2, $nt) = ($v1-$offset, $v2-$offset, $type); 
        		    push @$store_edges, [$nv1, $nv2, $nt];
                    # print "   Edge stored: $nv1, $nv2, $nt\n";
        		} else {
        		    # find the jcount for v2
        		    my $v2_jcount = $jhash->{$v2};
                    
                    # print "for vertex... v1: $v1, v2: $v2 | count: $jvert_hash->{$v2}->{'count'},  extracount: $jvert_hash->{$v2}->{'extracount'}\n";
        		    my $threshold_val = $jvert_hash->{$v2}->{'count'} - $jvert_hash->{$v2}->{'extracount'};
                    # print "jcount: $jvert_hash->{$v2}->{'count'} <= threshold: $threshold_val\n";
        		    
        		    if ($threshold_val <= $embedlevel) { 
            			#keep the edge
            			my ($nv1, $nv2, $nt) = ($v1-$offset, $v2-$offset, $type); 
            			push @$store_edges, [$nv1, $nv2, $nt];
                        # print "   Edge stored: $nv1, $nv2, $nt\n";
        		    } 
        		  
        		    # check to see if $v1 has any o or x relationship with the 
        		    # existing neighbor of $v2, if yes decrement the extra count 
        		    # and remove $v1 from neighborlist
        		    if ($jvert_hash->{$v2}->{'extracount'} > 0) {
            			my $nbor_hash = $jvert_hash->{$v2}->{'nbor'};
            			my $ve = $dfsgraph->vertex($v1);
            			# print "v1: $v1, v2:$v2\n";
            			foreach my $edge (@$ve) {
            			    my ( $ver1, $ver2, $vertype ) = @{$edge};
            			    if ($v1 == $ver2) {
                				# my $tmp = $ver1;
                				$ver2 = $ver1;
                				$ver1 = $v1;
                				if ($vertype == $included) {
                				    $vertype = 1;
                				}
            			    }
            			    next if ($v2 == $ver2);
            			    
            			    #if ver2 is present in nbor hash... 
            			    if ( exists $nbor_hash->{$ver2} ) {
                                # print "      neighbor present: $ver2 ||| $nbor_hash->{$ver2}\n";
                				if ( ($vertype == $exclusive) || ($vertype == $overlap)) {
                				    $jvert_hash->{$v2}->{'extracount'}--;
                				    delete $nbor_hash->{$v1}; #delete v1 from neghbor hash
                				    # $extracount=1;
                				    # print "---------------  EXIT vertex: $ver2\n";
                				    last;
                				}
            			    }
                        }
                    } # END of if ($jvert_hash->{$v2}->{'extracount'} > 0)
    
                    $jvert_hash->{$v2}->{'count'}--;
                } # END Included edges
    	    }     # END, else defined edge
    	}         # END foreach edge
    }             # END vertex sorted by i edge
    
    # What about other vertices that that had no i edges
    foreach my $num ($vertex_begin .. $vertex_end) {
    	my $ve = $dfsgraph->vertex( $num );
        # my $v = $num-$offset; # will not need to use it
    	next if (defined($ihash->{$num}));
    	
        # $new_vertex_loc->{$v} = $vertex_loc_hash->{$num};
    	foreach my $edge (@$ve) {
    	    if (defined($edge_hash->{$edge})) {
                next;
    	    } else {
        		$edge_hash->{$edge}++;
        		my ( $v1, $v2, $type ) = @{$edge};
        		if ( ($type==$exclusive) ||  ($type==$overlap) ) {
        		    my ($nv1, $nv2, $nt) = ($v1-$offset, $v2-$offset, $type); 
        		    push @$store_edges, [$nv1, $nv2, $nt];
        		}
    	    }
    	}
    }
    
    #### code for the printing the reduced edgelist
    # foreach my $edge (@$store_edges) {
    #     print "@$edge\n";
    # }
    #########################
    my $total_vertex = $vertex_end - $vertex_begin + 1;     
    my $newedgelist = [$store_edges, $total_vertex];
    my $description = $dfsgraph->description;
    my $descript = $description->[$struct];
    
    my $new_dfs = new DfsGraph( {'type' => 'edgedata', 
                                 'edgelist'=> $newedgelist, 
                                 'description'=> $descript, 
                                 'vertexloc' =>  $new_vertex_loc 
                                } );
    my $Ndiffgraph = $new_dfs->dfsDiffGraphPrint(0); 
  
    # ($allsize, $isize, $osize, $xsize, $vertex_num) = basicDfsGraphStatsOne($new_dfs,1);
    # print "\n\ntotal: $allsize, i=$isize, o=$osize, x=$xsize, vertex_num = $vertex_num\n";
    
    #return a dfsgraph with a reduced set of edges
    return $Ndiffgraph;
}

# end of reduceEdge

#-----------------------------------------------------------------------------
# dfsEdge
#
# returns the contents of the nth edge in the current DFS structure. Edges are 
# numbered starting at zero
# 
# USAGE
#   $edge_ref = $self->dfsEdge( $n );
#-----------------------------------------------------------------------------
sub dfsEdge{
    my ( $self, $n ) = @_;

    return $self->{dfs}[$n];
}

# end of dfsEdge

#-----------------------------------------------------------------------------
# dfsAddEdge
#
# adds the edge to the end of the DFS.  The edge should be a reference to a
# hash with the keys v1, v2, type. Status is true if the provided reference
# has the correct hash keys. Status is the number of edges in the DFS if
# true.
#
# USAGE
#   $status = $self->dfsAddEdge( $next_edge_ptr );
#-----------------------------------------------------------------------------
sub dfsAddEdge{
    my ( $self, $e ) = @_;

    my $status = 0;

    if ( defined $e and
         defined $e->{v1} and
         defined $e->{v2} and
         defined $e->{type} ) 
    {

        my $n = $self->{dfsLen};
        $self->{dfs}[$n] = $e;
        $n++;
        # $self->dfsLen( $n );
        $self->{'dfsLen'} = $n;
		$status = $n;
        my $vertex = $self->{dfs_max_vertex} || $e->{v1};
        if ( $e->{v1} > $vertex ) { $vertex = $e->{v1} };
        if ( $e->{v2} > $vertex ) { $vertex = $e->{v2} };
        $self->{dfs_max_vertex}  = $vertex;
    }

    return $status;
}

# end of dfsAddEdge

#-----------------------------------------------------------------------------
# dfsLastEdge
#
# Returns a reference to a hash with the last edge that was tried in the 
# search tree.  Status is $SUCCESS if the edge was successfully added, and
# $FAILURE if not.  Success and Failure are determined by the support.
#-----------------------------------------------------------------------------
sub dfsLastEdge{
    my ( $self ) = @_;

    my $n = $self->{dfsLen};
    my $e = $self->{dfs}[$n];

    return $e;
}

# end of dfsLastEdge

#-----------------------------------------------------------------------------
# AUTOLOAD
#
# This function will not create new hash keys.  This keeps typos from creating
# crazy new attributes for objects.  Anything that is not a hash key or 
# a function in this package fails. Many functions act on arrays where the 
# hash element in DfsGraph is a reference to the array.  These functions, 
# listed in the hash %array_function are dealt with separately.  For these
# functions the usage is arrayFunction( index, value )
#-----------------------------------------------------------------------------
sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;
    my %array_function = ( vertex_index=>1, vertexIndex=>1, 
                           edge_index=>1, edgeIndex=>1, 
                           edge=>1, vertex=>1 );
  
    my $unknown = $AUTOLOAD;                    # The name of the unknown function
    $unknown =~ s/.*:://;                       # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;        # Skip all-caps methods like DESTROY

    if ( defined $array_function{ $unknown } ) {

        # convert a function call such as $dfs->thisFunction to act on
        # $dfs->{ this_function )

        $unknown =~ s/([A-Z])/_$1/g;            
        $unknown =~ tr/A-Z/a-z/;

        if ( @_ ) {                     # if a structure number is given
            my $i = shift;                      # this only deals with the first structure
            my $j = shift;
            # cases
            # $i is int, $j is array ref - save of array element $i
            # $i is int, $j undef - get of array element $i
            # $i is array ref - this is a save of the entire array
            if ( defined $j ) {
                # print STDERR "DfsGraph:AUTOLOAD   setting $unknown (array function) - $i:$j\n";
                $self->{ $unknown }[ $i ] = $j;
            } else {
                # if no array element is given assume it is the ref to the whole array
                if ( ref $i eq 'ARRAY' ) {
                    # print STDERR "DfsGraph:AUTOLOAD   setting $unknown (array function) - $i\n";
                    $self->{ $unknown } = $i;
                } else {
                    # print STDERR "DfsGraph:AUTOLOAD   retrieving $unknown (array function) - $i (no array index)\n";
                    return  $self->{ $unknown }[ $i ];
                }
            }
            # return the array reference, consistent with the get call below
            return $self->{ $unknown }; 

        } else {                                # return whole array
            # print STDERR "DfsGraph:AUTOLOAD   array function - retrieve all\n";
            return $self->{ $unknown };
        }
    } 
  
    if ( exists $self->{ $unknown } ) {         # do not create new fields
        $self->{ $unknown } = shift if @_;      # Set new value if one was supplied
        return( $self->{ $unknown } );          # Return current value
    }
        
    die "DfsGraph:AUTOLOAD   unknown function $unknown\n\n";
}

# end of AUTOLOAD

#-----------------------------------------------------------------------------
# $Log: DfsGraph.pm,v $
# Revision 1.26.4.1  2015/10/13 17:52:55  huang147
# this is the working version.
#
# Revision 1.26.6.6  2013/02/08 16:35:19  gribskov
# Added our to lexsort function in addVE to eliminate warnings
#
# Revision 1.26.6.5  2012/10/14 15:12:47  gribskov
# Fixed sortVE and tested.  lexicographic sort looks OK.
#
# Revision 1.26.6.4  2012/10/12 18:06:48  gribskov
# Updated documentation.  No functional change.
#
# Revision 1.26.6.3  2012/10/12 16:08:18  gribskov
# Fixed error in %EDGENUM hash.  Compiles and appears to load XIOS file OK.
#
# Revision 1.26.6.2  2012/10/12 15:56:33  gribskov
# Added addXios function to add a structure in XIOS format.
# Modified addDiffgraph to use global definition if edge types for consistency.
# Add edgeTypeToNumber to convert consistently between edge letter and number.
#
# Revision 1.26.6.1  2012/10/12 00:19:29  gribskov
# Updated documentation and fixed indentation.  No functional changes.
#
# Revision 1.26  2012/04/29 06:31:54  rahmanr
#
# updated the reduce edge code to fix a bug
#
# Revision 1.25  2012/04/25 02:40:01  rahmanr
#
# Added subroutine addXiosText to added xios data to DfsGraph Object
#
# Revision 1.24  2012/02/07 18:47:34  rahmanr
#
# Changed the initialization of dfsgraph description attribute
# Added new subroutine reduceEdge to reduce i edges and modified dfsDiffGraphPrint
# Also fixed bugs regarding the use vertex_loc attribute
#
# Revision 1.23  2011/11/16 20:09:20  gribskov
# brought current version in this directory up to date with repository.
# should be no functional chnage.
#
# Revision 1.22  2011/08/05 00:14:33  rahmanr
# updated documentation
#
# Revision 1.21  2011/08/04 08:53:55  rahmanr
#
# Updated the way edgeIndex works, I believe it works properly
#
# Revision 1.20  2011/08/02 23:17:33  rahmanr
#
# Added "use warnings;" to the package.
#
# Revision 1.19  2011/08/02 22:54:15  rahmanr
#
# Fixed the part that broke my code, added documentations, and removed unnecessary text
#
# Revision 1.18  2011/07/29 16:56:48  gribskov
# Fixed the behavior of the autoload accessor function for array functions.
# Hard to see how it was working at all before.  Use -w people!
#
# Revision 1.17  2010/06/23 11:57:32  rahmanr
#
# Added a subroutine to print Dfsgraph as diffgraph adjaceny list
#
# Revision 1.16  2010/06/15 20:42:28  rahmanr
#
# Couple of subroutines added to map nodes to actual stem coordinates
#
# Revision 1.15  2010/05/12 14:45:58  gribskov
# Added $JNCLUDED as an exported edge type (meaning a reverse I edge)
#
# Revision 1.14  2010/04/17 18:10:00  gribskov
# Updataed documentation, no functional changes.
#
# Revision 1.13  2010/04/05 06:57:39  rahmanr
#
# fixing the bug that i introduced in the last update
#
# Revision 1.10  2009/05/03 09:39:03  gribskov
# Added input for diffgraph files
#
# Revision 1.9  2008/04/04 20:15:26  likejie
#
# 1.8 locally modified
#
# Revision 1.7  2007/12/04 20:49:07  rahmanr
#
# Minor bugs
#
# Revision 1.6  2007/10/29 16:40:20  gupta31
# Added a line to sort half-stems
#
# Revision 1.5  2007/10/24 18:10:38  gribskov
# merged from dfs-0-2mrg
#
# Revision 1.4.2.8.2.8  2007/10/21 20:23:32  gribskov
# VE list works now, sorting function for VE list works.
#
# Revision 1.4.2.8.2.7  2007/10/21 20:04:32  gribskov
# still don't have the right vertex-edge structure.  needs to be an array
# of edges.  committed this version before ripping it apart again.
#
# Revision 1.4.2.8.2.6  2007/10/21 18:33:07  gribskov
# fixed assignment of edges to vertices in addEdgeToVE.  Updated more calls
# to methods instead of direct variable access.
#
# Revision 1.4.2.8.2.5  2007/10/21 12:26:47  gribskov
# Corrected edgeList to always return array of edges instead of edge indices.
#
# Revision 1.4.2.8.2.4  2007/10/21 11:56:47  gribskov
# Converted addStructure to use class methods instead of direct access to
# class structure.
#
# Revision 1.4.2.8.2.3  2007/10/21 11:31:56  gribskov
# Removed commented out methods and functions that may be no longer needed.
#
# Revision 1.4.2.8.2.2  2007/10/21 11:27:32  gribskov
# Converted addEdge to use class methods instead of direct access to
# class structure.
#
# Revision 1.4.2.8.2.1  2007/10/21 03:35:58  gribskov
# Added vertex-edge list and appropriate dump.
# Commented out vertex function, now autoloaded
# Compiles, runs, and tests
#
# Revision 1.4.2.8  2007/10/20 13:48:03  gribskov
# Commented out duplicate function dfsLastEdge.  Notu sue which one to keep;
# I commented out what looked like the newer one.  Removed a few more
# redeclarations that give compile warnings.   Clean compile with perl -cw
#
# Revision 1.4.2.7  2007/10/20 13:42:19  gribskov
# Fixed a few minor compiler warnings.  Should make no difference.
#
# Revision 1.4.2.6  2007/10/20 13:30:07  gribskov
# Added export of global variable
#
# Revision 1.4.2.5  2007/10/20 12:47:52  gribskov
# Added edgeList function to get the edges for specified structures.  This
# version compiles.
#
# Revision 1.4.2.4  2007/10/20 12:39:46  gribskov
# Added printing string to dump.
#
# Revision 1.4.2.3  2007/10/19 14:31:37  gribskov
# Added AUTOLOAD function for accessors.  includes accessors for array
# elements.  Added Vertex function to get lists of vertices for all or
# specified structures.
#
# Revision 1.4.2.2  2007/10/19 13:47:56  gribskov
# Moved log to end of file, no functional changes.
#
# Revision 1.4.2.1  2007/10/19 13:37:56  gribskov
# Fixed bug in reading pair graphs as strings.
# Graph string is now stored in field 'description'.
# Updated dump to display description.
#
# Revision 1.4  2007/10/18 20:18:50  likejie
#         All my modifications
#
# Revision 1.3  2007/10/08 21:15:12  gribskov
# Partially debugged.  I doubt this works
#
# Revision 1.2  2007/08/09 04:10:01  gribskov
# Added missing functions (all of them??).  COmpiles with perl -c
#
# Revision 1.1  2007/08/09 03:06:39  gribskov
# Initial version.  This version has undefined methods and variables and
# generally does not work.  The gernaral outline should be right for
# generating the DFS search tree.
#-----------------------------------------------------------------------------

# end of dfsGraph

1;
