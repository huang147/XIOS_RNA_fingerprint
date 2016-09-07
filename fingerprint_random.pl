#!/usr/bin/perl -w
#------------------------------------------------------------------------------
# $Id: fingerprint_random.pl,v 1.1.4.9 2016/08/24 05:28:05 huang147 Exp $
# 
# Calculate a fingerprint using a random sampling approach.
#
# Michael Gribskov     Dec 19, 2012
#
# Revision log at end of file
#------------------------------------------------------------------------------
use strict;
use Data::Dumper;
use Getopt::Std;
use Time::HiRes qw( gettimeofday tv_interval );
use List::Util qw( min max );
use Cwd qw( abs_path );
                      
#use lib '/home/gribskov/RNA/src/perl_src/';
#use lib '/home/huang147/RNA/src/perl_src/';
#use lib '/home/huang147/mrg121003/RNA/src/perl_src/';
use lib './';
use DfsGraph;
use Dfsgenerator;
use Motif_db;
use Fingerprint;
use Topology;

my $DEFAULT_N             = 50000;
my $DEFAULT_SUBGRAPH_SIZE = 7;
my $DEFAULT_SAMPLE_LIMIT  = $DEFAULT_N;
my $DEFAULT_MOTIF_DB      = 'db.storable';
#my $DEFAULT_MOTIF_DB      = '/cbio/proj/rna/data/Motif_fingerprint/2_to_7_stems_topologies.removed_not_true.mini_dfs.txt.removed_redundant.with_label.motif.storable';
my $DEFAULT_QUIET         = 0;
my $DEFAULT_VERBOSE       = 0;
my $DEFAULT_WORKDIR       = $ENV{'PBS_O_WORKDIR'};

my $id = '$Id: fingerprint_random.pl,v 1.1.4.9 2016/08/24 05:28:05 huang147 Exp $';
my ( $program, $version ) = $id =~ /Id: ([^,]+),v ([^ ]+)/;

my $USAGE = qq{$program $version <options> <xios_file>
    -h     this usage message
    -d     <filename>, read in the storable motif database
    -g     <int>, graph size to sample (default=$DEFAULT_SUBGRAPH_SIZE)
    -l     <int>, stop after every motif has this many hits (default=number or random samples)
    -n     <int>, maximum number of random samples (default=$DEFAULT_N)
    -q     quiet (default=$DEFAULT_QUIET)
    -v     show DFS code for each motif (verbose, default=$DEFAULT_VERBOSE)
    -w     working directory
    -x     create decoy directory for holding sampled subgraphs as .xios files  
};

# command line options

my $option = {};
getopts( 'd:g:hl:n:qvw:x', $option );

# help
if ( $$option{h} ) {
    print "$USAGE\n";
    exit 1;
}

#motif database
my $motif_fname = $DEFAULT_MOTIF_DB;
if ( $$option{d} ) {
    $motif_fname = $$option{d};
}

my $motifdb = new Motif_db( { db=>$motif_fname } );
my $n_motif   = $motifdb->n;
if ( $n_motif < 1 ) {
    print STDERR "No motifs read from motif database ($motif_fname)\n\n";
    exit 1;
}

# subgraph size
my $subgraph_size = $DEFAULT_SUBGRAPH_SIZE;
if ( $$option{g} ) {
    $subgraph_size = $$option{g};
}

# number of iterations
my $iterations = $DEFAULT_N;
if ( $$option{n} ) {
    $iterations = $$option{n};
}

# sample limit
my $limit = $iterations;
if ( $$option{l} ) {
    $limit = $$option{l};
}

my $quiet = $DEFAULT_QUIET;
if ( $$option{q} ) {
    $quiet = 1;
}

my $verbose = $DEFAULT_VERBOSE;
if ( $$option{v} ) {
    $verbose = 1;
}

my $workdir = $DEFAULT_WORKDIR;
if ( $$option{w} ) {
    $workdir = $$option{w};
}
$workdir 	= abs_path( $workdir ).'/';


my $decoydir;
if ( $$option{x} ) {
    $decoydir = $workdir."/decoy/";
    system( "mkdir $decoydir" ) unless ( -e $decoydir );
}


#------------------------------------------------------------------------------
# main program
#------------------------------------------------------------------------------

my $xios_file = $ARGV[0];

unless ( $quiet ) {
    print STDERR "# v$version $program \n";
    print STDERR "# Graph file:     $xios_file\n";
    print STDERR "# Motif database: $motif_fname\n";
    print STDERR "#       n motifs: $n_motif\n";
    print STDERR "# Subgraph size:  $subgraph_size\n";
    print STDERR "# Iterations:     $iterations\n";
    print STDERR "# Sample limit:   $limit\n";
}

my $t0 = [gettimeofday];
my $graph = DfsGraph->new( {type=>'xios', file=>$xios_file} );

# setup fingerprint
my $fingerprint = Fingerprint->new;
$fingerprint->type( 'random' );
$fingerprint->program( "$program v$version" );
$fingerprint->queryId( $xios_file );
$fingerprint->databaseId( $motif_fname );
$fingerprint->queryEdge( $graph->n_edge );
$fingerprint->queryVertex( $graph->n_vertex );

my %count;
my %first;
my %seen_before;
my $mincount = 0;
my $minmapping = 0;
my %mapcount;
my $minindex;
my $iter = 0;
while ( $iter < $iterations ) {
    $iter++;
    my ($vstring, $subgraph ) = sampleByVertex( $graph, $subgraph_size );


    my $index;
    if ( defined $seen_before{$vstring} ){
        $index = $seen_before{$vstring};
    } else {
        my $dfs_gen_obj = new Dfsgenerator ( {'dfsgraph'=>$subgraph, 'struct_num' => 0} );
        my @dfscode = @{ $dfs_gen_obj->getDfscode };
        $index = $motifdb->encodeDfsRowArray( @dfscode ); 
        $seen_before{$vstring} = $index;       
	$mapcount{ $index }++;
    }
    #$minmapping = min( values %mapcount );


    unless ( defined $count{$index} ) {
#        print "   new: $iter     $index\n";
        $first{$index} = $iter;
        $mincount = 0;
        $minindex = $index;
    }
    $count{$index}++;
    if ( $index eq $minindex ) {
        $mincount++;
    } else {
        if ( $count{$index} < $mincount ) { 
            $mincount = $count{$index}; 
            $minindex = $index;
        }
    }
    
    # check to see that all subgraphs have reached the minimum count
    if ( $mincount >= $limit ) {
        my $m;
        foreach my $k ( keys %count ) {
            if ( $count{$k} < $limit ) {
                $m = $k;
                last;
            }
            
        }
        # if $m is defined, at least one count is low
        if ( $m ) {
            $mincount = $count{$m};
            $minindex = $m;  
        } else {
            last;
        }
    }
}

my $t1 = [gettimeofday];
my $elapsed = tv_interval ( $t0, $t1 );
$fingerprint->timeElapsed( $elapsed );

# TODO need to add to fingerprint class
my @mappings = keys %seen_before;
my $mappings = @mappings;
$fingerprint->totalMapping( $mappings );

#my %mapcount;
#foreach my $vstr ( keys %seen_before ) {
#    $mapcount{ $seen_before{$vstr} }++;
#}
my @countkeys = keys %count;
my $motifcount = @countkeys;
$fingerprint->iteration( $iter );

my $topology = new Topology;  
my $nstem = $topology->XIOSRead( $xios_file ); 
my $stemlist = $topology->stemList;
foreach my $k ( sort { $count{$b} <=> $count{$a} } keys %count ) {
    if ( $decoydir ) { 
    	my $id = $motifdb->lookupMotif( $k );
    	my $xios_motif = $decoydir.$id.".decoy";
    	my $vstr = pop @mappings; 
    	my @vstr = split /\_/, $vstr;
    	my $motif_topology = new Topology; 
    	foreach my $v ( @vstr ) {
	    $motif_topology->addStem( $$stemlist[$v] );
    	}
    	$motif_topology->XIOSWrite( $xios_motif ) unless ( -s $xios_motif );
    }
    my $enc = $k;
    $fingerprint->addMotif( {   id             => $motifdb->lookupMotif( $k ),
                                count          => $count{$k},
                                first_observed => $first{$k},
                                encoded_dfs    => $k,
                                mapping        => $mapcount{$k},
                            }
                           );
}

print $fingerprint->asXML;

exit 0;

#------------------------------------------------------------------------------
# Sampling Functions:
# sampling functions only differ by how the next vertex is selected.  the 
# sampling portion of the code is delimited by #-+-+...
#
# all sampling functions return a subgraph as a DfsGraph object and a string with
# the sampled vertices concatenated as v1_v2_v3 ... The returned subgraph has 
# vertices labeled 0 to n_vertex-1.
#
# obviously this could be handled using a subroutine, but since the sampling is
# called many times, its better to leave use separate functions
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# sampleByVertex
#
# sample a random subgraph from a graph.  Each vertex that has a neighbor in the
# currently selected graph has equal probability.
#
# The returned subgraph has vertices labeled 0 to n_vertex-1.
#
# USAGE
#    my ( $subgraph, $vertex_string )  = sampleByVertex( $graph, $n_vertex );
#------------------------------------------------------------------------------
sub sampleByVertex{
	my ( $graph, $n_target ) = @_;
	
	my @v  = @{ $graph->vertex };

	restart:
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	# select an initial random vertex, assumes there is only one graph present
	my $nv = @v;
	my $r  = int( rand($nv) );
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-	
		
	# build a hash of neighbor vertices
	# select an additional vertex and update neighbor hash until the required 
	# number of vertices are selected
	
    my %neighbor = ();
    my %selected = ();
    my %excluded = ();
    my $n_selected = 0;
    
    my @test = keys %selected;
#    if ( @test ) {
#        print "keys retained:@test\n";
#    }
    
    $neighbor{$r}++;
	while ( $n_selected < $n_target && $nv) {
	    # $nv checks to make sure there are more vertices available in the neighbor list
	    $selected{$r} = $n_selected;
	    $n_selected++;
	    # remove selected vertex from neighbor list
        delete $neighbor{$r};
	    
	    # add any new neighbors of selected edge
	    foreach my $edge ( @{$v[$r]} ) {
	        if ( $$edge[2] == 3 ) {
	            # x edges.  any neighbor of the new edge related by an x edge 
	            # cannot be added to the sampled graph
	            $excluded{$$edge[0]}++;
                $excluded{$$edge[1]}++;
                # also remove from neighbor list
                delete $neighbor{$$edge[0]};
                delete $neighbor{$$edge[1]};
	        }
	        foreach my $i ( 0 .. 1 ) {
	            # examine v1 and v2 of the edge
	            next if ( defined $selected{$$edge[$i]} );
	            next if ( defined $excluded{$$edge[$i]} );
	            $neighbor{$$edge[$i]}++;
	        }	        
	    }
	    
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	    # select new random vertex
	    my @nlist = keys %neighbor;
	    $nv = @nlist;
	    $r  = $nlist[ int(rand($nv)) ];	    
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	}
	my @vlist = sort {$a<=>$b} keys %selected;
	my $vstring = join "_", @vlist;
	
	# if there are vertices only connected by x edges, you can end up with a
	# one vertex, no edge graph, in this case do over
	if ( $n_selected == 1 ) { goto restart }
		
	# add all edges between selected vertices
		
	my $subgraph = new DfsGraph;
	$subgraph->addStructure();
	$subgraph->addVertices( $n_target );
	
    my @alledges = @{ $graph->edge };
    foreach my $edge ( @alledges ) {
       my ( $v1, $v2, $e ) = @{$edge};
       if ( defined $selected{$v1} && defined $selected{$v2} ) {
           my $this_edge = $subgraph->addEdge( $selected{$v1}, $selected{$v2}, $e ); 
           $subgraph->addEdgeToVE( $selected{$v1}, $selected{$v2}, $this_edge );
       }
   }
	    		
	return ( $vstring, $subgraph );
}

# End of sampleByVertex

#------------------------------------------------------------------------------
# sampleByHistory
#
# sample a random subgraph, weighted by the number of times a vertex has 
# previously been sampled. maximum weight on a vertex is 1 / (nv-1). Minimum 
# weight is zero (check this)
#
# The returned subgraph has vertices labeled 0 to n_vertex-1.
#
# USAGE
#    my ( $subgraph, $vertex_string ) = sampleByHistory( $graph, $n_vertex );
#------------------------------------------------------------------------------
{
	# define vcount here so it is retained between calls
	my @vcount;                    # number of times a vertex has been sampled
    
sub sampleByHistory{
	my ( $graph, $n_target ) = @_;

	my @v  = @{ $graph->vertex };
	
	restart:
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	# select a random vertex, assume there is only one graph present
	my $nv = @v;
	my $vsum = 0;                  # total vertex counts
	
    unless ( @vcount ) {
        # initialize vertex count if needed
        foreach my $v ( 0 .. $nv-1 ) {
            $vcount[$v] = 0;
        }
    }
    	
    foreach my $vc ( @vcount ) {
        $vsum += $vc + 1;         
    }   

	# sample based on an inverted vertex count. The inverted count is the sum of 
	# the counts minus the vertex count.  total * (nv-1) is the sum of the 
	# inverted counts
	my $c  = int( rand($vsum*($nv-1)) );
	my $r = 0;                     # index of selected vertex
	my $count = 0;
	foreach my $vertex ( 0 .. $nv-1 ) {
	    $count += $vsum - $vcount[$vertex] -1;
	    if ( $c <= $count ) {
	        last;
	    }
	    $r++;	    
	}
	
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-		
	# build a hash of neighbor vertices
	# select an additional vertex and update neighbor hash until the required 
	# number of vertices are selected
	
    my %neighbor = ();
    my %selected = ();
    my %excluded = ();
    my $n_selected = 0;
    
    $neighbor{$r}++;
	while ( $n_selected < $n_target && $nv) {
	    # $nv checks to make sure there are more vertices available in the neighbor list
	    $selected{$r} = $n_selected;
	    $n_selected++;
	    # remove selected vertex from neighbor list
        delete $neighbor{$r};
	    
	    # add any new neighbors of selected edge
	    foreach my $edge ( @{$v[$r]} ) {
	        if ( $$edge[2] == 3 ) {
	            # x edges.  any neighbor of the new edge related by an x edge 
	            # cannot be added to the sampled graph
	            $excluded{$$edge[0]}++;
                $excluded{$$edge[1]}++;
                # also remove from neighbor list
                delete $neighbor{$$edge[0]};
                delete $neighbor{$$edge[1]};
	        }
	        foreach my $i ( 0 .. 1 ) {
	            # examine v1 and v2 of the edge
	            next if ( defined $selected{$$edge[$i]} );
	            next if ( defined $excluded{$$edge[$i]} );
	            $neighbor{$$edge[$i]}++;
	        }	        
	    }	    
	    
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	    $vcount[$r]++;                 # update vertex count
	    # select new random vertex
	    my @nlist = keys %neighbor;
	    $nv = @nlist;
    	my $vsum = 0;                  # total vertex counts	
    	foreach my $v ( @nlist ) {
	       $vsum += $vcount[$v] + 1;         
        }		
	    
    	$c  = int( rand($vsum*($nv-1)) );
    	$count = 0;
    	foreach my $vertex ( @nlist ) {
    	    $count += $vsum - $vcount[$vertex] - 1;
    	    if ( $c <= $count ) {
    	        $r = $vertex;
    	        last;
    	    }
    	}
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	}
	my @vlist = sort {$a<=>$b} keys %selected;
	my $vstring = join "_", @vlist;
	
	# if there are vertices only connected by x edges, you can end up with a
	# one vertex, no edge graph, in this case do over
	if ( $n_selected == 1 ) { goto restart }
	
	# add all edges between selected vertices
	
	my $subgraph = new DfsGraph;
	$subgraph->addStructure();
	$subgraph->addVertices( $n_target );
	
    my @alledges = @{ $graph->edge };
    foreach my $edge ( @alledges ) {
       my ( $v1, $v2, $e ) = @{$edge};
       if ( defined $selected{$v1} && defined $selected{$v2} ) {
           #next unless ( $selected{$v1} < $selected{$v2} );
           
           my $this_edge = $subgraph->addEdge( $selected{$v1}, $selected{$v2}, $e );
           $subgraph->addEdgeToVE( $selected{$v1}, $selected{$v2}, $this_edge );
       }
   }
	    		
	return ( $vstring, $subgraph );
}

# End of sampleByHistory
}

#------------------------------------------------------------------------------
# sampleByDegree
#
# sample a random subgraph, weighted by the number of times a vertex has 
# previously been sampled. maximum weight on a vertex is 1 / (nv-1). Minimum 
# weight is zero (check this)
#
# The returned subgraph has vertices labeled 0 to n_vertex-1.
#
# USAGE
#    my ( $subgraph, $vertex_string ) = sampleByDegree( $graph, $n_vertex );
#------------------------------------------------------------------------------
{
	# define vcount here so it is retained between calls
	my @vcount;                    # number of times a vertex has been sampled
    
sub sampleByDegree{
	my ( $graph, $n_target ) = @_;

	my @v  = @{ $graph->vertex };
	
	restart:
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	# select a random vertex, assume there is only one graph present
	my $nv = @v;
	my $vsum = 0;                  # total vertex counts	
    unless ( @vcount ) {
        # initialize degree count, degree never changes
        foreach my $v ( 0 .. $nv-1 ) {
            $vcount[$v] = @{$v[$v]};
            $vsum += $vcount[$v];
        }
    }
    	
	# sample based on an inverted vertex count. The inverted count is the sum of 
	# the counts minus the vertex count.  total * (nv-1) is the sum of the 
	# inverted counts
	my $c  = int( rand($vsum) );
	my $r = 0;                     # index of selected vertex
	my $count = 0;
	foreach my $vertex ( 0 .. $nv-1 ) {
	    $count += $vcount[$vertex];
	    $r = $vertex;
	    if ( $c <= $count ) {
	        last;
	    }
	}	
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-		

	# build a hash of neighbor vertices
	# select an additional vertex and update neighbor hash until the required 
	# number of vertices are selected
	
    my %neighbor = ();
    my %selected = ();
    my %excluded = ();
    my $n_selected = 0;
    
    $neighbor{$r}++;
	while ( $n_selected < $n_target && $nv) {
	    # $nv checks to make sure there are more vertices available in the neighbor list
	    $selected{$r} = $n_selected;
	    $n_selected++;
	    # remove selected vertex from neighbor list
        delete $neighbor{$r};
	    
	    # add any new neighbors of selected edge
	    foreach my $edge ( @{$v[$r]} ) {
	        if ( $$edge[2] == 3 ) {
	            # x edges.  any neighbor of the new edge related by an x edge 
	            # cannot be added to the sampled graph
	            $excluded{$$edge[0]}++;
                $excluded{$$edge[1]}++;
                # also remove from neighbor list
                delete $neighbor{$$edge[0]};
                delete $neighbor{$$edge[1]};
	        }
	        foreach my $i ( 0 .. 1 ) {
	            # examine v1 and v2 of the edge
	            next if ( defined $selected{$$edge[$i]} );
	            next if ( defined $excluded{$$edge[$i]} );
	            $neighbor{$$edge[$i]}++;
	        }	        
	    }	    
	    
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	    # select new random vertex
        my @nlist = keys %neighbor;
        $nv = @nlist;
        $vsum = 0;                      # total vertex counts    
        foreach my $v ( @nlist ) {
           $vsum += $vcount[$v];   
           unless ( defined $vcount[$v] ) {
               print "undef\n";
           }      
        }       

        my $c  = int( rand($vsum) );
        my $count = 0;
        foreach my $v ( @nlist ) {
            $count += $vcount[$v];
            $r = $v;
            if ( $c <= $count ) {
                last;
            }
        }
	    #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
	}
	my @vlist = sort {$a<=>$b} keys %selected;
	my $vstring = join "_", @vlist;
	
	# if there are vertices only connected by x edges, you can end up with a
	# one vertex, no edge graph, in this case do over
	if ( $n_selected == 1 ) { goto restart }
	
	# add all edges between selected vertices
	
	my $subgraph = new DfsGraph;
	$subgraph->addStructure();
	$subgraph->addVertices( $n_target );
	
    my @alledges = @{ $graph->edge };
    foreach my $edge ( @alledges ) {
       my ( $v1, $v2, $e ) = @{$edge};
       if ( defined $selected{$v1} && defined $selected{$v2} ) {
           #next unless ( $selected{$v1} < $selected{$v2} );
           
           my $this_edge = $subgraph->addEdge( $selected{$v1}, $selected{$v2}, $e );
           $subgraph->addEdgeToVE( $selected{$v1}, $selected{$v2}, $this_edge );
       }
   }
	    		
	return ( $vstring, $subgraph );
}

# End of sampleByDegree
}

#------------------------------------------------------------------------------
# $Log: fingerprint_random.pl,v $
# Revision 1.1.4.9  2016/08/24 05:28:05  huang147
# added the library.
#
# Revision 1.1.4.8  2015/10/12 19:37:05  huang147
# added the option of working directory and decoy directory.
#
# Revision 1.1.4.6  2015/10/12 17:36:21  huang147
# improved the code in calculating motif mappings.
#
# Revision 1.1.4.5  2013/03/21 19:25:16  huang147
# changed the order in the few headlines.
#
# Revision 1.1.4.4  2013/03/21 19:24:13  huang147
# Got this from mrg121003 branch.
#
# Revision 1.1.2.16  2013/03/14 18:13:38  gribskov
# Added option -q for quiet.
#
# Revision 1.1.2.15  2013/03/13 19:49:50  gribskov
# Converted to use Fingerprint.pm.  Works.
# sampleByvertex seems to have a bug.
#
# Revision 1.1.2.14  2013/01/07 17:53:14  gribskov
# Added and debugged selection by degree.  Fixed error with total number of
# vertices (last is not null)
#
# Revision 1.1.2.13  2013/01/07 16:07:53  gribskov
# Updated documentation and cleaned up and rearranged a little.  Not functional
# changes.
#
# Revision 1.1.2.12  2013/01/07 13:14:44  gribskov
# working version.  Doesn't improve finding unique motifs.  double check if the code is right.
#
# Revision 1.1.2.11  2013/01/05 14:05:10  gribskov
# Adds check to skip DFS generation for vertex sets that have already been seen.
# Adds reports for number of mappings and motifs found.
#
# Revision 1.1.2.10  2013/01/05 12:45:26  gribskov
# No functional changes.  This version used for timings before adding check for vertex sets.
#
# Revision 1.1.2.9  2012/12/27 14:26:19  gribskov
# Fixed timing to report elapsed seconds instead of epoch seconds.
#
# Revision 1.1.2.8  2012/12/26 22:32:30  gribskov
# Added timing and report of vertices and edges in query graph.
#
# Revision 1.1.2.7  2012/12/24 11:23:14  gribskov
# Added gribskov library.
# Changed default motif file.
#
# Revision 1.1.2.6  2012/12/20 15:45:04  gribskov
# Added option -l to stop sampling when all subgraphs have a minimum number
# of hits.
# Added verbose option to control printing of DFS
#
# Revision 1.1.2.5  2012/12/20 14:17:13  gribskov
# Added printout of motif DFS.  maybe this should be optional.
#
# Revision 1.1.2.4  2012/12/20 13:54:54  gribskov
# This version works.  Not exhaustively tested, but result seems OK.
#
# Revision 1.1.2.3  2012/12/19 19:55:55  gribskov
# Includes motif hex indexing and lookup.  Something not quite right since the strings do not match
#
# Revision 1.1.2.2  2012/12/19 18:52:02  gribskov
# Excludes x edges, works in tests.
#
# Revision 1.1.2.1  2012/12/19 17:02:26  gribskov
# Initial version.  Samples graphs but does not omit x edges.
#
#------------------------------------------------------------------------------
