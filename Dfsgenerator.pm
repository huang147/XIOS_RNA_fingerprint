package Dfsgenerator; 
#-----------------------------------------------------------------------------
# $Id: Dfsgenerator.pm,v 1.17 2011/09/01 13:35:44 rahmanr Exp $
# $Author: rahmanr $     $Date: 2011/09/01 13:35:44 $
#
# this package generates the dfscode for a given dfsgraph object and the 
# structure number. This code is optimized to generate the dfscode as quickly as possible. 
#
# USAGE
#
#     my ($dfs_gen_obj) = new Dfsgenerator ({'dfsgraph'=>$dfsgraph, 'struct_num' => 0});
#     my @dfscode =@{$dfs_gen_obj->getDfscode()}; 
#
# @dfscode has the minimum dfscode for the graph represented by structure 0;
# use this with Maplite package
#
# Revision log at end of file.
#-----------------------------------------------------------------------------
use strict;
use warnings;
use DfsGraph;
use Maplite;
my $DEBUG=0;
my $TODO = 0;

#-----------------------------------------------------------------------------
# new:
#     construcutor for dfsGenerator
#
#----------------------------------------------------------------------------


sub new {
    my ( $class, $param) = @_;
    my $self = {};
    bless $self, $class;
    $self->{dfscode} = [];
    if (defined $param->{'dfsgraph'} && defined $param->{'struct_num'}) {
	$self->{dfscode} = &dfsGenerator($param);
    }
    return $self;
}

sub setDfscode{
    my ($self, $dfscode )= @_;

    if ( defined $dfscode ) {
	$self->{dfscode} = $dfscode;
    }
    return;
}

sub getDfscode{
    my $self= shift;
    return $self->{dfscode};
}

#----------------------------------------------------------------------------
# dfsGenerator
#
# Generates the minimum dfscode for a given graph. This subroutine is called by the constructor, 
# When the graph is loaded into a DfsGraph object and the strucutre number is known
#
#----------------------------------------------------------------------------

sub dfsGenerator  {
    my ($param) = @_;
    my $dfsgraph = $param->{'dfsgraph'};
    my $struct_num = $param->{'struct_num'};
    my ($map_list, $vlist) = &initialMap($dfsgraph, $struct_num);
#    foreach my $mp (@$map_list) {
#	print $mp;
#	$mp->dump("TEST");
#    }  
#    die;
    my $dfscode_main = [];
    my $rmp = [0];   
    while ( @$map_list) {
#	my @elist;
	my $graph = $struct_num;
	#reconstruct rmp path here
	if (@$dfscode_main) { #true except the first round
	    $rmp = findRmp($dfscode_main);
	}
	my %ehash;
	my ($ehash_ref, $mindfs, $count_hash) = findExtensions($graph, $map_list, $vlist, $rmp, \%ehash);
	#$ehash_ref is just the reference to the 

#also need backward todo and count of numbers of maps for each of them;


	my ( $m1, $m2, $mt );

	# if there is an extension with support, $mindfs is true.  Update @all_map and go to next cycle
	
	if ( $mindfs ) {
	    ( $m1, $m2, $mt ) = split ",", $mindfs;
	    my $row = [$m1, $m2, $mt];
	    #	    print "     $mindfs |||||||||||\n";
	    my $curr_map_size =  @{$ehash{$mindfs}->{dfs}};
	   $TODO && print "$mindfs     size=$curr_map_size\n";
	    my $prev_map_num = @$map_list;
	    
	    foreach my $map_obj (@{$map_list}) {
		undef %$map_obj;
	    }
	    
	    @$map_list = ();
	    my $ee = 0;
	    while ( $ee < @{$ehash{$mindfs}->{dfs}} ) {
		#print "adding to allmap $ehash{$mindfs}->{dfs}->[$ee], @{$ehash{$mindfs}->{dfs}->[$ee+1]} \n";
#		my $graph_num = $ehash{$mindfs}->{dfs}->[$ee+2];
		my $new_map_obj = new Maplite({'dfsrow'=>$row, 'dtg'=>$ehash{$mindfs}->{dfs}->[$ee] });
		push (@{$map_list}, $new_map_obj);
		$ee++;
		
#		$new_map_obj->dump("ADDED");
	    }
	    push @$dfscode_main, [$m1, $m2, $mt];
	    %ehash = ();
#what about the backedges that have been found, lets add them in one go
	    my @backward_todo=();
	    foreach my $key (keys %$count_hash) {
		next if ($mindfs eq $key);
	$TODO && 	print "     backward todo: $key, size= $count_hash->{$key}\n";
		my ($d1, $d2, $dt) = split ",", $key;
		push @backward_todo, [$d1, $d2, $dt];
	    }
	    my $record_map_num;

	    foreach my $e ( 
			    sort {  if ( defined $$a[1] && defined $$b[1] ) {
				return $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]; #sort them by order of valuse in dfscode
			    } else {
				return $$a[2] <=> $$b[2];
			    } } @backward_todo			
			    ) {
		#    print "Sorted backward edges: @$e|$e\n";
		my $dfs_str = "$e->[0],$e->[1],$e->[2]";
		my $dfs_row_ref = [$e->[0],$e->[1],$e->[2]];
		my $next_map_num = $count_hash->{$dfs_str}; 
	$TODO && 	print "  BACKWARD edge under consideration: $dfs_str | size = $next_map_num\n";
		if (($curr_map_size == $next_map_num) && ($prev_map_num == $curr_map_size)) {
		    foreach my $map_obj (@{$map_list}) {
			#set dfs row, since this is backward edge, no need to update dtg maps
			$map_obj->setDfsrow($dfs_row_ref);
		    }
		    push @$dfscode_main, [$e->[0],$e->[1],$e->[2]];
		   $TODO &&  print "backward edges added: @$e ||| maps: [$next_map_num == $curr_map_size && $prev_map_num == $curr_map_size]\n";
		} elsif (($curr_map_size >= $next_map_num) && ($prev_map_num >= $curr_map_size)) { 
		    # All other cases, when the first condition is not met
		    my ($min_maparr, $non_minmap_arr) = &find_known_backward_extension($graph, $map_list, $vlist, $rmp, $e);
		    if (@$min_maparr) {
			$map_list = $min_maparr;
		$TODO && 	print "  backward edges added: @$e ||| maps: [$curr_map_size >= $next_map_num) && ($prev_map_num >= $curr_map_size]\n";
			$record_map_num=  $curr_map_size;
			$curr_map_size = @$map_list;
                        #	print "backward edges added: @$e ||| maps: [$new_map_num == $maplist_size]\n";
			#the mapobjects have already been updated, now update the master dfscode
			push @$dfscode_main, [$e->[0],$e->[1],$e->[2]];
			#delete map objects in @$non_minmap_arr;
			if (@$non_minmap_arr) {
			    foreach my $map_obj (@{$non_minmap_arr}) {
				undef %$map_obj;
			    }
			}

			if ($record_map_num==  $curr_map_size) {
			    last; #not sure if this last, makes a difference in speed
			}
		    } else {
			$map_list = $non_minmap_arr;
			last; #go back to findextension and try again
		    }
		} elsif (($curr_map_size < $next_map_num) && ($prev_map_num > $curr_map_size)) {
		    last; #prev = 10, curr = 2 but next_map = 5
		    
		} else {
		    print "should not reach here, ERROR!\n";
		    print "$curr_map_size ? $next_map_num) && ($prev_map_num ? $curr_map_size\n";
		    die;
		    #last;
		}
		
	    }
	    %$count_hash = ();
	} else {
	    #need more work here...
#	    my @dfscode =  @$dfscode_main; #@{$first_map_obj->getDfs()}; #need the master dfscode her
	    my $ref = \@$dfscode_main;
	    foreach my $map_obj (@{$map_list}) {
		undef %$map_obj;
	    }
	    @$map_list = ();
	    
#	    foreach my $rr (@dfscode) {
#		print "@$rr\n";
#	    }
	    %ehash = ();
	    return $ref;
	}
    }   
}


#-----------------------------------------------------------------------------
# findRmp
#
# Returns a reference to a list of vertices in the rightmost path of $dfs.
# Start from the last row of the DFS and work backwards.  The last last
# forward edge is the rightmost vertex and its parent.  Trace the chain of
# parents back to the first row to get the rightmost path.
#
# USAGE
#   $rightmost_path = findRmp( $dfs );
#-----------------------------------------------------------------------------
sub findRmp{
   my ($dfs) = @_;
   my @rp = ();                            # rightmost path

   my $i = $#{$dfs};
   my $last;

   # skip terminal backward edges
   my $p = $dfs->[$i];                     # pointer to last row of DFS
   while ( ($last=$$p[0]) > $$p[1] ) {     # v1 > v2; skip backward edges
       $p = $$dfs[--$i];
   }
   push ( @rp, ($last,$$p[1]) );           # push v1 and v2 on rightmost path
   # $$p[1] is rightmost vertex
   # $last ($$p[0]) is its parent
   while ( $i >= 0) {
       $p = $$dfs[--$i];                   # decrement row in DFS
       if ($last == $$p[1] && $$p[1] > $$p[0] ) {   # most recent vertex is child (v2), skip backwards edges
           unshift (@rp, ($last=$$p[0]) );          # add parent to rightmost path, and save as most recent
       }
   }
   return \@rp;
} 


#----------------------------------------------------------------------------
# findExtensions
#
# given the current DFS and dtg maps, find possible extensions for a single
# graph.
#
# USAGE
#    my @elist = findExtensions( $map[$graph] );
#-----------------------------------------------------------------------------
sub findExtensions{
    my ($graph, $p_map, $p_vlist, $rmp_arr, $ehash ) = @_;
    my @vlist = @{$p_vlist};
    my @elist;
    my $map = $#{$p_map};
    my %rplist;
#    print "\nENTERING FINDEXTENSION\n";

    my @rmp = @$rmp_arr;
    my $minimum_dfs;
    my ( $m1, $m2, $mt );
    my $count_hash={};

    my $skipped = 0;
    while ( $map >= 0 ) {
        $DEBUG && print "\n---------------------------\n* NEXT MAP STEP - GRAPH: ($map) *\n---------------------------\n\n";
        #print "   map=$map\n";
		my $map_obj = $p_map->[$map];
	#	print "MAP OBJ: $map_obj   ||| $map\n";
	#	$map_obj->dump("TEST");
		$map--;
		my @dtg = @{$map_obj->getDtg()}; #@{ $$p_map[$map] };
		my @code = @{$map_obj->getDfsrow()}; #@{ $$p_map[$map] };
	#	my @rmp = @$rmp_arr; #@{$map_obj->getRmp()};  
		#may be i don't have to instantiate it as well
		my @rp = ();
	#	print "dtg: |@dtg|, code: |@code|, rmp: |@rmp|\n"; 
		
		# construct gtd
		my @gtd;
		foreach my $i ( 0 .. $#dtg ) {
			if ( defined $dtg[$i] ) {
				$gtd[ $dtg[$i] ] = $i;
				#print "dtg[$i] = $dtg[$i]     gtd[$dtg[$i]] = $gtd[$dtg[$i]]\n";
			}
		}
		$DEBUG && print "   gtd: @gtd\n";
		$DEBUG && print "   dtg: @dtg\n";
		$DEBUG && print "    rmp: @rmp\n";
		# intiialize rightmost path (@rp)
	#	print "dtg: |@dtg|, code: |@code|, rmp: |@rmp|, gtd: |@gtd|\n";
		
		#get the dfs code of the last row
		my ( $lc1, $lc2, $lct );  
		if (@code) {
			($lc1, $lc2, $lct )=  @code;#@{ $code[$#code] };
			$DEBUG &&     print "Last code: $lc1, $lc2, $lct\n";
		}
		my $right_v = $#dtg;
		
		#Now loop over the rightmost path
		# special case: righmost vertex goes through a special case
		# special case 2: and when rightmost vertex is 0 and when no dfs rows have been made
		my $index = $#rmp;
		for (my $i=$index; $i>=0; $i--) {
		#	foreach my $node (@rmp) { #$node represents the dfs value
			my $node = $rmp[$i];
			my @edge_list = @{ $vlist[$dtg[$node]] };
			my $nbor = [];
			# following sort extended nodes:
			# 1. if both nodes are used, put the smaller dfs index one in the front
			# 2. if only one of them is used, put the used node to the front
			# 3. if both nodes were NOT used, put the smaller edge type to the front
			foreach my $e (
				   sort {  if ( defined $gtd[$$a[1]] && defined $gtd[$$b[1]] ) {
					   return $gtd[$$a[1]] <=> $gtd[$$b[1]]; #sort them by order of valuse in dfscode
				   } elsif ( defined $gtd[$$a[1]] ) {
					   return -1;
				   } elsif ( defined $gtd[$$b[1]] ) {
					   return 1;
				   } else {
					   return $$a[2] <=> $$b[2];
				   } } @edge_list 
				   ) {
				if ($node == $right_v) { #rightmost node can add both forward and back edges
					if (($node ==0) && ($#code == 0)) { #0 is the right most vertex, and no dfs rows have been added
						next if ( defined $gtd[$$e[1]] );       #skip backward edges or used edges
						#push @{$nbor},  [$$e[1],$$e[2]];        # neighbor vertex and edge type in graph numbering
						$DEBUG && print "   saving neighbor (rightmost begin at: 0)  $node: $dtg[$node] -> $$e[1]   type:$$e[2]\n"; 	
					} else {
						my $incoming_vertex = $rmp[$#rmp-1]; #the second last element form the @rmp
						next if ( defined $gtd[$$e[1]] && $incoming_vertex == $gtd[$$e[1]] ); # skip the incoming edge 
						if (@code && $lc1> $lc2) { #the last added dfs is backward
							$DEBUG && print "gtd[e1]: $gtd[$$e[1]] <= lc2: $lc2\n";
							next if ( (defined $gtd[$$e[1]]) && ($gtd[$$e[1]] <= $lc2)); 
							# skip all backward edges smaller than equal to the last backward edge
							#don't skip if the vertex has not been defined
						} 
						$DEBUG && print "   saving neighbor $node (rightmost): $dtg[$node] -> $$e[1]   type:$$e[2]\n";
						push @{$nbor},  [$$e[1],$$e[2]]; 
					}
				} else {
					# other nodes, including node 0, can add only forward edges
					next if ( defined $gtd[$$e[1]] );       #skip backward edges or used edges
					push @{$nbor},  [$$e[1],$$e[2]];        # neighbor vertex and edge type in graph numbering
					$DEBUG && print "   saving neighbor $node: $dtg[$node] -> $$e[1]   type:$$e[2]\n"; 
				}
				#this could be a heuristic to use when, number of edges are high...
				####if some backward edges have already been added to the @$nbor, then I can stop looking through forward edges...
				if (@{$nbor}) { #interesting to see how much effect it will have.
					my ($f2, $ft) = @{$nbor->[0]};
					#if the first element of the neighbor list is backward edges, and the current is forward, then last
					if ((defined $gtd[$f2]) && !(defined $gtd[$$e[1]])) { 
						last; #forward edges are being considered now, no need to consider them
					}
				}
	# I can last if the current entry is a forward egde, where gtd map is undef
			}
			
	#	if there are neighbors, save on rightmost path
			if ( @{$nbor} ) {
				#push @rp, [0, $nbor];
				push @rp, [$node, $nbor];
		#		print "NODE: $node, Neighbors: ";
		#		foreach my $nb (@$nbor) {
		#		    my ($g2, $gt) = @$nb;
		#		    if (defined $gtd[$g2]) {
		#			print "B $gtd[$g2], $gt|";
		#		    } else {
		#			print "F $g2, $gt|";
		#		    }
		#		} 
		#		print "\n";
				last; #no need to look for more neighbors
			} else {
				$DEBUG && print "   no neighbors for NODE: $node\n";
			}
		}
		
	 #---------------------------------------------------------------------
			# end of loop over DFS finding all neighbors
			#---------------------------------------------------------------------
		
			#---------------------------------------------------------------------
			# check along rightmost path for extensions
			#---------------------------------------------------------------------
		
		my $rpstring = "";
		foreach my $v ( @rmp ) {
			$rpstring .= "$dtg[$v],";
		}
		#print "   rightmost path: $rpstring\n";
		if ( defined $rplist{$rpstring} ) {
			$rplist{$rpstring}++;
			#print "Path skipped: $rpstring\n";
			$skipped++;
	#	    print "     skipping rmp, gvertex: $rpstring||| dvertex: ";
	#	    foreach my $v (@rmp) {
	#		print "$v,";
	#	    }
	#	    print "||| $lc1, $lc2, $lct ||| $skipped\n";
			next;
		} else {
			$rplist{$rpstring} = 1;
		}
		
		# find all possible extensions. the minimum extension for this 
		# mapping is the first one.  this depends in the
		# nbor list being already sorted in lexicographic order
		

		my $new_right_v = $right_v + 1;         # update to next available dfs index 
		my $first = 1;
		my $minedge;
		while ( @rp ) {
			my ($c, $nbor ) = @{ pop @rp };
			next unless ( $nbor );
			#print "  popping from rp:  c=$c   nbor=@{$nbor}   rp=@rp\n";
			my $c2;
			foreach my $e ( @{$nbor} ) {
				#print "   path: $c -> $$e[0]  t:$$e[1]     ";
		#		my @temp_rmp= ();
				my $new_d = undef;
				if ( defined $gtd[$$e[0]] ) {
					$c2 = $gtd[$$e[0]];
					#rmp does not need to be updated, if it is a backward edge
				} else {
		#           @temp_rmp = @rmp;
					$c2 = $new_right_v;
					$new_d = $$e[0];
				}
					
				#print "   possible extension dfs: $c, $c2, $$e[1]    dtg: @dtg   new vertex: $new_d   graph:$graph\n";
				if ( $first ) {
					$first = 0;
					$minedge = $$e[1];
				}
				my $ismin = 0;
				if ( $$e[1] == $minedge ) {
					$ismin = 1;
				}
				#########################
				my $dfs_str = "$c,$c2,$$e[1]";
				#i should do some kinda minimum test here as well
				# edges are not minimum, really does not need to be in the ehash...
				if ($minimum_dfs) {
					my ( $t1, $t2, $tt ) = ($c,$c2,$$e[1]);    
					if ($t1> $t2) {
						$count_hash->{$dfs_str}++;
					}
					if ( ( $t2 < $m2 ) or
					 ( $t2 == $m2 && $t1 > $m1 ) or
					 ( $t2 == $m2 && $t1 == $m1 && $tt< $mt ) ) {
						# found a new minimum
						$m1 = $t1;  $m2 = $t2; $mt = $tt; 
						$minimum_dfs = $dfs_str;  #"$c,$c2,$$e[1]";
					} 
				} else {
					#first round
					$minimum_dfs = $dfs_str; #"$c,$c2,$$e[1]";
					( $m1, $m2, $mt ) = ($c,$c2,$$e[1]);
					if ($m1> $m2) {
						$count_hash->{$dfs_str}++;
					}
					
				}
					

				#put stuff in the ehash directly && only when it is minimum		
				if ( $m1==$c && $m2==$c2 && $mt==$$e[1]) {
					if (@{$nbor} > 1) {
						#multiple edges using the same dtg map
						my $dtg_arr =  [@dtg];
						if (defined $new_d) { 
							push @$dtg_arr, $new_d; 
						}
						push @{$ehash->{$dfs_str}->{dfs}}, $dtg_arr;
						$ehash->{$dfs_str}->{ngraph}->[$graph] = 1;
						$ehash->{$dfs_str}->{ismin} = $ehash->{$dfs_str}->{ismin} || $ismin;	
					} else { #this will be risky for the mapping project
						my $dtg_arr =  \@dtg;
						if (defined $new_d) { 
							push @$dtg_arr, $new_d;
						}
						push @{$ehash->{$dfs_str}->{dfs}}, $dtg_arr;
						$ehash->{$dfs_str}->{ngraph}->[$graph] = 1;
						$ehash->{$dfs_str}->{ismin} = $ehash->{$dfs_str}->{ismin} || $ismin;
					}
				}
				#put stuff in the ehash directly					
				#	}
			}       # end of each nbor
			# make sure only the minimum extension matches to minedge
			if ( !$first && @{$nbor} ) {
				$minedge = 1000;
			}
		}           # end of finding all the possible extension along the rightmost path
    }               # end of each mapping
	#   print "    Number of maps skipped: $skipped ($minimum_dfs)\n";
    return ($ehash, $minimum_dfs, $count_hash);
}

#----------------------------------------------------------------------------
# initialMap
#
# Given DfsGraph object, structure num
# Returns a list of dtg maps for the first vertex, and the edge list associated with each vertex 
#
# USAGE
#      my ($map_list, $vlist) = &initialMap($dfsgraph, $struct_num);
#
#----------------------------------------------------------------------------


sub initialMap {
    my ($dfs, $struct_num) = @_;
    #print "\nInitializing vertex lists and graph maps\n\n";
    my $graph= $struct_num;
    # make a list of edges by vertex, and initialize list of starting vertices in @map
    my @used;
    my @vlist;
    my @map;
    my $edge = $dfs->edgeList( $graph );
    my $minedge;
    if ($#$edge > -1) { 
    	# there were errors "Can't use an undefined value as an ARRAY reference at Dfsgenerator.pm line xx"
    	# added this if statement to avoid the empty edge case
	    foreach my $e ( sort { $a->[2]<=>$b->[2] } @{$edge} ) {
	        my ( $c1, $c2, $type ) = @{$e};
		#    print "      (c1,c2,type) = ($c1,$c2,$type)\n";
	        unless ( defined $minedge ) { $minedge = $type; }
	        unless ( defined $used[$c1] ) {
	            if ( $type == $minedge ) {
			my $map_obj = new Maplite({'dfsrow'=>[], 'dtg'=>[$c1]});
			unshift (@map, $map_obj);
	#		unshift @map, ( [], [$c1] );
	#                print "         pushing $c1 on map\n";
	            } 
	            $used[$c1] = 1;
	        }
		
	        if ( $type == 0 || $type == 1) { 
	            if ( $c1 < $c2 ) {
	                push @{$vlist[$c1]}, [$c1,$c2,$type];
	                push @{$vlist[$c2]}, [$c2,$c1,1-$type];
			
	            } elsif ($c1 > $c2) {
			        push @{$vlist[$c1]}, [$c1,$c2,$type];
	                push @{$vlist[$c2]}, [$c2,$c1,1-$type];	
		        }
            } else {
	            push @{$vlist[$c1]}, [$c1,$c2,$type];
	            push @{$vlist[$c2]}, [$c2,$c1,$type];
	        }
	    }
    }
    
    return (\@map, \@vlist);
}

#----------------------------------------------------------------------------
# find_known_backward_extension
#
# Given structure num, list of dtg maps, vertex list, rightmost path, edge
# Returns array ref to dtg maps for minimum extension and non-minimum extension for a known backward edge.
#
# USAGE
#     my ($min_maparr, $non_minmap_arr) = &find_known_backward_extension($graph, $map_list, $vlist, $rmp, $e);
#
#----------------------------------------------------------------------------

sub find_known_backward_extension {
    my ($graph, $map_arr, $p_vlist, $rmp_arr, $e) = @_;
    # $e is the edge in the backward_todo_edgelist
    my @vlist = @{$p_vlist}; 
    my $SHOULD_ADD=0;
    my $min_maparr =[];
    my $non_minmap_arr=[];
    my @rmp = @$rmp_arr;
    while (@{$map_arr}) {
	#		    foreach my $map_obj (@{$map_arr}) {
	my $mp_size = @{$map_arr};
	my $map_obj = shift @{$map_arr};
	#		print "   Map obj shifted: $map_obj |||    map index # $graph:$mp_size; \n";
	#add the edge $dfs-str to the map object... need to check to make sure the edge exist!
# the backward edge exists for sure
#pass in map and vertex list...
	my @dtg = @{$map_obj->getDtg()}; 
	my @code = @{$map_obj->getDfsrow()};
	#print "dtg: |@dtg|, code: |@code|, rmp: |@rmp|\n"; 
	# construct gtd
	my @gtd;
	foreach my $i ( 0 .. $#dtg ) {
	    if ( defined $dtg[$i] ) {
		$gtd[ $dtg[$i] ] = $i;
	    }
	}
#print "dtg: |@dtg|, code: |@code|, rmp: |@rmp|, gtd: |@gtd|\n"; 
	my ( $lc1, $lc2, $lct );  
	if (@code) {
	    ($lc1, $lc2, $lct )=  @code;
	    $DEBUG &&     print "Last code: $lc1, $lc2, $lct\n";
	    #		    print "      Last code: $lc1, $lc2, $lct\n";
	}
	my $right_v = $#dtg;
	my $rm_node = $rmp[$#rmp];
	my @edge_list = @{ $vlist[$dtg[$rm_node]] };
	foreach my $eg (
			sort {  if ( defined $gtd[$$a[1]] && defined $gtd[$$b[1]] ) {
			    return $gtd[$$a[1]] <=> $gtd[$$b[1]]; #sort them by order of valuse in dfscode
			} elsif ( defined $gtd[$$a[1]] ) {
			    return -1;
			} elsif ( defined $gtd[$$b[1]] ) {
			    return 1;
			} else {
			    return $$a[2] <=> $$b[2];
			} } @edge_list 
			) {
	    my $incoming_vertex = $rmp[$#rmp-1]; #the second last element form the @rmp
	    next if ( defined($gtd[$$eg[1]]) && ($incoming_vertex == $gtd[$$eg[1]]) ); # skip the incoming edge 
	    
	    if (@code && $lc1> $lc2) { #the last added dfs is backward
		# $DEBUG && print "gtd[e1]: $gtd[$$e[1]] <= lc2: $lc2\n";
		next if ( (defined $gtd[$$eg[1]]) && ($gtd[$$eg[1]] <= $lc2)); 
		# skip all backward edges smaller than equal to the last backward edge
	    } 
	    if ($e->[0] != $rm_node) { 
		print "$e->[0] != $rm_node this case should not happen\n"; 
		die; #error check, should never happen...
	    }
#	    print "    edge: @$eg|||  $e->[1] == $gtd[$$eg[1]]) && ($e->[2] == $$eg[2]) || dtg: |@dtg|; gtd: |@gtd|\n";
	    if (!defined($gtd[$$eg[1]])) { #this is a forward edge
		last;
	    } elsif (($e->[1] == $gtd[$$eg[1]]) && ($e->[2] == $$eg[2]) ) {
		$SHOULD_ADD = 1;
		last;
		#the desired backward found, and ready to to be added
	    } else { #other backward edges
		if ($gtd[$$eg[1]] > $e->[1]) { #back edge is bigger than the one we are looking for
		    last;
		}
	    }
	}
	if ($SHOULD_ADD) {
#			  my $dfs_str = "$e->[0],$e->[1],$e->[2]";
#	    print "SHOULD_ADD: $SHOULD_ADD| @$e\n";
	    $map_obj->setDfsrow($e);

	#    push @{$map_obj->getDfs()}, $e;
	    #		    print "backward edges added: @$e ||| map number varied $map_obj\n\n";
	    push @$min_maparr, $map_obj; 
	    $SHOULD_ADD=0;
	} else {
	    push @$non_minmap_arr, $map_obj; #all the maps that did not the backward edge
	}		
    } #end while loop
    

    return ($min_maparr, $non_minmap_arr);
}



##############################################################################
# AUTOLOAD
#
# This function will not create new hash keys.  This keeps typos from creating
# crazy new attributes for objects.  Anything that is not a hash key or 
# a function in this package fails.
#-----------------------------------------------------------------------------

sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;
    
    my $unknown = $AUTOLOAD;                    # The name of the unknown function
    $unknown =~ s/.*:://;                       # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;        # Skip all-caps methods like DESTROY
  
    if ( exists $self->{ $unknown } ) {         # do not create new fields
        $self->{ $unknown } = shift if @_;      # Set new value if one was supplied
        return( $self->{ $unknown } );          # Return current value
    }
    return;        
}


#-----------------------------------------------------------------------------
# $Log: Dfsgenerator.pm,v $
# Revision 1.17  2011/09/01 13:35:44  rahmanr
#
# updated the documentation, hopefully it is easier to read now
#
# Revision 1.16  2011/07/29 16:22:55  gribskov
# Reformatted slightly.  No functional changes.
#
# Revision 1.15  2011/07/29 16:19:01  gribskov
# Latest version modified by kejie.  it was not checked in so i don't
# know if it works
#
# Revision 1.14  2011/01/05 02:35:39  likejie
# bug found from line 546-553 $type could be 0, 1, 2 and 3
#
# type 0 and 1 should be put together (not 1,2,3 in a group)
#
# Revision 1.13  2010/04/22 00:54:52  rahmanr
#
# Updated with a faster version findRmp
#
# Revision 1.12  2010/03/10 16:05:06  rahmanr
#
# removed print statements
#
# Revision 1.11  2010/03/10 11:58:21  rahmanr
#
# Minor update
#
# Revision 1.10  2010/03/09 11:03:36  rahmanr
#
# Fixed some bugs regarding addition of backward edges, should be faster than before
#
# Revision 1.9  2010/03/04 19:02:07  rahmanr
#
# removed the print statements
#
# Revision 1.8  2010/03/04 18:49:30  rahmanr
#
# Fixed the bug that was creating the error... will try to make it run faster for the next update
#
# Revision 1.7  2010/03/04 10:25:33  rahmanr
#
# Made changes to reduce computation in calculating backward edges All backward needs to searched only once... even if the number maps are reduced after each addition of backward edges
#
# Revision 1.6  2010/03/03 09:14:56  rahmanr
#
# made major changes to Dfsgenerator.pm to make it more memory efficient without losing much speed. Removed the use of @elist completely and use ehash as efficiently as possible!
#
# Revision 1.5  2010/03/02 06:03:23   rahmanr
#
# the package now uses maplite.pm and thus uses significantly less memory
# 
# Revision 1.4  2010/03/02 04:44:17   rahmanr
# 
# Revision 1.3  2010/02/25 19:47:14   likejie
#
# added one if statement to avoid empty edge list case
#
# Revision 1.2  2010/02/24 15:36:46   rahmanr
# 
# Commenting out some print statement
#
# Revision 1.1  2010/02/24 09:55:28   rahmanr
#
# This package can be used to generate the minimum dfscode for a graph give the dfsgraph object and the structure number.
# 
#-----------------------------------------------------------------------------

# end of Dfsgenerator package

1;
