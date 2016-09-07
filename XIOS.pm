#-----------------------------------------------------------------------------
# $Id: XIOS.pm,v 1.5 2008/12/27 13:44:43 gribskov Exp $
#
# This class stores a XIOS graph, namely the vertices and their edges.
# Exclusive (X), inside (J), outside (I), overlap (O) and serial (S)  edges 
# are supported
#
# revision log at end of file
#-----------------------------------------------------------------------------

package XIOS;
use strict;
#use vars qw(@ISA);
#@ISA = qw { AutoLoader };   

#---------------------------------------------------
# new
#
# this is the constructor for this class. 
# the class structure is a hash within hash:
# xios = { $order_index1 => {$index1  => {'Label' => $label, 'X' => [], 'J' => [], 'I' => [], 'O' => []},
#                            $index2  => {'Label' => $label, 'X' => [], 'J' => [], 'I' => [], 'O' => []},
#                            $index3  => {'Label' => $label, 'X' => [], 'J' => [], 'I' => [], 'O' => []},
#                            $index4  => {'Label' => $label, 'X' => [], 'J' => [], 'I' => [], 'O' => []},
#                            ...
#                           }
#          $order_index2 => ...
#        }
#
# USAGE
#      $xios_obj = XIOS->new();
#----------------------------------------------

sub new {
    my ($class, $topology_ptr) = @_;
	my $xios = {};
	bless $xios, $class;
	if (defined($topology_ptr)) {
        $xios->setXIOS($topology_ptr);
    } else {
		print STDERR "WARNING:\nWhen you were using the XIOS class xios->new() function, you failed to provide the topology information\n";
		print STDERR "The DEFAULT set of the xios_obj is 2 stems with \"SS\" expression\n";
        $xios->initialize();
    }
    return $xios;   
}

#---------------------------------------------------
# initialize
#
# This is called by the constructor "new" to initialize the different attributes 
# of the class.It initializes the hash and creates the necessary attributes for the XIOS Object
#
# USAGE
#      $xios_obj->initialize(); 
# TODO: This is used locally within the class, AND THE DEFAULT VALUE WOULD BE 2 STEMS??
#----------------------------------------------

sub initialize {
    my ($self)=@_;
	$$self{'1'}{'1'}={'X' => [], 'J' => [], 'I' => [], 'O' => []};
	$$self{'1'}{'2'}={'X' => [], 'J' => [], 'I' => [], 'O' => []};
    return $self;    
}
#---------------------------------------------------
# setXIOS
#
# updates all the values of the according to a given hash
#
# USAGE
#      $xios_obj->setXIOS($hash_ref)
# TODO: modify the code so that it deals with a list of array
#----------------------------------------------

sub setXIOS {
	my ($self, $topology_ptr)=@_;
	my $topology = ${$topology_ptr}[0];
	my $index = 0;
	foreach my $order (@$topology_ptr) {
		my $nstem = @{$order}/2;
		$index++;
		my $porder = printOrderAsPairs($order);
		my $sign = convertXIOS($porder);
		foreach my $n (1..$nstem) {
			$$self{$index}{$n} = $$sign{$n};
		}
	}
    return;
}


#-----------------------------------------------------------------------------
# printOrderAsPairs
#
# print the topology to stdout.  The paired stems are printed in parentheses,
# i.e., (1,2) (3,4) is a topology with two serial stems
#
# USAGE
#   printOrderAsPairs( $ptr_to_order );
#
# 8 June 2007   Gribskov
#-----------------------------------------------------------------------------
sub printOrderAsPairs{
    my ( $order ) = @_;
    my ( @porder );

    my $nstems = @{$order};

    # reorder into print order

    my $i = 0;
    for my $stem ( @$order ) {
        $i++;
        $porder[$stem-1] = $i;
    }

    #print "porder:@porder:\n";
    #print "order: @$order\n";
    my $new_pair = 1;
    my $pos = 1;
    for my $i ( @porder ) {

        if ( $new_pair ) {
            #print "($i,";
            $new_pair = 0;
        } else {
            #print "$i)";
            $new_pair = 1;
            unless ( $pos == $nstems ) {
                #print " ";
            } else {
                #print "\n";
            }
        }

        $pos++;

    }

    return (\@porder);
}

#-----------------------------------------------------------------------------
# convertXIOS
#
# convert all the topologies to XIOS expression
#
# USAGE
#   convertXIOS( $ptr_to_porder )
#   porder is the order returned by printOrderAsPairs sub
#
# 7/12/2007   Kejie Li
#-----------------------------------------------------------------------------

sub convertXIOS{
	my ( $porder ) = @_;
	
	my %sign;
	my %new_sign;

	my $nstem = @{$porder}/2;
	my ($forward_string, $reverse_string);
	$reverse_string ="";
	foreach my $pos (@$porder) {
		$forward_string .= $pos." ";
		my $comp_pos = 2*$nstem+1-$pos;
		$reverse_string = $comp_pos." ".$reverse_string;
	}
	#print "for: $forward_string, rev: $reverse_string;\n";
	my @split = split " ",$forward_string,2*$nstem;
	foreach my $n1 (1..$nstem) {		
		foreach my $n2 ($n1+1..$nstem) {
			my $x1 = 2*$n1-2;
			my $x2 = 2*$n1-1;
			my $x3 = 2*$n2-2;
			my $x4 = 2*$n2-1;
			if ( $split[$x3] < $split[$x2] ) {
				$sign{$n1}{$n2} = "O";
				$sign{$n2}{$n1} = "O";
				if ( $split[$x4] < $split[$x2] ) {
					$sign{$n1}{$n2} = "I";
					$sign{$n2}{$n1} = "J";
				}
			} else {
				$sign{$n1}{$n2} = "S";
				$sign{$n2}{$n1} = "S";
			}
			#print "$n1 $n2 $sign{$n1}{$n2}\n";
		}
	}

	
	foreach my $n1 (1..$nstem) {
		my @x_array;
		my @i_plus_array;
		my @i_minus_array;
		my @o_array;
		my $label ="";
		foreach my $n2 (1..$nstem) {
			if ($n1 == $n2) {
				next;
			}
			#print "***$sign{$n1}{$n2}***\n";
			unless ( $sign{$n1}{$n2} eq "S" ) {
				$label .= $sign{$n1}{$n2};
			}
			#print "now label is $label\n";
			if ($sign{$n1}{$n2} eq "J") {
				#print "$n1 $n2 $sign{$n1}{$n2}\n";
				push @i_plus_array,$n2;
			}
			if ($sign{$n1}{$n2} eq "I") {
				#print "$n1 $n2 $sign{$n1}{$n2}\n";
				push @i_minus_array,$n2;
			}
			if ($sign{$n1}{$n2} eq "X") {
				#print "$n1 $n2 $sign{$n1}{$n2}\n";
				push @x_array,$n2;
			}
			if ($sign{$n1}{$n2} eq "O") {
				#print "$n1 $n2 $sign{$n1}{$n2}\n";
				push @o_array,$n2;
			}
		}
		#print "Label is $label\n";
		$new_sign{$n1}{'Label'} = $label;
		$new_sign{$n1}{'X'} = \@x_array;
		$new_sign{$n1}{'O'} = \@o_array;
		$new_sign{$n1}{'J'} = \@i_plus_array;
		$new_sign{$n1}{'I'} = \@i_minus_array;
		#foreach my $key (keys %{$new_sign{$n1}}) {
		#	print "key: $key\tvalue:@{$new_sign{$n1}{$key}}\n";
		#}
	}
#	print "\t";
#	foreach my $n1 (1..$nstem) {
#		foreach my $n2 ($n1+1..$nstem) {
#			print "$sign{$index}{$n2}{$n1}";
#		}
#	}
#	print "\n";
	
    return \%new_sign;
}


#---------------------------------------------------
# printXIOS
#
# print the xios object
#
# USAGE
#      $xios_obj->printXIOS($xios_obj_ptr)
#	OR $xios_obj->printXIOS($xios_obj_ptr, $ptr_to_porder)
#
# TODO: need to change the code to support the list topologies with different stem numbers
#----------------------------------------------


#print the xios object
sub printXIOS {
    my ($self, $porder)=@_;

	my @topology_keys = keys %{$self};
	my @keys = keys %{$$self{'1'}}; #THIS ONLY WORKS FOR THE SAME STEM NUMBER
	my $nstem = @keys;
	if ( defined $porder ) {
		foreach my $topo_key ( @topology_keys ) {
			print "For topology @{$$porder[$topo_key-1]}\n";
			foreach my $n2 (1..$nstem) {
				print "$n2\nLabel $$self{$topo_key}{$n2}{'Label'}\n X @{$$self{$topo_key}{$n2}{'X'}}\n O @{$$self{$topo_key}{$n2}{'O'}}\nJ @{$$self{$topo_key}{$n2}{'J'}}\nI @{$$self{$topo_key}{$n2}{'I'}}\n\n";
			}
		}
	} else {
		foreach my $topo_key ( @topology_keys ) {
			print "For topology @{$$porder[$topo_key-1]}\n";
			foreach my $n2 (1..$nstem) {
				print "$n2\nLabel $$self{$topo_key}{$n2}{'Label'}\n X @{$$self{$topo_key}{$n2}{'X'}}\n O @{$$self{$topo_key}{$n2}{'O'}}\nJ @{$$self{$topo_key}{$n2}{'J'}}\nI @{$$self{$topo_key}{$n2}{'I'}}\n\n";
			}
		}
	}
}


#---------------------------------------------------
# printNautyFormat
#
# print the xios object as Nauty format
#
# USAGE
#      $xios_obj->printXIOS($xios_obj_ptr)
#	OR $xios_obj->printXIOS($xios_obj_ptr, $ptr_to_porder)
#----------------------------------------------


#print the xios object as Nauty format
sub printNautyFormat {
    my ($self, $porder)=@_;

	my @topology_keys = keys %{$self};
	my @keys = keys %{$$self{'1'}};
	my $nstem = @keys;
	if ( defined $porder ) {
		foreach my $topo_key ( @topology_keys ) {
			print "For topology @{$$porder[$topo_key-1]}\n";
			foreach my $n2 (1..$nstem) {
				my @neighbor_array;
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'X'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'J'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'I'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'O'}};
				my @sorted_neighbor_array = sort numerically @neighbor_array;				
				print "$n2: @sorted_neighbor_array";
				if ($n2 == $nstem) {
					print ".\n";
				} else {
					print ";\n";
				}
			}
		}
	} else {
		foreach my $topo_key ( @topology_keys ) {
			print "For topology @{$$porder[$topo_key-1]}\n";
			foreach my $n2 (1..$nstem) {
				my @neighbor_array;
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'X'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'J'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'I'}};
				push @neighbor_array,@{$$self{$topo_key}{$n2}{'O'}};
				my @sorted_neighbor_array = sort numerically @neighbor_array;
				print "$n2: @sorted_neighbor_array";
				if ($n2 == $nstem) {
					print ".\n";
				} else {
					print ";\n";
				}
			}
		}
	}
}



# AUTOLOAD
# When methods refer unknown functions to their parent, they all end up 
# here.  A last attempt is made to treat unknown functions as accessor 
# functions
# 
# $subclass->unknown_function( $value )
#   sets $subclass->{unknown_function} = $value and returns $value
#
#  $subclass->unknown_function
#   returns $subclass->{unknown_function}
#
# This function will not create new hash keys.  This keeps typos from
# crazy new attributes for objects.  Anything that is not a hash key or 
# a function in this package fails.

sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;
    
    my $unknown = $AUTOLOAD;                # The name of the unknown function
    $unknown =~ s/.*:://;                   # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;    # Skip all-caps methods like DESTROY
    
    if ( exists $self->{ $unknown } ) {         # do not create new fields
        $self->{ $unknown } = shift if @_;      # Set new value if one was supplied
        return( $self->{ $unknown } );          # Return current value
    }
    
    # the buck stops here, there are no more places to check
    print STDERR "UNKNOWN Function called\n";
   # $self->dieStackTrace( "Object::AUTOLOAD: Undefined function $unknown: $!\n" );
}

sub numerically    { $a <=> $b };

#-----------------------------------------------------------------------------
# $Log: XIOS.pm,v $
# Revision 1.5  2008/12/27 13:44:43  gribskov
# Updated documentation.  Added cvs tags.
#
#-----------------------------------------------------------------------------
1;
