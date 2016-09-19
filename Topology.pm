package Topology;

#-----------------------------------------------------------------------------
# $Id: Topology.pm,v 1.56.4.29 2016/08/24 05:28:32 huang147 Exp $
# A topology is a set of stems that may form in a single sequence.  It may
# contain overlapping and exclusive structures.
#
# Stems are contained in the stem_list element of the class hash.  They can be
# accessed with
#   my @stemlist = @{ $topology->stemList };
# usage
#
#   create a  new topology
#
#       $rna = new Topology;
#       my $nstem = $rna->XIOSRead( $xios_file );           # read XIOS
#       my $nstem = $rna->mfoldMultiple( $ct, 1, $ddg );    # read unafold subopt
#
#   read a structure in stem coordinate format, e.g., 1 7  66 72
#
#    my $rna = Topology->new;
#    while ( my $line = <> ) {
#        next unless ( $line =~ /\w/ );
#        my @value = split " ", $line;
#        my $stem = Stem->new( \@value );
#        $rna->addStem( $stem );
#    }
#
#   read a structure in pair array format, e.g. 1,5,2,6,3,4
#
#    my $rna = Topology->new;
#    while ( my $line = <> ) {
#        next unless ( $line =~ /\w/ );
#
#        $line =~ s/[\(\),]/ /g;    # replace parentheses and commas with spaces
#        my @value = split " ", $line;
#        $rna->pairArray( \@value );
#    }
#
#   process the list of stems
#
#       foreach my $stem ( @{$rna->stemlist} ) {...}
#
# Revision log at end of file
#-----------------------------------------------------------------------------
use strict;
use Data::Dumper;
use XML::Simple;

#use lib '/cbio/proj/rna/gribskov/RNA/src/perl_src/';
#use lib "/home/huang147/RNA/src/perl_src/";
#use lib '/home/huang147/code/';
use lib './';
use Stem;
use Exporter;
use Statistics::Descriptive;

use vars qw( @ISA @EXPORT @EXPORT_OK $VERSION );
@ISA       = qw( Exporter );
@EXPORT    = qw( );
@EXPORT_OK = qw(basepairToStem splitStems);

my $REVISION = '$Revision: 1.56.4.29 $';
($VERSION) = $REVISION =~ /Revision: ([\d.]*) /;

# define the allowable types of parentheses and whether they are
# opening (level_type=1) or closing (level_type=-1)

my $nbracket_types = 4;

my @bracket = ( '.', '(', '[', '{', '<', '-' );

my %index = (
    '(' => 1,
    ')' => 1,
    '[' => 2,
    ']' => 2,
    '{' => 3,
    '}' => 3,
    '<' => 4,
    '>' => 4,
    '.' => 0,
    '-' => 0
);

my %value = (
    '(' => 1,
    ')' => -1,
    '[' => 1,
    ']' => -1,
    '{' => 1,
    '}' => -1,
    '<' => 1,
    '>' => -1,
    '.' => 0,
    '-' => 0
);

#-----------------------------------------------------------------------------
# new
#
# this is the constructor for this class. If some attributes of the class
# are provided in $attributes hash it calls setTopology($attribute_hash),
# otherwise it calls initialize() to set up the base attributes
# of the class.
#
# USAGE
#      $topology_obj = Topology->new();
#
# NEED TO Change code to add different attributes
#-----------------------------------------------------------------------------
sub new {
    my ( $class, $attribute_hash ) = @_;

    my $topology = {};
    bless $topology, $class;

    if ( defined($attribute_hash) ) {
        $topology->setTopology($attribute_hash);
    }
    else {
        $topology->initialize();
    }

    return $topology;
}

#-----------------------------------------------------------------------------
# initialize
#
# This is called by the constructor "new" to set up empty attributes of the
# class.
#
# USAGE
#      $Topology_obj->initialize();
#
#-----------------------------------------------------------------------------

sub initialize {
    my ($self) = @_;
    $self->{stem_list}    = [];
    $self->{vienna}       = '';
    $self->{level_array}  = [];
    $self->{brackets}     = ();
    $self->{sequence}     = '';
    $self->{sequence_id}  = '';
    $self->{sequence_doc} = '';
    $self->{comment}      = [];    # a generic place to store other information
    $self->{length}       = 0;
    return;
}

#-----------------------------------------------------------------------------
# setTopology
#
# If certain attribute of the topology object are already specified, then when
# the topology object is created those attributes are already set.
#
# options:
#   unpaired => <int> : set maximum number of unpaired bases in a stem to <int>
#   data     => <string>  : get structure from Vienna parenthesis string
#   vienna   => <string>  : same as data
#   mfold    => <file>    : get structure in the file in Zuker mfold/ct format
#   ct       => <file>    : same as mfold
#   pair    => [ [base1,base1'], [base2,base2'],...]
#
# USAGE
#      $Topology_obj->setTopology($attribute_hash);
#
#-----------------------------------------------------------------------------
sub setTopology {
    my ( $self, $attribute_hash ) = @_;

    $self->initialize;

    my $max_unpaired = 2;
    my $findstems    = 0;
    my $file         = "unknown";

    foreach my $attr ( keys %$attribute_hash ) {
        if ( $attr =~ /unpaired/i ) {
            $max_unpaired = $$attribute_hash{$attr};

        }
        elsif ( $attr =~ /data|vienna/i ) {
            $self->vienna( $$attribute_hash{$attr} );
            my $length = $self->findStems($max_unpaired);

        }
        elsif ( $attr =~ /mfold|ct/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->mfold2( $file, $max_unpaired );

        }
        elsif ( $attr =~ /rnaml/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->rnaml( $file, $max_unpaired );

        }
        elsif ( $attr =~ /plot/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->energyplot($file);

        }
        elsif ( $attr =~ /pair/ ) {
            my $pairs = $$attribute_hash{pair};
            my $nstem = $self->pairs( $pairs, $max_unpaired );

        }
        elsif ( $attr =~ /parray/ ) {
            my $parray = $$attribute_hash{parray};
            my $nstem  = $self->pairArray($parray);

        }
        else {
            $self->{$attr} = $attribute_hash->{$attr};
        }
    }

    return;
}

#-----------------------------------------------------------------------------
# findStems
#
# extract the stems from the vienna structure. left stem positions are pushed
# on a stack until the matching right parenthesis is found.  Stems are split
# in two when $max_unpaired unpaired bases are found in either side of the
# stem.
#
# USAGE
#   $nstems = $topology->findStems( $max_unpaired );
#
#-----------------------------------------------------------------------------
sub findStems {
    my ( $self, $max_unpaired ) = @_;

    my @stack;

    $max_unpaired++; # threshold is one larger than max number of unpaired bases

    # trim the vienna string and convert to an array

    my $vienna = $self->vienna;

    # print "vienna $vienna\n";
    $vienna =~ s/\s+//g;
    my @varray = split( '', $vienna );
    my $size = @varray;
    if ( $size == 0 ) { return undef; }

    my $nb = $self->setBracketTypes;

    # print "$nb brackets found\n";

    my $first_left  = -1;
    my $first_right = 0;
    my $old_left    = 0;
    my $old_right   = 0;
    my $old_type    = 0;
    my $right       = 0;
    my $left        = 0;

    my $pos   = 0;
    my $nstem = 0;

    foreach my $char (@varray) {
        $pos++;
        my $type = $index{$char};
        if ($type) {    # skip spaces (type=0)

            if ( $value{$char} > 0 ) {    # left parentheses have value>0
                push @{ $stack[$type] }, $pos;
            }
            elsif ( $value{$char} < 0 ) {    # right parentheses have value<0
                $right = $pos;
                $left  = pop @{ $stack[$type] };
                last
                  unless defined $left
                ; # this should only happen if parentheses are unbalanced, possibly because of un-annotated pseudoknot
                if (   $right - $old_right > $max_unpaired
                    or $old_left - $left > $max_unpaired
                    or $type != $old_type )
                {
                    unless ( $first_left < 0 ) {
                        $nstem++;
                        my $leftv = substr(
                            $vienna,
                            $old_left - 1,
                            $first_left - $old_left + 1
                        );
                        my $rightv = substr(
                            $vienna,
                            $first_right - 1,
                            $old_right - $first_right + 1
                        );

#print "stem $nstem: $left,$old_left,$first_left,$leftv to $right,$first_right,$old_right,$rightv\n";
                        my $stem = Stem->new(
                            [
                                $old_left,  $first_left, $first_right,
                                $old_right, $leftv,      $rightv
                            ]
                        );
                        $self->addStem($stem);
                    }
                    $first_left  = $left;
                    $first_right = $right;
                }
                $old_right = $right;
                $old_left  = $left;
                $old_type  = $type;
            }
        }
    }

    # get the last stem

    if ( defined $left ) {
        $nstem++;
        my $leftv =
          substr( $vienna, $old_left - 1, $first_left - $old_left + 1 );
        my $rightv =
          substr( $vienna, $first_right - 1, $old_right - $first_right + 1 );

#print "final stem $nstem: $left,$old_left,$first_left to $right,$first_right, $old_right\n";
        my $stem = Stem->new(
            [
                $old_left,  $first_left, $first_right,
                $old_right, $leftv,      $rightv
            ]
        );
        $self->addStem($stem);
    }

    return ($nstem);
}

# end of findStems

#-----------------------------------------------------------------------------
# setBrackettypes
#
# checks the vienna string and fill the brackets attribute with a list of
# the indices of the brackets found.  '.' and '-' are ignored.
#
# USAGE
#   $ntypes = setBracketTypes;
#-----------------------------------------------------------------------------
sub setBracketTypes {
    my ($self) = @_;

    my $found = 0;
    foreach my $i ( 1 .. $nbracket_types ) {    # remember: type 0 is '.'
        my $search = qq{\\$bracket[$i]};
        if ( $self->vienna =~ /$search/ ) {
            push @{ $self->{brackets} }, $i;
            $found++;
        }
    }

    return $found;
}

# end of setBracketTypes

#-----------------------------------------------------------------------------
# addStem
#
# This accessor method adds stem object to the stem list array.
# Output: return;
#
# USAGE
#     $Topology_obj->addStem($stem_obj);
#-----------------------------------------------------------------------------

sub addStem {
    my ( $self, $stem_obj ) = @_;

    push( @{ $self->stemList }, $stem_obj );
    return;
}

#-----------------------------------------------------------------------------
# addTopology
#
# This accessor method adds all the stems in another topology object to the stem list 
# array of the existing topology.
#
# USAGE
#     $Topology_obj->addTopology($another_topology_obj);
#-----------------------------------------------------------------------------
sub addTopology {
    my ( $self, $another_topo ) = @_;
    foreach my $stem_obj ( @{ $another_topo->stemList } )  {
        $self->addStem( $stem_obj );
    }
    return;
}

#-----------------------------------------------------------------------------
# printStemListOld
#
# This prints all the stem objects present in the stem List array of the
# Topology object This is the one originally written by Reazur
#
# Output: return;
#
# USAGE
#      $Topology_obj->printStemList();
#-----------------------------------------------------------------------------

sub printStemListOld {
    my ($self) = @_;
    foreach my $stem_obj ( @{ $self->stemList } ) {
        $stem_obj->printStem();    #not sure about this line
    }
    return;
}

#-----------------------------------------------------------------------------
# printStemList
#
# Prints all the stem objects present in the stem List array of the Topology
# object
#
# Output: returns an array of stems;
#
# USAGE
#      $Topology_obj->printStemList();
#-----------------------------------------------------------------------------

sub printStemList {
    my ($self) = @_;
    my @stems = ();
    foreach my $stem_obj ( @{ $self->stemList } ) {
        my $stem = $stem_obj->printStem();
        push @stems, $stem;
    }
    return (@stems);
}

#-----------------------------------------------------------------------------
# getStemList
#
# This accessor method returns the reference to the StemList Array
#
# USAGE
#      $arr_ref = $Topology_obj->getStemList();
#-----------------------------------------------------------------------------

sub getStemList {
    my ($self) = @_;
    return $self->stemList;
}

#--------------------------------------------------------------------------------
# pruneStems
#
# Given predicted and curated stems, find out the predicted stems that match the
# curated ones, and prune the non-matching predicted stems; return the predicted
# stems after pruning
# my $n_stem = $topology->pruneStems( $test_stems, $target_stems );
#--------------------------------------------------------------------------------
sub pruneStems {

    my ( $self, $test_stem, $target_stem ) = @_;
    print "target stem: $target_stem\n";
    my %matchlist;

    my @target_stem = @{$target_stem};
    my @test_stem   = @{$test_stem};

    foreach my $i ( 0 ... $#test_stem ) {
        my @a = $test_stem[$i]->getStemAsArray;

        print "a: @a\n";

#my ( $a1, $a2, $a3, $a4 ) = $test_stem[$i]->getStemAsArray; # lmin, lmax, rmin, rmax
        my ( $a1, $a2, $a3, $a4 ) = @a;

        foreach my $j ( 0 ... $#target_stem ) {

            my @b = $target_stem[$j]->getStemAsArray;
            my ( $b1, $b2, $b3, $b4 ) = @b;

            if (
                (
                    max( $a1, $a2, $b1, $b2 ) - min( $a1, $a2, $b1, $b2 ) <
                    abs( $a1 - $a2 ) + abs( $b1 - $b2 )
                )
                && ( max( $a3, $a4, $b3, $b4 ) - min( $a3, $a4, $b3, $b4 ) <
                    abs( $a3 - $a4 ) + abs( $b3 - $b4 ) )
              )
            {

                $matchlist{$i}++;
            }
        }
    }
    foreach my $i ( sort { $a <=> $b } keys %matchlist ) {
        $self->addStem( $test_stem[$i] );
    }

    return $self->nStem;
}

# End of pruneStems

#-------------------------------------------------------------------------------
# randomTopol
#
# generates random topologies for a given number of stems.
# Takes in the number of stems, and the number of random topologies to be
# generated. Returns a 2D array where each row is a new random topology.
#
# USAGE
#   @random_topologies = randomTopol( $n_stem, $n_rand_topol );
#
# 28 June 2007     Aditi
#-------------------------------------------------------------------------------
sub randomTopol {
    my ( $self, $n_stem, $n_rand ) = @_;
    my ( @n_list, @result ) = ();
    for ( my $i = 1 ; $i <= ( $n_stem * 2 ) ; $i++ ) {

        # generate the list of half-stems
        push @n_list, $i;
    }

    # call random with the generated list
    @result = random( $n_rand, @n_list );
    return @result;
}

#-------------------------------------------------------------------------------
# randomPartialBlock
#
# generates partially random topologies when provided with a partially specified
# graph. The non-random part of graph stays as a block on the left side of
# topology. Takes in number of stems, number of partially random topologies to
# be generated, and the partially specified topology. Returns a 2D array where
# each row is a new partially random topology.
#
# USAGE
#   @part_random_topol = randomPartialBlock( $n_stem, $n_rand, @fixed_topol );
#
# 29 June 2007     Aditi
#-------------------------------------------------------------------------------
sub randomPartialBlock {
    my ( $self, $n_stem, $n_rand, @part_fixed ) = @_;
    my ( @n_list, @rand_list, @result ) = ();
    for ( my $i = 1 ; $i <= ( $n_stem * 2 ) ; $i++ ) {

        # generate list of half stems
        push @n_list, $i;
    }

    # remove those half stems that are in the partially specified list of half
    # stems.
    for ( my $j = 0 ; $j < @part_fixed ; $j++ ) {
        for ( my $k = 0 ; $k < @n_list ; $k++ ) {
            if ( $n_list[$k] == $part_fixed[$j] ) {

                splice( @n_list, $k, 1 );
            }
        }
    }
    @rand_list = random( $n_rand, @n_list );
    for my $i ( 0 .. $#rand_list ) {

        # non-random part is at the begining of topology
        # concatenate the part specified and part random lists to give a
        # partially random topology
        @{ $result[$i] } = ( @part_fixed, @{ $rand_list[$i] } );

        #       print "part random topol: @{$result[$i]}\n";
    }

    return @result;
}

#-------------------------------------------------------------------------------
# random
#
# generates randomized versions of a given list(topology) - used by randomTopol,
# randomPartialBlock and randomPartialInterspersed.
# Takes in the number of randomized versions to be generated, and the list to be
# randomized.
# Returns a 2D array where each row is a new randomized version of the list.
#
# USAGE
#   @random_lists = random( $n_rand, @list );
#
# 29 June 2007     Aditi
#-------------------------------------------------------------------------------
sub random {
    my ( $self, $n_rand, @n_list ) = @_;

    my ( $pick, $range, $pos, $n );
    my ( @rand_list, @result ) = ();
    $n = scalar(@n_list);
    my @temp = @n_list;

    for ( my $j = 0 ; $j < $n_rand ; $j++ ) {

        # for each new random topology
        my @sort_list;
        @n_list    = @temp;
        @rand_list = ();

        # randomly pick an index from array of half-stems(@n_list), push the
        # corresponding number(half-stem) into @rand_list, delete that number
        # from array of half-stems, repeat till no numbers are left in @n_list
        for ( my $i = 0 ; $i < $n ; $i++ ) {
            $range = @n_list;
            $pick  = int( rand($range) );
            push @rand_list, $n_list[$pick];
            splice( @n_list, $pick, 1 );
        }

        # in @rand_list, sort within each pair of half-stems
        for ( my $z = 0 ; $z < @rand_list ; ) {
            ( $rand_list[$z], $rand_list[ $z + 1 ] ) =
              sort { $a <=> $b } ( $rand_list[$z], $rand_list[ $z + 1 ] );
            $z = $z + 2;
        }

        # push this sorted random topology as a new row in 2D array @result
        @sort_list = sortPairs(@rand_list);
        push @result, \@sort_list;
    }
    return @result;
}

#-------------------------------------------------------------------------------
# randomPartialInterspersed
#
# generates partially random graph where the non-random part is not contained as
# a block.  Takes in number of stems, number of partially random topologies to
# be generated, and the partially specified topology. Returns a reference for a
# 2D array- @result, where each row is a new partially random topology, and a
# reference to another 2D array- @mapping, which contains the half stems in the
# new random topology, that corresponds to the fixed graph specified initially.
#
# USAGE:
#   ($r_partial_topo, $ref_mapping) = randomPartialInterspersed ($num_stems,
#                                                                $num_random_str,
#                                                                @fixed_topol );
#
# 19 July 2007      Aditi
#-------------------------------------------------------------------------------
sub randomPartialInterspersed {
    my ( $self, $total_stem, $n_rand, @part_fixed ) = @_;

    my ( @result, @mapping ) = ();
    my $total_half_stem = $total_stem * 2;    # length of final topology
    my $fixed_half_stem = @part_fixed;        # length of fixed part of topology
    my $range = $fixed_half_stem + 1;
    for ( my $num = 0 ; $num < $n_rand ; $num++ ) {

        # for each new random topology
        my ( @n_list, @new_fixed, @rand_list, @rem_list, @one_result ) = ();

        # generate random numbers between 0 and $fixed_half_stem + 1
        for ( my $r = 0 ; $r < $total_half_stem - $fixed_half_stem ; $r++ ) {
            my $rand = rand($range);
            push @rand_list, $rand;
        }

        # add these random numbers to the list of stems and sort
        @n_list = ( @part_fixed, @rand_list );
        @n_list = sort { $a <=> $b } @n_list;

     # in the list, get the index of list number that is present in fixed list
     # push that index into @new_fixed - a random arrangement for the fixed list
        foreach my $fixed (@part_fixed) {
            for ( my $pos = 0 ; $pos < @n_list ; $pos++ ) {
                if ( $fixed == $n_list[$pos] ) {
                    push @new_fixed, $pos + 1;
                }
            }
        }
        push @mapping, \@new_fixed;

        # get the indices of remaining numbers, which are floating points
        # and push them into @rem_list (the remaining half stems)
        for ( my $pos = 0 ; $pos < @n_list ; $pos++ ) {
            if ( $n_list[$pos] =~ m/\d+\.\d+/ ) {
                push @rem_list, $pos + 1;
            }
        }

        # get random arrangement of numbers in @rem_list
        my @rem_random = random( 1, @rem_list );
        for my $i ( 0 .. $#rem_random ) {

            # combine @new_fixed and random arrangement of @rem_list
            # and sort the pairs in this final list
            my @temp = ( @new_fixed, @{ $rem_random[$i] } );
            @one_result = sortPairs(@temp);
        }
        push @result, \@one_result;

    }
    return ( \@result, \@mapping );

}

#------------------------------------------------------------------------------
# randomAddPair
#
# Randomly add a single new pair to an existing pair array.  this is used to add
# variable stems to a common core graph.
#
# for a pair array with n elements, there are n+1 intervals where the endpoints
# of a new graph could fall, for instance a to e for a 2 stem array
#
#   ___0___1___2___3___
#    a   b   c   d   e
#
# this function simple chooses two intervals for the location of the new stem
# and updates the previous numbers
#
# USAGE
#    my @parray_plus = $topology->randomAddPair( @parray );
#------------------------------------------------------------------------------
sub randomAddPair {
    my ( $self, @parray ) = @_;

    my $intervals = @parray;
    my @new;
    push @new, int( rand( $intervals + 1 ) );
    push @new, int( rand( $intervals + 1 ) );

    # make sure the smaller endpoint comes first in the pair
    if ( $new[0] > $new[1] ) { ( $new[0], $new[1] ) = ( $new[1], $new[0] ); }

    unshift @parray, @new;
    foreach my $i ( 2 .. $#parray ) {
        foreach my $newpos (@new) {
            if ( $parray[$i] >= $newpos ) {
                $parray[$i]++;
            }
        }
    }

    # if the new stem has both ends in the same interval, increment the value of
    # the right half
    if ( $parray[0] == $parray[1] ) { $parray[1]++; }

    return sortPairs(@parray);
}

# End of randomAddPair

#------------------------------------------------------------------------------
# randomAddPairWeighted
#
# Randomly add a single new pair to an existing pair array, but with a
# decreasing weight for each successively further away interval.  this is used
# to add variable stems to a common core graph.
#
# for a pair array with n elements, there are n+1 intervals where the endpoints
# of a new graph could fall, for instance a to e for a 2 stem array
#
#   ___0___1___2___3___
#    a   b   c   d   e
#
# this function simple chooses two intervals for the location of the new stem
# and updates the previous numbers
#
# USAGE
#    my @parray_plus = $topology->randomAddPairWeighted( $weight, @parray );
#------------------------------------------------------------------------------
sub randomAddPairWeighted {
    my ( $self, $factor, @parray ) = @_;

    my $intervals = @parray + 1;

    # first end point
    my $p1 = int( rand($intervals) );

    # weight the intervals so the second end point is more likely to be in the
    # same interval as the first (tries to be more like biological graph)
    my @weight;
    $weight[$p1] = 1.0;
    #print STDERR "p1 $p1 $weight[$p1] ints $intervals\n";
    foreach my $i ( $p1 + 1 .. $intervals - 1 ) {
        $weight[$i] = $weight[ $i - 1 ] * $factor;
    #    print STDERR "begin w $i $weight[$i]\n";
    }
    foreach my $i ( 1 .. $p1 ) {
        $weight[ $p1 - $i ] = $weight[ $p1 - $i + 1 ] * $factor;
        my $j = $p1 - $i;
    #    print STDERR "end w $j $weight[$j]\n";
    }

    # convert wieght to thresholds
    foreach my $i ( 1 .. $intervals - 1 ) {
        $weight[$i] += $weight[ $i - 1 ];
    }
    #foreach my $i ( 0 ... $intervals -1 ) {
    #    print STDERR "summary w $i $weight[$i]\n";
    #}

    # second end point
    my $r2 = rand( $weight[ $intervals - 1 ] );

    my $p2 = 0;
    while ( $weight[$p2] < $r2 ) {
        $p2++;
    }
    #print STDERR "pick thres $r2 $p2 $weight[$p2]\n";

    my @new = ( $p1, $p2 );

    # make sure the smaller endpoint comes first in the pair
    if ( $new[0] > $new[1] ) { ( $new[0], $new[1] ) = ( $new[1], $new[0] ); }
    #print STDERR "parray @parray\n";

    unshift @parray, @new;
    #print STDERR "new @new parray @parray\n";

    foreach my $i ( 2 .. $#parray ) {
        foreach my $newpos (@new) {
            if ( $parray[$i] >= $newpos ) {
                $parray[$i]++;
            }
        }
    }

    #print STDERR "newnew @parray\n";
    # if the new stem has both ends in the same interval, increment the value of
    # the right half
    if ( $parray[0] == $parray[1] ) { $parray[1]++; }
    #print STDERR "newnewnew @parray\n\n";

    return sortPairs(@parray);
}

# End of randomAddPairWeighted

#-----------------------------------------------------------------------------
# sortPairs
#
# for sorting pairs in an array- used by random and randomPartialInterspersed.
#
# USAGE:
#   @sorted_list = sortPairs(@unsorted_list);
#
# 17 July 2007        Gribskov
#-----------------------------------------------------------------------------
sub sortPairs {
    my @rand_list = @_;

    my ( @firstpos, @sort_list );
    foreach my $k ( 0 .. $#rand_list ) {
        next if ( $k % 2 );
        push @firstpos, $k;
    }

   # sort the pairs so that they are ordered by the first position of the paired
   # halfstem
    foreach my $pos ( sort { $rand_list[$a] <=> $rand_list[$b] } @firstpos ) {
        push @sort_list, @rand_list[ $pos, $pos + 1 ];
    }

    return @sort_list;
}

# end of sortPairs


#------------------------------------------------------------------------------
# stemsToPairArray
#
#   Given a list of stems, transfer into pair array format.
#
# USAGE
#   my $pairs = stemsToPairArray( \@stem_list );
#------------------------------------------------------------------------------
sub stemsToPairArray {
    my ($stems) = @_;
    my %coord;
    my @coord;
    my %pairs;
    my @pair = ();
    foreach my $stem ( @{$stems} ) {
        my ( $a1, $a2, $a3, $a4 ) =
          ( $stem->[0], $stem->[1], $stem->[2], $stem->[3] );
        push @coord, ([$a1, $a2], [$a3, $a4] );
        my $ll   = join ".", ( $a1, $a2 );
        my $rr   = join ".", ( $a3, $a4 );
        my $llrr = join ".", ( $ll, $rr );
        $coord{$ll}   = undef;
        $coord{$rr}   = undef;
        $pairs{$llrr} = undef;
    }

    @coord = sort { $a->[0] <=> $b->[0] } @coord;
    #print STDERR "coord: @coord\n";
    my $pos = 0;
    while ( $pos < @coord ) {
        my $cc = join ".", @{$coord[$pos]};
        $coord{$cc} = $pos;
        $pos++;
    }

    foreach my $llrr ( keys %pairs ) {

        my @cc = split /\./, $llrr;
        my $ll = join ".", @cc[ 0 ... 1 ];
        my $rr = join ".", @cc[ 2 ... 3 ];
        my $lc = $coord{$ll};
        my $rc = $coord{$rr};
        #print STDERR "ll $ll rr $rr cc @cc lc $lc rc $rc\n";
        push @pair, ( $lc, $rc );
    }

    return Topology::sortPairs(@pair);

}

# End of stemsToPairArray


#-----------------------------------------------------------------------------
# topologyToVienna
#
# generates dot-bracket (Vienna) representation of RNA structure from
# topology representation (list of stem pairs).
# Takes in number of dots and number of brackets desired in the Vienna
# representation, and the topology array. Returns a string containing
# Vienna representation.
#
# this is used to create a dummy vienna string for a given set of half-stem
# pairs
#
# USAGE:
#   $vienna_rep = topologyToVienna( $ndots, $nbrackets, @topology );
#
# 10 July, 2007
#-----------------------------------------------------------------------------
sub topologyToVienna {
    my ( $self, $ndot, $nbrackets, @topol ) = @_;

    my @out;
    my %bracket = (
        '(' => 0,
        ')' => 0,
        ']' => 1,
        '[' => 1,
        '{' => 2,
        '}' => 2,
        '>' => 3,
        '<' => 3,
    );

    my @pairs = ( [ "(", ")" ], [ "[", "]" ], [ "{", "}" ], [ "<", ">" ], );

    while (@topol) {

        # get the begin and end of a stem
        my $begin = shift @topol;
        my $end   = shift @topol;
        my @used;
        foreach my $i ( ( $begin - 1 ) .. ( $end - 1 ) ) {
            my $level;
            if ( defined $out[$i] ) {

                # if a bracket is present between begin and end of the stem
                # get the level of that bracket from %bracket
                $level = $bracket{ $out[$i] };
                $used[$level]++;
            }
        }

        my $found = -1;
        foreach my $level ( 0 .. $#pairs ) {

            # for all 4 levels of brackets in @pairs get the level in use
            if ( !defined( $used[$level] ) || !$used[$level] % 2 ) {
                $found = $level;
                last;
            }
        }
        if ( $found == -1 ) {
            print STDERR "begin=$begin   end=$end   found=$found\n";
            print STDERR "used=@used\n";
            print STDERR "out=@out\n";
            print STDERR "too many bracket levels\n";
            die;
        }
        else {

            # fill in the brackets according to the level in use
            $out[ $begin - 1 ] = $pairs[$found][0];
            $out[ $end - 1 ]   = $pairs[$found][1];
        }
    }

    # join individual brackets to get vienna representation string
    # expand structure by inserting dots and multiples of each bracket
    my $vienna = join( "", @out );
    foreach my $onebrace ( keys %bracket ) {
        my ( $from, $to ) = expandBrackets( $onebrace, $ndot, $nbrackets );
        $vienna =~ s/$from/$to/g;
    }
    return $vienna;
}

#-----------------------------------------------------------------------------
# expandBrackets
#
# create the strings needed to change a single brace to a series of dots and
# braces - used by topologyToVienna
#
# USAGE
#   ($to,$from) = expandBrackets( "{", $ndots, $nbrackets );
#
# 10 July 2007     Gribskov
#-----------------------------------------------------------------------------
sub expandBrackets {
    my ( $bracket, $ndots, $nbrackets ) = @_;
    my $in  = "\\$bracket";
    my $out = "";
    for ( 0 .. $ndots ) {
        $out .= ".";
    }
    for ( 0 .. $nbrackets ) {
        $out .= $bracket;
    }
    for ( 0 .. $ndots ) {
        $out .= ".";
    }

    return $in, $out;
}

#-----------------------------------------------------------------------------
# data
#
# retrieve the Vienna RNA representation of this structure.  This function
# does the same thing as $self->vienna but is less mnemonic. mainly included
# for back-compatability.
#
# USAGE
#   $vienna_str = $topology->vienna;
#   $topology->vienna( $new_vienna_str );
#-----------------------------------------------------------------------------
sub data {
    my $self = shift;
    return $self->data(@_);
}

#-----------------------------------------------------------------------------
# mfold
#
# DEPRECATED: Conversion of CT to topology by this function is not reliable
# for pseudoknots
#
# open a file and read the structure from a mfold/ct formatted file.  The ct
# file is converted to a vienna string and then parsed as a vienna string.
# this will not work if the structure contains pseduoknots.
#
# SAMPLE FILE:
#
#369     dG = -138.2     sacce S.cerevisiae RNase P RNA
#1       G       0       2       0       1       0       0
#2       U       1       3       367     2       0       3
#3       G       2       4       366     3       2       4
#4       G       3       5       365     4       3       5
#5       A       4       6       364     5       4       6
#6       A       5       7       363     6       5       7
#
# USAGE
#   $ n_stems  = $topology->mfold( $filename );
#-----------------------------------------------------------------------------
sub mfold {
    my ( $self, $filename ) = @_;

    my $n_stem = 0;

    open( CT, "<$filename" ) || die "unable to open ct file $filename\n\n";

    my $line = <CT>;
    my ( $length, $energy, $id, $doc ) = split " ", $line, 4;
    my @vienna_array;
    foreach my $i ( 0 .. $length - 1 ) {
        $vienna_array[$i] = ".";
    }

    my $seq = "";
    while ( $line = <CT> ) {

        my ( $pos, $base, $left, $x1, $right, $x2 ) = split " ", $line;

        $seq .= $base;

        # in CT file, pairs are listed from both directions
        next unless ( $right && $left < $right );
        $vienna_array[ $left - 1 ]  = "(";
        $vienna_array[ $right - 1 ] = ")";

    }
    my $vienna = join( "", @vienna_array );
    return ( $seq, $vienna );

    #print "vienna:$vienna_str\n";
}

#-----------------------------------------------------------------------------
# mfold
#
# open a file and read the structure from a mfold/ct formatted file.  The ct
# file is converted to a vienna string and then parsed as a vienna string.
# this will not work if the structure contains pseduoknots.
#
# SAMPLE FILE:
#
#369     dG = -138.2     sacce S.cerevisiae RNase P RNA
#1       G       0       2       0       1       0       0
#2       U       1       3       367     2       0       3
#3       G       2       4       366     3       2       4
#4       G       3       5       365     4       3       5
#5       A       4       6       364     5       4       6
#6       A       5       7       363     6       5       7
#
# USAGE
#   $ n_stems  = $topology->mfold( $filename );
#-----------------------------------------------------------------------------
sub mfold2 {
    my ( $self, $filename, $max_unpaired ) = @_;

    open( CT, "<$filename" ) || die "unable to open ct file $filename\n\n";

    my $line = <CT>;
    my ( $length, $dg, $equal, $energy, $id, $doc ) = split " ", $line, 6;
    $self->{energy} = $energy;

    #my @vienna_array;
    #foreach my $i ( 0 .. $length-1 ) {
    #    $vienna_array[ $i ] = ".";
    #}

    # get the basepairs from the CT file

    my @pair;

    my $seq   = "";
    my $npair = 0;
    while ( $line = <CT> ) {

        next if ( $line =~ /#/ );
        last if ( $line =~ /dG/ );    #stop after first structure

        my ( $left, $base, $pos, $x1, $right, $x2 ) = split " ", $line;
        next if ( ( !( defined($base) ) ) && ( !( defined($pos) ) ) );
        $seq .= $base;

        # in CT file, pairs are listed from both directions
        next unless ( $right && $left < $right );
        $pair[$npair] = [ $left, $right ];
        $npair++;
    }

    #foreach my $ii ( 0..$#pair ) {
    #    print " pair[$ii] = $pair[$ii]->[0]:$pair[$ii]->[1]\n";
    #}
    $self->sequence($seq);
    $self->sequenceId($filename);
    my @stemlist = basepairToStem(@pair);
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }

    return ($n_stem);
}

#------------------------------------------------------------------------------
# mfoldMultiple
#
# Read the all the suboptimal CT sructures in the specified file down to the
# delta deltaG limit.  If no limit is specified, all structures are read.  The
# list of stems is filtered to remove duplicates.
#
# SAMPLE FILE:
#
#369     dG = -138.2     sacce S.cerevisiae RNase P RNA
#1       G       0       2       0       1       0       0
#2       U       1       3       367     2       0       3
#3       G       2       4       366     3       2       4
#4       G       3       5       365     4       3       5
#5       A       4       6       364     5       4       6
#6       A       5       7       363     6       5       7
#
# USAGE
#    my $n_stem = $topology->mfoldMultiple( $filename, $max_unpaired, $ddG );
#------------------------------------------------------------------------------
sub mfoldMultiple {
    my ( $self, $filename, $max_unpaired, $ddG ) = @_;
    my $n_stem;

    open( my $in, "<", $filename )
      or die
      "Topology::mfoldMultiple, unable to open input CT file ($filename)\n\n";

    $self->sequenceId($filename);

    # read header for first structure

    my $line = <$in>;
    my ( $length, $dg, $equal, $mfe, $id, $doc ) = split " ", $line, 6;
    $self->length($length);

    my @stemlist;
    my ( $seq, @pair );
    while ( my $line = <$in> ) {
        next if ( $line =~ /#/ );    # skip comments
        next unless ( $line =~ /\w/ );    # skip blank lines

        if ( $line =~ /dG/ or $line =~ /ENERGY/ or eof ) {

            # process the current structure
            foreach my $p ( @pair ) {
    #            print STDERR qq{$p->[0]\t$p->[1]\n};
            }
            push @stemlist, basepairToStem(@pair);
            #print STDERR "$stemlist[0]\n";
            my $dnb = basepairToDotsAndBrackets( \@pair, $length );
            push my @comment_list, ($dnb);
            $self->comment( \@comment_list );

            my ( $length, $dg, $equal, $energy, $id, $doc ) = split " ", $line,
              6;
            #last if ( $energy > $mfe + $ddG );    # done if energy limit reached

            $seq = uc $seq;
            $self->sequence($seq) unless ($self->sequence);    # sequence for all structures should be the same
            $seq  = "";
            @pair = ();
        }
        else {

            # add to current structure.  Initially structures are stored as an
            # array of base-paired positions
            my ( $left, $base, $pos, $x1, $right, $x2 ) = split " ", $line;

            # not sure if the following test is necessary
            # next if ( (!(defined($base))) && (!(defined($pos) )) );

            $seq .= $base;

            # in CT file, pairs are listed from both directions
            next unless ( $right && $left < $right );
            push @pair, [ $left, $right ];
        }
    }

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    #print Dumper(\%unique ), "\n";

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of mfoldMultiple

#------------------------------------------------------------------------------
# bpseqToTopology
#
# Read the .bpseq file and return a stem list and a dots-and-bracket file
# with sequence
#
# SAMPLE FILE:
#
# Filename: a.I1.b.Bacteriophage.SP01.A2.DP.g31.bpseq
# Organism: Bacillus phage SPO1
# Accession Number: M37686
# Citation and related information available at http://www.rna.ccbb.utexas.edu
# 1 C 0
# 2 G 0
# 3 U 0
# 4 U 0
# 5 U 0
# 6 G 0
# 7 A 0
# 8 G 0
# 9 U 0
# 10 A 0
#
#
# USAGE
#    my $n_stem = $topology->bpseqToTopology( $filename, $max_unpaired );
#------------------------------------------------------------------------------
sub bpseqToTopology {
    my ( $self, $filename, $max_unpaired ) = @_;
    my $n_stem;

    open( my $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input BPSEQ file ($filename)\n\n";

    my ( $organism, $access_number );
    my ( $seq, $head, @pair, @stemlist );
    while ( my $line = <$in> ) {
        next if ( $line =~ /#/ );    # skip comments
        next unless ( $line =~ /\w/ );    # skip blank lines
        chomp $line;
        ($organism) = $line =~ /^Organism: (.*)/ if ( $line =~ /^Organism/ );
        ($access_number) = $line =~ /^Accession Number: (.*)/
          if ( $line =~ /^Accession Number: / );
        if (eof) {

            # process the current structure
            if ($head) {
                my $headn = length($head);
                print STDERR "   $headn leading N/n's removed from $filename\n";
            }

            #my ( $npair, $nseq ) = basepairPrune( \@pair, $seq );
	    my $npair = \@pair; 
	    my $nseq = $seq; 

            $self->sequence($nseq);
            my $length = length($nseq);
            $self->length($length);

            push @stemlist, basepairToStem( @{$npair} );
            my $dnb = basepairToDotsAndBrackets( $npair, $length );
            push my @comment_list, ($dnb);
            $self->comment( \@comment_list );
        }
        if ( $line =~ /^\d/ ) {

            # add to current structure.  Initially structures are stored as an
            # array of base-paired positions

            my ( $left, $base, $right ) = split " ", $line;
            if ( $base eq 'n' or $base eq 'N' ) {

                # remove leading N/n's from the sequence
                $head .= $base unless ($seq);
            }
            else {
                $seq .= $base;

                my $headn = 0;
                $headn = length($head) if ($head);
                $left  = $left - $headn;
                $right = $right - $headn;

                # in BPSEQ file, pairs are listed from both directions
                next unless ( $right && $left < $right );
                push @pair, [ $left, $right ];

            }

        }
    }
    my $id = ">$filename Organism: $organism Accession Number: $access_number";
    $self->sequenceId($id);

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of bpseqToTopology




#------------------------------------------------------------------------------
# ViennaToTopology
#
#  Reads in a .va file, convert it to RNA topology. 
#
# USAGE
#    $nstems = $rna->ViennaToTopology( $filename, $max_unpaired );
#------------------------------------------------------------------------------
sub ViennaToTopology {
    my ( $self, $filename, $max_unpaired ) = @_;
    open ( IN, "<$filename" ) || die "can not open $filename!\n";
    while ( my $line = <IN> ) {
	next unless ( $line =~ /\(/ or $line =~ /\)/ );
	chomp $line; 
	my $vienna = $line; 
	my $hash = {
	    vienna => $vienna, 
	};
	my $xrna = new Topology( $hash ); 
        $self->addTopology( $xrna );
    }
    close IN; 


    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );
    #print $self->XIOSWrite;
    return $self->nStem;
}
# End of VieenaToTopology 


#------------------------------------------------------------------------------
# HNumToTopology
#
# Read the .h-num file and return a stem list
# The stems being read in should have a h-num value lower than $h_num_threshold
#
# SAMPLE FILE:
#
# level   length  istart  jstart  h-num
# 1       5       29      43      3.4
# 1       4       42      53      3.25
# 1       4       9       63      3.25
# 1       5       51      65      3.2
#
#
# USAGE
#    my $n_stem = $topology->HnumToTopology( $filename, $max_unpaired,
#       $h_num_threshold, $h_num_precent );
#------------------------------------------------------------------------------
sub HNumToTopology {
    my ( $self, $filename, $max_unpaired, $h_num_max, $h_num_percent ) = @_;
    my $h_num_threshold;

    my $in;
    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    my @h_num;
    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;
        my ( $level, $length, $istart, $jstart, $h_num ) = split " ", $line;
        push @h_num, $h_num;
    }

    @h_num = sort { $a <=> $b } @h_num;
    my $nstem = scalar @h_num;
    print STDERR "$nstem stems read!\n";
    if ($h_num_percent) {
        # keep stems with the lowest 10% h-num values
        $h_num_threshold = $h_num[ int( @h_num * $h_num_percent / 100 ) ];
        print STDERR "Keeping $h_num_percent% stems!\n";
    }
    else {
        $h_num_threshold = $h_num_max;
    }
    print STDERR "Pruning stems using h-num threshold $h_num_threshold\n";
    close $in;

    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;

        my ( $level, $length, $istart, $jstart, $h_num ) = split " ", $line;
        my $left1        = $istart;
        my $left2        = $istart + $length - 1;
        my $right1       = $jstart - $length + 1;
        my $right2       = $jstart;
        my $vienna_left  = "(" x $length;
        my $vienna_right = ")" x $length;

        if ( $h_num > $h_num_threshold ) {
            print STDERR "Skip stem [$left1 $left2 $right1 $right2] with h-num value $h_num\n";
            next;
        }

        my $stem = new Stem;
        $stem->setStem(
            [ $left1, $left2, $right1, $right2, $vienna_left, $vienna_right ] );
        $self->addStem($stem);
    }

    close $in;

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of HNumToTopology

#------------------------------------------------------------------------------
# EnergyToTopology
#
# Read the .plot energy file and return a stem list
# The stems being read in should have an energy value lower than $dG_threshold
#
# SAMPLE FILE:
#
# level   length  i       j       energy
# 1       2       321     326     -1688
# 1       2       317     325     -1679
# 1       2       315     326     -1689
# 1       2       316     323     -1689
#
# USAGE
#    my $n_stem = $topology->EnergyToTopology( $filename, $max_unpaired, $dG_threshold );
#------------------------------------------------------------------------------
sub EnergyToTopology {
    my ( $self, $filename, $max_unpaired, $dG_threshold ) = @_;
    my %dG;

    my $in;

#open( $in, "<", $filename ) or
#   die "Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    #my @dG;
    #while ( my $line = <$in> ) {
    #    next if ( $line =~ /level/ );       # skip comments
    #    next unless ( $line =~ /\w/ );  # skip blank lines
    #    chomp $line;
    #    my ( $level, $length, $istart, $jstart, $dG ) = split " ", $line;
    #    push @dG, $dG;
    #}

    #@dG = sort { $a <=> $b } @dG;
    #unless ( $dG_threshold ) {
    #    $dG_threshold = $h_num[int(@h_num * 0.65)];
    #    print STDERR "Pruning stems using delta-G threshold $dG_threshold\n"
    #}
    #close $in;

    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;

        my ( $level, $length, $istart, $jstart, $dG ) = split " ", $line;
        my $left1        = $istart;
        my $left2        = $istart + $length - 1;
        my $right1       = $jstart - $length + 1;
        my $right2       = $jstart;
        my $vienna_left  = "(" x $length;
        my $vienna_right = ")" x $length;

        my $id = $left1 . '.' . $left2 . '.' . $right1 . '.' . $right2;
        $dG{$id} = $dG;

#if ( $dG > $dG_threshold ) {
#    print STDERR "Skip stem [ $left1 $left2 $right1 $right2 ] with delta-G value $dG\n";
#    next;
#}

        my $stem = new Stem;
        $stem->setStem(
            [ $left1, $left2, $right1, $right2, $vienna_left, $vienna_right ] );
        $self->addStem($stem);
    }

    close $in;

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    #return $self->nStem;
    my $nstem = $self->nStem;

    return ( $nstem, \%dG );

}

# End of EnergyToTopology

#-----------------------------------------------------------------------------
# pairs
#
# read in a structure from an array of paired bases.  basically the same as
# except the paired bases are provided in a data structure of
#
#   [ [base1, base1'], [base2, base2'] ...]
#
# USAGE
#   $ n_stems  = $topology->pairs( $filename );
#-----------------------------------------------------------------------------
sub pairs {
    my ( $self, $pairs, $max_unpaired ) = @_;

    my @stemlist = basepairToStem( @{$pairs} );
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }

    return ($n_stem);
}

#-----------------------------------------------------------------------------
# pairArray
#
# read in a structure from an the paired array format.  this allows an abstract
# topology to be read.  A paired array might be something like
# @p = ( 1, 5, 2, 4, 3, 6 ), indicating three stems: stem 1 pairing position 1
# and 5, stem 2 pairing postions 2 and 4, and stem 3 pairing positions 3 and 6.
#
# USAGE
#   $ n_stems  = $topology->pairArray( $ref_pair_array );
#-----------------------------------------------------------------------------
sub pairArray {
    my ( $self, $parray ) = @_;

    # stem = { 'left1'  => left start,
    #          'left2'  => left end,
    #          'right1' => right start,
    #          'right2' => right end
    #          'vienna_left' => vienna string for left side
    #          'vienna_right => vienna string for right side
    #         }
    my $n_stem = 0;

    my $pos = 0;
    while ( $pos < @$parray ) {
        my $stem = Stem->new(
            {
                left1        => $$parray[$pos],
                left2        => $$parray[$pos],
                right1       => $$parray[ $pos + 1 ],
                right2       => $$parray[ $pos + 1 ],
                vienna_left  => '(',
                vienna_right => ')'
            }
        );
        $pos += 2;
        $self->addStem($stem);
        $n_stem++;
    }

    return ($n_stem);
}

#----------------------------------------------------------------------------
# energyplot
#
# open energy dot plot file output by unafold(mfold) and read out stems, sorts
# them by energy. Stems that are in conflict (both substems) with lower-energy
# stems are dropped.
#
# USAGE
#      $n_stems = $topology->energyplot( $file );
#---------------------------------------------------------------------------
sub energyplot {
    my ( $self, $file ) = @_;

    my @stems      = ();
    my @uniq_stems = ();
    open( PLOT, "$file" ) or die "can't open plotfile $file: $!\n";
    while (<PLOT>) {
        if ( $_ =~ /\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(.+)/ ) {
            if ( $1 > 1 ) {    # skip one base stems
                    # stem: left start, right end, stem length and energy
                push @stems, [ $2, $3, $1, $4 ];
            }
        }
    }

    # sort stems by energy, left start, right end
    my @sorted_stems =
      sort { $a->[3] <=> $b->[3] || $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
      @stems;
    my $best_e = $sorted_stems[0][3];

    for my $newstem (@sorted_stems) {
        my $ls_new = ${$newstem}[0];
        my $le_new = ${$newstem}[0] + ${$newstem}[2] - 1;
        my $re_new = ${$newstem}[1];
        my $rs_new = ${$newstem}[1] - ${$newstem}[2] + 1;

        if ( ${$newstem}[3] == $best_e ) {
            push @uniq_stems, [ $ls_new, $le_new, $rs_new, $re_new, "", "" ];
            next;
        }
        else {
            my $flag = 0;    # flag for whether the new stem is unique or not
            for my $uniqstem (@uniq_stems) {
                if ( $ls_new <= ${$uniqstem}[1] && $le_new >= ${$uniqstem}[0] )
                {

                    # left substems overlap
                    if (   $rs_new <= ${$uniqstem}[3]
                        && $re_new >= ${$uniqstem}[2] )
                    {

                        # right substems overlap
                        $flag = 1;
                        last;
                    }
                }
            }
            if ( $flag != 1 ) { # if there is no overlap with prev stems, add it
                push @uniq_stems,
                  [ $ls_new, $le_new, $rs_new, $re_new, "", "" ];
                $flag = 0;
            }
        }
    }
    my $n_stem = 0;
    for my $stem_coor (@uniq_stems) {
        my $stem = Stem->new($stem_coor);
        $self->addStem($stem);

        #      print "@$stem_coor\n";
        #      print "stem obj: $stem\n";
        #      $stem->dump;
        #      print " stem $n_stem =>  @$stem_coor \n";
        $n_stem++;
    }

    #  print "Num stems: $n_stem\n";
    return ($n_stem);
}

#-----------------------------------------------------------------------------
# rnaml
#
# read a structure in rnaml format. Will not corrrectly read multiple
# structures in one file.
#
# USAGE
#   $n_stems = $topology->rnaml( $rnaml_file );
#-----------------------------------------------------------------------------
sub rnaml {
    my ( $self, $filename, $max_unpaired ) = @_;

    open( RNAML, "<$filename" )
      || die "Unable to open RNAML file $filename\n\n";

    my $in_sequence  = 0;
    my $in_data      = 0;
    my $in_structure = 0;

    my $sequence = "";
    my $seqlen   = 0;
    my @pairs;

    while ( my $line = <RNAML> ) {
        next if ( $line =~ /^\s*$/ );
        chomp $line;

        if ( $line =~ /<\s*rna/ ) {
            my ($title) = $line =~ /<rna\s+name\s*=\s*\"([^"]*)\"/;
            my ( $id, $doc ) = split " ", $title, 2;
            $self->sequenceId($id);
            $self->sequenceDoc($doc);
        }
        if ( $line =~ /<\s*sequence/ ) { $in_sequence = 1; }
        if ( $line =~ /<\s*data/ ) {
            if ($in_sequence) { $in_data = 1; next; }
        }
        if ( $line =~ /<\s*structure/ ) { $in_structure = 1; next; }

        # closing tags

        # read only the first structure
        if ( $line =~ /<\s*\/\s*rnaml/ )     { last; }
        if ( $line =~ /<\s*\/\s*sequence/ )  { $in_sequence = 0; }
        if ( $line =~ /<\s*\/\s*data/ )      { $in_data = 0; }
        if ( $line =~ /<\s*\/\s*structure/ ) { $in_structure = 0; }

        # data element (raw string)

        if ($in_data) {

            $line =~ s/\W//g;   # remove non-alphabetic characters from sequence
            $sequence .= $line;
            $seqlen = length($sequence);

  # structure elements, e.g.,  <base-pair _5p-base-num="2"  _3p-base-num="403"/>

        }
        elsif ($in_structure) {

            my ( $l, $r ) = $line =~
/base-pair\s*_5p-base-num\s*=\s*\"(\d+)\"\s*_3p-base-num\s*=\s*\"(\d+)\"/;
            unless ( defined($l) && defined($r) ) {
                print "undefined l or r:$line\n";
            }

            if ( $l < $r ) {
                push @pairs, [ $l, $r ];
            }

        }
    }
    close RNAML;

    $self->sequence($sequence);
    my @stemlist = basepairToStem(@pairs);
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }
    $self->sequence($sequence);

    return ($n_stem);
}

# end of rnaml

#-----------------------------------------------------------------------------
# basepairToStem
#
# Convert a list of paired bases to an array of stems.  Basepairs are provided
# as an array of pairs, each pair an array of the positions of the 5' and 3'
# bases.
#
# internally, the pairs are entered into a position index where each element
# of the index array is either undefined (unpaired) or has the value to the
# index of that base in the pair array.
#
# USAGE
#   @stems = basepairToStem( @basepairs );
#-----------------------------------------------------------------------------
sub basepairToStem {
    my (@pair) = @_;

    # make the index

    my $n_pair = 0;
    my @index;

    foreach my $pair (@pair) {
        my ( $l, $r ) = @{$pair};
        $index[$l] = $n_pair;
        $index[$r] = $n_pair;
        $n_pair++;
    }

    # find the stems

    my @stemlist;
    my $n_stem = 1;
    my $stem   = new Stem;
    my ( $left, $right ) = @{ shift @pair };
    $stem->setStem( [ $left, $left, $right, $right, "(", ")" ] );

    # get each pair from stack

    while ( my $pair = shift @pair ) {

        my ( $left_new, $right_new ) = @{$pair};

 # check if any pair has an index between old and new position on left and right

        my $need_new_stem = 0;
        my $left_vienna   = "";
        for ( my $l = $left + 1 ; $l < $left_new ; $l++ ) {
            $left_vienna .= ".";
            if ( defined( $index[$l] ) ) {
                $need_new_stem = 1;
                last;
            }
        }
        $left_vienna .= "(";

        # in case of pseudoknot, the new r may be larger than the old.
        # don't need to check this for left side because pair list is in order
        # of left base (I think)

        my $right_vienna = "";
        my $r0           = $right;
        my $r1           = $right_new;
        if ( $right_new > $right ) {
            $r0 = $right_new;
            $r1 = $right;
        }
        for ( my $r = $r0 - 1 ; $r > $r1 ; $r-- ) {
            $right_vienna .= ".";
            if ( defined( $index[$r] ) ) {
                $need_new_stem = 1;
                last;
            }
        }
        $right_vienna = ")" . $right_vienna;

        # if either search found an indexed position before finding its pair,
        # start a new stem

        if ($need_new_stem) {

            push @stemlist, $stem;
            $stem = new Stem;
            $stem->setStem(
                [ $left_new, $left_new, $right_new, $right_new, "(", ")" ] );
            $n_stem++;

        }
        else {

            # if no other index was found this is a continuing stem

            $stem->left2($left_new);
            $stem->right1($right_new);

            my $lv = $stem->viennaLeft . $left_vienna;
            $stem->viennaLeft($lv);

            my $rv = $right_vienna . $stem->viennaRight;
            $stem->viennaRight($rv);

        }

        $left  = $left_new;
        $right = $right_new;

    }
    push @stemlist, $stem;
    #my @values = values %$stem;
    #print STDERR "@values\n";

    return (@stemlist);
}

# end of basepairToStem

#-----------------------------------------------------------------------------
# basepairToDotsAndBrackets
#
# Convert a list of paired bases to a sequence of dots-and-brackets.
#
# USAGE
#   $dnb = basepairToDotsAndBrackets( @basepairs, $length );
#-----------------------------------------------------------------------------
sub basepairToDotsAndBrackets {
    my ( $pair, $length ) = @_;
    my @p = @{$pair};

    # make the index
    my %index;
    foreach my $p (@p) {
        my ( $l, $r ) = @{$p};
        $index{$l} = "(";
        $index{$r} = ")";
    }

    my $db;
    for my $i ( 1 ... $length ) {
        if ( $index{$i} ) {
            $db .= $index{$i};
        }
        else {
            $db .= ".";
        }
    }

    return ($db);
}

# end of basepairToDotsAndBrackets



#-----------------------------------------------------------------------------
# stemsToDotsAndBrackets
#
# Convert a list of stems to a sequence of dots-and-brackets.
#
# USAGE
#   $dnb = stemsToDotsAndBrackets( @stems, $length );
#-----------------------------------------------------------------------------
sub stemsToDotsAndBrackets {
    my ( $pair, $length ) = @_;
    my @p = @{$pair};

    # make the index
    my %index;
    foreach my $p (@p) {
        my ( $l, $r ) = @{$p};
        $index{$l} = "(";
        $index{$r} = ")";
    }

    my $db;
    for my $i ( 1 ... $length ) {
        if ( $index{$i} ) {
            $db .= $index{$i};
        }
        else {
            $db .= ".";
        }
    }

    return ($db);
}

# end of stemsToDotsAndBrackets



#-----------------------------------------------------------------------------
# basepairPrune
#
# Given a list of paired bases, prune the redundant bases
#
# USAGE
#   $basepairs = basepairPrune( @basepairs, $sequence );
#-----------------------------------------------------------------------------
sub basepairPrune {
    my ( $pair, $seq ) = @_;
    my $length = length($seq);
    my @seq = split "", $seq;
    unshift @seq, "n";

    $pair = [ sort { $a->[0] <=> $b->[0] } @{$pair} ];
    my $l = $pair->[0]->[0];
    my ( $left, $leftmv );
    if ( $l > 5 ) {
        $left = $l - 5;
    }
    else {
        $left = 1;
    }
    $leftmv = $left - 1;

    $pair = [ sort { $b->[1] <=> $a->[1] } @{$pair} ];
    my $r = $pair->[0]->[1];
    my ( $right, $rightmv );
    if ( $length - $r > 5 ) {
        $right = $r + 5;
    }
    else {
        $right = $length;
    }
    $rightmv = $length - $right;
    print STDERR
      "   $leftmv (left) and $rightmv (right) redundant bases removed\n";

    my $nseq;
    for my $s ( $left ... $right ) {
        $nseq .= $seq[$s];
    }

    for my $p ( @{$pair} ) {
        $p->[0] = $p->[0] - $leftmv;
        $p->[1] = $p->[1] - $leftmv;
    }

    return ( $pair, $nseq );

}

# end of basepairPrune

#-----------------------------------------------------------------------------
# splitStems
#
# check each stem in a list of stems and see if it has more than $max_unpaired
# unpaired bases.  If yes, break it into two stems.
#
# USAGE
#   @newstems = splitStems( $max_unpaired, @stems ); # not an object function
#-----------------------------------------------------------------------------
sub splitStems {
    my ( $max_unpaired, @stems ) = @_;

    my @new_stems;

    # target is a string with $max_unpaired dots

    my $target = qq{\\.\{$max_unpaired\}};

    while ( my $stem = shift @stems ) {

        # if stem is OK, remove from list and put in output

        my $left  = $stem->viennaLeft;
        my $right = $stem->viennaRight;

        unless ( $left =~ /$target/ || $right =~ /$target/ ) {

            #neither half stem has >= max unpaired contiguous unpaired bases
            push @new_stems, $stem;
            next;
        }

       # use whichever run of dots is longer to split on, the template is the
       # split vienna string from the left or right with the longer run of dots.

        #print "\nsplitting $left     $right\n";
        $left =~ /\.{$max_unpaired,}/;
        my $ll       = $`;
        my $lc       = $&;
        my $lr       = $';
        my $leftdots = length($&);

        $right =~ /\.{$max_unpaired,}/;
        my $rl        = $`;
        my $rc        = $&;
        my $rr        = $';
        my $rightdots = length($&);

        my $template_left  = $ll;
        my $template_dots  = $lc;
        my $template_right = $lr;
        if ( $rightdots > $leftdots ) {
            $template_left  = $rl;
            $template_dots  = $rc;
            $template_right = $rr;
        }

    # split the stem.  the regex extract the portions with the right number
    # of parens skipping padding dots.  left and right cases are not symmetrical
    # becasue of greedy matching by regex engine.

        my $nparenl = $template_left =~ tr/()/()/;

        #print " left template:$template_left   n:$nparenl\n";
        my ($newl1) = $left =~ /((\.*[()]){$nparenl})/;

        #print "      new left:$newl1    from $ll\n";
        my ($newr1) = $right =~ /(([()]\.*){$nparenl})$/;

        #print "     new right:$newr1    from $rl\n";

        my $nparenr = $template_right =~ tr/()/()/;

        #print "right template:$template_right   n:$nparenr\n";
        my ($newl2) = $left =~ /(([()]\.*){$nparenr})$/;

        #print "      new left:$newl2\n";
        my ($newr2) = $right =~ /((\.*[()]){$nparenr})/;

        #print "     new right:$newr2\n";

        # get coordinates of split stems.  push new stems back on stack in case
        # they need further splitting

        my ( $ls, $le, $rs, $re, $vl, $vr ) = $stem->getStemAsArray;
        my $stem1 = new Stem;
        $stem1->left1($ls);
        $stem1->left2( $ls + length($newl1) - 1 );
        $stem1->right1( $re - length($newr1) + 1 );
        $stem1->right2($re);
        $stem1->viennaLeft($newl1);
        $stem1->viennaRight($newr1);
        push @stems, $stem1;

        my $stem2 = new Stem;
        $stem2->left1( $le - length($newl2) + 1 );
        $stem2->left2($le);
        $stem2->right1($rs);
        $stem2->right2( $rs + length($newr2) - 1 );
        $stem2->viennaLeft($newl2);
        $stem2->viennaRight($newr2);
        push @stems, $stem2;

    }

    return @new_stems;
}

# end of splitStems

#------------------------------------------------------------------------------
# nStem
#
# Returns the number of stems stored in the topology.
#
# USAGE
#    $nstem = $topology->nStem;
#------------------------------------------------------------------------------
sub nStem {
    my ($topology) = @_;

    my $stemlist = $topology->stemList;

    return @{$stemlist};
}

# End of nStem

#-----------------------------------------------------------------------------
# dump
#
# Print the contents of the topology object.
#
# USAGE
#   $topology->dump;
#-----------------------------------------------------------------------------
sub dump {
    my ($self) = @_;

    print "vienna string: ", $self->vienna, "\n";

    my $i = 0;
    foreach my $stem ( @{ $self->stemList } ) {
        $i++;
        print "stem $i => $stem\n";
        $stem->dump;
    }
    return;
}

#-----------------------------------------------------------------------------
# dumpToFile
#
# Print the contents of the topology object to the lexical filehandle passed
# as an argument.
#
# USAGE
#   $topology->dumpToFile( $FH );
#-----------------------------------------------------------------------------
sub dumpToFile {
    my ( $self, $FH ) = @_;

    print $FH "vienna string: ", $self->vienna, "\n";

    my $i = 0;
    foreach my $stem ( @{ $self->stemList } ) {
        $i++;
        print $FH "stem $i => $stem\n";
        for my $key ( keys %$stem ) {
            print $FH "$key => $stem->{$key}\n";
        }
        print $FH "\n";
    }
    return;
}

# End of dumpToFile

#------------------------------------------------------------------------------
# XIOSWrite
#
# write a XIOS file describing the topology.  This functions is the
# authoritative definition of the format of a XIOS file.  Make sure that any
# changes made here are reflected in readXIOS.
#
# A XIOS file has up to three sections.  Each section is delimited by an XML
# like tag so that it is easier to read them independently.
#
# USAGE
#   $xios_string = $topology->XIOSWrite( );
#   $xios_string = $topology->XIOSWrite( $filename );
#------------------------------------------------------------------------------
sub XIOSWrite {
    my ( $self, $filename ) = @_;
    my $str = "";

    $self->adjacencyMatrix;

    $str .= qq{<XIOS version='$VERSION'>\n};

    $str .= $self->informationFormat;
    $str .= $self->stemlistFormat;
    $str .= $self->edgelistFormat;
    $str .= $self->adjacencyFormat;

    $str .= qq{</XIOS>\n};

    if ( defined $filename ) {

        # if no filname is provided just return the string
        open( my $out, '>', $filename )
          or die "Topology::XIOSWrite unable to open $filename for output\n\n";
        print $out $str;
    }

    return $str;
}

# End of XIOSWrite

#------------------------------------------------------------------------------
# writeXIOS
#
# forwards to XIOSwrite.  For back compatibility, do not use.
#
# USAGE
#   $xios_string = $topology->writeXIOS( );
#   $xios_string = $topology->writeXIOS( $filename );

#------------------------------------------------------------------------------
sub writeXIOS {
    my ( $self, $filename ) = @_;

    return $self->XIOSWrite($filename);
}

# End of writeXIOS

#------------------------------------------------------------------------------
# XIOSRead
#
# Read a XIOS file into a topology.  The XIOS file comprises three sections, a
# list of stems (stemlistRead), a list of edges (not needed), and an adjacency
# matrix (adjacencyRead).
#
# USAGE
#    $topology->XIOSRead( $filename );
#------------------------------------------------------------------------------
sub XIOSRead {
    my ( $topology, $filename ) = @_;

    open( my $in, "<", $filename )
      or die "Topology::XIOSRead, unable to open XIOS file ($filename)\n\n";

    my $file = new XML::Simple;
    my $xios = $file->XMLin($filename);

    # TODO add reading for sequence information (<information> block)
    my $seq = $xios->{information}->{sequence}->{unformatted};
    if ($seq) {
    $topology->sequence($seq);
    $topology->sequenceId($filename);
    }

    my $nstem = $topology->stemlistRead( $xios->{stem_list} );

    #print "$nstem stems read\n";

    my $nrow = $topology->adjacencyRead( $xios->{adjacency} );

    return $nstem;
}

# End of XIOSRead

#-------------------------------------------------------------------------------
# informationFormat
#
# Format a section with other information which may include the sequence,
# provenence, or other metadata
#
# usage
#	$formatted = $topology->informationFormat();
#-------------------------------------------------------------------------------
sub informationFormat {
    my ($topology) = @_;
    my $str = qq{  <information>\n};

    # comments
    my @commentlist = @{ $topology->comment };
    if (@commentlist) {
        $str .= qq{    <comment_list>\n};
        foreach my $comment (@commentlist) {
            $str .= qq{      <comment>$comment</comment>\n};
        }
        $str .= qq{    </comment_list>\n};
    }

    # sequence
    my $seq    = $topology->sequence;
    my $seqid  = $topology->sequenceId;
    my $seqdoc = $topology->sequenceDoc;
    if ( $seq || $seqid || $seqdoc ) {
        $str .= qq{    <sequence>\n};
        if ($seqid) {
            $str .= qq{      <id>$seqid</id>\n};
        }
        if ($seqdoc) {
            $str .= qq{      <doc>$seqdoc</doc>\n};
        }
        if ($seq) {
            $seq =~ s/\W//g;
            $str .= qq{      <unformatted>$seq</unformatted>\n};
        }
        $str .= qq{    </sequence>\n};
    }

    $str .= qq{  </information>\n};
    return $str;
}

# end of informationFormat

#------------------------------------------------------------------------------
# stemlistFormat
#
# Return a string with a formatted version of the stem_list component of the
# XIOS file.  The column widths are determined from the data so that all the
# information aligns nicely in columns.
#
# 0 184.5 [     2     9   360   367 ]    ((((((((    ))))))))
#
# USAGE
#    $formatted = $topology->stemlistFormat( );
#------------------------------------------------------------------------------
sub stemlistFormat {
    my ($self) = @_;
    my $str = "  <stem_list>\n";

    # find the sizes of each column
    my @size = ( 0, 0, 0, 0, 0, 0 );
    foreach my $stem ( @{ $self->stemList } ) {

        # $stem is a stem object
        if ( $size[0] < length $stem->{left1} ) {
            $size[0] = length $stem->{left1};
        }
        if ( $size[1] < length $stem->{left2} ) {
            $size[1] = length $stem->{left2};
        }
        if ( $size[2] < length $stem->{right1} ) {
            $size[2] = length $stem->{right1};
        }
        if ( $size[3] < length $stem->{right2} ) {
            $size[3] = length $stem->{right2};
        }
        if ( $size[4] < length $stem->{vienna_left} ) {
            $size[4] = length $stem->{vienna_left};
        }
        if ( $size[5] < length $stem->{vienna_right} ) {
            $size[5] = length $stem->{vienna_right};
        }
    }

    my $nstem = @{ $self->stemList };
    $size[6] = length $nstem;
    my $max_stem = $size[1];
    if ( $size[3] > $max_stem ) { $max_stem = $size[3] }
    my $dsize = $max_stem + 2;    # field width for stem center

    # make up the format string
    my $fmt =
"    \%$size[6]d \%$dsize.1f [ \%$size[0]d \%$size[1]d   \%$size[2]d \%$size[3]d ]  \%$size[4]s   \%-$size[5]s\n";

    # print "format:$fmt\n";

    my $n = 0;
    foreach my $s ( sort { $a->{left1} <=> $b->{left1} } @{ $self->stemList } )
    {

        # stems are sorted by left1
        my $center = ( $s->{left2} + $s->{right1} ) / 2.0;
        $str .= sprintf $fmt, $n, $center,
          $s->{left1},       $s->{left2},
          $s->{right1},      $s->{right2},
          $s->{vienna_left}, $s->{vienna_right};
        $n++;
    }
    $str .= "  </stem_list>\n";
    return $str;
}

# End of stemlistFormat

#-------------------------------------------------------------------------------
# edgelistFormat
#
# # Return a string with a formatted version of the edges of the topology for
# edgelist component of the XIOS file.
#
# usage
#	$str = $topology->edgelistFormat;
#-------------------------------------------------------------------------------
sub edgelistFormat {
    my ($self) = @_;
    my $str = "  <edge_list>\n";

    my @adj  = @{ $self->adjacency };
    my $wid  = length $#adj;
    my $fmtd = "    %$wid" . "d:";
    my $fmts = " %$wid" . "d%s";
    foreach my $i ( 0 .. $#adj ) {
        $str .= sprintf $fmtd, $i;
        foreach my $j ( $i + 1 .. $#adj ) {
            next if ( $adj[$i][$j] eq 's' );

            $str .= sprintf $fmts, $j, $adj[$i][$j];
        }
        $str .= "\n";
    }
    $str .= "  </edge_list>\n";

    return $str;
}

# end of edgelistFormat





#------------------------------------------------------------------------------
# stemlistRead
#
# Read a set of stems from the from the <stem_list> section of the XIOS format.
# The XML element <stem_list </stem_list> should already have been read and
# removed.
#
# format of each line is
#     0 203.0 [   1  16   390 404 ]  (((((((((((..(((   ))).)))))))))))
#
# USAGE
#    my $n_stem = $topology->stemlistRead( $text );
#------------------------------------------------------------------------------
sub stemlistRead {
    my ( $topology, $text ) = @_;

    foreach my $line ( split "\n", $text ) {
        my @field = split " ", $line;
        next unless ( $line =~ /\d/ );    # skip blank lines

        my $stemin = {
            left1        => $field[3],
            left2        => $field[4],
            right1       => $field[5],
            right2       => $field[6],
            vienna_left  => $field[8],
            vienna_right => $field[9]
        };
        unless ( defined $$stemin{vienna_left} ) { $$stemin{vienna_left} = "" }
        unless ( defined $$stemin{vienna_right} ) {
            $$stemin{vienna_right} = "";
        }
        if ( $$stemin{right2} > $topology->length ) {
            $topology->length( $$stemin{right2} );
        }

        my $stem = new Stem($stemin);
        $topology->addStem($stem);
    }

    return $topology->nStem;
}

# End of stemlistRead

#------------------------------------------------------------------------------
# stemListAsArray
#
# Returns the stemlist as an array of arrays.  the fields of the Stem objects,
# left1, left2, right1, right2 are returned as the elements of an array.
#
# USAGE
#    $stemlist_array = $topology->stemListAsArray;
#------------------------------------------------------------------------------
sub stemListAsArray {
    my ($topology) = @_;
    my $stemarray;

    foreach my $stem ( @{ $topology->stemList } ) {
        push @{$stemarray},
          [ $stem->left1, $stem->left2, $stem->right1, $stem->right2 ];
    }

    return $stemarray;
}

# End of stemListAsArray

#------------------------------------------------------------------------------
# adjacencyMatrix
#
# Calculates an adjacency matrix for the topology.  Each pair of stems are
# compared and classified as
#
# x exclusive, the all or part of the base-paired regions use the same bases
# o overlapping, the stems are pseudoknotted with respect to each other
# s serial, the beginning of the second stem is after the end of the first
# i included, one stem is nested within the other.  for the adjacency matrix,
#   if the stem represented by the row is the outer stem, the stem will be
#   will be labeled i, if it is the inner stem, it is labeled j.
#
# Example
# (((....(((..[[....))).]]..)))
# AAA....BBB..CC....bbb.cc..aaa
# there are three stems, with the basepaired regions shown in capital/lower case
# pairs.  The adjacency matrix is
#
#       A   B   C
#   A   -   i   i
#   B   j   -   o
#   C   j   o   -
#
# USAGE
#    $adjacency_matrix = $topology->adjacencyMatrix;
#------------------------------------------------------------------------------
sub adjacencyMatrix {
    my ($self) = @_;
    my $adj;

    my $relation;
    my @stem = sort { $a->{left1} <=> $b->{left1} } @{ $self->stemList };
    foreach my $i ( 0 .. $#stem ) {
        my $s = $stem[$i];
        my ( $a1, $a2, $a3, $a4 ) =
          ( $s->left1, $s->left2, $s->right1, $s->right2 );
        my @row;
        foreach my $j ( 0 .. $#stem ) {
            if ( $i == $j ) {
                push @row, '-';
                next;
            }

            my $s = $stem[$j];
            my ( $b1, $b2, $b3, $b4 ) =
              ( $s->left1, $s->left2, $s->right1, $s->right2 );

            if ( $b1 > $a4 || $a1 > $b4 ) {

                # s edge
                $relation = 's';
            }
            elsif ( $b1 > $a2 && $b4 < $a3 ) {

                # i edge
                $relation = 'i';
            }
            elsif ( $a1 > $b2 && $a4 < $b3 ) {

                # j edge
                $relation = 'j';
            }
            elsif (( $a2 < $b1 && $a3 > $b2 && $a4 < $b3 )
                || ( $b2 < $a1 && $b3 > $a2 && $b4 < $a3 ) )
            {

                # o edge
                $relation = 'o';
            }
            else {

                # x edge
                $relation = 'x';
            }
            push @row, $relation;
        }
        push @{$adj}, \@row;
    }
    $self->{adjacency} = $adj;

    return $adj;
}

# End of adjacencyMatrix

#------------------------------------------------------------------------------
# adjacencyRead
#
# Read the adjacency matrix from the <adjacency> section of the XIOS format.  If
# absent, the input will be empty and nothing will be stored ($nrows will be
# zero).  The text should have the <adjacency </adjacency> tags already removed.
#
# USAGE
#    $nrows = $topology->adjacencyRead( $text );
#------------------------------------------------------------------------------
sub adjacencyRead {
    my ( $topology, $text ) = @_;
    my $nrow = 0;

    my $adj;
    my @lines = split "\n", $text;
    foreach my $l (@lines) {
        next unless ( $l =~ /[-ijoxs]/i );    # skip blank rows

        my $ncol = 0;
        foreach my $c ( split " ", $l ) {
            if ( $ncol == 0 ) {

                # skip the stem number at the beginning of the line
                $ncol++;
                next;
            }

            $adj->[$nrow]->[ $ncol - 1 ] = $c;
            $ncol++;
        }
        $nrow++;
    }

    if ($nrow) {
        $topology->{adjacency} = $adj;
    }

    return $nrow;
}

# End of adjacencyRead

#-------------------------------------------------------------------------------
# adjacencyFormat
#
# Produce a formatted version of the adjacency matrix for use as a component of
# the XIOS file.
#
# usage
#	$str = $topology->adjacencyFormat();
#-------------------------------------------------------------------------------
sub adjacencyFormat {
    my ($self) = @_;
    my $str = "  <adjacency>\n";

    my @adj = @{ $self->adjacency };

   # the column width is the length of the largest index of the adjacency matrix
    my $wid  = length $#adj;
    my $fmts = "%$wid" . "s ";
    my $fmtd = "%$wid" . "d ";

    $str .= "    ";
    $str .= sprintf $fmts, " ";
    foreach my $col ( 0 .. $#{ $adj[0] } ) {
        $str .= sprintf $fmtd, $col;
    }
    $str .= "\n";

    foreach my $row ( 0 .. $#adj ) {
        $str .= "    ";
        $str .= sprintf $fmtd, $row;

        foreach my $col ( 0 .. $#{ $adj[$row] } ) {
            $str .= sprintf $fmts, $adj[$row]->[$col];
        }
        $str .= "\n";
    }
    $str .= "  </adjacency>\n";

    return $str;
}

# end of adjacencyFormat

#------------------------------------------------------------------------------
# overlapType
#
# Returns an enumerated overlap type from 0 - 69 based on the numer of two stems
# in the stemlist
#
# USAGE
#    my $info = $Topology->overlapType( $stem_a, $stem_b );
#------------------------------------------------------------------------------
sub overlapType {
    my ( $self, $n_a, $n_b ) = @_;

    my $stemlist = $self->stemList;
    my $sa       = $$stemlist[$n_a];
    my $sb       = $$stemlist[$n_b];
    my @coord    = (
        $sa->left1, $sa->left2, $sa->right1, $sa->right2,
        $sb->left1, $sb->left2, $sb->right1, $sb->right2
    );

    # get the left to right order of the indices
    my @order;
    my $order_str = "";
    foreach my $i ( sort { $coord[$a] <=> $coord[$b] } 0 .. 7 ) {
        $order_str .= $i;
        push @order, $i;
    }

    # this list of overlap types was generated by a script and hence should be
    # correct the three values are the enumerated stem number, the reduced
    # (grouped) stem number, and the symmetry operation that converts the stem
    # to the canonical reduced form

    my %type = (
        "01234567" => [ 0,  0,  0 ],
        "45670123" => [ 69, 0,  1 ],
        "01243567" => [ 1,  1,  0 ],
        "45607123" => [ 68, 1,  1 ],
        "01245367" => [ 2,  2,  0 ],
        "45601723" => [ 67, 2,  1 ],
        "45067123" => [ 64, 2,  2 ],
        "01423567" => [ 5,  2,  3 ],
        "01245637" => [ 3,  3,  0 ],
        "45601273" => [ 66, 3,  1 ],
        "40567123" => [ 54, 3,  2 ],
        "04123567" => [ 15, 3,  3 ],
        "01245673" => [ 4,  4,  0 ],
        "45601237" => [ 65, 4,  1 ],
        "04567123" => [ 34, 4,  2 ],
        "40123567" => [ 35, 4,  3 ],
        "01425367" => [ 6,  5,  0 ],
        "45061723" => [ 63, 5,  1 ],
        "01425637" => [ 7,  6,  0 ],
        "45061273" => [ 62, 6,  1 ],
        "40561723" => [ 53, 6,  2 ],
        "04125367" => [ 16, 6,  3 ],
        "01425673" => [ 8,  7,  0 ],
        "45061237" => [ 61, 7,  1 ],
        "04561723" => [ 33, 7,  2 ],
        "40125367" => [ 36, 7,  3 ],
        "01452367" => [ 9,  8,  0 ],
        "45016723" => [ 60, 8,  1 ],
        "01452637" => [ 10, 9,  0 ],
        "45016273" => [ 59, 9,  1 ],
        "40516723" => [ 50, 9,  2 ],
        "04152367" => [ 19, 9,  3 ],
        "01452673" => [ 11, 10, 0 ],
        "45016237" => [ 58, 10, 1 ],
        "04516723" => [ 30, 10, 2 ],
        "40152367" => [ 39, 10, 3 ],
        "01456237" => [ 12, 11, 0 ],
        "45012673" => [ 57, 11, 1 ],
        "40156723" => [ 44, 11, 2 ],
        "04512367" => [ 25, 11, 3 ],
        "01456273" => [ 13, 12, 0 ],
        "45012637" => [ 56, 12, 1 ],
        "04156723" => [ 24, 12, 2 ],
        "40512367" => [ 45, 12, 3 ],
        "01456723" => [ 14, 13, 0 ],
        "45012367" => [ 55, 13, 1 ],
        "04125637" => [ 17, 14, 0 ],
        "40561273" => [ 52, 14, 1 ],
        "04125673" => [ 18, 15, 0 ],
        "40561237" => [ 51, 15, 1 ],
        "04561273" => [ 32, 15, 2 ],
        "40125637" => [ 37, 15, 3 ],
        "04152637" => [ 20, 16, 0 ],
        "40516273" => [ 49, 16, 1 ],
        "04152673" => [ 21, 17, 0 ],
        "40516237" => [ 48, 17, 1 ],
        "04516273" => [ 29, 17, 2 ],
        "40152637" => [ 40, 17, 3 ],
        "04156237" => [ 22, 18, 0 ],
        "40512673" => [ 47, 18, 1 ],
        "40156273" => [ 43, 18, 2 ],
        "04512637" => [ 26, 18, 3 ],
        "04156273" => [ 23, 19, 0 ],
        "40512637" => [ 46, 19, 1 ],
        "04512673" => [ 27, 20, 0 ],
        "40156237" => [ 42, 20, 1 ],
        "04516237" => [ 28, 21, 0 ],
        "40152673" => [ 41, 21, 1 ],
        "04561237" => [ 31, 22, 0 ],
        "40125673" => [ 38, 22, 1 ],
    );

    return ( $order_str, @{ $type{$order_str} } );
}

# End of overlapType

#------------------------------------------------------------------------------
# overlapToCanonical
#
# Given the numbers of a pair of stems  and a symmetry opearation (0..3) return
# an array of eight canonically ordered coordinates.  Returns the order string
# needed to return the coordinates to their original AAaaBBbb order, and the
# coordinates in canonical order.
#
#
# USAGE
#    ($order,@coordinates) = $Topology->overlapToCanonical( $a, $b, $order, $symm );
#------------------------------------------------------------------------------
sub overlapToCanonical {
    my ( $self, $n_a, $n_b, $order, $symm ) = @_;
    my @coor;

    my $stemlist = $self->stemList;
    my $sa       = $$stemlist[$n_a];
    my $sb       = $$stemlist[$n_b];

    @coor = (
        $sa->left1, $sa->left2, $sa->right1, $sa->right2,
        $sb->left1, $sb->left2, $sb->right1, $sb->right2
    );

    # put coodinates into sorted order as indicated by $order
    my @o_index = split "", $order;
    my $i = 0;
    my @sorted;
    foreach my $c (@o_index) {
        $sorted[$i] = $coor[$c];
        $i++;
    }

    # symmetry types 0 and 1 only need to be sorted according to the $order to
    # be in canonical form.  symmetry types 1 and 2 need to have the order
    # string and coordinates inverted
    my @symm_coor;
    if ( $symm > 1 ) {

        # reverse the order string
        $order        = reverse $order;
        $symm_coor[0] = $sorted[7];
        $symm_coor[1] = $sorted[6];
        $symm_coor[2] = $sorted[5];
        $symm_coor[3] = $sorted[4];
        $symm_coor[4] = $sorted[3];
        $symm_coor[5] = $sorted[2];
        $symm_coor[6] = $sorted[1];
        $symm_coor[7] = $sorted[0];
        @sorted       = @symm_coor;

    }
    return $order, @sorted;
}

# End of overlapToCanonical

#------------------------------------------------------------------------------
# overlapFromCanonical
#
# Returns a canonicalized set of coordinates for two strings to the original
# AAaaBBbb order.  The coordinates are then stores as the coordinates of stem
# n_a and n_b
#
# USAGE
#    @coordinates = $Topology->overlapFromCanonical( $order, $n_a, $n_b, @canonical_coor );
#------------------------------------------------------------------------------
sub overlapFromCanonical {
    my ( $self, $order, $n_a, $n_b, @coor ) = @_;

    # put coodinates into sorted order as indicated by $order
    my @o_index = split "", $order;
    my $i = 0;
    my @orig;
    foreach my $c (@o_index) {
        $orig[$c] = $coor[$i];
        $i++;
    }

    my $stemlist = $self->stemList;
    my $sa       = $$stemlist[$n_a];
    my $sb       = $$stemlist[$n_b];

    $sa->left1( $orig[0] );
    $sa->left2( $orig[1] );
    $sa->right1( $orig[2] );
    $sa->right2( $orig[3] );
    $sb->left1( $orig[4] );
    $sb->left2( $orig[5] );
    $sb->right1( $orig[6] );
    $sb->right2( $orig[7] );

    return;
}

# End of overlapFromCanonical

#-----------------------------------------------------------------------------
# graphmlHeader
#
# Returns a string with a header for a graphml file suitable for display with
# yEd.
#
# USAGE
#   $string = graphmlHeader;
#-----------------------------------------------------------------------------
sub graphmlHeader {
    my ( $self, $param ) = @_;

    my $title = "XIOS";    # default title
    foreach my $option ( keys %$param ) {
        if ( $option =~ /title/ ) {
            $title = $$param{title};
        }
    }

    return qq{<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <graphml xmlns="http://graphml.graphdrawing.org/xmlns/graphml" 
             xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
             xmlns:y="http://www.yworks.com/xml/graphml" 
             xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/graphml 
                                 http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd">
    <key for="node" id="d0" yfiles.type="nodegraphics"/>
    <key attr.name="description" attr.type="string" for="node" id="d1"/>
    <key for="edge" id="d2" yfiles.type="edgegraphics"/>
    <key attr.name="description" attr.type="string" for="edge" id="d3"/>
    <key for="graphml" id="d4" yfiles.type="resources"/>
    <graph id="$title" edgedefault="directed">

    };

}

# end of graphmlHeader

#-----------------------------------------------------------------------------
# graphmlNode
#
# Return a string with the definition for one node (vertex) in the graphml
# file. A hash of parameters is used to size and position the node.
#
# available parameters:
#   xpos => <real>
#   ypos => <real>
#   fill_color => <RGB_hexadecimal> e.g., 'FFCC00'
#
# USAGE
#   $string = graphmlNode( $id );
#-----------------------------------------------------------------------------
sub graphmlNode {
    my ( $self, $id, $param ) = @_;

    my $xpos       = 128;
    my $ypos       = 128;
    my $height     = 30;
    my $width      = 30;
    my $fill_color = 'FFCC00';
    foreach my $option ( keys %$param ) {
        if ( $option =~ /xpos/i ) {
            $xpos = $$param{$option};
        }
        if ( $option =~ /ypos/i ) {
            $ypos = $$param{$option};
        }
        if ( $option =~ /fill/i ) {
            $fill_color = $$param{$option};
        }
        if ( $option =~ /height/i ) {
            $height = $$param{$option};
        }
        if ( $option =~ /width/i ) {
            $width = $$param{$option};
        }

    }

    return qq{ 
    <node id="$id">
        <data key="d0">
            <y:ShapeNode>
            <y:Geometry height="$height" width="$width" x="$xpos" y="$ypos"/>
            <y:Fill color="#$fill_color" transparent="false"/>
            <y:BorderStyle color="#000000" type="line" width="1.0"/>
            <y:NodeLabel alignment="center" autoSizePolicy="content" 
               fontFamily="Dialog" fontSize="12" fontStyle="plain" 
               hasBackgroundColor="false" hasLineColor="false" 
               height="18.701171875" modelName="internal" modelPosition="c" 
               textColor="#000000" visible="true" 
               width="19.0" x="5.5" y="5.6494140625">$id</y:NodeLabel>
            <y:Shape type="roundrectangle"/>
            </y:ShapeNode>
        </data>
        <data key="d1"/>
    </node>
    };
}

# end of graphmlNode

#-----------------------------------------------------------------------------
# graphmlEdge
#
# Return a string with the definition for one edge in the graphml
# file.
#
# USAGE
#   $string = graphmlEdge( $source, $target );
#-----------------------------------------------------------------------------
sub graphmlEdge {
    my ( $self, $source, $id, $target, $param ) = @_;

    my $type         = "";
    my $line_color   = '0000FF';
    my $line_type    = 'line';
    my $line_width   = 2.0;
    my $directed     = 'false';
    my $arrow_source = 'none';
    my $arrow_target = 'none';
    foreach my $option ( keys %$param ) {
        if ( $option =~ /type/i ) {
            $type = $$param{$option};
        }
        elsif ( $option =~ /color/i ) {
            $line_color = $$param{$option};
        }
        elsif ( $option =~ /directed/i ) {
            if ( $$param{$option} =~ /true/i || $$param{$option} eq '1' ) {
                $directed     = 'true';
                $arrow_source = 'none';
                $arrow_target = 'standard';
            }
        }
        elsif ( $option =~ /type/i ) {
            $line_type = $$param{$option};
        }
        elsif ( $option =~ /width/i ) {
            $line_width = $$param{$option};
        }
    }

    return qq{
    <edge id="$type$id" source="$source" target="$target">
        <data key="d2">
            <y:PolyLineEdge>
            <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0"/>
            <y:LineStyle color="#$line_color" type="$line_type" width="$line_width"/>
            <y:Arrows source="$arrow_source" target="$arrow_target"/>
            <y:EdgeLabel alignment="center" distance="2.0" 
                         fontFamily="Dialog" fontSize="12" fontStyle="plain" 
                         hasBackgroundColor="false" hasLineColor="false" 
                         height="4.0" modelName="six_pos" modelPosition="tail" 
                         preferredPlacement="anywhere" ratio="0.5" textColor="#000000" visible="true" 
                         width="4.0" x="61.0" y="23.0"/>
            <y:BendStyle smoothed="false"/>
            </y:PolyLineEdge>
        </data>
        <data key="d3"/>
    </edge>
    };

}

# end of graphmlEdge

#-----------------------------------------------------------------------------
# graphmlFooter
#
# Return a string that closes the graph and graphml elements
#
# USAGE
#   print $topology->graphmlFooter;
#-----------------------------------------------------------------------------
sub graphmlFooter {
    my ($self) = @_;

    return qq{
    </graph>\n</graphml>\n
    };

}

#-----------------------------------------------------------------------------
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
#-----------------------------------------------------------------------------
sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;

    my $unknown = $AUTOLOAD;    # The name of the unknown function
    $unknown =~ s/.*:://;       # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;    # Skip all-caps methods like DESTROY

    # convert capitals to _ followed by lowercase

    $unknown =~ s/([A-Z])/_$1/;
    $unknown =~ tr/A-Z/a-z/;

    # if this hash key exists, treat as accessor

    if ( exists $self->{$unknown} ) {       # do not create new fields
        $self->{$unknown} = shift if @_;    # Set new value if one was supplied
        return ( $self->{$unknown} );       # Return current value
    }

    print STDERR qq{Topology:AUTOLOAD   unknown function "$unknown" called\n};
}

#-----------------------------------------------------------------------------
# $Log: Topology.pm,v $
# Revision 1.56.4.29  2016/08/24 05:28:32  huang147
# revised the library for Perl.
#
# Revision 1.56.4.28  2015/10/30 21:23:48  huang147
# removed some comments.
#
# Revision 1.56.4.27  2015/10/30 20:13:45  huang147
# added subroutine stemsToPairArray.
#
# Revision 1.56.4.26  2015/09/07 21:32:44  huang147
# Fixed some bug in bpseqToTopology by commenting the basepairPrune subroutine.
#
# Revision 1.56.4.25  2015/09/07 16:17:36  huang147
# Updated XIOSRead so it reads in the sequence of RNA now.
#
# Revision 1.56.4.24  2014/11/06 17:38:21  huang147
# Fixed a bug in ViennaToTopology.
#
# Revision 1.56.4.23  2014/11/05 16:53:10  huang147
# Added function ViennaToTopology.
#
# Revision 1.56.4.22  2014/10/21 18:55:42  huang147
# Changed the file header for the subroutine mfoldMultiple to take. Now it also takes .ct files created by the program Fold at RNAstructure.
#
# Revision 1.56.4.21  2013/10/24 02:46:42  huang147
# Fixed the comments in the function RandomAddPairsWeighted.
#
# Revision 1.56.4.20  2013/08/19 02:23:43  huang147
# Added subroutine addTopology, which adds a topology to an existing topology by adding all of its stems to the stem list of the existing topology.
#
# Revision 1.56.4.19  2013/08/14 16:06:40  huang147
# Fixed some bugs.
#
# Revision 1.56.4.18  2013/08/14 15:30:11  huang147
# Added $h_num_percent into the subroutine HNumToTopology.
#
# Revision 1.56.4.17  2013/08/08 19:30:50  huang147
# Added calculation of 65% quantile h-num value to do the stem pruning.
# Changed the CVS: ----------------------------------------------------------------------
#
# Revision 1.56.4.15  2013/08/08 15:36:44  huang147
# Added the option of $max_unpaired in HNumToTopology.
#
# Revision 1.56.4.14  2013/08/08 04:11:51  huang147
# Added a subroutine HNumToTopology, which converts .h-num files to .xios files by reading in a list of stems.
#
# Revision 1.56.4.13  2013/07/23 03:11:52  huang147
# Added basepairPrune to remove redundant bases in a given sequence.
#
# Revision 1.56.4.12  2013/07/19 22:09:44  huang147
# Changed the printout. No functional change.
#
# Revision 1.56.4.11  2013/07/19 22:06:56  huang147
# Added sub basepairToDotsAndBrackets.
#
# Revision 1.56.4.10  2013/07/19 20:46:01  huang147
# In sub bpseqToTopology, removed leading N/n's from the sequence when reading in .bpseq.
#
# Revision 1.56.4.9  2013/07/19 16:29:23  huang147
# Added subroutine bpseqToToplogy, which takes a bpseq file and writes it into an RNA topology.
#
# Revision 1.56.4.8  2013/05/29 18:29:31  huang147
# Fixed the bug in removing duplicate stems.
#
# Revision 1.59.2.16  2013/03/31 17:25:05  gribskov
# Added a weighted random strategy for creating random graphs.  This should make graphs
# that are more like biological RNAs.
#
# Revision 1.59.2.15  2013/03/21 11:53:05  gribskov
# Updated examples.  No functional changes.
#
# Revision 1.59.2.14  2013/01/11 21:27:46  gribskov
# No functional changes
#
# Revision 1.59.2.13  2013/01/11 17:12:06  gribskov
# Added function randomAddPair to randomly add one new stem to a pair Array.
# this can be used to augment a core graph.
#
# Revision 1.59.2.12  2013/01/11 16:10:46  gribskov
# Added pairArray function to read in graphs as pair arrays where an integer
# array indicates the relative positions of each pair of stems in the graph.
#
# Revision 1.59.2.11  2013/01/11 15:48:40  gribskov
# Simplified setTopology.
#
# Revision 1.59.2.10  2013/01/11 15:17:02  gribskov
# Converted random to be an object function.  Some minor format improvements.
#
# Revision 1.59.2.9  2013/01/10 12:54:40  gribskov
# Added some documentation and reformatted to remove long lines.  No functional changes.
#
# Revision 1.59.2.8  2012/12/24 11:21:31  gribskov
# Small changes to make sequence block work better.
#
# Revision 1.59.2.7  2012/12/21 18:37:31  gribskov
# Added <information> section to the XML output for metadata.  Currently holds
# sequence and comments.
#
# Revision 1.59.2.6  2012/12/19 14:53:12  gribskov
# Saved current version.  Updates to GraphML i think.
#
# Revision 1.59.2.5  2012/10/07 22:08:45  gribskov
# Added a length field.  not sure it is needed.  should there be a sequence block in xios
# format?  seems likely.
#
# Revision 1.59.2.4  2012/10/07 12:13:03  gribskov
# Added length to class.
#
# Revision 1.59.2.3  2012/10/05 21:20:35  gribskov
# more updates to graphml code.  proably has bugs.
#
# Revision 1.59.2.2  2012/10/04 21:24:57  gribskov
# Added documentation on overall usage.
# Minor updates to some function documentation.
#
# Corrected bug in splitStems.  Only stems with greater than max_unpaired unpaired
# bases in left and right half stems were being split.
# Added parameters for xpos, ypos, and fill color to graphmlNode.
#
# Revision 1.59.2.1  2012/10/03 16:01:52  gribskov
# bug in splitStems.  Stems were only split if both half stems had >= max_unpaired
# bases.
#
# Revision 1.59  2012/09/21 15:01:38  gribskov
# Fixed incorrect syntax in OverlapType.  Works in tests.
#
# Revision 1.58  2012/09/21 12:54:13  gribskov
# Added functions to determine stem overlap, and to convert to and from canonical
# form.  These functions have not been tested but should not have any side effects on
# any existing code.
#
# Revision 1.57  2012/09/14 21:36:16  gribskov
# Added overlapType to determine overlap type between stems.  This function is not
# finished.  Need to add automatic determination of whether the overlap should be
# determined by an array of coordinates or stem indices.  Also need to add conversion
# to and from reduced overlap type 0-22.
#
# Revision 1.56  2012/06/27 20:38:52  gribskov
# Removed commented out lines.  tested, seems ok.
#
# Revision 1.55  2012/06/26 21:06:37  huang147
# Fixed the bug of identifying "i" and "j" edges.
#
# Revision 1.54  2012/06/22 18:57:30  huang147
# Deleted the debugging code so the xios file can print out correctly.
#
# Revision 1.53  2012/06/22 15:12:16  gribskov
# Fixed calculation of adjacency matrix so that lower triangle is correct.
# Fixed apparent error in identifying o edges in adjacency matrix.
#
# Revision 1.52  2012/04/05 00:33:56  gribskov
# Added a function to return the coordinates of the stems as an array of arrays
# rather than an array of Stem objects: stemListAsArray.  This is useful when
# updating older code that does not use Stem objects, probably should bite
# the bullet and make it use objects.  Too tired.
#
# Revision 1.51  2012/04/04 23:06:40  gribskov
# Removed duplicate functions stemlistNoViennaFormat, writeXIOSNoVienna,
# XIOSwriteNoVienna.  stemlistFormat, and XIOSWrite work fine for topologies with
# no Vienna format information. the  writeXIOS function is only needed for back
# compatibility and should not be propagated to new code.
#
# Revision 1.50  2012/04/04 22:56:31  gribskov
# Added mfoldMultiple to read CT files with multiple structures into a single topology.
#
# Revision 1.49  2012/04/04 20:10:34  huang147
# changed writeXIOSNoVienna
#
# Revision 1.48  2012/04/04 19:57:03  huang147
# Added 3 new functions: stemlistNoViennaFormat, writeXIOSNoVienna, XIOSwriteNoVienna
# CiVS: ----------------------------------------------------------------------
#
# Revision 1.47  2012/04/02 21:47:03  gribskov
# Fixed adjacencyRead to ignore the label line, and the stem label at the
# beginning of each row.
#
# Revision 1.46  2012/04/02 20:27:20  gribskov
# renamed writeXIOS -> XIOSwrite
# renamed formatAdjacency ->adjacencyFormat
# renamed formatEdgelist -> edgelistFormat
# ranmaed formatStemlist -> stemlistFormat
# retained writeXIOS for backcompatability (deprecated) but did not do anything
# for the others because they should only be in use within the package.
#
# Revision 1.45  2012/04/02 20:19:30  gribskov
# Added XIOS read to read a topology in the new XIOS format.
#
# Revision 1.44  2012/03/31 13:34:59  gribskov
# Finished implementing XIOS format in writeXIOS.  Works with test data.
#
# Revision 1.43  2012/03/29 16:52:04  gribskov
# Added adjacency matrix calculation.
#
# Revision 1.42  2012/03/29 13:13:44  gribskov
# Implemented and tested code for printing stem_list section of XIOS format.
#
# Revision 1.41  2012/03/28 15:30:45  gribskov
# Added use for standard package variables.
# Added VERSION and REVISION
# Added writeXIOS function.  this function relies on several additional functions that
# have not yet been written so it will not wrk yet.
#
# Revision 1.40  2012/03/28 12:14:12  gribskov
# Added support for pair format.
#
# Revision 1.39  2011/09/04 05:42:48  gupta31
# added subroutine dumpToFile to print stem object to a file
#
# Revision 1.38  2010/01/14 13:45:47  likejie
# fixed the undef $newl1 problem.
# It was a  regex mistake on line 1103, the paren put to a wrong place
# Also removed the my's from 1092 to 1094
#
# Revision 1.37  2009/11/02 16:43:42  gupta31
# Edited energyplot function to add stems correctly. Tested out dump function on it.
#
# Revision 1.36  2009/10/19 17:38:28  gupta31
# Added funtion 'energyplot' to read stems from a plot file. See test_plotfile.pl on an example usage.
#
# Revision 1.35  2009/08/19 14:04:49  gribskov
# Fixed bug in findstems that caused pseudoknots to be incorrectly identified when there
# were multiple closing right brackets in series, e.g., ))]]}}
#
# Revision 1.34  2009/08/12 14:08:57  gribskov
# Modified mfold2 to stop after reading a single structure.  I don't know why
# my files have two complete copies of the structure, probably an error in
# how i ran unafold.
#
# Revision 1.33  2009/08/12 13:09:03  gribskov
# Modified mfold2 to return energy as part of object.  simplified code to skip comments.
#
# Revision 1.32  2009/07/28 18:54:54  gribskov
# Added Exporter to allow export of basepairToStem and splitStems functions.
#
# Revision 1.31  2009/04/23 14:25:50  rahmanr
#
# Added some error check code
#
# Revision 1.30  2009/04/23 11:44:05  gribskov
# Fixed error in finding pseudoknots from ct files.  original ct files did
# not have any since MFOLD/unafold cannot find them.  Changes in function
# basepairtostem.
#
# Revision 1.29  2009/01/06 19:06:18  gribskov
# Fixed small bug in stem splitting function.
#
# Revision 1.28  2009/01/02 16:51:45  gribskov
# Added Sequence ID and documentation to class structure.
# Added code to set sequence ID and doc in rnaml and mfold2.
# Removed debugging printout.
#
# Revision 1.27  2009/01/02 16:40:30  gribskov
# Added function to split stems with unpaired bases.  debugged RNAML reading
# Agrees with my hand translation for sacce.rnasepdb.rnaml
#
# Revision 1.26  2009/01/02 15:04:14  gribskov
# Added rnaml function to read RNAML file.
#
# Revision 1.25  2009/01/01 14:39:43  gribskov
# Created new function, basepairToStem, from mfold2 so that it can also be
# used for rnaml.  works for ct files.
#
# Revision 1.24  2009/01/01 14:09:35  gribskov
# Added rnaml to setTopology.  No function to read rnaml yet.
#
# Revision 1.23  2009/01/01 13:55:57  gribskov
# Modified setTopology for new mfold2.
#
# Revision 1.22  2008/12/31 17:16:08  gribskov
# Added new function, mfold2, which does a better job getting stems from
# CT file.  Does not rely on conversion to Vienna string, so should do a
# better job on pseudoknots.  Did not add maximum bulge size for breaking
# up stems yet.
#
# Revision 1.21  2008/12/27 13:42:09  gribskov
# Added sequence attribute to topology.  When reding ct (mfold) format, the
# sequence of the last strcuture read is stored.
#
# Revision 1.20  2008/12/19 12:39:34  gribskov
# Added mfold function to read Zuker/mfold/ct structures.
#
# Revision 1.19  2008/12/18 18:52:32  gribskov
# Updated documentation.  Added stub function for reading mfold format but
# did not change any functional code.
#
# Revision 1.18  2008/05/18 12:13:09  gribskov
# Added graphml writing functions.  Still need function to write entire
# topology in GraphML.
# Changed vienna function to data.  Calling by name "vienna" is now the
# default (autoloaded)
#
# Revision 1.17  2008/05/15 21:46:42  gribskov
# changed a few function calls to OO form, fixed typo in call to
# expandBrackets.
#
# Revision 1.16  2008/05/15 21:29:59  gribskov
# Complete replaced vienna stem finding with much more straightforward
# method (findStems).  Stems can now be split after a specified number of
# unpaired bases.
#
# Revision 1.15  2008/04/22 03:04:50  gribskov
# Added dump function
# Added vienna function to get structure string.  Should rename internal
# variable from data to vienna.
# Many small modifications to increase uniformity and enforce structure of
# Topology as an array of Stem objects.
#
# Revision 1.14  2008/04/01 17:05:02  gribskov
# Added indentation to findstem.  Cleaned up comments and added some
# documentation
#
# Revision 1.13  2008/04/01 16:05:58  gribskov
# Converted storage of stems in findStem to use stem objects.
#
# Revision 1.12  2008/03/28 02:59:42  gribskov
# updated documentation and formatting.  Should be no functional changes
#
# Revision 1.11  2008/03/21 23:11:39  gupta31
# Renamed original printStemList to printStemListOld, and added a new
# printStemList that returns the stems in an array instead of printing them
#
# Revision 1.10  2007/07/31 22:29:21  rahmanr
# Documentation chages and some otehre major changes. Now, vienna strings with
# pseudknots can be parsed. In addition, addData($data) has been uspdateD so
# that it calls find stem when data is added to topology object. The
# constructor is modified so that if any atrribute is passed ot it then it can
# put it in the Topology object by itself
#
# Revision 1.9  2007/07/31 22:12:00  rahmanr
# changes to documentation
# - revision comments
# Revision 1.6  2007/07/27 17:14:13  gupta31
# Added the mapping feature in randomPartialInterspersed. Now can also tell
# the half stems in new random topology corresponding to the specified topology.
# Made some revisions in subroutine explainations.
#
# Revision 1.5  2007/07/23 19:31:46  gupta31
# Added subroutine to generate Vienna format of RNA structure from topology
# format. Subroutine name : topologyToVienna
#
# Revision 1.4  2007/07/23 17:20:48  gupta31
# Added subroutines for random generation of Topologies.
#
# Revision 1.3  2007/07/18 17:32:57  rahmanr
#-----------------------------------------------------------------------------

