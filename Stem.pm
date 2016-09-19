package Stem;
#-----------------------------------------------------------------------------
# $Id: Stem.pm,v 1.12 2009/01/02 16:41:48 gribskov Exp $
#
# This class stores all information that is necessary for stems, namely
# left1 (left start), left2 (left end), right1 (right start) and right2 
# (right end)# 
# 
# 
#-------------------------------------------------------
use strict;
#use vars qw(@ISA);
#@ISA = qw { AutoLoader };   

#---------------------------------------------------
# new
#
# this is the constructor for this class. This calls the initialize subroutine
# to intialize different attributes of the class. 
# the class structure is a hash:
# stem = { 'left1'  => left start,
#          'left2'  => left end,
#          'right1' => right start,
#          'right2' => right end
#          'vienna_left' => vienna string for left side
#          'vienna_right => vienna string for right side
#         }
#
# USAGE
#      $stem_obj = Stem->new();
#      $stem_obj = Stem->new( $coord_hash );
#      $stem_obj = Stem->new( $coord_array );
#----------------------------------------------
sub new {
    my ($class, $coord) = @_;
    
    my $stem = {};	
    bless $stem, $class;

	$stem->initialize();	
    if ( defined( $coord ) ) {
	    $stem->setStem($coord);
    }
    
    return $stem;   
}

#---------------------------------------------------
# initialize
#
# This is called by the constructor "new" to initialize the different attributes 
# of the class.It initializes the hash and creates the necessary attributes for the Stem Object
#
# USAGE
#      $stem_obj->initialize(); 
# NOTE: This is used locally within the class 
#----------------------------------------------

sub initialize {
    my ($self)=@_;
    $self->{'left1'} = 0; 
    $self->{'left2'} = 0; 
    $self->{'right1'} = 0; 
    $self->{'right2'} = 0; 
    $self->{'vienna_left'}  = "";
    $self->{'vienna_right'} = "";
    
    return $self;    
}
#-----------------------------------------------------------------------------
# setStem
#
# sets the attributes of a stem based on a hash or an array.  If the argument 
# is a hash, the keys are:
#   left1, left2, right1, right2, vienna_left, vienna_right
# if the argument is a hash, the four elements of the array should be the 
# coordinates in the order given above.
#
# USAGE
#      $stem_obj->setStem($hash-ref)
#      $stem_obj->setStem($array-ref)
#----------------------------------------------

sub setStem {
    my ($self, $argument) = @_;

    if ( ref($argument)=~/hash/i ) {
        $self->left1(  $argument->{'left1'} ); 
        $self->left2(  $argument->{'left2'} ); 
        $self->right1( $argument->{'right1'} ); 
        $self->right2( $argument->{'right2'} ); 
        $self->viennaLeft( $argument->{'vienna_left'} );
        $self->viennaRight( $argument->{'vienna_right'} );
    } else {
        $self->left1(  $argument->[0] ); 
        $self->left2(  $argument->[1] ); 
        $self->right1( $argument->[2] ); 
        $self->right2( $argument->[3] ); 
        $self->viennaLeft( $argument->[4] );
        $self->viennaRight( $argument->[5] );
    }

    return ;
}

#-----------------------------------------------------------------------------
# getStemAsArray
#
# gets the information in the stem as an array;
#
# USAGE
#   @stem = $stem->getStemAsArray;
#-----------------------------------------------------------------------------
sub getStemAsArray{
    my ( $self ) = @_;

    return ( $self->left1,      $self->left2, 
             $self->right1,     $self->right2, 
             $self->viennaLeft, $self->viennaRight );
}

# end of getStemAsArray

#-----------------------------------------------------------------------------
# dump
#
# print out the stem structure
#
# USAGE
#   $stem->dump;
#-----------------------------------------------------------------------------
sub dump{
    my ( $self ) = @_;

    foreach my $key ( sort keys %{$self} ) {
        print "$key  =>  ", $self->{$key}, "\n";
    }
    print "\n";

    return;
}

#---------------------------------------------------
# printStemOld
#
# print the stem object
# original version written by reazur
#
# USAGE
#      $Stem_obj->printStem;
#----------------------------------------------


#print the stem object
sub printStemOld {
    my ($self)=@_;
    print "left1: $self->{'left1'}; left2: $self->{'left2'}; right1: $self->{'right1'}; right2: $self->{'right2'}\n";
}


#---------------------------------------------------
# printStem
#
# returns the stem object as <left start>-<left end>:<right start>-<right end>
#
# USAGE
#      $stem = $Stem_obj->printStem;
#----------------------------------------------


#print the stem object
sub printStem {
    my ($self)=@_;
    my $stem = sprintf "%s-%s  %s       %s-%s  %s\n",
        $self->left1,$self->left2,$self->vienna_left,$self->right1,$self->right2,$self->vienna_right;
    print "$stem\n";
    return $stem;
}

#############################################################################
# AUTOLOAD
# When methods refer unknown functions to their parent, they all end up 
# here.  A last attempt is made to treat unknown functions as accessor 
# functions
# 
# $class->unknown_function( $value )
#   sets $subclass->{unknown_function} = $value and returns $value
#
#  $class->unknown_function
#   returns $subclass->{unknown_function}
# 
# capitals in function names get translated to _lowercase, e.g.,
# viennaString changes to vienna_string
#
# This function will not create new hash keys.  This keeps typos from
# crazy new attributes for objects.  Anything that is not a hash key or 
# a function in this package fails.
#############################################################################
sub AUTOLOAD {
    my $self = shift;
    our $AUTOLOAD;
    
    my $unknown = $AUTOLOAD;                # The name of the unknown function
    $unknown =~ s/.*:://;                   # Object::name becomes name
    return unless $unknown =~ m/[^A-Z]/;    # Skip all-caps methods like DESTROY
    
    # convert capitals to _ followed by lowercase

    $unknown =~ s/([A-Z])/_$1/;
    $unknown =~ tr/A-Z/a-z/;

    # if this hash key exists, treat as accessor


    if ( exists $self->{ $unknown } ) {         # do not create new fields
        $self->{ $unknown } = shift if @_;      # Set new value if one was supplied
        return( $self->{ $unknown } );          # Return current value
    }
    
    # the buck stops here, there are no more places to check

    print STDERR "Stem.pm::UNKNOWN Function $unknown called\n";
   
}

#-----------------------------------------------------------------------------
# $Log: Stem.pm,v $
# Revision 1.12  2009/01/02 16:41:48  gribskov
# Added getStemAsArray to return whole stem information in an array.  Should
# probably have getStemAsHash function #TODO
#
# Revision 1.11  2008/12/31 17:07:50  gribskov
# Added sort to dump function.
#
# Revision 1.10  2008/05/15 21:37:18  gribskov
# Updated annotation.  No functional change.
#
# Revision 1.9  2008/05/15 21:28:44  gribskov
# Made proinStem really print.  Added vienna structure to printStem.
#
# Revision 1.8  2008/04/22 03:01:18  gribskov
# Added dump function.
# Added printout of function name when not found in AUTOLOAD.
#
# Revision 1.7  2008/04/01 17:18:55  gribskov
# Added storage for left and right Vienna strings.  Autoload now converts
# capitals to corresponding _ lowercase when auto creating accessors.
#
# Revision 1.6  2008/04/01 17:03:58  gribskov
# Added indentation to findstem.  Added documentation.  No functional
# changes.
#
# Revision 1.5  2008/04/01 13:44:18  gribskov
# Added ability to enter stem coordinates from array in }setStem.  Should
# be completely back-compatible.  Did not test.
#
# Revision 1.4  2008/04/01 13:00:24  gribskov
# Removed spurious carriage return characters (^M). Moved revision log to
# end of file.  No functional changes.
#
# Revision 1.3  2008/03/21 23:21:47  gupta31
# Renamed original printStemList to printStemListOld, and  to printStemOld, 
# and added new printStem which returns a stem as a scalar instead of printing 
# it
#
# Revision 1.2  2007/07/06 01:15:11  rahmanr
# update documentation   VS: Modified Files:
#-----------------------------------------------------------------------------

1;
