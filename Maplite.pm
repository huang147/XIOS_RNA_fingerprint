##############################################################################
# Maplite package
#
# data structure to hold 
# 1. Current dtg map
# 2. Current dfs row (last dfsrow) 
##############################################################################
package Maplite;

#-----------------------------------------------------------------------------
# new
#
# class structure is an array of dfs code rows, and an array of dtg mappings.
# each mapping corresponds to the single dfs code.
#
# USAGE
#   my $map = new Map;
#-----------------------------------------------------------------------------
sub new{
    my ( $class, $param ) = @_;

    my $self = {};
    bless $self, $class;
    $self->_init; 

    if ((defined $param->{'dfsrow'}) && defined ($param->{'dtg'})) {
	$self->{dfsrow} = $param->{'dfsrow'};
	$self->{dtg} = $param->{'dtg'};
    }  
   
    return $self;
}


#-----------------------------------------------------------------------------
# _init:
#
# create the attributes for the DfsGraph object.  Try to keep all initializ-
# ation here.
#
# USAGE:
#   $dfs->_init;
#-----------------------------------------------------------------------------
sub _init{
    my ( $self ) = @_;
    $self->{dfsrow} = [];
    $self->{dtg} = [];
    return;
}


#-----------------------------------------------------------------------------
# setDfs
#
# set the dfs code array.
#
#
# USAGE
#   $map->setDfs($dfs);
# where $dfs represents the array of dfs code
#-----------------------------------------------------------------------------
sub setDfsrow{
    my ($self, $dfs )= @_;
    
    if ( defined $dfs ) {
	$self->{dfsrow} = $dfs;
    }
    #return $self->{dfs};
    return;
}

#-----------------------------------------------------------------------------
# getDfs
#
# return the dfs code array.
#
#
# USAGE
#   $dfs = $map->getDfs();
#   
#-----------------------------------------------------------------------------
sub getDfsrow{
    my $self= shift;
    return $self->{dfsrow};
}


#-----------------------------------------------------------------------------
# setDtg
#
# set dtg map
#
# USAGE
#   $map->setDtg( $dtg );
#----------------------------------------------------------------------------
sub setDtg {
    my ($self, $dtg) = @_;
    
    if ( defined $dtg ) {
	$self->{dtg} = $dtg;
    }
    return;
}



#-----------------------------------------------------------------------------
# getDtg
#
# get dtg map
#
# USAGE
#   $dtg = $map->getDtg();
#-----------------------------------------------------------------------------
sub getDtg {
    my ($self) = @_;
    return $self->{dtg};
}



#-----------------------------------------------------------------------------
# dump
#
# print out the content of the map structure
#
# USAGE
#   $map->dump;
#-----------------------------------------------------------------------------
sub dump{
    my ( $self, $text ) = @_;

    if ( defined $text ) {
        print "\n$text\n";
    }
    print "DFSrow\n";
    my @dfs = @{$self->getDfsrow()};
 #   my ( $v1, $v2, $edge ) = @{$row};

    my $n = 0;
    foreach my $row ( @dfs ) {
        my ( $v1, $v2, $edge ) = @{$row};
        printf "   %4d %4d %4d\n", $v1, $v2, $edge;
        $n++;
    }

    print "\nDTG maps\n";
    my @dtg = @{$self->getDtg()};
  
    $n = 0;
    foreach my $row ( @dtg ) {
        printf "   %-3d ", $n;
        print $row,"\n";
        $n++;
    }

    return;
}

##############################################################################


##############################################################################
#-----------------------------------------------------------------------------
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
# $Log: Maplite.pm,v $
# Revision 1.3  2011/08/05 12:52:34  gribskov
# Committed existing version.  I is the same as the older one except for
# formatting
#
# Revision 1.2  2010/03/02 05:54:43  rahmanr
#
# This is a lite version of the map package to be used with DFsgenerator package. The maplite object only stores
# the dtg map and the current dfsrow
#
# Revision 1  2010/01/20 09:39:03  rahmanr
# 
#-----------------------------------------------------------------------------

# end of Map package

1;
