#!/apps/group/bioinformatics/apps/perl-5.16.1/bin/perl -w 
#------------------------------------------------------------------------------
# $Id: fingerprint_distance.pl,v 1.1.2.69 2016/08/24 05:23:31 huang147 Exp $
#
# This program is to calculate the distance of every query structure towards
#  all other structures, based on the features selected for each structure.
#
# usage
#   fingerprint_distance.pl *.xfp
#------------------------------------------------------------------------------
use strict;
use lib './';
use Data::Dumper;
use Storable;
use Getopt::Std;
use List::Util qw(sum max min shuffle);
use Fingerprint;
use Cwd;
use GD;
use GD::Graph;

my $revision = '$Revision: 1.1.2.69 $';
my $USAGE .= qq{$revision [-h] [-i|em] [-s|abcdehmnrtvx] <.fpout>

    -d      input directory of fingerprint files (*.xpt format)
   	-h  	display this help message
	-i|m	original fingerprint
	-i|e	extended fingerprint
   	-s|b	BinarySimilarity
	-s|c	CosineSimilarity
	-s|d	DiceSimilarity
	-s|t	TanimotoSimilarity
	-s|e	EdgeNumberSimilarity
	-s|v	VertexNumberSimilarity
    -w      output directory for clustering results

};


my $cwd       = getcwd;
my @xfp       = ();
my @xios      = ();
my @mfp       = ();
my %xfp2motif = ();
my %xfp2extended = ();
my %xfp2dp = (); # direct parent: only one vertex less
my %xfp2vertex = ();
my %xfp2edge  = ();
my %mfp2motif = ();
my %xfpfam    = ();
my %mfpfam    = ();
my %extended_sp = ();
my %sp = ();

my %option;
getopts( "d:hi:s:w:", \%option );

if ( defined $option{h} ) {
    print STDERR $USAGE;
    exit 1;
}

my $output_dir = $cwd."/outputs/curated_graphs/clustering/";
if ( defined $option{w} ) {
    $output_dir = $option{w}."/";
}
unless (-e $output_dir) {
    system("mkdir -p $output_dir");
}
my $report = $output_dir."AUC.txt";


my $extendedDB = $cwd."/Motif_fingerprint/2_to_7_parent_child.storable";
my $pc = retrieve( $extendedDB );
my $parent = $$pc{parent};
my $child  = $$pc{child};
my %dp;

my $input_dir = $cwd."inputs/curated_graphs/fingerprint/";
if ( defined $option{d} ) {
    $input_dir = $option{d}."/";
}

my $DEFAULT_INDEX = "me";
my $current_index = $DEFAULT_INDEX;
if ( defined $option{i} ) {
    $current_index = $option{i};	
}

my $DEFAULT_SIMILARITY = "bcdht";
my $current_similarity = $DEFAULT_SIMILARITY;
if ( defined $option{s} ) {
    $current_similarity = $option{s};
}

opendir( DH, $input_dir ) || die "can not open directory $input_dir!\n";
while ( my $file = readdir DH ) {
    if ( $file =~ /\.fpout$/ || $file =~ /\.xpt$/ ) {
        push @xfp, $file;
    }
    elsif ( $file =~ /manual\.structure_out\.result$/ ) {
        push @mfp, $file;
    }
}
close DH;


my $graph_stats_string = qq{fam\tname\tvnum\tenum\tmonum\tmenum\n};
print STDERR "Graph Stats Calculation...\n";
my $nxfp = scalar @xfp;
my $ixfp = 1;
foreach my $xfp (@xfp) {
    print STDERR "$xfp read in ($ixfp / $nxfp)\n";
    $ixfp++;
    my ($name) = $xfp =~ /(.*)\.xios/;
    my ($fam) = split /[\.\_]+/, $xfp;

    # set up fingerprint
    my %motif_prob_index = ();
    my %motif_index = ();
    my %extended_index = ();
    my $fingerprint = Fingerprint->new;
    my $n_motif     = $fingerprint->readXML("$input_dir/$xfp");
    my $motiflist   = $fingerprint->motifList;
    my $monum = scalar @{$motiflist};
    my $query_vertex = $fingerprint->queryVertex;
    my $query_edge = $fingerprint->queryEdge;
    my $iteration = $fingerprint->iteration;


    foreach my $motif ( @{$motiflist} ) {

        my $id    = $$motif{id};
        my $count = $$motif{count};

        my ( $msize ) = $id =~ /(.*)\_/;
        $motif_index{$id} = 1;
        $extended_index{$id} = 1;
        $dp{$id}{$id} = 1;

        foreach my $p ( @{$$parent{$id}} ) {
            my ( $psize ) = $p =~ /(.*)\_/;
            $extended_index{$p} = 1;
            $dp{$id}{$id} = 1;
        }
    }

    my $menum = scalar keys %extended_index;
    $graph_stats_string .= qq{$fam\t$name\t$query_vertex\t$query_edge\t$monum\t$menum\n};

    foreach my $p ( keys %extended_index ) {
        $extended_sp{$p}++;
    }


    $xfp2vertex{$name} = $query_vertex;
    $xfp2edge{$name} = $query_edge;
    $xfp2motif{$name} = \%motif_index;
    $xfp2extended{$name} = \%extended_index;
    $xfpfam{$fam}++;
}
print STDERR "Graph Stats Calculation Done!\n";


foreach my $mfp (@mfp) {
    my $motif_index = readmanout($mfp);
    my ($fam, $species) = split /\./, $mfp;
    my $name = $fam.'.'.$species;
    $mfp2motif{$name} = $motif_index;
    $mfpfam{$fam} = 1;
}

my @index = split "", $current_index;
while (@index) {

    my $index = shift @index;

    if ( $index eq "m" ) {

        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\n\nSimple Fingerprint\n";
	print STDERR "\n\nSimple Fingerprint\n";
        close OUT;
        Similarity2ROC( \%xfp2motif, \%xfp2vertex, \%xfp2edge, $current_similarity );

    } elsif ( $index eq "e" ) {

        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\n\nExtended Fingerprint\n";
	print STDERR "\n\nExtended Fingerprint\n";
        close OUT;
        Similarity2ROC( \%xfp2extended, \%xfp2vertex, \%xfp2edge,  $current_similarity );

    } 
}





sub Similarity2ROC {

    my ( $xfp2motif, $xfp2vertex, $xfp2edge, $current_similarity ) = @_;

    my @simi_func = split "", $current_similarity;

    while (@simi_func) {

        my $simi_func = shift @simi_func;

        my ( $similarity, $distance );

        my $similarity_string;

        if ( $simi_func eq "b" ) {

            $similarity_string = "BinarySimilarity";
            ( $similarity, $distance ) = BinarySimilarity($xfp2motif);

        }
        elsif ( $simi_func eq "c" ) {

            $similarity_string = "CosineSimilarity";
            ( $similarity, $distance ) = CosineSimilarity($xfp2motif);

        }

        elsif ( $simi_func eq "d" ) {

            $similarity_string = "DiceSimilarity";
            ( $similarity ) = DiceSimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "h" ) {

            $similarity_string = "HammingSimilarity";
            ( $similarity ) = HammingSimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "t" ) {

            $similarity_string = "TanimotoSimilarity";
            ( $similarity, $distance ) = TanimotoSimilarity($xfp2motif);

        }


        elsif ( $simi_func eq "v" ) {

            $similarity_string = "VertexNumberSimilarity";
            ( $similarity, $distance ) = VertexNumberSimilarity($xfp2motif);

        }


        elsif ( $simi_func eq "e" ) {

            $similarity_string = "EdgeNumberSimilarity";
            ( $similarity, $distance ) = EdgeNumberSimilarity($xfp2motif);

        }


	#my $mega = $output_dir.$similarity_string.".meg";
	#if ( defined $distance ) {
	#    DistanceToMega( $distance,  $mega );
	#} else {
	#    SimilarityToMega( $similarity, $mega );
	#}
	#my $mtx = $output_dir.$similarity_string.".mtx";
	#SimilarityToMatrix( $similarity, $mtx );

        open( OUT, ">>$report" ) || die "can not open $report!\n";
        print OUT "Similarity function: $similarity_string\n";
        print STDERR "Current similarity function is $similarity_string!\n";
        close OUT;

        my %similarity = %{$similarity};

        my $plot_total = similarity4ROCAllFam( \%similarity );

        my $roc_total = $output_dir."$similarity_string" . "_total_ROC.png";
        my $fam       = "total";
        my ( $im_total, $coordinates_total )  = ROC( $plot_total, $report, $fam );
        ROCplot( $roc_total, $im_total );
	
	#my $coordinates_total_file = $output_dir.$similarity_string."_roc_total_coordinates.txt";
	#open ( my $out, ">$coordinates_total_file" ) || die "can not open $coordinates_total_file!\n";
	#print $out $coordinates_total;
	#close $out;

        for my $fam ( keys %xfpfam ) {
            my $plot = similarity4ROCbyFam( $fam, \%similarity );

            my $roc = $output_dir."$similarity_string" . "_" . "$fam" . "_ROC.png";
            my ( $im, $coordinates) = ROC( $plot, $report, $fam );
            ROCplot( $roc, $im );

	    #my $coordinates_file = $output_dir.$similarity_string."_roc_".$fam."_coordinates.txt";
	    #open ( my $out, ">$coordinates_file" ) || die "can not open $coordinates_file!\n";
	    #print $out $coordinates;
	    #close $out;
        }
	
	my $similarity_sorted = sortSimilarity( \%similarity, \%xfp2motif, \%xfp2extended );
        open( OUT, ">>$report" ) || die "can not open $report!\n";
	#print OUT "\nrnaX(motif_n/extended_motif_n)\trnaY(motif_n/extended_motif_n)\t(common_motif_n/extended_common_motif_n)\tSimilarity\n";
	#for my $pair ( sort { $similarity_sorted->{$b} <=> $similarity_sorted->{$a} } keys %$similarity_sorted ) {
	#    print OUT "$pair\t$similarity_sorted->{$pair}\n";
	#} 
        print OUT "End of similarity function: $similarity_string\n\n";
        close OUT;
    }
}


sub sortSimilarity {
    my ( $s, $xfp2motif, $xfp2extended ) = @_;
    my %s = %{$s};
    my @rna = sort keys %s;
    my %s_sorted;
    my ( $s_basic, $d_basic ) = BinarySimilarity($xfp2motif);
    my ( $se_basic, $de_basic ) = BinarySimilarity($xfp2extended);

    for my $x ( 0 ... $#rna ) {
	for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my ( $famx ) = split /\./, $rnax; 
	    my $nmx = scalar keys %{$xfp2motif->{$rnax}};
	    my $nmex = scalar keys %{$xfp2extended->{$rnax}};
	    my $rnay = $rna[$y];
	    my ( $famy ) = split /\./, $rnay;
	    my $nmy = scalar keys %{$xfp2motif->{$rnay}};
	    my $nmey = scalar keys %{$xfp2extended->{$rnay}};
	    my $nmxy = $s_basic->{$rnax}{$rnay};
	    my $nmexy = $se_basic->{$rnax}{$rnay};
	    my $pair;
	    if ( $famx eq $famy ) {
	    	$pair = "$rnax($nmx/$nmex)\t$rnay($nmy/$nmey)\t($nmxy/$nmexy)";
	    } else {
		$pair = "##FALSE##  $rnax($nmx/$nmex)\t$rnay($nmy/$nmey)\t($nmxy/$nmexy)";
		next unless ( $nmx && $nmy );
		if ( (( $nmx eq $nmxy ) || ( $nmy eq $nmxy )) && ( $nmx/$nmy + $nmy/$nmx < 3 ) ) {
		    $pair = "##FALSE EXTREME##  $rnax($nmx/$nmex)\t$rnay($nmy/$nmey)\t($nmxy/$nmexy)";
		}
		if ( ( $nmx eq $nmxy ) && ( $nmy eq $nmxy ) ) {
		    $pair = "##FALSE IDENTICAL##  $rnax($nmx/$nmex)\t$rnay($nmy/$nmey)\t($nmxy/$nmexy)";
		}
	    }
	    $s_sorted{$pair}=$s{$rnax}{$rnay}; 
	} 
    }
    return (\%s_sorted);
}


sub readmanout {
    my ($manout) = @_;
    my %motif;

    open( MAN, "<$manout" ) || die "can not open $manout!\n";
    while ( my $line = <MAN> ) {
        next unless ( $line =~ /\_stem_motif\_/ );
        my @line1 = split " ", $line;
        my @line2 = split ":", $line1[0];
        my ($motif)       = $line2[1];
        my ($motif_type)  = $motif =~ /(\d)\_stem_motif/;
        my ($motif_index) = $motif =~ /stem_motif\_(.*)/;
        $motif = $motif_type . "_" . $motif_index;
        $motif{$motif} = 1;
    }
    close MAN;
    return ( \%motif );
}


sub DistanceToMega {

    my ( $distance, $mega ) = @_;
    my %distance = %{$distance};
    my @rna = sort keys %distance;
    my $total_n = scalar @rna;

    my $out_string = qq{
#mega
!Format DataType=Distance DataFormat=LowerLeft NTaxa=$total_n;

    };

    foreach my $n ( 1 ... $total_n ) {
	my $nn = $n - 1;
	my $rna = $rna[$nn];
	$out_string .= qq{[$n] #$rna\n};
    }
    $out_string .= qq{

    };
    $out_string .= qq{[\t};
    foreach my $n ( 1 ... $total_n ) {
	$out_string .= qq{$n\t};
    }
    $out_string .= qq{]\n};
    foreach my $n ( 1... $total_n ) {
	$out_string .= qq{[$n]\t};
	if ( $n > 1 ) {
	    foreach my $m ( 1 ... $n-1 ) { 
		my $d = $distance{$rna[$n-1]}{$rna[$m-1]};
		$out_string .= qq{$d\t};
	    }
	}
        $out_string .= qq{\n};
    }

    open ( my $out, ">$mega" ) or die "can not open $mega!\n";
    print $out $out_string;
    close $out;

}



sub SimilarityToMega {

    my ( $similarity, $mega ) = @_;
    my %similarity = %$similarity; 
    my %distance;
    my $s_max = 0;
    my @rna = sort keys %similarity;
    my $total_n = scalar @rna;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {

	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

	    my $s;
	    if ( defined $similarity{$rnax}{$rnay} ) {
	        $s = $similarity{$rnax}{$rnay}; 
 	    } elsif ( defined $similarity{$rnay}{$rnax} ) {
		$s = $similarity{$rnay}{$rnax}; 
	    }
	    $s_max = $s  if ( $s > $s_max );
	}
    }

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {

	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

	    my $s = $similarity{$rnax}{$rnay};
	    my $d;
	    if ( $rnax eq $rnay ) {
		$d = 0;
	    } else {
		$d = -log($s + 0.000000000001)/log(2);
		$d = sprintf( "%.15f", $d );
	    }
	    $distance{$rnax}{$rnay} = $d;
	    $distance{$rnay}{$rnax} = $d;
	}	    
    }

    my $out_string = qq{
#mega
!Format DataType=Distance DataFormat=LowerLeft NTaxa=$total_n;

    };

    foreach my $n ( 1 ... $total_n ) {
	my $nn = $n - 1;
	my $rna = $rna[$nn];
	$out_string .= qq{[$n] #$rna\n};
    }
    $out_string .= qq{

    };
    $out_string .= qq{[\t};
    foreach my $n ( 1 ... $total_n ) {
	$out_string .= qq{$n\t};
    }
    $out_string .= qq{]\n};
    foreach my $n ( 1... $total_n ) {
	$out_string .= qq{[$n]\t};
	if ( $n > 1 ) {
	    foreach my $m ( 1 ... $n-1 ) { 
		my $d = $distance{$rna[$n-1]}{$rna[$m-1]};
		$out_string .= qq{$d\t};
	    }
	}
        $out_string .= qq{\n};
    }

    open ( my $out, ">$mega" ) or die "can not open $mega!\n";
    print $out $out_string;
    close $out;

}





sub SimilarityToMatrix {

    my ( $similarity, $mtrx ) = @_;
    my %similarity = %$similarity; 
    my @rna = sort { $a cmp $b } keys %similarity;

    my $out_string;

    for my $rna ( @rna ) {
	my ( $fam ) = split /[\.\_]+/, $rna; 
	$out_string .= "$fam\t";
    }
    $out_string .= "\n";

    for my $rna ( @rna ) {
	$out_string .= "$rna\t";
    }
    $out_string .= "\n";


    for my $x ( @rna ) {
        for my $y ( @rna ) {
	    $out_string .= "$similarity{$x}{$y}\t";
	}
	$out_string .= "\n";
    }




    open ( my $out, ">$mtrx" ) or die "can not open $mtrx!\n";
    print $out $out_string;
    close $out;

}


#--------------------------------------------------------------------------------
# BinarySimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity: 
# 
# S(X,Y) = |X, Y intersection| 
#
# USAGE
#   ( $similarity ) = BinarySimilarity( $xfp2motif ); 
#--------------------------------------------------------------------------------
sub BinarySimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;
    my %distance;
    my $max_similarity = 0;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

            my %mx = %{$xfp2motif->{$rnax}};
            my %my = %{$xfp2motif->{$rnay}};

            my $similarity = 0;

            for my $m ( keys %mx ) {
                $similarity ++ if ( defined $my{$m} );
            }
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
	    $max_similarity = $similarity if ( $max_similarity < $similarity );

        }
    }

    for my $x ( 0... $#rna ) {
        for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];
	    if ( $x == $y ) {
		$distance{$rnax}{$rnay} = 0; 
	    } else {
		$distance{$rnax}{$rnay} = $max_similarity - $similarity{$rnax}{$rnay}; 
		$distance{$rnay}{$rnax} = $max_similarity - $similarity{$rnay}{$rnax}; 
	    }

    	}	 
    }


    return ( \%similarity, \%distance );
}
# End of BinarySimilarity



#--------------------------------------------------------------------------------
# CosineSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their distance and similarity: 
#
# S(X, Y) = |X && Y|/sqrt(|X||Y|)
#
# USAGE
#   ( $similarity ) = CosineSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub CosineSimilarity {

    my ($xfp2motif) = @_;
    my %xfp2motif = %{$xfp2motif};
    my @rna = sort keys %$xfp2motif;

    my %similarity;
    my %distance;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {

	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

            my %motif_x = %{$xfp2motif{$rnax}};
            my %motif_y = %{$xfp2motif{$rnay}};

	    my $similarity = 0;
	    my $xsize = scalar keys %motif_x;
	    my $ysize = scalar keys %motif_y;

            for my $motif ( keys %motif_x ) {
		$similarity++ if ( defined $motif_y{$motif} );
            }
	    $similarity = $similarity / ( $xsize * $ysize );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
	    if ( $x == $y ) {
		$distance{$rnax}{$rnay} = 0;
	    } else {
		$distance{$rnax}{$rnay} = 1 - $similarity;
		$distance{$rnay}{$rnax} = 1 - $similarity;
		
  	    }
        }
    }


    return ( \%similarity, \%distance );
}
# End of CosineSimilarity



#--------------------------------------------------------------------------------
# DiceSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity: 
# 
# S(X,Y) = 2 * |intersection of X & Y| / (|X| + |Y|)
#
# USAGE
#   ( $similarity ) = DiceSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub DiceSimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
	my $rnax = $rna[$x];
	my $mx = $xfp2motif->{$rnax};

        for my $y ( $x ... $#rna ) {

	    my $rnay = $rna[$y];
            my $my = $xfp2motif->{$rnay};

	    my ( $nmx, $nmy, $nmxy ) = ( 0, 0, 0 );
            for my $m ( keys %$mx ) {
                $nmxy+= ( $mx->{$m} )*( $my->{$m} ) if ( defined $my->{$m} );
		$nmx += ( $mx->{$m} )*( $mx->{$m} ); 
            }
	    $nmx = $nmx**0.5;
            for my $m ( keys %$my ) {
		$nmy += ( $my->{$m} )*( $my->{$m} ); 
            }
	    $nmy = $nmy**0.5;
            my $similarity =
		2 * $nmxy / ( $nmx + $nmy + 0.00000000001 );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of DiceSimilarity



#--------------------------------------------------------------------------------
# HammingSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their distance and similarity: 
# 
# D(X,Y) = # of positions at which the corresponding motifs are different 
#
# S(X,Y) = total # of motifs in database ( default: N = 55728) - D(X,Y) 
#        = N - |X| - |Y| + 2 * |X & Y intersection|
# USAGE
#   ( $similarity ) = HammingSimilarity($xfp2motif, $ndb); 
#--------------------------------------------------------------------------------
sub HammingSimilarity {

    my ($xfp2motif, $ndb) = @_;

    my $DEFAULT_N = 55728;
    unless ( defined $ndb ) {
	$ndb = $DEFAULT_N;
    }

    my @rna = sort keys %$xfp2motif;

    my %similarity;

    for my $x ( 0 ... $#rna ) {
	my $rnax = $rna[$x];
	my $mx = $xfp2motif->{$rnax};
	my $nmx = scalar keys %$mx;

        for my $y ( $x ... $#rna ) {

	    my $rnay = $rna[$y];
            my $my = $xfp2motif->{$rnay};
	    my $nmy = scalar keys %$my;

	    my $nmxy = 0;
            for my $m ( keys %$mx ) {
                $nmxy++ if ( defined $my->{$m} );
            }
            my $similarity =
		$ndb - $nmx - $nmy + 2 * $nmxy;
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }
    return ( \%similarity );
}

# End of HammingSimilarity



#--------------------------------------------------------------------------------
# TanimotoSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their distance and similarity: 
#
# S(X, Y) = |X && Y|/|X or Y|
#
# USAGE
#   ( $similarity ) = TanimotoSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub TanimotoSimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
	my $rnax = $rna[$x];
	my $mx = $xfp2motif->{$rnax};

        for my $y ( $x ... $#rna ) {

	    my $rnay = $rna[$y];
            my $my = $xfp2motif->{$rnay};

	    my ( $nmx, $nmy, $nmxo, $nmyo, $nmxy ) = (0, 0, 0, 0, 0);
            for my $m ( keys %$mx ) {
		$nmx += ( $mx->{$m} )**2;
                if ( defined $my->{$m} ) {
		    $nmxy+= ( $mx->{$m} )* ( $my->{$m} );
		} else {
		    $nmxo+= ( $mx->{$m} )**2;
		}
            }
	    for my $m ( keys %$my ) {
		$nmy += ( $my->{$m} )**2;
		$nmyo+= ( $my->{$m} )**2 unless ( defined $mx->{$m} );
	    }
            my $similarity =
		$nmxy /( $nmx + $nmy - $nmxy + 0.000000001 ); 
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of TanimotoSimilarity



#--------------------------------------------------------------------------------
# VertexNumberSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity 
# according to their vertex numbers: 
# 
# S(X,Y) = N_verticesX * N_verticesY / (N_verticesX + N_verticesY)**2
#
# USAGE
#   ( $similarity ) = VertexNumberSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub VertexNumberSimilarity {

    my ($xfp2vertex) = @_;
    my @rna = sort keys %$xfp2vertex;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

	    my $vnx = $xfp2vertex{$rnax};
	    my $vny = $xfp2vertex{$rnay};

            my $similarity = min( $vnx, $vny ) / max( $vnx, $vny );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }
    return ( \%similarity );
}
# End of VertexNumberSimilarity


#--------------------------------------------------------------------------------
# EdgeNumberSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity 
# according to their edge numbers: 
# 
# S(X,Y) = N_edgeX * N_edgeY / (N_edgeX + N_edgeY)**2
#
# USAGE
#   ( $similarity ) = EdgeNumberSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub EdgeNumberSimilarity {

    my ($xfp2edge) = @_;
    my @rna = sort keys %$xfp2edge;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

	    my $enx = $xfp2edge{$rnax};
	    my $eny = $xfp2edge{$rnay};

	    my $similarity = min( $enx, $eny ) / max( $enx, $eny );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}
# End of EdgeNumberSimilarity








sub similarity4ROCAllFam {

    my ($similarity) = @_;
    my %s = %{$similarity};

    my $fam = "total";

    my ( @true, @false, @total );

    my @rna = sort keys %s;

    foreach my $x ( 0 ... $#rna ) {
        foreach my $y ( $x ... $#rna ) {

            my $rnax = $rna[$x];
            my $rnay = $rna[$y];

	    my ( $famx ) = split /[\_,\.]+/, $rnax; 
	    my ( $famy ) = split /[\_,\.]+/, $rnay; 

            # TODO test on the undef similarity!! 
            my $similarity = $s{$rnax}{$rnay};

            if ( "$famx" eq "$famy" ) {
                push @true, $similarity;
                push @total, [ "true", $similarity ];
            }
            else {
                push @false, $similarity;
                push @total, [ "false", $similarity ];
            }

        }
    }

    @true  = sort @true;
    @false = sort @false;
    @total = sort { $b->[1] <=> $a->[1]; } @total;

    my $plot = checkIdentical4ROC( \@total );

    my @plots;
    push @plots, $plot;
    return ( \@plots );
}

sub similarity4ROCbyFam {

    my ( $fam, $similarity ) = @_;
    my %s = %{$similarity};

    my ( @true, @false, @total );

    my @rna = sort keys %s;

    foreach my $x ( 0 ... $#rna ) {

        foreach my $y ( $x ... $#rna ) {

            my $rnax = $rna[$x];
            my $rnay = $rna[$y];

	    my ( $famx ) = split /[\.,\_]+/, $rnax; 
	    my ( $famy ) = split /[\.,\_]+/, $rnay; 


            next if ( ( "$famx" ne "$fam" ) && ( "$famy" ne "$fam" ) );

            if ( "$famx" eq "$famy" ) {
                push @true, $s{$rnax}{$rnay};
                my $similarity = $s{$rnax}{$rnay};
                push @total, [ "true", $similarity ];
            }
            else {
                push @false, $s{$rnax}{$rnay};
                my $similarity = $s{$rnax}{$rnay};
                push @total, [ "false", $similarity ];
            }

        }
    }

    @true  = sort @true;
    @false = sort @false;
    @total = sort { $b->[1] <=> $a->[1]; } @total;

    my $plot = checkIdentical4ROC( \@total );

    my @plots;
    push @plots, $plot;
    return ( \@plots );
}

sub checkIdentical4ROC {

    my ($total) = @_;
    my @total = @{$total};

    my @total_idxed;

    my @true;
    my @false;

    my $idx = 0;
    for my $pair (@total) {

        # add index of each pair so the format now is
        # [index, true/false, distance]
        push @total_idxed, [ $idx, $pair->[0], $pair->[1] ];
        if ( $pair->[0] eq "true" ) {
            push @true, [ $idx, $pair->[0], $pair->[1] ];
        }
        else {
            push @false, [ $idx, $pair->[0], $pair->[1] ];
        }
        $idx++;
    }

    my ( $true_percent, $false_percent );
    if ( @true ) { $true_percent = 1/($#true + 1); } 
    else { $true_percent = 0; }
    if ( @false ) { $false_percent = 1/($#false + 1); }

    @total = @total_idxed;

    my %plot;

    my $dist  = -1;
    my $index = -1;

    my %identicalidx;
    my %identicaltf;

    while (@total) {
        my $pair = shift @total;

        if ( $pair->[2] != $dist ) {
            push @{ $identicalidx{ $pair->[0] } }, $pair->[0];
            push @{ $identicaltf{ $pair->[0] } },  $pair->[1];
            $dist  = $pair->[2];
            $index = $pair->[0];
        }
        else {
            if ( defined $index ) {
                push @{ $identicalidx{$index} }, $pair->[0];
                push @{ $identicaltf{$index} },  $pair->[1];
            }
        }
    }

    for my $i ( keys %identicalidx ) {
        my @identical = @{ $identicalidx{$i} };
        my @TF        = @{ $identicaltf{$i} };
        my $size      = @identical;
        my ( $true_size, $false_size );
        for my $tf (@TF) {
            $true_size++  if ( $tf eq "true" );
            $false_size++ if ( $tf eq "false" );
        }
        if ( $size <= 1 ) {
            push @{ $plot{$i} }, ( $false_percent, 0 )
              if ( $total_idxed[$i]->[1] eq "false" );
            push @{ $plot{$i} }, ( 0, $true_percent )
              if ( $total_idxed[$i]->[1] eq "true" );

        }
        elsif ( $true_size && $false_size ) {

            for my $ident (@identical) {
                my $false_move = $false_percent * $false_size / $size;
                my $true_move  = $true_percent * $true_size / $size;
                push @{ $plot{$ident} }, ( $false_move, $true_move );
            }

        }
        else {

            for my $ident (@identical) {
                if ( $total_idxed[$ident]->[1] eq "false" ) {
                    push @{ $plot{$ident} }, ( $false_percent, 0 );
                }
                elsif ( $total_idxed[$ident]->[1] eq "true" ) {
                    push @{ $plot{$ident} }, ( 0, $true_percent );
                }
            }
        }
    }

    return ( \%plot );

}

sub ROC {

    my ( $plots, $report, $fam ) = @_;
    my @plots = @{$plots};
    my $coordinates = qq{FPR\tTPR\n};

    # ROC

    # create new image object
    my $im = new GD::Image( 800, 800 );

    # color
    my $white    = $im->colorAllocate( 255, 255, 255 );
    my $black    = $im->colorAllocate( 0,   0,   0 );
    my $red      = $im->colorAllocate( 255, 0,   0 );
    my $blue     = $im->colorAllocate( 0,   0,   255 );
    my $green    = $im->colorAllocate( 0,   128, 0 );
    my $yellow   = $im->colorAllocate( 255, 255, 0 );
    my $macgenta = $im->colorAllocate( 255, 0,   255 );

    my @colors = ( $yellow, $green, $macgenta, $blue, $red );

    my $pixels_per_unit = 700;

    # axes
    my $x_axis_origin = 50;
    my $x_axis_end    = 750;
    my $y_axis_origin = 750;
    my $y_axis_end    = 50;

    # x-axis
    $im->line( $x_axis_origin, $y_axis_origin, $x_axis_end, $y_axis_origin,
        $black );

    # line of (0,1) to (1,1)
    $im->line( $x_axis_origin, $y_axis_end, $x_axis_end, $y_axis_end, $black );

    # x-axis label
    $im->string( gdGiantFont, 300,
        $y_axis_origin + 25,
        "False Positive Rate", $black
    );

    # y-axis
    $im->line( $x_axis_origin, $y_axis_origin, $x_axis_origin, $y_axis_end,
        $black );

    # line of (1,0) to (1,1)
    $im->line( $x_axis_end, $y_axis_origin, $x_axis_end, $y_axis_end, $black );

    # y-axis label
    $im->stringUp( gdGiantFont, $x_axis_origin - 50,
        500, "True Positive Rate", $black );

    my @x_ticks = ( 0, 0.20, 0.40, 0.60, 0.80, 1.00 );
    foreach my $tick (@x_ticks) {
        my $tick_pos = $tick * $pixels_per_unit + $x_axis_origin;
        $im->line(
            $tick_pos, $y_axis_origin - 3, $tick_pos, $y_axis_origin + 3,
            $black
        );
        $im->string( gdSmallFont,
            $tick_pos - 5,
            $y_axis_origin + 10,
            $tick, $black
        );
    }

    my @y_ticks = ( 0, 0.20, 0.40, 0.60, 0.80, 1.00 );
    foreach my $tick (@y_ticks) {
        my $tick_pos = $y_axis_origin - $tick * $pixels_per_unit;
        $im->line(
            $x_axis_origin - 3, $tick_pos, $x_axis_origin + 3, $tick_pos,
            $black
        );
        $im->stringUp( gdSmallFont,
            $x_axis_origin - 20,
            $tick_pos + 10,
            $tick, $black
        );
    }

    # line indicating random
    $im->dashedLine( $x_axis_origin, $y_axis_origin, $x_axis_end, $y_axis_end,
        $black );

    for my $plot (@plots) {

        my %plot = %{$plot};

        my $color = pop @colors;

        my %pos;
        my ( $false_pos, $true_pos ) = ( 0, 0 );
        push @{ $pos{0} }, ( $false_pos, $true_pos );

        my @idx = sort keys %plot;
        foreach my $i ( 0 ... $#idx ) {

            my ( $false_move, $true_move ) = @{ $plot{$i} };
            $false_pos = $false_pos + $false_move;
            $true_pos  = $true_pos + $true_move;

            push @{ $pos{ $i + 1 } }, ( $false_pos, $true_pos );
	    $coordinates .= qq{$false_pos $true_pos\n};
        }

        my $AUC  = 0;
        my @idxx = sort keys %pos;

        foreach my $i ( 0 ... $#idxx - 1 ) {

            my ( $false_pos,      $true_pos )      = @{ $pos{$i} };
            my ( $false_pos_next, $true_pos_next ) = @{ $pos{ $i + 1 } };
            my $area =
              0.5 *
              ( $false_pos_next - $false_pos ) *
              ( $true_pos_next + $true_pos );
            $AUC = $AUC + $area;

            my $x_pos->[0] = $false_pos * $pixels_per_unit + $x_axis_origin;
            my $y_pos->[0] = $y_axis_origin - $true_pos * $pixels_per_unit;
            $x_pos->[1] = $false_pos_next * $pixels_per_unit + $x_axis_origin;
            $y_pos->[1] = $y_axis_origin - $true_pos_next * $pixels_per_unit;
            $im->line( $x_pos->[0], $y_pos->[0], $x_pos->[1], $y_pos->[1],
                $color );
        }

        my $legend = sprintf "$fam AUC = %1.4f", $AUC;

        $im->string( gdSmallFont, 10, 20, $legend, $black );

        print STDERR "$legend\n";

        open( OUT, ">>$report" ) || die "can not open $report!\n";
        print OUT "$legend\n";
        close OUT;

    }

    return ($im, $coordinates);

}

sub ROCplot {

    my ( $rocplot, $im ) = @_;

    open( GRAPH, ">$rocplot" ) || die "can not open graph file: $!\n";
    binmode GRAPH;
    print GRAPH $im->png;
    close GRAPH;
}



sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}


#------------------------------------------------------------------------------
# $Log: fingerprint_distance.pl,v $
#------------------------------------------------------------------------------
