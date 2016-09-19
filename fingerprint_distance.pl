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
#use lib '/home/huang147/RNA/src/perl_src';
#use lib '/home/huang147/perl_src_huang_120827/';
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
    -f      weight file
   	-h  	display this help message
	-i|m	original fingerprint
	-i|e	extended fingerprint
	-r|a    remove any random motifs in each fingerprint
	-r|f    remove frequent motifs in each fingerprint
	-r|r    remove rare motifs in each fingerprint
	-r|m    remove middle motifs in each fingerprint
	-r|n    no removal 
	-p	    fraction of motifs (in xx%) to remove
  	-s|a 	AntiProbabilitySimilarity
   	-s|b	BinarySimilarity
	-s|c	CosineSimilarity
	-s|d	DiceSimilarity
	-s|e	EdgeNumberSimilarity
	-s|h	HammingSimilarity
	-s|l	MarylandBridgeSimilarity
	-s|m	MotifNumberSimilarity
	-s|n	ParentTanimotiSimilarity
	-s|r	ParentCountSimilarity
	-s|t	TanimotoSimilarity
	-s|v	VertexNumberSimilarity
	-s|x	ExtendedTanimotoSimilarity
    -w      output directory for clustering results
    -z      percentage of stems to remove

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
my %xfp2prob  = ();
my %xfp2tfidf = ();
my %xfp2idf = ();
my %xfp2idfextended = ();
my %extended_sp = ();
my %sp = ();
my %prob_index = ();


my %option;
getopts( "d:f:hi:s:w:z:", \%option );

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


my $weight_file; 
if ( defined $option{f} ) {
    $weight_file = $option{f};
}



my $DEFAULT_REMOVE_PERCENTAGES = ""; 
my $current_remove_percentages = $DEFAULT_REMOVE_PERCENTAGES;
if ( defined $option{z} ) {
    if ( length($option{z}) % 4   ) {
	print STDERR "Please input a string that contains pairs of percentages!\n";
	exit 4;
    } else {
    	$current_remove_percentages = $option{z};
    }
}  
if ( $current_remove_percentages =~ /\d/ ) {
    my @remove_percentages = $current_remove_percentages =~ /\d{4}/g;
    open( OUT, ">$report" ) || die "can not open $report!\n";
    print OUT "\nPercentages of Motifs (Sorted by Probability) to Remove:\n";
    print STDERR "\nPercentages (start::end) of Motifs (Sorted by Probability) to Remove:\n";
    print OUT "@remove_percentages\n";
    print STDERR "@remove_percentages\n";
    close OUT;
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
    #my ( $fam, $subfam ); 
    #my ($fam_subfam, $species) = split /\./, $xfp;
    #if ( $fam_subfam =~ /\_/ ) {
#	my @fam_subfam = split /\_/, $fam_subfam;
#	$fam = $fam_subfam[0];
#	$subfam = join "_", @fam_subfam[1...$#fam_subfam];
#    } else {
#	$fam = $fam_subfam; 
#    }
#    my $name; 
#    if ( $subfam ) {
#        $name = $fam.".".$subfam."_".$species; 
#    } else {
#	$name = $fam.".".$species; 
#    }
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
	my $prob = $count/$iteration;

    	unless ( $motif_prob_index{$id} ) {
	    $motif_prob_index{$id} = {
		id			=> $id,
		indicator		=> 1,
		prob			=> $prob,
		prob_aver_all		=> undef,
		sp_frac_all		=> undef,
		tf			=> undef,
		idf			=> undef, 
		tfidf 			=> undef,
	    };
	} 



    	unless ( $xfp2prob{$fam}{$id} ) {
	    $xfp2prob{$fam}{$id} = {
		id			=> $id,
		prob			=> undef,
		prob_aver	 	=> undef,
		support 		=> 0,
		support_frac		=> undef,
		support_ratio   	=> undef,
	    };
	} 
	push @{$xfp2prob{$fam}{$id}{prob}}, $prob;
	$xfp2prob{$fam}{$id}{support}++;


    	unless ( $xfp2prob{all}{$id} ) {
	    $xfp2prob{all}{$id} = {
		id			=> $id,
		prob			=> undef,
		prob_aver 		=> undef,
		support		 	=> 0,
		support_frac		=> undef,
		support_ratio		=> undef,
	    };
	} 
	push @{$xfp2prob{all}{$id}{prob}}, $prob;
	$xfp2prob{all}{$id}{support}++;
	$sp{$id}++;


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
    $prob_index{$name} = \%motif_prob_index;
    $xfpfam{all}++;
    $xfpfam{$fam}++;
}
print STDERR "Graph Stats Calculation Done!\n";
my $graph_stats = $output_dir."graph_stats.txt";
open ( G, ">$graph_stats" ) or die "can not open $graph_stats!\n";
print G $graph_stats_string;
close G;


foreach my $fam ( sort keys %xfp2prob ) {

    my $prob_stats_string = qq{id\tprob_aver\tprob_all\tsupport_frac\tsupport_all\tsupport_ratio\n};

    my $pf = $xfp2prob{$fam};
    my $nfam = $xfpfam{$fam};
    my $nall = $xfpfam{all};

    foreach my $id ( keys %$pf ) {
	my $support_frac = $$pf{$id}{support} / $nfam; 
	my $support_all = $xfp2prob{all}{$id}{support}/$nall;
	my $support_ratio = $support_frac/$support_all;
	my $prob_all = sum( @{$xfp2prob{all}{$id}{prob}} )/$nall; 
	my $prob_aver = sum( @{$$pf{$id}{prob}} )/$nfam; 
	$$pf{$id}{support_frac} = $support_frac;
	$$pf{$id}{support_ratio} = $support_ratio;
	$$pf{all}{$id}{support_ratio} = $support_ratio;
	$$pf{$id}{prob_aver} = $prob_aver; 
	$xfp2prob{all}{$id}{prob_aver} = $prob_all; 
	$prob_stats_string .= qq{$id\t$prob_aver\t$prob_all\t$support_frac\t$support_all\t$support_ratio\n};
    }

    #my $prob_stats = $output_dir.$fam."_motif_prob.txt";
    #open ( my $ps, ">$prob_stats" ) or die "can not open $prob_stats!\n";
    #print $ps $prob_stats_string;
    #close $ps;
}




foreach my $name ( sort keys %prob_index  ) { 

    my %motif_prob_index = %{$prob_index{$name}};
    my %tfidf = ();
    my %idf = ();
    my %extended_idf = ();
    #my @sps = sort { $a <=> $b } values %sp;
    #my @probs = sort { $a<=>$b } values %probs;
    #my $prob_max = max( values %probs ); 
    #my $sp_lo = $sps[int(@sps*0.25)];
    #my $sp_up = $sps[-int(@sps*0.25)];
    #my $prob_lo = $probs[int(@probs*0.25)];
    #my $prob_up = $probs[-1];
    #my $prob_lo = $probs[0];
    #my $prob_up = $probs[-int(@probs*0.25)];

    foreach my $id ( sort keys %motif_prob_index ) {

	$motif_prob_index{$id}{sp_frac_all} = $xfp2prob{all}{$id}{support_frac};
	$motif_prob_index{$id}{prob_aver_all} = $xfp2prob{all}{$id}{prob_aver};
	$motif_prob_index{$id}{extended_sp_frac_all} = $xfp2prob{all}{$id}{extended_sp_frac};


	#my $tf = log( $motif_prob_index{$id}{prob} )/log(10);
	#my $idf = log( 1 / $motif_prob_index{$id}{sp_frac_all} ) / log(10); 
	my $tf = $motif_prob_index{$id}{prob}; 
	my $idf = 1- $motif_prob_index{$id}{sp_frac_all}; 
	my $tfidf = $tf * $idf;
	#print STDERR "$name\t$id\t$tf\t$idf\t$tfidf\t$extended_idf\n";
	$idf{$id} = $idf;
	$tfidf{$id} = $tfidf;
	#print STDERR "$probs{$id}\t$prob_lo\t$prob_up\n";
	#my $tf = 1;
	#my $idf = 1;
	#$idf = 0 if ( $probs{$id} < $prob_lo || $probs{$id} > $prob_up );
	#$idf = 0 if ( $probs{$id} > $prob_lo && $probs{$id} < $prob_up );
	#$idf = 0 if ( $sp{$id} < $sp_lo || $sp{$id} < $sp_up );

    }

    my $total_support = scalar @xfp; 
    foreach my $p ( sort keys %{$xfp2extended{$name}} ) {
	my $extended_idf = 1 - $extended_sp{$p}/$total_support; 
	$extended_idf{$p} = $extended_idf;
    }


    $xfp2idfextended{$name} = \%extended_idf;
    $xfp2tfidf{$name} = \%tfidf;
    $xfp2idf{$name} = \%idf;

}



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

    } elsif ( $index eq "t" ) {

        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\n\nSimple Fingerprint with TF-IDF Weighting\n";
	print STDERR "\n\nSimple Fingerprint with TF-IDF Weighting\n";
        close OUT;
        Similarity2ROC( \%xfp2tfidf, \%xfp2vertex, \%xfp2edge,  $current_similarity );

    } elsif ( $index eq "i" ) {

        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\n\nSimple Fingerprint with IDF Weighting\n";
	print STDERR "\n\nSimple Fingerprint with IDF Weighting\n";
        close OUT;
        Similarity2ROC( \%xfp2idf, \%xfp2vertex, \%xfp2edge,  $current_similarity );

    } elsif ( $index eq "x" ) {

        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\n\nExtended Fingerprint with IDF Weighting\n";
	print STDERR "\n\nExtended Fingerprint with IDF Weighting\n";
        close OUT;
        Similarity2ROC( \%xfp2idfextended, \%xfp2vertex, \%xfp2edge,  $current_similarity );
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
        elsif ( $simi_func eq "t" ) {

            $similarity_string = "TanimotoSimilarity";
            ( $similarity, $distance ) = TanimotoSimilarity($xfp2motif);

        }
        elsif ( $simi_func eq "c" ) {

            $similarity_string = "CosineSimilarity";
            ( $similarity, $distance ) = CosineSimilarity($xfp2motif);

        }

        elsif ( $simi_func eq "n" ) {

            $similarity_string = "ParentTanimotoSimilarity";
            ( $similarity ) = ParentTanimotoSimilarity($xfp2motif);
	    #( $similarity ) = TanimotoSimilarity( $xfp2motif );
        }

        elsif ( $simi_func eq "r" ) {

            $similarity_string = "ParentCountSimilarity";
            ( $similarity ) = ParentCountSimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "d" ) {

            $similarity_string = "DiceSimilarity";
            ( $similarity ) = DiceSimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "h" ) {

            $similarity_string = "HammingSimilarity";
            ( $similarity ) = HammingSimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "a" ) {

            $similarity_string = "AntiProbabilitySimilarity";
            ( $similarity ) = AntiProbabilitySimilarity($xfp2motif);
        }

        elsif ( $simi_func eq "x" ) {

            $similarity_string = "ExtendedTanimotoSimilarity";
            ( $similarity, $distance ) = ExtendedTanimotoSimilarity($xfp2motif);
            #( $similarity, $distance ) = TanimotoSimilarity(\%xfp2extended);
        }

        elsif ( $simi_func eq "m" ) {

            $similarity_string = "MotifNumberSimilarity";
            ( $similarity ) = MotifNumberSimilarity($xfp2motif);
        }


        elsif ( $simi_func eq "q" ) {

            $similarity_string = "SquareSimilarity";
            ( $similarity ) = SquareSimilarity($xfp2motif);
        }

	elsif ( $simi_func eq "w" ) {
	    $similarity_string = "WeightedTanimotoSimilarity";
	    ( $similarity ) = WeightedTanimotoSimilarity( $xfp2motif, $weight_file );
	}

	elsif ( $simi_func eq "v" ) {
	    $similarity_string = "VertexNumberSimilarity";
	    ( $similarity ) = VertexNumberSimilarity( $xfp2vertex );
	}

	elsif ( $simi_func eq "e" ) {
	    $similarity_string = "EdgeNumberSimilarity";
	    ( $similarity ) = EdgeNumberSimilarity( $xfp2edge );
	}

	elsif ( $simi_func eq "l" ) {
	    $similarity_string = "MBSimilarity";
	    ( $similarity ) = MBSimilarity( $xfp2edge );
	}



	my $mega = $output_dir.$similarity_string.".meg";
	if ( defined $distance ) {
	    DistanceToMega( $distance,  $mega );
	} else {
	    SimilarityToMega( $similarity, $mega );
	}
	my $mtx = $output_dir.$similarity_string.".mtx";
	SimilarityToMatrix( $similarity, $mtx );

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
	
	my $coordinates_total_file = $output_dir.$similarity_string."_roc_total_coordinates.txt";
	open ( my $out, ">$coordinates_total_file" ) || die "can not open $coordinates_total_file!\n";
	print $out $coordinates_total;
	close $out;

        for my $fam ( keys %xfpfam ) {
            my $plot = similarity4ROCbyFam( $fam, \%similarity );

            my $roc = $output_dir."$similarity_string" . "_" . "$fam" . "_ROC.png";
            my ( $im, $coordinates) = ROC( $plot, $report, $fam );
            ROCplot( $roc, $im );

	    my $coordinates_file = $output_dir.$similarity_string."_roc_".$fam."_coordinates.txt";
	    open ( my $out, ">$coordinates_file" ) || die "can not open $coordinates_file!\n";
	    print $out $coordinates;
	    close $out;
        }
	
	my $similarity_sorted = sortSimilarity( \%similarity, \%xfp2motif, \%xfp2extended );
        open( OUT, ">>$report" ) || die "can not open $report!\n";
	print OUT "\nrnaX(motif_n/extended_motif_n)\trnaY(motif_n/extended_motif_n)\t(common_motif_n/extended_common_motif_n)\tSimilarity\n";
	for my $pair ( sort { $similarity_sorted->{$b} <=> $similarity_sorted->{$a} } keys %$similarity_sorted ) {
	    print OUT "$pair\t$similarity_sorted->{$pair}\n";
	} 
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

#    for my $x ( 0 ... $#rna ) {
#        for my $y ( $x ... $#rna ) {

#	    my $rnax = $rna[$x];
#	    my $rnay = $rna[$y];

#	    my $s;
#	    if ( defined $similarity{$rnax}{$rnay} ) {
#	        $s = $similarity{$rnax}{$rnay}; 
# 	    } elsif ( defined $similarity{$rnay}{$rnax} ) {
#		$s = $similarity{$rnay}{$rnax}; 
#	    }
#	    $s_max = $s  if ( $s > $s_max );
#	}
#    }

#    for my $x ( 0 ... $#rna ) {
#        for my $y ( $x ... $#rna ) {

#	    my $rnax = $rna[$x];
#	    my $rnay = $rna[$y];

	    #my $s = $similarity{$rnax}{$rnay};
#	    my $d = $similarity{$rnax}{$rnay};
#	    if ( $rnax eq $rnay ) {
#		$d = 0;
#	    } else {
#	        #$d = $s_max - $s;
#		$d = -log($s + 0.000000000001)/log(2);
#		$d = sprintf( "%.15f", $d );
#	    }
#	    $distance{$rnax}{$rnay} = $d;
#	    $distance{$rnay}{$rnax} = $d;
#	}	    
#    }

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
	        #$d = $s_max - $s;
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
	#my ( $fam ) = $rna =~ /(.*)\_/;
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

            #my $similarity = $vnx * $vny / ($vnx + $vny)**2;
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

            #my $similarity = $enx * $eny / ($enx + $eny)**2;
	    my $similarity = min( $enx, $eny ) / max( $enx, $eny );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}
# End of EdgeNumberSimilarity


#--------------------------------------------------------------------------------
# MotifNumberSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity 
# according to their motif numbers: 
# 
# S(X,Y) = N_motifX * N_motifY / (N_motifX + N_motifY)**2
#
# USAGE
#   ( $similarity ) = MotifNumberSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub MotifNumberSimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {

	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

            my $nmx = scalar keys %{$xfp2motif{$rnax}};
            my $nmy = scalar keys %{$xfp2motif{$rnay}};

            #my $similarity = $nmx * $nmy / ($nmx+$nmy)**2;
	    my $similarity = min( $nmx, $nmy ) / max( $nmx, $nmy );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }
    return ( \%similarity );
}

# End of MotifNumberSimilarity


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
	#my $nmx = scalar keys %$mx;

        for my $y ( $x ... $#rna ) {

	    my $rnay = $rna[$y];
            my $my = $xfp2motif->{$rnay};
	    #my $nmy = scalar keys %$my;

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
# RemoveMotifs
# 
# Given a list of motifs in one fingerprint, remove some of them according to 
#	user-defined criteria.
#  
# USAGE
#   ( $motif_index, $n_motif ) = RemoveMotifs(\%motif_index, $string_remove); 
#--------------------------------------------------------------------------------
sub RemoveMotifs {

    my ( $mr, $remove_percentages ) = @_;

    my @m = sort { $mr->{$a} <=> $mr->{$b} } keys %$mr;
    my $nm = scalar @m;

    my @remove_percentages = $remove_percentages =~/\d{4}/g;
    my @remove;
    foreach my $p ( @remove_percentages ) {
	my ( $p1, $p2 ) = $p =~ /\d{2}/g;
	my $r1 = int($p1*0.1*$nm);
	my $r2 = int($p2*0.1*$nm) - 1;
	push @remove, @m[$r1...$r2];
    }

    for my $m ( @remove ) {
        delete $mr->{$m};           
    }

    my $n_motif = scalar keys %$mr;

    return ( $mr, $n_motif );
}

# End of RemoveMotifs    



#--------------------------------------------------------------------------------
# RemoveRareMotifs
# 
# Given a list of motifs in one fingerprint, remove the rare motifs (with low 
#   probabilities to be found in random sampling) with a user-defined fraction. 
#
# USAGE
#   ( $motif_index, $n_motif ) = RemoveRareMotifs(\%motif_index, $frac_remove); 
#--------------------------------------------------------------------------------
sub RemoveRareMotifs {

    my ( $mr, $frac_rm ) = @_;

    my @counts = sort { $b <=> $a } values %$mr;
    my $count_keep = $counts[int(($#counts +1)*( 1 - $frac_rm ))]; 

    for my $m ( sort { $mr->{$a} <=> $mr->{$b} } keys %$mr  ) {
        delete $mr->{$m} if ( $mr->{$m} < $count_keep );           
    }

    my $n_motif = scalar keys %$mr;

    return ( $mr, $n_motif );
}

# End of RemoveRareMotifs    



#--------------------------------------------------------------------------------
# RemoveFreqMotifs
# 
# Given a list of motifs in one fingerprint, remove the frequent motifs (with high
#   probabilities to be found in random sampling) with a user-defined fraction. 
#
# USAGE
#   ( $motif_index, $n_motif ) = RemoveFreqMotifs(\%motif_index, $frac_remove); 
#--------------------------------------------------------------------------------
sub RemoveFreqMotifs {

    my ( $mr, $frac_rm ) = @_;

    my @counts = sort { $a <=> $b } values %$mr;
    my $count_keep = $counts[int(($#counts +1)*( 1 - $frac_rm ))]; 

    for my $m ( sort { $mr->{$a} <=> $mr->{$b} } keys %$mr  ) {
        delete $mr->{$m} if ( $mr->{$m} > $count_keep );           
    }

    my $n_motif = scalar keys %$mr;

    return ( $mr, $n_motif );
}
# End of RemoveFreqMotifs




#--------------------------------------------------------------------------------
# RemoveAnyMotifs
# 
# Given a list of motifs in one fingerprint, remove the frequent motifs (with high
#   probabilities to be found in random sampling) with a user-defined fraction. 
#
# USAGE
#   ( $motif_index, $n_motif ) = RemoveFreqMotifs(\%motif_index, $frac_remove); 
#--------------------------------------------------------------------------------
sub RemoveAnyMotifs {

    my ( $mr, $frac_rm ) = @_;

    my @m = keys %$mr; 
    @m = shuffle @m;
    my @rm = @m[ 0 ... int(($#m + 1)*($frac_rm))];
    if ( $#rm == 0 ) {@rm = ();}
    #my @counts = sort { $a <=> $b } values %$mr;
    #my $count_keep = $counts[int(($#counts +1)*( 1 - $frac_rm ))]; 

    for my $m ( @rm ) {
        delete $mr->{$m};           
    }

    my $n_motif = scalar keys %$mr;

    return ( $mr, $n_motif );
}
# End of RemoveFreqMotifs



#--------------------------------------------------------------------------------
# RemoveMidMotifs
# 
# Given a list of motifs in one fingerprint, remove the middle motifs (with neither 
# low or high probabilities to be found in random sampling) with a user-defined fraction. 
#
# USAGE
#   ( $motif_index, $n_motif ) = RemoveMidMotifs(\%motif_index, $frac_remove); 
#--------------------------------------------------------------------------------
sub RemoveMidMotifs {

    my ( $mr, $frac_rm ) = @_;

    my $m_motif = scalar keys %$mr;

    my @counts = sort { $a <=> $b } values %$mr;
    my $count_keep_lo = $counts[int(($#counts +1)*( 0.5 - $frac_rm/2 ))]; 
    my @counts_rev = sort { $b <=> $a } values %$mr;
    my $count_keep_hi = $counts_rev[int(($#counts +1)*( 0.5 - $frac_rm/2 ))]; 


    for my $m ( sort { $mr->{$a} <=> $mr->{$b} } keys %$mr  ) {
        delete $mr->{$m} unless ( $mr->{$m} >= $count_keep_hi || $mr->{$m} <= $count_keep_lo );           
    }

    my $n_motif = scalar keys %$mr;
    #print STDERR "bfr $m_motif aft $n_motif\n";
    #if ( $n_motif == 0 ) {
    #    print STDERR "lo $count_keep_lo hi $count_keep_hi  bfr $m_motif aft $n_motif\n";
    #}
    return ( $mr, $n_motif );
}
# End of RemoveMidMotifs




#--------------------------------------------------------------------------------
# RemoveRareAndFreqMotifs
# 
# Given a list of motifs in one fingerprint, remove the middle motifs (with neither 
# low or high probabilities to be found in random sampling) with a user-defined fraction. 
#
# USAGE
#   ( $motif_index, $n_motif ) = RemoveMidMotifs(\%motif_index, $frac_remove); 
#--------------------------------------------------------------------------------
sub RemoveRareAndFreqMotifs {

    my ( $mr, $frac_rm ) = @_;

    my $m_motif = scalar keys %$mr;

    my @counts = sort { $a <=> $b } values %$mr;
    my $count_keep_lo = $counts[int(($#counts +1)*( 0.5 - $frac_rm/2 ))]; 
    my @counts_rev = sort { $b <=> $a } values %$mr;
    my $count_keep_hi = $counts_rev[int(($#counts +1)*( 0.5 - $frac_rm/2 ))]; 


    for my $m ( sort { $mr->{$a} <=> $mr->{$b} } keys %$mr  ) {
        delete $mr->{$m} if ( $mr->{$m} > $count_keep_hi || $mr->{$m} < $count_keep_lo );           
    }

    my $n_motif = scalar keys %$mr;
    #print STDERR "bfr $m_motif aft $n_motif\n";
    #if ( $n_motif == 0 ) {
    #    print STDERR "lo $count_keep_lo hi $count_keep_hi  bfr $m_motif aft $n_motif\n";
    #}
    return ( $mr, $n_motif );
}
# End of RemoveRareAndFreqMotifs




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


sub MBSimilarity {

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
		    #my $mmm = ( $mx->{$m} ); 
		    #my $nnn = ( $my->{$m} );
		    #my $qqq = $mmm * $nnn;
		    #print STDERR "$mmm\t$nnn\t$qqq\n" unless ( $qqq );
		    $nmxy+= ( $mx->{$m} )* ( $my->{$m} );
		} else {
		    $nmxo+= ( $mx->{$m} )**2;
		}
            }
	    for my $m ( keys %$my ) {
		$nmy += ( $my->{$m} )**2;
		$nmyo+= ( $my->{$m} )**2 unless ( defined $mx->{$m} );
	    }
	    #print STDERR "$nmx $nmy nmxy $nmxy \n" unless ($nmxy);
	    #next unless ( $nmx && $nmy );
            my $similarity =
		#$nmxy/( $nmxo + $nmyo + $nmxy + 0.000000001 );
		#$nmxy /( $nmx + $nmy - $nmxy + 0.000000001 ); 
		0.5 * ($nmxy / ($nmx+0.00000000001) + $nmxy / ($nmy+0.00000000001) );
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of MarylandBridgeSimilarity



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
		    #my $mmm = ( $mx->{$m} ); 
		    #my $nnn = ( $my->{$m} );
		    #my $qqq = $mmm * $nnn;
		    #print STDERR "$mmm\t$nnn\t$qqq\n" unless ( $qqq );
		    $nmxy+= ( $mx->{$m} )* ( $my->{$m} );
		} else {
		    $nmxo+= ( $mx->{$m} )**2;
		}
            }
	    for my $m ( keys %$my ) {
		$nmy += ( $my->{$m} )**2;
		$nmyo+= ( $my->{$m} )**2 unless ( defined $mx->{$m} );
	    }
	    #print STDERR "$nmx $nmy nmxy $nmxy \n" unless ($nmxy);
	    #next unless ( $nmx && $nmy );
            my $similarity =
		#$nmxy/( $nmxo + $nmyo + $nmxy + 0.000000001 );
		$nmxy /( $nmx + $nmy - $nmxy + 0.000000001 ); 
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of TanimotoSimilarity


#--------------------------------------------------------------------------------
# WeightedTanimotoSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their distance and similarity: 
#
# S(X, Y) = |X && Y|/|X or Y|
#
# USAGE
#   ( $similarity ) = WeightedTanimotoSimilarity($xfp2motif, $weight_file); 
#--------------------------------------------------------------------------------
sub WeightedTanimotoSimilarity {

    my ($xfp2motif, $weight_file) = (@_);

    my %w;
    open ( my $in, "<$weight_file" ) || die "can not open $weight_file\n";
    while ( my $line = <$in> ) {
        chomp $line;
        my ( $motif, $weight ) = split " ", $line;
        $w{$motif} = $weight;
        #print STDERR "$motif $weight\n";
    }
    close $in;


    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
        my $rnax = $rna[$x];
        my $mx = $xfp2motif->{$rnax};
#        print STDERR "mx:$mx\n";
#        foreach my $id ( keys %$mx ) {
#            print STDERR "$id $$mx{$id}\n";
#        }

        for my $y ( $x ... $#rna ) {
	        my $rnay = $rna[$y];
            my $my = $xfp2motif->{$rnay};

            my ( $x_and_y, $x_or_y ) = (0, 0);
	    #my ( $nmx, $nmy, $nmxo, $nmyo, $nmxy ) = (0, 0, 0, 0, 0);
        #    for my $m ( keys %$mx ) {
		#$nmx += ( $mx->{$m} )**2;
        #        if ( defined $my->{$m} ) {
		#    $nmxy+= ( $mx->{$m} )* ( $my->{$m} );
		#} else {
		#    $nmxo+= ( $mx->{$m} )**2;
		#}
        #    }
	    #for my $m ( keys %$my ) {
		#$nmy += ( $my->{$m} )**2;
		#$nmyo+= ( $my->{$m} )**2 unless ( defined $mx->{$m} );
	    #}
            for my $m ( keys %$mx ) {
#                print STDERR "$m\n";
#                $x_or_y += 1;
                $x_or_y += $w{$m};
                if ( defined $my->{$m} ) {
                #    $x_and_y += 1;
                    $x_and_y += $w{$m};
                } 
            }
            for my $m ( keys %$my ) {
                unless ( defined $mx->{$m} ) {
                #    $x_or_y += 1;
                    $x_or_y += $w{$m};
                }
            }


            my $similarity = $x_and_y / $x_or_y; 
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of WeightedTanimotoSimilarity






sub TanimotoSimilarityOld {

    my ($xfp2motif) = @_;
    my %xfp2motif = %{$xfp2motif};

    my %similarity;
    my %distance;

    my $max_similarity = 0;

    for my $rnax ( sort keys %xfp2motif ) {
        for my $rnay ( sort keys %xfp2motif ) {

            my $motif_x = $xfp2motif{$rnax};
            my $motif_y = $xfp2motif{$rnay};

            my %motif_x = %{$motif_x};
            my %motif_y = %{$motif_y};

            my $denominator_x = scalar keys %motif_x;
            my $denominator_y = scalar keys %motif_y;

            my $numerator = 0;

            for my $motif ( keys %motif_x ) {
                $numerator++ if ( defined $motif_y{$motif} );
            }
            my $similarity =
              $numerator / ( $denominator_x + $denominator_y - $numerator );
            $similarity{$rnax}{$rnay} = $similarity;
            $max_similarity = $similarity if ( $max_similarity < $similarity );
        }
    }

    for my $rnax ( keys %similarity ) {
        for my $rnay ( keys %similarity ) {
            if ( $rnax eq $rnay ) {
                my $distance = 0;
                $distance{$rnax}{$rnay} = $distance;
            }
            else {
                my $distance =
                  -log( $similarity{$rnax}{$rnay} + 0.00000000001 ) / log(2);
                $distance{$rnax}{$rnay} = $distance;
            }
        }
    }

    return ( \%similarity, \%distance );
}




sub SquareSimilarity {

    my ($xfp2motif) = @_;
    my %xfp2motif = %{$xfp2motif};
    my @rna = sort keys %$xfp2motif;

    my %similarity;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {

	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

            my %motif_x = %{$xfp2motif{$rnax}};
            my %motif_y = %{$xfp2motif{$rnay}};

	    my $xsize = scalar keys %motif_x;
	    my $ysize = scalar keys %motif_y;
	    my $similarity = 0;
            my %sim;
            for my $motif ( keys %motif_x ) {
                if ( defined $motif_y{$motif} ) {
                    my ( $msize ) = split /\_/, $motif;
                    $sim{$msize}++;                    
                }
            }
            foreach my $msize ( keys %sim ) {
                $similarity += ( $sim{$msize} )**2;
            }
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}







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



sub ExtendedTanimotoSimilarity {

    my ($xfp2motif) = @_;
    my %xfp2motif = %{$xfp2motif};
    my @rna = sort keys %xfp2motif;

    my %similarity;
    my %distance;

    my $max_similarity = 0;

    for my $index_x ( 0 ... $#rna ) {

	my $rnax = $rna[$index_x];
        my $motif_x = $xfp2motif{$rnax};
        my %motif_x = %{$motif_x};
	my %extended_x;

	# get the extended set of x
        for my $mx ( keys %motif_x ) {
	my ( $mxsize, $mxindex ) = split /\_/, $mx;

	    $extended_x{$mx}++;
	    my @px = @{$$parent{$mx}};

	    foreach my $px ( @px ) {
		my ( $pxsize, $pxindex ) = split /\_/, $px;
		next unless ( $pxsize == $mxsize - 1 );
		$extended_x{$px}++;
	    }
	}
	next unless ( %extended_x );

        for my $index_y ( $index_x ... $#rna ) {

	    my $rnay = $rna[$index_y];
            my $motif_y = $xfp2motif{$rnay};
            my %motif_y = %{$motif_y};
	    my %extended_y;
	    my %extended_xy;

	    # get the extended set of y
	    for my $my ( keys %motif_y ) {
	    	my ( $mysize, $myindex ) = split /\_/, $my;

          	$extended_y{$my}++;
		$extended_xy{$my}++ if (defined $extended_x{$my});

		my @py = @{$$parent{$my}};

		foreach my $py ( @py ) {
		    my ( $pysize, $pyindex ) = split /\_/, $py;
		    next unless ( $pysize == $mysize - 1 );
		    $extended_y{$py}++;
		    $extended_xy{$py}++ if (defined $extended_x{$py});
		}
	    }

	    my $exy = scalar ( keys %extended_xy );
	    my $ex = scalar ( keys %extended_x );
	    my $ey = scalar ( keys %extended_y );
	    my $similarity = $exy / ( $ex + $ey - $exy );

            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
            #$max_similarity = $similarity if ( $max_similarity < $similarity );

	 

	    my $distance;
	    if ( $rnax eq $rnay ) {
		$distance = 0;
	    } else {
		if ( $similarity == 0 ) { 
		    $similarity = $similarity + 0.00000000000001;
	        }
		$distance = -log($similarity )/log(2);
		$distance = sprintf( "%.15f", $distance );
	    }
	    $distance{$rnax}{$rnay} = $distance;
	    $distance{$rnay}{$rnax} = $distance;

	    
	}

    }
    return ( \%similarity, \%distance );
}


#--------------------------------------------------------------------------------
# ParentTanimotoSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity: 
# 
# S(X,Y) = sum(direct parent similarity of all pairs of X motifs vs Y motifs) 
#
# USAGE
#   ( $similarity ) = ParentTanimotoSimilarity($xfp2motif); 
#--------------------------------------------------------------------------------
sub ParentTanimotoSimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
	my $rnax = $rna[$x];
	my %mx = %{$xfp2motif->{$rnax}};

        for my $y ( $x ... $#rna ) {

	    my $rnay = $rna[$y];
            my %my = %{$xfp2motif->{$rnay}};
	    my $similarity = 0;
	    
	    my ( $enum, $denom ) = ( 0, 0 );
 
	    for my $mx ( sort keys %mx ) {
		for my $my ( sort keys %my ) {
		    my %pxy;
		    my ( $npxy, $npx, $npy ) = ( 0, 0, 0 );
		    my @px = sort keys %{$dp{$mx}};
		    my @py = sort keys %{$dp{$my}};
		    next unless ( @px && @py );

		    for my $px ( @px ) {
			$pxy{$px}++; $npx++;
		    }   
		    for my $py ( @py ) {
			$pxy{$py}++; $npy++;
			$npxy++ if ( $pxy{$py} eq 2 );
		    }   
		    my $spxy = $npxy / ( $npx + $npy - $npxy ); 
		    $enum = $enum + $npxy;
		    $denom = $denom + ($npx + $npy - $npxy);
		}
	    }
	    $similarity = $enum / $denom;
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
}

# End of ParentTanimotoSimilarity



#--------------------------------------------------------------------------------
# ParentCountSimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity: 
# 
# S(X,Y) = sum(direct parent similarity of all pairs of X motifs vs Y motifs) 
#
# USAGE
#   ( $similarity ) = ParentTanimotoSimilarity($xfp2motif, $parent); 
#--------------------------------------------------------------------------------
sub ParentCountSimilarity {

    my ($xfp2motif) = @_;
    my @rna = sort keys %$xfp2motif;
    my %similarity;

    for my $x ( 0 ... $#rna ) {
        for my $y ( $x ... $#rna ) {
	    my $rnax = $rna[$x];
	    my $rnay = $rna[$y];

            my %mx = %{$xfp2motif->{$rnax}};
            my %my = %{$xfp2motif->{$rnay}};

	    my $similarity = 0;
	    my $numerator = 0;
	    my $denominator = 0;
            for my $mx ( keys %mx ) {
		for my $my ( keys %my ) {
		    my @pxy;
		    my @px = @{$$parent{$mx}};
		    my @py = @{$$parent{$my}};

		    next unless ( @px && @py );
		    my %pxy;
		    $pxy{$mx}++; $pxy{$my}++;
		    my $npxy = 0;
		    my $npx = scalar @px;
		    my $npy = scalar @py;

		    foreach my $px ( @px ) {
			$pxy{$px}++;
		    }
		    foreach my $py ( @py ) {
			$pxy{$py}++;
			$npxy++ if ( $pxy{$py} eq 2 );
		    }
		    $numerator = $numerator + $npxy;
		    $denominator = $denominator + $npx + $npy - $npxy;  
		}
            }
	    $similarity = $numerator / $denominator;
            $similarity{$rnax}{$rnay} = $similarity;
            $similarity{$rnay}{$rnax} = $similarity;
        }
    }

    return ( \%similarity );
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
            #    $similarity += $mx{$m} * $my{$m}  if ( defined $my{$m} );
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
# AntiProbabilitySimilarity
# 
# Given 2 fingerprints of RNA X & Y, calculate their similarity: 
# 
# S(X,Y) = |X, Y intersection| 
#
# USAGE
#   ( $similarity ) = AntiProbabilitySimilarity( $xfp2motif ); 
#--------------------------------------------------------------------------------
sub AntiProbabilitySimilarity {

    my ($xfp2motif) = @_;
    my %xfp2motif = %{$xfp2motif};

    my %similarity;
    my %distance;

    my $max_similarity = 0;

    for my $rnax ( keys %xfp2motif ) {
        for my $rnay ( keys %xfp2motif ) {

            my $motif_x = $xfp2motif{$rnax};
            my $motif_y = $xfp2motif{$rnay};

            my %motif_x = %{$motif_x};
            my %motif_y = %{$motif_y};

            my ( $totalcount_x, $totalcount_y ) = ( 0, 0 );
            for my $mx ( keys %motif_x ) {
                $totalcount_x += $motif_x{$mx};
            }

            for my $my ( keys %motif_y ) {
                $totalcount_y += $motif_y{$my};
            }

            my $similarity = 0;

            for my $motif ( keys %motif_x ) {
                next unless ( defined $motif_y{$motif} );
                my $probx = $motif_x{$motif}/$totalcount_x;
		my $proby = $motif_y{$motif}/$totalcount_y;
		$similarity += 
			log ($probx * $proby * ( 1 - $probx ) * ( 1 - $proby )+0.00000001);
            }
            $similarity{$rnax}{$rnay} = $similarity;
            $max_similarity = $similarity if ( $max_similarity < $similarity );
        }
    }

    for my $rnax ( keys %similarity ) {
        for my $rnay ( keys %similarity ) {
            if ( $rnax eq $rnay ) {
                my $distance = 0;
                $distance{$rnax}{$rnay} = $distance;
            }
            else {
                my $distance = $max_similarity - $similarity{$rnax}{$rnay};
                $distance{$rnax}{$rnay} = $distance;
            }

        }
    }

    return ( \%similarity, \%distance );
}
# End of AntiProbabilitySimilarity




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

            #my ($famx) = $rnax =~ /^([a-zA-Z0-9]*)\./;
            #my ($famy) = $rnay =~ /^([a-zA-Z0-9]*)\./;
	    my ( $famx ) = split /[\_,\.]+/, $rnax; 
	    my ( $famy ) = split /[\_,\.]+/, $rnay; 

            # TODO test on the undef similarity!! 
            my $similarity = $s{$rnax}{$rnay};
	    #print STDERR "undef s[$rnax][$rnay] = $s{$rnax}{$rnay}\n" unless ( $s{$rnax}{$rnay} );
	    #print STDERR "def s[$rnax][$rnay] = $s{$rnax}{$rnay}\n" if ( $s{$rnax}{$rnay} );

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

            #my ($famx) = $rnax =~ /^([a-zA-Z0-9]*)\./;
            #my ($famy) = $rnay =~ /^([a-zA-Z0-9]*)\./;
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

            for my $_ (@identical) {
                my $false_move = $false_percent * $false_size / $size;
                my $true_move  = $true_percent * $true_size / $size;
                push @{ $plot{$_} }, ( $false_move, $true_move );
            }

        }
        else {

            for my $_ (@identical) {
                if ( $total_idxed[$_]->[1] eq "false" ) {
                    push @{ $plot{$_} }, ( $false_percent, 0 );
                }
                elsif ( $total_idxed[$_]->[1] eq "true" ) {
                    push @{ $plot{$_} }, ( 0, $true_percent );
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

    # TODO tick marks (entend 3 pixels on either side of the axis)

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
# Revision 1.1.2.69  2016/08/24 05:23:31  huang147
# revised the input and output diretories.
#
# Revision 1.1.2.68  2016/08/24 04:07:28  huang147
# revised the input directory.
#
# Revision 1.1.2.67  2016/08/23 05:56:36  huang147
# removed unnecessary output files.
#
# Revision 1.1.2.66  2016/08/23 05:51:33  huang147
# adeed the output and input directory options.
#
# Revision 1.1.2.65  2016/03/08 15:28:30  huang147
# added WeightedTanimotoSimilarity.
#
# Revision 1.1.2.64  2015/12/04 03:40:04  huang147
# added SquareSimilarity, but it doesn't work very well. :(
#
# Revision 1.1.2.63  2014/10/15 15:01:25  huang147
# Added a similarity function MBSimilarity, using Maryland-Bridge coefficient.
#
# Revision 1.1.2.62  2014/09/05 15:38:29  huang147
# Fixed the RNA name problem:
# Now it works for the formats below:
# (1) fam_subfam.Genus_species.xios
# (2) fam.Genus_species.xios
#
# Revision 1.1.2.61  2014/05/06 18:45:28  huang147
# Added a subroutine SimilarityToMatrix.
#
# Revision 1.1.2.60  2014/04/29 03:36:48  huang147
# Modified so it fits the data with subclasses labeled.
#
# Revision 1.1.2.59  2013/12/13 03:57:40  huang147
# Fixed a bug in ExtendedTanimotoSimilarity calculation.
#
# Revision 1.1.2.58  2013/11/12 03:48:56  huang147
# Added DistanceToMega function.
#
# Revision 1.1.2.57  2013/11/11 22:23:48  huang147
# Changed the calculation of distances.
# TODO: should change the way of calculating mega files.
#
# Revision 1.1.2.56  2013/10/30 15:33:48  huang147
# Updated the calculation of VertexNumberSimilarity, EdgeNumberSimilarity, and MotifNumberSimilarity.
#
# Revision 1.1.2.55  2013/10/28 15:45:21  huang147
# Added graph stats calculation before similarity calculation. Prints out as an output.
#
# Revision 1.1.2.54  2013/10/07 05:00:23  huang147
# Added a different calculation for distance in building mega files.
#
# Revision 1.1.2.53  2013/10/03 02:48:49  huang147
# Addd extended tfidf calculation. Fixed the bug. Works!
#
# Revision 1.1.2.52  2013/10/02 16:45:18  huang147
# Fixed the artificial inflation of probability calculation in the previous versions.
#
# Revision 1.1.2.51  2013/10/01 21:01:35  huang147
# Added the functional to calculate the support fraction, support ratio, of each motif.
#
# Revision 1.1.2.50  2013/09/30 16:58:20  huang147
# Added a subroutine RemoveMotifs, which can remove motifs according to a user-defined string. Much easier to use than the other removal subroutines before.
#
# Revision 1.1.2.49  2013/09/26 05:01:30  huang147
# Fixed some bug.
#
# Revision 1.1.2.48  2013/09/26 02:34:23  huang147
# Added 3 functions to remove motifs according to their frequencies in the random sampling.
#
# Revision 1.1.2.47  2013/09/23 03:23:00  huang147
# Added the printout for mega 5 input ( .meg files ).
#
# Revision 1.1.2.46  2013/09/11 15:46:13  huang147
# Changed the output of AUC to be 4 digits after the decimal point.
#
# Revision 1.1.2.45  2013/09/11 14:35:53  huang147
# Added the prefix of similarity functions to the roc_coordinates files.
#
# Revision 1.1.2.44  2013/09/11 02:52:41  huang147
# Added the print-outs for TPR and FPR coordintes for ROC curve plotting in an outside software (like R).
#
# Revision 1.1.2.43  2013/07/31 02:32:19  huang147
# Changed the comments in the log part in AUC.txt.
#
# Revision 1.1.2.41  2013/07/17 04:03:27  huang147
# Removed subroutine ProbabilitySimilarity.
#
# Revision 1.1.2.40  2013/07/16 19:42:30  huang147
# Improved the speed of calculation by changing the similarity calculation.
# Removed some unnecessary comments.
#
# Revision 1.1.2.39  2013/07/16 18:26:26  huang147
# Added more printouts in the AUC.txt output.
#
# Revision 1.1.2.38  2013/07/16 16:18:33  huang147
# Added subroutine sortSimilarity to print output of similarities in sorted order.
#
# Revision 1.1.2.37  2013/07/15 02:03:55  huang147
# Changed the calculation of ParentSimilarity.
#
# Revision 1.1.2.35  2013/07/08 20:29:57  huang147
# Added some more printings to the output file.
#
# Revision 1.1.2.34  2013/07/05 15:57:33  huang147
# Added HammingSimilarity.
#
# Revision 1.1.2.33  2013/07/04 20:49:27  huang147
# The previous version is wrong!!!!
# This one is the one that added DiceSimilarity.
#
# Revision 1.1.2.32  2013/07/04 16:20:32  huang147
# Added Dice Similarity function.
#
# Revision 1.1.2.30  2013/06/04 18:16:55  huang147
# Fixed some bug in calculation.
# No functional change.
#
# Revision 1.1.2.28  2013/06/03 20:19:54  huang147
# Fixed some bug in calculation.
#
# Revision 1.1.2.26  2013/04/25 03:20:55  huang147
# Chnaged the code in ExtendedSimilarity so it calculates much faster.
#
# Revision 1.1.2.25  2013/04/25 01:46:21  huang147
# Added ExtendedSimiarity. Tested on Rossmann RCAC server. Works.
#
# Revision 1.1.2.23  2013/04/09 03:01:12  huang147
# Added ParentCountSimilarity. Works perfect in ROC curve.
#
# Revision 1.1.2.22  2013/04/08 16:36:42  huang147
# Added the ParentSimilarity, which only looks at the common parents of each motif pair. (No longer look at the motif itself!!)
#
# Revision 1.1.2.21  2013/04/08 03:57:52  huang147
# Changed the print out of output file AUC.txt
#
# Revision 1.1.2.20  2013/04/08 03:42:15  huang147
# Added extended fingerprint and corresponding distance function. Works on rossmann.rcac cluster.
#
# Revision 1.1.2.19  2013/04/03 01:16:19  huang147
# Used perltidy to make the code clean.
#
# Revision 1.1.2.18  2013/04/02 20:12:04  huang147
# This is a working version.
# Removed some unnecessary comments.
#
# Revision 1.1.2.17  2013/04/02 19:47:00  huang147
# Put in the calculation of AUC. But need to make the code tidy.
#
# Revision 1.1.2.16  2013/04/02 14:20:36  huang147
# Added a working ROC curve.
# Still need to add the calculatioin of AUC.
#
# Revision 1.1.2.15  2013/03/30 04:10:38  huang147
# Trying to fix the ROC plot bug in identical distance. Not done yet. Not debugged yet.
#
# Revision 1.1.2.14  2013/03/27 02:19:49  huang147
# Fixed the ranking of distance/similarity pairs in sub ROC.
#
# Revision 1.1.2.13  2013/03/26 02:21:38  huang147
# Changed the output format.
# Added a way to output all the similarity functions together.
#
# Revision 1.1.2.12  2013/03/26 02:01:33  huang147
# Added the output of AUC.txt.
#
# Revision 1.1.2.11  2013/03/26 01:13:47  huang147
# Changed some printout formats. No function change.
#
# Revision 1.1.2.10  2013/03/26 01:11:07  huang147
# Took out some unnecessary comments.
#
# Revision 1.1.2.9  2013/03/26 00:46:09  huang147
# Fixed the bug in reading motifs from XML format.
#
# Revision 1.1.2.8  2013/03/25 16:34:51  huang147
# Added the choice of different similarity functions as getopts options.
#
# Revision 1.1.2.7  2013/03/22 16:40:55  huang147
# Removed some comments.
#
# Revision 1.1.2.6  2013/03/22 16:10:16  huang147
# Added more than one ROC curve in a single plot.
# -------------------------------------------------------------------
#
# Revision 1.1.2.4  2013/03/22 02:28:23  huang147
# Added the calculation of ROC according to distance.
#
# Revision 1.1.2.3  2013/03/22 02:02:03  huang147
# Changed the subroutines of similarity4ROC and ROC.
# Removed unnecessary codes and comments.
#
# Revision 1.1.2.2  2013/03/21 22:08:44  huang147
# Just added ROC curve. Not debugged yet!!
#
# Revision 1.1.2.1  2013/03/21 21:24:11  huang147
# Initial version. Use new Fingerprint.pm to read in .xfp files.
#
#
#------------------------------------------------------------------------------

