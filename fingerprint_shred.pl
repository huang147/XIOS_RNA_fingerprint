#!/apps/group/bioinformatics/apps/perl-5.16.1/bin/perl
#PBS -N fingerprint
#PBS -q mgribsko
#PBS -l nodes=1:ppn=16
################################################################################
# $Id: fingerprint_shred.pl,v 1.1.4.9 2016/08/24 05:24:52 huang147 Exp $
#
# Run many random fingerprint jobs.
# Requires
#   fingerprint_random.pl
#   DfsGraph.pm
#   Dfsgenerator.pm
#   Motif_db.pm
################################################################################
use strict;
#use lib '/home/huang147/perl_src_huang_120827';
use lib './';
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Std;
use Time::HiRes qw( gettimeofday tv_interval );
use Cwd qw( abs_path );
use Fingerprint;

my $DEFAULT_EXEDIR        = $ENV{'PBS_O_WORKDIR'};
my $DEFAULT_INPUTDIR      = $ENV{'PBS_O_WORKDIR'};
my $DEFAULT_NJOBS         = 8;
my $DEFAULT_OUTPUT        = 'fingerprint';
my $DEFAULT_WORKDIR       = $ENV{'PBS_O_WORKDIR'};

my $DEFAULT_MOTIFDB       = 'Motif_fingerprint/2_to_7_stems_topologies.removed_not_true.mini_dfs.txt.removed_redundant.with_label.motif.storable';
my $DEFAULT_SUBGRAPH_SIZE = 7;
my $DEFAULT_LIMIT         = 10;
my $DEFAULT_MAPPING       = 1;
my $DEFAULT_N             = 100000;
my $DEFAULT_VERBOSE       = 0;
my $DEFAULT_WALLTIME_HOUR = 4;


my $id = '$Id: fingerprint_shred.pl,v 1.1.4.9 2016/08/24 05:24:52 huang147 Exp $';
my ( $program, $version ) = $id =~ /Id: ([^,]+),v ([^ ]+)/;

my $USAGE = qq{v$version $program  
    -h            this usage message
    -e <dir>      directory for executable code (default=$DEFAULT_EXEDIR)
    -i <dir>      input directory for xios files (default=$DEFAULT_INPUTDIR)
    -j <int>      number of executables to run (default=$DEFAULT_NJOBS)
    -w <dir>      working directory (default=$DEFAULT_WORKDIR)
    
    fingerprint_random.pl parameters
    -d <motif_db> motif database file (default=$DEFAULT_MOTIFDB)
    -g <int>      graph size to sample (default=$DEFAULT_SUBGRAPH_SIZE)
    -l <int>      minimum count required for each motif (default=number of random samples)
    -m <int>      minimum mapping required for each motif (default=1)
    -n <int>      maximum number of random samples (default=$DEFAULT_N)
    -v            show DFS code for each motif (verbose, default=$DEFAULT_VERBOSE)
    -t            walltime for queue on cluster (in hours, default=$DEFAULT_WALLTIME_HOUR) 
    -x            create decoy directory for holding sampled subgraphs as .xios files  
};

# command line options

my $option = {};
getopts( 'd:e:g:hi:j:l:m:n:vw:t:x', $option );

# help
if ( $$option{h} ) {
    print "$USAGE\n";
    exit 1;
}

# parallel job parameters
my $exedir   = $DEFAULT_EXEDIR;
my $inputdir = $DEFAULT_INPUTDIR;
my $njobs    = $DEFAULT_NJOBS;
my $workdir  = $DEFAULT_WORKDIR;
my $decoydir;

if ( $$option{e} ) { $exedir   = $$option{e} }  # directory for executable files
if ( $$option{i} ) { $inputdir = $$option{i} }  # input directory for xios files
if ( $$option{j} ) { $njobs    = $$option{j} }  # number of concurrent jobs
if ( $$option{w} ) { $workdir  = $$option{w} }  # working directory

$exedir 	= abs_path( $exedir ).'/';
$inputdir   	= abs_path( $inputdir ).'/';
unless (-e $workdir) {
    print STDERR "$workdir doesn't exist, creating...";
    system("mkdir -p $workdir");
}
$workdir 	= abs_path( $workdir ).'/';

# fingerprint_random.pl parameters
my $motif_db        = $DEFAULT_MOTIFDB;
my $subgraph_size   = $DEFAULT_SUBGRAPH_SIZE;
my $motif_limit     = $DEFAULT_LIMIT;
my $mapping_limit   = $DEFAULT_MAPPING;
my $motif_iteration = $DEFAULT_N;
my $verbose         = $DEFAULT_VERBOSE;
my $hangtime        = 3600 * $DEFAULT_WALLTIME_HOUR - 1000;

if ( $$option{d} ) { $motif_db        = $$option{d} }
if ( $$option{g} ) { $subgraph_size   = $$option{g} }
if ( $$option{l} ) { $motif_limit     = $$option{l} }
if ( $$option{m} ) { $mapping_limit   = $$option{m}; $motif_limit = $DEFAULT_N; }
if ( $$option{n} ) { $motif_iteration = $$option{n} }
if ( $$option{t} ) { $hangtime        = 3600 * $$option{t} - 1000 }
if ( $$option{x} ) {$decoydir = $$option{x} }
# verbose is processed in child command section

my $hostname = `hostname`;
chomp $hostname;

#-------------------------------------------------------------------------------
# set up parallel manager and callbacks
#-------------------------------------------------------------------------------
my %query;
my $manager = new Parallel::ForkManager( $njobs );

$manager->run_on_finish( sub { 
        my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_;
        my ( $jobnum, $jobname ) = split "_", $ident, 2;
        postProcess( $jobnum, $jobname );
        
        my $job = $query{$jobname};
        my $status = "finished";
        my $fpt = $$job{fingerprint};
        my $min = 0;
        if ( defined $fpt ) { $min = $fpt->motifMin; } # Minimum count for all motifs
        
        if ( $$job{done} ) {  $status = "done";} 
    
        printf "%24s  %6d  %6d  %8s  %-30s  %d\n", 
               scalar localtime(), $pid, $jobnum, $status, $jobname,$min;
    }
);

$manager->run_on_start( sub { 
    my ( $pid, $ident ) = @_;
    my ( $jobnum, $jobname ) = split "_", $ident, 2;
    my $job = $query{$jobname};
    $$job{count}++;
    my $fpt = $job->{fingerprint};
    my $iteration = 0;
    if ( defined $fpt ) {
        $iteration = $fpt->iteration || 0;
    }
    printf "%24s  %6d  %6d  started   %-30s  %d\n", 
           scalar localtime(), $pid, $jobnum, $jobname, $iteration;        
    }
);

$manager->run_on_wait(
    sub {
#        print "Waiting for children\n";
#        print "Finish time: ".localtime(time)."\n";
    },
);

#-------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------
print "$program $version\n";
print "    Host: $hostname\n";
print "    Executable directory: $exedir\n";
print "    Input directory: $inputdir\n";
print "    Working directory: $workdir\n";
print "    Concurrent jobs: $njobs\n";
print "    Start time: ".localtime(time)."\n";

print "    Motif database: $motif_db\n";
print "    Subgraph size: $subgraph_size\n";
print "    Motif iterations: $motif_iteration\n";
print "    Motif limit: $motif_limit\n\n";

# read in a list of xios files.  for each one create an entry in the query hash
# to store information about this structure 

my $t0 = [gettimeofday];

opendir( my $dh, $inputdir ) || die "unable to open input directory ($inputdir)\n\n";

my $nquery = 0;
my $nfpt = 0;
while( my $file = readdir $dh) {
    if ( $file =~ /^\./ ) { next; }
    chomp $file;
    
    if ( $file =~ /\.xios$/ ) {
	my $file_done = $file.'.done';
	next if ( -e $file_done );
	print "    Reading xios file $file\n";
        $nquery++;
    	my $original_fingerprint_file = $workdir.$file.".xpt";
    	my $original_fpt = Fingerprint->new;
	# read the fingerprint from original fingerprint file if it exists
	if ( -e $original_fingerprint_file && -s $original_fingerprint_file ) {
	    print "    Reading fingerprints from the original file $original_fingerprint_file\n";
    	    $original_fpt->readXML( $original_fingerprint_file );
	    $nfpt++;
	}  else {
	    $original_fpt = undef;
	}
        $query{$file} = {  jobname     => $file,
                           infile      => $inputdir.$file,
                           summaryfile => $workdir.$file.".summary",
                           count       => 0,
                           done        => 0,
#                           iteration   => 0,
#                           time        => 0,
#                           ve          => "",
#                           motiflist   => {},
#                           nmotif      => 0,
#                           min         => 0,
                           fingerprint => $original_fpt,
        };
    }
}
closedir $dh;
print "\n    $nquery query xios files read from $inputdir\n";
print "    $nfpt original fingerprint files read from $inputdir\n\n"; 

# run commands in parallel

my $jobnum = 0;
while ( my $job = getNextQuery(\%query,$manager) ) {
    my $t1 = [gettimeofday];
    my $elapsed = tv_interval ( $t0, $t1 );
    print STDERR "The program has elapsed $elapsed seconds\n";
    if ( $elapsed > $hangtime ) {
	print STDERR "Exceed the hanging time $hangtime seconds! Stop the program\n";
	last;
    }
    last if ( $job eq 'finished');

    $jobnum++;

    # Fork and return the pid for the child process   
    my $ident = $jobnum."_$$job{jobname}";
    my $pid = $manager->start( $ident ) and next; 

    #---------------------------------------------------------------------------
    # begin child command
    #---------------------------------------------------------------------------

    my $outfile = $workdir.$$job{jobname}.".fingerprint_$jobnum";
    
    my $command = "";
    $command .= "$exedir"."fingerprint_random.pl -q";
    $command .= " -d $motif_db";
    $command .= " -g $subgraph_size";
    $command .= " -n $motif_iteration";
    $command .= " -w $workdir";         
    $command .= " -x" if ( $decoydir );
    if ( $$option{v} ) { $command .= " -v"; }
    # turned limit check on individual jobs off, doesn't save much and makes
    # short graphs terminate very rapidly  (possibly too rapidly)
    #if ( $$option{l} ) { $command .= " -l $motif_limit" }
    $command .= " $$job{infile}";
    $command .= " > $outfile";
   
    #my $xios_done = $workdir.$$job{jobname}.".done"; 
    exec( $command );    

    #---------------------------------------------------------------------------
    # end child command
    #---------------------------------------------------------------------------
  
}

$manager->wait_all_children;

#print out complete results for all queries
#foreach my $r ( sort keys %query ) {
    
#    my $fingerprint_file = $workdir.$r.".xpt";
#    my $fingerprint = $query{$r}->{fingerprint};
#    $fingerprint->asXML( $fingerprint_file );
#}

exit 0;

#-------------------------------------------------------------------------------
# getNextQuery
#
# examine the queries and fine a query that is 1) not currently running and 
# 2) not done.
#
# if the total number of queries available is less than $njobs, reduce $njobs
#
# usage
#   while ( my $job = getNextQuery(\%query,$manager) ) { ...
#-------------------------------------------------------------------------------
sub getNextQuery{
    my ( $query, $manager ) = @_;

    my $target;
    my $nquery = 0;
    my $navail = 0;
    foreach my $file ( keys %$query ) {
        my $job = $query{$file};
        next if $$job{done};        # job is done
        $nquery++;
        $navail++;
        #the target job is the one that has been run the fewest times
        if ( $target ) {
            if ( $$target{count} > $$job{count} ) {
                $target = $job;
            }
        } else {
            $target = $job;
        }
    }
    
    if ( $nquery == 0 ) { return 'finished'; }

    # this section decreases the number of jobs when there are fewer queries 
    # than the number of jobs.  this ensures only one of each kind of query
    # runs at the same time
    
    # very bad, reaching into the manager object here, but there is no
    # accessor function provided
#    my $njobs = $manager->{max_proc};
#    if ( $nquery < $njobs ) {
#        my $old_njobs = $njobs;
#        $manager->set_max_procs( $nquery );
#        print STDERR "setting njobs to $nquery, previous value was $old_njobs\n";
#        print "setting njobs to $nquery, previous value was $old_njobs\n";
#    }

    return $target;
}

# end of getNextQuery

#-------------------------------------------------------------------------------
# postProcess
#
# Read the output file and extract the information about timing and motifs
# found.
#
# uses global variable %query, $motif_limit
# usage
#   postProcess( $info );
#-------------------------------------------------------------------------------
sub postProcess {
    my ( $jobnum, $jobname ) = @_;
        
    my $outfile = $workdir.$jobname.".fingerprint_$jobnum";
    my $partial_fpt = Fingerprint->new;
    $partial_fpt->readXML( $outfile );    
    
    my $parent_fpt;
    my $info = $query{$jobname};
    if ( defined $$info{fingerprint} ) {
           $parent_fpt = $$info{fingerprint};  
           $parent_fpt->sum( $partial_fpt );      
    } else {
        $parent_fpt = $partial_fpt;
        $$info{fingerprint} = $partial_fpt;
    }
        
    # print out the fingerprint file as soon as one job finishes 
    # in case the program crashes, we still have the intermediate files
    my $fingerprint_file = $workdir.$jobname.".xpt";
    $parent_fpt->asXML( $fingerprint_file );
    
    my $xios_file = $inputdir.$jobname;
    my $xios_done = $inputdir.$jobname.".done";


    if ( $parent_fpt->motifMin >= $motif_limit ) { 
	$$info{done} = 1; 
	unless ( -e $xios_done ) {
	    print STDERR "\nQuery $jobname has reached the minimun count of motifs: $motif_limit\n";
	    print STDERR "cp $xios_file $xios_done\n\n";
	    system( "cp $xios_file $xios_done" );
	}
    } 



    open my $summary, ">>", $$info{summaryfile} or warn "unable to open $$info{summaryfile}\n";
    printf $summary "%30s %6d %8d %9.2f %6d %6d\n",
           $jobname, $jobnum, $parent_fpt->iteration, $parent_fpt->timeElapsed, 
           $parent_fpt->motifN, $parent_fpt->motifMin;
    close $summary;
    `rm $outfile`;

    return;
}

# end of postProcess

################################################################################
# $Log: fingerprint_shred.pl,v $
# Revision 1.1.4.9  2016/08/24 05:24:52  huang147
# revised the input and output directories.
#
# Revision 1.1.4.8  2015/10/12 19:38:46  huang147
# added the option for inputting decoy directory into fingerprint_random.pl
#
# Revision 1.1.4.6  2013/08/17 04:24:44  huang147
# Revised the way of writing out xios.done files: keep both .xios and .xios.done files.
# Fixed the bug that the code will die if one of multiple jobs for the same file is done.
# TODO: fingerprint_random.pl doesn't run on graphs with a size smaller than the expected size.
# Should add a check to make sure the code doesn't get killed by this reason.
#
# Revision 1.1.4.5  2013/07/02 22:44:11  huang147
# Added time out option.
# Rename the xios files that reach the running criteria so they won't be bothered to check into when starting new jobs.
# Added some comments and printing-outs.
#
# Revision 1.1.4.3  2013/06/26 19:13:34  huang147
# Changed from Michael's original version.
# Added the reading of original fingerpirnt files if they exists.
# Also added writing of fingerprint files as soon as a job finishes in case the program crashes we still have the intermediate files and can resume at that point.
#
# Revision 1.1.2.8  2013/03/14 18:09:50  gribskov
# Debugged.
#
# Revision 1.1.2.7  2013/03/13 20:51:02  gribskov
# Started converting to use Fingerprint.pm.
#
# Revision 1.1.2.6  2013/03/10 12:17:46  gribskov
# Tried to resolve TODOs.
# Decided default file is OK; PBS_O_WORKDIR is root dir outside of PBS.
# Added vertex/edge numbers to data structure ($$data{ve}) and to report.
# Added individual files for fingerprint.
#
# Revision 1.1.2.5  2013/03/10 11:58:46  gribskov
# Added exedir for executables.
# Switched to using number of jobs as the cri9terion for which job to run next in getNextQuery.
# Turned on limit checking (-l).
#
# Revision 1.1.2.4  2013/03/09 10:28:54  gribskov
# Finally made this work.  need to clean it up.  some issues:
# 1. workdir doesn't work quite right
# 2. too many digits in start output sometimes
# 3. job choosing strategy needs improvement
#
# Revision 1.1.2.2  2013/03/03 21:16:22  gribskov
# Improved testing for completion, hopefully eliminating race condition
#
# Revision 1.1.2.1  2013/03/03 20:31:02  gribskov
# Initial version.  Runs but does not properly enforce end condition due
# to race situation between job loop and completion test.
#
################################################################################

__END__
