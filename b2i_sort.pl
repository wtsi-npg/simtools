#!/software/bin/perl
#
# b2i_sort.pl
#
# $Id: b2i_sort.pl 1281 2010-08-27 14:25:53Z mg10 $
#
# By default, b2i sorts its output files into ascending order of probe position. It does
# this by using system commands. Sometimes there does not appear to be a system processor
# available, so this final sorting does not take place.
#
# This script is really just a lift of the b2i system() commands converted into Perl
# and given a command-line wrapper. It can be used to sort the b2i intensity output
# files into ascending probe position, assuming b2i wase otherwise successful and said,
# "Sorry, unable to order file. Aside from that, b2i did what it was supposed to."
#

use strict;
use warnings;
use File::Temp ();
use Getopt::Long;

$| = 1; # Flush buffer.

my ( $unordered, $ordered, $sort );

#$unordered = "19.unordered.txt";
#$ordered = "temp.txt";

&process_command_line_arguments();
&order_output_by_pos();

print "\nDone.\n\n";



###################################################
#
# sub create_temp_file()
#
# Create a File::Temp object.
# Return this object and its filename.
#
###################################################

sub create_temp_file {

  $File::Temp::KEEP_ALL = 1; # Keep the file created.
  my $tempfile = new File::Temp;
  $tempfile->autoflush(1);
  my $filename = $tempfile->filename;

  return ( $tempfile, $filename );

}



###################################################
#
# sub order_output_by_pos()
#
# Use various system commands to re-order the file
# into ascending probe co-ordinate order.
#
###################################################

sub order_output_by_pos {

  my ( $tmp, $tmp_fn ) = &create_temp_file();

  if ( system("sort -c -k 2,2 -n -s $unordered > /dev/null 2>&1") == 0 ) {
    print STDERR "\nInput appears to be already sorted!\n";
	if ( $sort ) { print "\nSorting regardless...\n";}
    else { exit (1); }
  };


  #
  # Copy first line of the file (lists sampleIDs):
  #
  my $cmd = "head -1 $unordered > $tmp_fn";
  print "\nExecuting cmd: $cmd\n";
  #system ($cmd) or die "\nCommand 1 ($cmd) failed : $! ($?)\n";
  qx /$cmd/;
  #
  # Now copy the rest of the file, sorting it by SNP position (col. 2) in the process.
  #
  $cmd = "tail -n +2 $unordered | sort -g --key=2 >> $tmp_fn";
  print "\nExecuting cmd: $cmd\n";
  #system($cmd) or die "\nCommand 2 ($cmd) failed : $! ($?)\n";
  qx /$cmd/;

  #
  # Finally, rename the temp file to the output file required.
  # (This may mean copying the new ordered file over the original (unsorted) file)
  #
  $cmd ="mv $tmp_fn $ordered";
  print "\nExecuting cmd: $cmd\n";
  #system($cmd) or die "\nCommand 3 ($cmd) failed: $! ($?)\n";
  qx /$cmd/;

}



##################################################
#
# sub process_command_line_arguments()
#
##################################################

sub process_command_line_arguments {

  my $force;

  GetOptions (
             'unordered|input=s' => \$unordered, # -u. Input (unordered) file.
             'ordered|output=s'  => \$ordered,   # -o. Output (ordered) file.
             'force'             => \$force,     # -f. Force overwrite
             'sort'              => \$sort,      # -s. Force sort.
             'help'              => sub {&usage()},

             );

  # print "\nDBG: $unordered, $ordered\n";

  if ( ! $unordered ) {
	print "\nPlease supply a file to process.\n\n";
    exit (1);
   }

  if ( ! ( -e $unordered ) ) {
    print "\nFile $unordered does not appear to exist.\n\n";
    exit (1);
  }

  if ( ! $ordered ) {
	if ( $force ) {
      print "\nI shall overwrite the original (unordered) file with the nice new one.\n";
	  $ordered = $unordered;
	}
	else {
	  print "\nPlease either supply an output filename, or use the -f option.\n\n";
	  exit (1);
	}
  }

}



##################################################
#
# sub usage()
#
##################################################

sub usage {

 my $explain = <<END;

 b2i_sort.pl : order "failed" b2i output files into ascending probe position order.

 Usage:
 perl b2i_sort.pl -input <unordered_file> -output <new_ordered_file> [-sort]
 perl b2i_sort.pl -input <unordered_file> -force [-sort]
 perl b2i_sort.pl -help

 where:
 -input,     -i : These switches are synonymous. Use one of these to
 -unordered, -u : specify the name of the file which is to be ordered.

 -output,    -o : These are also synonymous. Use one of these to
 -ordered,   -o : specify the name of the new, ordered file.

 -force,     -f : Use this to force the new file to overwrite the old one.
                : If you wish to do this then do not specify a new file name.

 -sort,      -s : The script will check the input file is not already sorted
                : before processing it. Use this option to sort regardless.

 -help,      -h : Display this help message.

END
  print $explain;
  exit();

} # End usage()

