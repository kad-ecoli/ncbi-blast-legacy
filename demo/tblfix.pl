#!/usr/bin/perl
use strict;
use warnings;

# tblfix.pl


#  NCBI is supplying this Perl script for the convenience of GenBank
#  submitters who need to edit 5-column feature tables prior to processing
#  with tbl2asn. The tblfix.pl script can perform various data conversions on
#  feature table files, with the conversion function specified by a command
#  line argument. Additional conversions can be supported by modification of
#  the script. Any suggestions for improvements or corrections are welcome.
#  However, given the number of possible variations of local data files and
#  conditions, NCBI cannot promise to implement them all and cannot be
#  expected to troubleshoot local installations. This script is being
#  supplied as is, but you are always welcome to make any changes or
#  improvements yourself.
#
#  For example, the script can create a codon_recognized qualifier from
#  either the tRNA product or note qualifier, and it can reverse complement
#  the triplet codon for the codon_recognized value. Since the script reads
#  from stdin and writes to stdout, individual conversion steps can be piped
#  together to accomplish multiple operations on a single command line. This
#  avoids requiring submitters to change the storage conventions of their
#  underlying data in order to prepare their submission.
#
#  Existing conversion arguments are described below:
#
#  -product_to_codon_recognized
#
#    Parses the product qualifier in a tRNA feature looking for a value that
#    fits the pattern "tRNA-Met-CAT". It will create separate product
#    "tRNA-Met" and codon_recognized "CAT" qualifiers.
#
#  -note_to_codon_recognized
#
#    Looks for a tRNA note qualifier that fits the pattern "codon recognized:
#    CAT". It will create a codon_recognized "CAT" qualifier and remove the
#    note.
#
#  -revcomp_codon_recognized
#
#    This will adjust the codon recognized if the original data provided the
#    reverse complement sequence. For codon_recognized "CAT" it will result
#    in a codon_recognized "ATG" qualifier, the correct codon for methionine.
#
#  -product_to_ec_number
#
#    This function parses an enzyme commission number at the end of a
#    /product qualifier in the form (EC 1.2.3.4). It will create a product
#    qualifier without the EC number and an EC_number qualifier of form
#    "1.2.3.4".
#
#  -remove_predicted_evidence
#
#    If the value of an evidence qualifier is the uninformative value
#    "predicted", it will be removed.
#
#  -clean_rrna_product
#
#    This moves secondary rRNA product names to a note qualifier.
#
#  -create_locus_tag
#
#    This takes a locus_tag prefix as an additional command-line argument
#    and creates gene features with locus_tag of the form PFX_0001, assigning
#    sequential numbers, on CDS, tRNA, rRNA, ncRNA, and tmRNA features.
#
#  -create_protein_id
#
#    This takes a general database identifier as an additional command-line
#    argument and creates a protein_id qualifier of the form gnl|DBX|PFX_0001,
#    using the existing locus_tag qualifier on a CDS feature.
#
#  The regular expressions used to detect these patterns can easily be
#  changed to handle different variants.


# Script to edit contents of 5-column feature table files.
# Uses STDIN and STDOUT so multiple instances can be piped.

# define global variables

my $line_number = 0;
my $locus_tag_prefix = "";
my $protein_id_prefix = "";
my $locus_tag_count = 0;

my $flag = shift or die "Must indicate desired command (e.g., -revcomp_codon_recognized)\n";

if ($flag =~ /^-create_locus_tag$/) {
  $locus_tag_prefix = shift or die "Must supply locus_tag prefix\n";
}

if ($flag =~ /^-create_protein_id$/) {
  $protein_id_prefix = shift or die "Must supply protein_id prefix\n";
}

# state variables for tracking current position in table file

my $thisline;
my $in_key;
my $in_qual;
my $current_key;
my $current_loc;
my $current_qual;
my $current_val;
my @master_loc;

sub clearflags {
  $thisline = "";
  $in_key = 0;
  $in_qual = 0;
  $current_key = "";
  $current_loc = "";
  $current_qual = "";
  $current_val = "";
  @master_loc = ();
}

# subroutine for printing key and location

sub printkey {
  my $thiskey = shift (@_);

  my $is_first = 1;
  foreach my $thisloc (@master_loc) {
    if ($is_first == 1) {
      print STDOUT "$thisloc\t$thiskey\n";
    } else {
      print STDOUT "$thisloc\n";
    }
    $is_first = 0;
  }
}

# subroutine for filtering by feature key

sub printfeat {
  my $thiskey = shift (@_);

  if ($flag =~ /^-create_locus_tag$/) {
    if ($thiskey eq "CDS" ||
        $thiskey eq "rRNA" || $thiskey eq "tRNA" ||
        $thiskey eq "ncRNA" || $thiskey eq "tmRNA") {
      printkey ("gene");
      my $locus_tag_string = "0000";
      $locus_tag_count++;
      $locus_tag_string = "0000";
      if ($locus_tag_count > 999) {
        $locus_tag_string = $locus_tag_prefix . "_" . $locus_tag_count;
      } elsif ($locus_tag_count > 99) {
        $locus_tag_string = $locus_tag_prefix . "_0" . $locus_tag_count;
      } elsif ($locus_tag_count > 9) {
        $locus_tag_string = $locus_tag_prefix . "_00" . $locus_tag_count;
      } else {
        $locus_tag_string = $locus_tag_prefix . "_000" . $locus_tag_count;
      }
      print STDOUT "\t\t\tlocus_tag\t$locus_tag_string\n";
      printkey ($thiskey);
      print STDOUT "\t\t\tlocus_tag\t$locus_tag_string\n";
      return;
    }
  }

  printkey ($thiskey);
}

# subroutine for filtering by qualifier and value

sub printqual {
  my $thisqual = shift (@_);
  my $thisval = shift (@_);
  my $thiskey = shift (@_);

  if ($flag =~ /^-product_to_codon_recognized$/) {
    if ($thiskey eq "tRNA" && $thisqual eq "product") {
      if ($thisval =~ /^tRNA-([A-Za-z]+)-([A-Za-z]{3})$/ || $thisval =~ /^([A-Za-z]+)-([A-Za-z]{3})$/) {
        print STDOUT "\t\t\tproduct\ttRNA-$1\n";
        print STDOUT "\t\t\tcodon_recognized\t$2\n";
        return;
      } elsif ($thisval =~ /^tRNA-Undet-\?\?\?$/) {
        print STDOUT "\t\t\tproduct\ttRNA-OTHER\n";
        return;
      }
    }
  }

  if ($flag =~ /^-revcomp_codon_recognized$/) {
    if ($thiskey eq "tRNA" && $thisqual eq "codon_recognized") {
      if ($thisval =~ /^([A-Za-z]{3})$/) {
        my $revcomp = reverse $1;
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        print STDOUT "\t\t\tcodon_recognized\t$revcomp\n";
        return;
      }
    }
  }

  if ($flag =~ /^-note_to_codon_recognized$/) {
    if ($thiskey eq "tRNA" && $thisqual eq "note") {
      if ($thisval =~ /^codon recognized: ([A-Za-z]{3})$/) {
        print STDOUT "\t\t\tcodon_recognized\t$1\n";
        return;
      }
    }
  }

  if ($flag =~ /^-product_to_ec_number$/) {
    if ($thiskey eq "CDS" && $thisqual eq "product") {
      if ($thisval =~ /^(.+) \(EC (.+)\)$/) {
        print STDOUT "\t\t\tproduct\t$1\n";
        print STDOUT "\t\t\tEC_number\t$2\n";
        return;
      }
    }
  }

  if ($flag =~ /^-clean_rrna_product$/) {
    if ($thiskey eq "rRNA" && $thisqual eq "product") {
      if ($thisval =~ /^([^;]+)\s*;\s*(.+)$/) {
        print STDOUT "\t\t\tproduct\t$1\n";
        print STDOUT "\t\t\tnote\t$2\n";
        return;
      }
    }
  }

  if ($flag =~ /^-remove_predicted_evidence$/) {
    if ($thisqual eq "evidence" && $thisval eq "predicted") {
      return;
    }
  }

  if ($flag =~ /^-create_protein_id$/) {
    if ($thiskey eq "CDS" && $thisqual eq "locus_tag") {
      print STDOUT "\t\t\tprotein_id\tgnl|$protein_id_prefix|$thisval\n";
      print STDOUT "\t\t\t$thisqual\t$thisval\n";
      return;
    }
  }

  print STDOUT "\t\t\t$thisqual\t$thisval\n";
}

#subroutine to print next feature key / location / qualifier line

sub flushline {
  if ($in_key == 1) {
    # call filter subroutine for feature keys
    printfeat ($current_key);
  } elsif ($in_qual == 1) {
    if ($current_val eq "") {
      print STDOUT "\t\t\t$current_qual\n";
    } else {
      # call filter subroutine for qualifiers
      printqual ($current_qual, $current_val, $current_key);
    }
  }
}

# initialize flags at start of program

clearflags ();

# main loop reads one line at a time

while ($thisline = <STDIN>) {
  $thisline =~ s/\r//;
  $thisline =~ s/\n//;
  $line_number++;

  if ($thisline =~ /^>Feature.+$/) {
    # new table
    flushline ();
    $in_key = 0;
    $in_qual = 0;
    $current_key = "";
    $current_loc = "";
    $current_qual = "";
    $current_val = "";
    @master_loc = ();
    print STDOUT "$thisline\n";

  } elsif ($thisline =~ /^(\d+)\t(\d+)\t(\S+)$/) {
    # new feature
    flushline ();
    $in_key = 1;
    $in_qual = 0;
    $current_key = $3;
    $current_loc = "$1\t$2";
    @master_loc = ();
    push (@master_loc, $current_loc);

  } elsif ($thisline =~ /^(\d+)\t(\d+)$/) {
    # multiple intervals
    $in_key = 1;
    $in_qual = 0;
    $current_loc = "$1\t$2";
    push (@master_loc, $current_loc);

  } elsif ($thisline =~ /^\t\t\t(\S+)\s*$/) {
    # singleton qualifier
    flushline ();
    $in_key = 0;
    $in_qual = 1;
    $current_qual = $1;
    $current_val = "";

  } elsif ($thisline =~ /^\t\t\t(\S+)\t(.+)$/) {
    flushline ();
    # qualifier
    $in_key = 0;
    $in_qual = 1;
    $current_qual = $1;
    $current_val = $2;
  }
}

flushline ();

# close input and output files

close (STDIN);
close (STDOUT);

