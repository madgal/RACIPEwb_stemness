#!/usr/bin/perl
# * $Author: alokito $
# * $RCSfile: splitFasta.pl,v $
# * $Revision: 1.2 $
# * $Date: 2004/04/06 04:41:42 $
# * $Name:  $
# *
# * This file is part of Java TreeView
# * Copyright (C) 2001-2003 Alok Saldanha, All Rights Reserved.
# *
# * This software is provided under the GNU GPL Version 2. In particular,
# *
# * 1) If you modify a source file, make a comment in it containing your name and the date.
# * 2) If you distribute a modified version, you must do it under the GPL 2.
# * 3) Developers are encouraged but not required to notify the Java TreeView maintainers at alokito@users.sourceforge.net when they make a useful addition. It would be nice if significant contributions could be merged into the main distribution.
# *
# * A full copy of the license can be found in gpl.txt or online at
# * http://www.gnu.org/licenses/gpl.txt



# splits fasta file into seq files
BEGIN {
    push(@INC, $ENV{'HOME'}.'/perl');
}
use strict;
use vars qw($opt_l $opt_o);
use Getopt::Std;
use ListManipulation;

sub usage() {
    print STDERR qq!Usage:

$0 [-l <listfile>] [-o <outputfile>] blah.fasta

by default, splits fasta file into seq files

Options:

-l only split out ids from list
-o output to single file instead of splitting into several
!;


}
getopts('l:o:'); # column number of pcl,cdt files to join on.

unless (@ARGV) {
    usage();
    exit(1);
}

my (%keep, @keep);
if ($opt_l) {
    my $list = ListManipulation::loadFromTdt($opt_l, 0);
    %keep = %{ListManipulation::counts($list)};
    @keep = @{$list};
}
my @missing = @keep;

if ($opt_o) {
    open (OUT, ">$opt_o");
}
my $i;
my @seen;
my $id = <>;
chomp($id);
while ($id =~ s/^\>//) {
    print STDERR '.';
    $i++;
    if ($i % 100 == 0) {
	print STDERR "$i\n";
    } elsif ($i % 10 == 0) {
	print STDERR " ";
    }
    my $seq = <>;
    while (( $_ = <>) !~ /^\>/) {
	chomp();
	last unless ($_);
	$seq .= $_;
    }
    chomp();
    my $next  = $_;
    # process current...
    if ($opt_l) {
	if (grep {$id =~ /$_/} @missing) {
	    print STDERR "Processing $id\n";
	    @missing = removeEl($id, @missing);
	    push (@seen, $id);
	    print OUT ">$id\n";
	    print OUT $seq, "\n";
	}
    } else {
	print STDERR "Processing $id\n";
	open (OUT, ">$id.fasta") || die " could not open file:$!";
	print OUT $seq, "\n";
	close(OUT);
    }
    $id = $next;
}
if ($opt_l) {
    close (OUT);
    for (@missing) {
	print STDERR "Missing $_\n";
    }
# check to see if we missed any...
}

sub removeEl {
    my $el = shift();
    my @in = @_;
    my @out;
    for (@in) {
	if ($el =~ /$_/) {
	    next;
	} else {
	    push (@out, $_);
	}
    }
    return @out;
}
