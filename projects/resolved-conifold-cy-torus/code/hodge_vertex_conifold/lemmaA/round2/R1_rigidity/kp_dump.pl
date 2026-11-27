#!/usr/bin/perl
# kp_dump.pl -- print the MV-normalised multi-leaf values Y(T) (KP propagation of kp_point.pl, exact) as linear forms in y_D..y_Lmax.
# Usage: perl kp_dump.pl D Lmax tau
use strict; use warnings; use FindBin;
my ($D, $Lmax, $tau) = @ARGV;
$ENV{FAST} = 1; $ENV{DUMPY} = 1;
exec("perl", "$FindBin::Bin/kp_point.pl", $D, $Lmax, $tau, "none", 0);
