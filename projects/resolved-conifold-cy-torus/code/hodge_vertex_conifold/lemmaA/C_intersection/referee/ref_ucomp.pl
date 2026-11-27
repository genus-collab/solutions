#!/usr/bin/perl
# Referee re-implementation (independent algorithm) of the Hodge-component extraction.
# Input: an exact rational function Xi^(c-norm)(x; j) = N(x, s, j)/(K x^a (1+x)^b) (Mathematica InputForm).
# Steps (all exact, integer arithmetic): rescale to c = 1 by (1+x)^E, E = 3h-3+2r (brief s.1); substitute s = -1 if s occurs;
# the result must be a Laurent polynomial in x that is a polynomial of degree <= h in u = -1/sigma = (1+x)^2/x
# (decomposition (C)). It is decomposed by peeling off the lowest power x^-k (u^k = x^-k (1+x)^2k).
# No sampling in x, no linear solve (the author's components.pl samples h+2 values of x).
# Output: components X_0..X_h as polynomials in j1..jr, and a comparison with the author's comp file if given.
# Usage: perl ref_ucomp.pl h r infile [field] [authorcompfile] [outfile]
#   field = index of the " | "-separated field holding the expression (default: last field); for multiprobe lines give 'mp'.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rpoly; use Math::BigInt; use Math::BigRat;
my ($h, $r, $in, $field, $author, $outf) = @ARGV;
open(my $fh, '<', $in) or die "cannot open $in"; my $txt = join('', <$fh>); close $fh;
my $expr;
if (defined $field && $field eq 'mp') { # multiprobe output: line "genus h symbolic r=..: EXPR" or "genus h r=1 symbolic: EXPR"
  for my $line (split /\n/, $txt) { if ($line =~ /^genus $h (?:symbolic r=$r|r=$r symbolic):\s*(.*)$/) { $expr = $1; } } die "no symbolic line for genus $h r $r" unless defined $expr; }
else { my @f = split /\s*\|\s*/, $txt; $expr = $f[(defined $field && $field ne '') ? $field : $#f]; }
my ($ns, $ds) = split_frac($expr);
my $N = parse_num($ns); my $Dn = parse_num($ds);
# denominator = K x^a (1+x)^b ?
my $a = ip_lowx($Dn); my $Dq = ip_shiftx($Dn, -$a); my $b = 0;
while (1) { my @k = keys %$Dq; last if (@k == 1 && $k[0] eq rpoly::zkey()); my $q = eval { ip_divx1($Dq) }; die "denominator not of the form K x^a (1+x)^b" unless defined $q; $Dq = $q; $b++; }
my $K = $Dq->{rpoly::zkey()};
my $E = 3 * $h - 3 + 2 * $r;
my $Pp = $N;
if ($E >= $b) { $Pp = ip_mul($Pp, ip_pow(ip_add(ip_const(1), ip_var('x')), $E - $b)); } else { for (1 .. $b - $E) { $Pp = ip_divx1($Pp); } }
$Pp = ip_subst_s($Pp, -1);           # s = -1 (no-op if s absent)
$Pp = ip_shiftx($Pp, -$a);           # Laurent polynomial in x
my $low = ip_lowx($Pp); die "power x^$low below x^-$h: not a polynomial of degree <= h in 1/sigma" if $low < -$h;
my @X;
for (my $k = $h; $k >= 0; $k--) {
  # coefficient of x^-k (polynomial in the j's), must be the lowest power left
  my %ck; for my $key (keys %$Pp) { my @e = split /,/, $key; die "unexpected lower power x^$e[0] at step k = $k" if $e[0] < -$k; next unless $e[0] == -$k; $e[0] = 0; $ck{join(',', @e)} = $Pp->{$key}->copy; }
  $X[$k] = rp_from_ip(\%ck, $K->bstr);
  next unless %ck;
  # subtract ck * x^-k (1+x)^(2k)
  my $sub = ip_shiftx(ip_mul(\%ck, ip_pow(ip_add(ip_const(1), ip_var('x')), 2 * $k)), -$k);
  $Pp = ip_add($Pp, ip_scale($sub, -1)); }
my $rem = scalar(keys %$Pp);
print "genus $h r $r ($in" . (defined $field ? ", $field" : '') . "): denominator K = $K, x^$a, (1+x)^$b; decomposition in u = (1+x)^2/x ", ($rem ? "FAILED (remainder with $rem terms)" : "exact"), "\n";
if (defined $outf && $outf ne '') { open(my $o, '>', $outf) or die; for my $k (0 .. $h) { print $o "$h | $r | $k | ", rp_str($X[$k]), "\n"; } close $o; }
if (defined $author && $author ne '' && -e $author) {
  open(my $af, '<', $author) or die; my @A; while (<$af>) { chomp; s/\r//; next unless /\|/; my @f = split /\s*\|\s*/; $A[$f[2]] = parse_rp($f[3]); } close $af;
  for my $k (0 .. $h) { my $d = rp_add($X[$k], rp_scale($A[$k] // {}, -1)); print "  k = $k: ", (rp_iszero($d) ? "agrees with $author" : "DIFFERS from $author: " . rp_str($d)), "\n"; } }
