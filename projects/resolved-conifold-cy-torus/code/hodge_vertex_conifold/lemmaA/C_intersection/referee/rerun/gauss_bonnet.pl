#!/usr/bin/perl
# Gauss-Bonnet check of Lemma 2 (c(T(-log)) = Lambda^vee(1) K): the constant term of the Hurwitz component
# Xi^(0)_h(0,..,0) = Int c_top(T_{Mbar_{h,r}}(-log)) must equal chi(M_{h,r}) (Harer-Zagier:
# chi(M_{h,1}) = -B_{2h}/(2h), chi(M_{h,r+1}) = (2-2h-r) chi(M_{h,r})).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use polylib; use Math::BigRat;
my @B = (1, "-1/2", "1/6", 0, "-1/30", 0, "1/42", 0, "-1/30");
opendir(my $dh, $FindBin::Bin) or die; my @files = sort grep { /^comp_g\d+_r\d+(_psi)?\.txt$/ } readdir($dh); closedir $dh;
for my $name (@files) { open(my $fh, '<', "$FindBin::Bin/$name") or die; while (<$fh>) { chomp; my @x = split /\s*\|\s*/; next unless $x[2] == 0;
  my ($h, $r) = @x[0, 1]; my $v = pcoef(parse_expr($x[3], {}), '');
  my $chi = Math::BigRat->new($B[2 * $h]) / (-2 * $h); for my $n (1 .. $r - 1) { $chi *= (2 - 2 * $h - $n); }
  printf "%-22s h=%d r=%d  Xi^H(0^r) = %-10s chi(M_{h,r}) = %-10s %s\n", $name, $h, $r, $v->bstr, $chi->bstr, ($v == $chi ? 'OK' : 'MISMATCH'); } close $fh; }
