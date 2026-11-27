#!/usr/bin/perl
# t_lambdag_laguerre.pl -- the Laguerre transform applied to the lambda_g input m_h(h, R+1^a) = b_h (J+a)^{2h-3+r+a}
# must reproduce Xi^(h)_h(R) = b_h (J-1)(J-2)...(J-N'), N' = 2h-3+r (proof_status 1.3).  Here b_h cancels: test the ratio.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my $bad = 0; my $n = 0;
for my $h (1..6) { for my $r (1..4) { for my $J ($r..($r+8)) {
  my $N = 3*$h - 3 + $r; my $s = R(0);
  for my $a (0..$N) { $s += R((-1)**$a) / fact($a) * lag($N - $a, 2*$h - 2 + $r + $a, R(-$J - $a)) * R($J + $a) ** (2*$h - 3 + $r + $a); }
  my $p = R(1); $p *= ($J - $_) for 1..(2*$h - 3 + $r); $n++;
  if ($s != $p) { $bad++; print "h=$h r=$r J=$J: laguerre=$s  product=$p\n"; } } } }
print "lambda_g Laguerre identity: $n cases, $bad failures\n";
