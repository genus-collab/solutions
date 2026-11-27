#!/usr/bin/perl
# Direct check, on the Schur side of the Marino-Vafa formula (no Hodge data), of the first modified-KP equation of the
# 2D Toda hierarchy for the pair (tau_0, tau_1) = (Z(p), Z(e^beta o p)), beta = tau x:
#     2 d_{p2}(F - F1) + d_{p1}^2 (F + F1) + (d_{p1}(F - F1))^2 = 0,   F = log Z,  F1 = log Z(e^beta o p),
#     Z = Sum_nu e^{(tau+1/2) kappa_nu x/2} prod_box (2 sinh(h x/2))^{-1} s_nu(p)   (Route A, (1.1)).
# For every monomial p_lambda, |lambda| <= d, the coefficient (times x^{|lambda|+2}) must vanish as a power series in x
# through x^P.  [p_mu] F = x^{-|mu|} s_mu(x) with s_mu from mvlib::mv_conn (set-partition cumulants of Schur sums).
# Usage: perl mv_mkp_direct.pl tau d P   (also tests the WRONG sign of the d_{p2} term as a control)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use mvlib;
$| = 1;
my ($taus, $d, $P) = @ARGV; $taus //= "2/7"; $d //= 4; $P //= 16;
my $tau = Math::BigRat->new($taus);
my %S;
sub s_of { my @mu = sort { $b <=> $a } @_; my $k = join(',', @mu); $S{$k} //= mv_conn(\@mu, $tau, $P); return $S{$k}; }
sub expser { my ($c) = @_; my @e = map { R(0) } 0..$P; my $t = R(1); for my $n (0..$P) { $e[$n] = $t; $t = $t * $c / ($n + 1); } return \@e; }   # e^{c x}
sub one_minus { my ($e) = @_; my @r = map { -$_ } @$e; $r[0] += 1; return \@r; }
sub one_plus { my ($e) = @_; my @r = map { $_->copy } @$e; $r[0] += 1; return \@r; }
sub mult { my ($m, $v) = @_; my %c; $c{$_}++ for @$v; return $c{$m} // 0; }
# all partitions (parts >= 1) of size <= d
my @lams; for my $n (0..$d) { push @lams, partitions($n); }
my ($n, $bad, $badw) = (0, 0, 0);
for my $lam (@lams) { my @l = @$lam; my $L = 0; $L += $_ for @l;
  my $eL2 = expser($tau * ($L + 2));
  # (i) 2 (m2+1) (1 - e^{beta(L+2)}) s_{lam+2}
  my $t1 = ser_scal(R(2 * (mult(2, \@l) + 1)), ser_mul(one_minus($eL2), s_of(@l, 2), $P));
  # (ii) (m1+2)(m1+1) (1 + e^{beta(L+2)}) s_{lam+1+1}
  my $m1 = mult(1, \@l);
  my $t2 = ser_scal(R(($m1 + 2) * ($m1 + 1)), ser_mul(one_plus($eL2), s_of(@l, 1, 1), $P));
  # (iii) sum over ordered splittings lam = al + be of c_al c_be,  c_al = (m1(al)+1)(1 - e^{beta(|al|+1)}) s_{al+1}
  my %cnt; $cnt{$_}++ for @l; my @vals = sort { $a <=> $b } keys %cnt; my @subs = ([]);
  for my $v (@vals) { my @new; for my $s (@subs) { for my $i (0..$cnt{$v}) { push @new, [@$s, ($v) x $i]; } } @subs = @new; }
  my $t3 = ser_zero($P);
  for my $al (@subs) { my %ca; $ca{$_}++ for @$al; my @be; for my $v (@vals) { push @be, ($v) x ($cnt{$v} - ($ca{$v} // 0)); }
    my $A = 0; $A += $_ for @$al; my $B = 0; $B += $_ for @be;
    my $ca_ = ser_scal(R(mult(1, $al) + 1), ser_mul(one_minus(expser($tau * ($A + 1))), s_of(@$al, 1), $P));
    my $cb_ = ser_scal(R(mult(1, \@be) + 1), ser_mul(one_minus(expser($tau * ($B + 1))), s_of(@be, 1), $P));
    $t3 = ser_add($t3, ser_mul($ca_, $cb_, $P)); }
  my $tot = ser_add(ser_add($t1, $t2), $t3);
  my $wrong = ser_add(ser_add(ser_scal(R(-1), $t1), $t2), $t3);
  my @nz = grep { !$tot->[$_]->is_zero } 0..$P; my @nzw = grep { !$wrong->[$_]->is_zero } 0..$P;
  $n++; $bad++ if @nz; $badw++ if @nzw;
  printf "lambda = (%s): %s ; wrong-sign control: %s\n", join(',', @l), (@nz ? "NONZERO at x^(@nz)" : "0 through x^$P"), (@nzw ? "nonzero (first x^$nzw[0])" : "zero");
}
printf "SUMMARY tau=%s d=%d P=%d: %d monomials, %d failures; wrong sign: %d of %d nonzero\n", $taus, $d, $P, $n, $bad, $badw, $n;
