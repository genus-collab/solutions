#!/usr/bin/perl
# Second referee: direct check, in BILINEAR (Hirota) form on the tau functions themselves (no logarithm, no cumulants), of the
# u_k-linear members P_k(D) tau_0 . tau_1 = 0,  P_k = D_{t_k} + p_k(-Dtilde),  of the 2D Toda hierarchy at equal second times,
# for tau_0 = Z(p), tau_1 = Z(e^{tau x} o p), Z = Marino-Vafa generating function (Schur side, referee's own R2mv.pm), mod p.
#   D^alpha f.g = sum_{beta subset alpha} (-1)^{|alpha - beta|} (d^beta f)(d^{alpha-beta} g),  d_{t_i} = i d_{p_i}.
# For every monomial p_lambda with |lambda| <= d:  x^{|lambda|+k} [p_lambda] P_k(D) tau_0.tau_1 must vanish through x^P.
# Controls: (b) partner tau_{-1} = Z(e^{-tau x} o p) instead of tau_1;  (c) reversed order tau_1 . tau_0.
# Usage: perl r2_mkp_bilinear.pl prime tau_mod_p k d P
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2mv; use R2rel;
$| = 1;
my ($prime, $tau, $k, $d, $P) = @ARGV; R2mv::setp($prime); my $p = $prime;
R2mv::prep($P + 4);
my $Pk = R2rel::hirota_Pk($k);
print "# P_$k = ", R2rel::show_P($Pk), "\n";
# zt_mu = x^{|mu|} [p_mu] Z  for |mu| <= d + k
my %zt; my $t0 = time;
for my $n (0..$d + $k) {
  my @mus; R2mv::each_partition($n, sub { push @mus, [ @{$_[0]} ]; });
  my %acc; $acc{join(',', @$_)} = [ (0) x ($P + 1) ] for @mus;
  for my $nu (@mus) { my $c = R2mv::conjugate($nu); my @h = R2mv::hooks($nu, $c); my %hm; $hm{$_}++ for @h; my $H = 1; $H = $H * $_ % $p for @h;
    my $w = R2mv::wser(\%hm, $H, R2mv::contsum($nu), $tau, $P);
    for my $mu (@mus) { my @np = grep { $_ >= 2 } @$mu; my $ch = R2mv::chi($nu, \@np); next unless $ch;
      my $cf = $ch * R2mv::inv(R2mv::zS(@$mu)) % $p; my $row = $acc{join(',', @$mu)}; for my $m (0..$P) { $row->[$m] = ($row->[$m] + $cf * $w->[$m]) % $p if $w->[$m]; } } }
  R2mv::clear_chi(); $zt{$_} = $acc{$_} for keys %acc; }
printf "# Schur side done (%d s)\n", time - $t0;
sub key { join(',', sort { $b <=> $a } @_) }
sub cnt { my %c; $c{$_}++ for @_; \%c }
sub expser { my ($c) = @_; my @e = (1); for my $n (1..$P) { $e[$n] = $e[$n-1] * $c % $p * R2mv::inv($n) % $p; } \@e }   # e^{c x}
my %EX; sub ex { my $c = R2mv::md($_[0]); $EX{$c} //= expser($c) }
# [p_lam] d^{ins} of Z(e^{s tau x} o p) times x^{|lam| + |ins|}:  mult * e^{s tau |mu|} zt_mu,  mu = lam + ins
sub dcoef { my ($lam, $ins, $s) = @_; my $cl = cnt(@$lam); my $ci = cnt(@$ins); my $m = 1;
  for my $v (keys %$ci) { my $a = $cl->{$v} // 0; for my $t (1..$ci->{$v}) { $m = $m * ($a + $t) % $p; } }
  my @mu = (@$lam, @$ins); my $z = $zt{key(@mu)} or die "no zt for " . key(@mu); my $N = 0; $N += $_ for @mu;
  my $e = ex($s * $tau % $p * $N % $p); my $q = R2mv::smul($z, $e, $P); [ map { $_ * $m % $p } @$q ] }
sub submultisets { my @l = @_; my $c = cnt(@l); my @vals = sort { $a <=> $b } keys %$c; my @s = ([]);
  for my $v (@vals) { my @n; for my $x (@s) { for my $i (0..$c->{$v}) { push @n, [@$x, ($v) x $i]; } } @s = @n; } @s }
sub minus { my ($l, $s) = @_; my $c = cnt(@$s); my @r; for my $v (@$l) { if (($c->{$v} // 0) > 0) { $c->{$v}--; } else { push @r, $v; } } \@r }
# x^{|lam|+k} [p_lam] P(D) f.g  with f = Z(e^{sf tau x} o p), g = Z(e^{sg tau x} o p)
sub hirota { my ($lam, $sf, $sg) = @_; my $tot = [ (0) x ($P + 1) ];
  for my $mono (@$Pk) { my ($c, $alpha) = @$mono; my $cR = R2mv::rat($c); my $q = @$alpha;
    for my $mask (0..(1 << $q) - 1) { my @bf = map { $alpha->[$_] } grep { ($mask >> $_) & 1 } 0..$q-1; my @bg = map { $alpha->[$_] } grep { !(($mask >> $_) & 1) } 0..$q-1;
      my $sgn = (@bg % 2) ? -1 : 1; my $pf = 1; $pf = $pf * $_ % $p for @bf, @bg;   # d_{t_i} = i d_{p_i}
      my $cf = R2mv::md($sgn * $cR % $p * $pf % $p);
      for my $l1 (submultisets(@$lam)) { my $l2 = minus($lam, $l1);
        my $a = dcoef($l1, \@bf, $sf); my $b = dcoef($l2, \@bg, $sg); my $pr = R2mv::smul($a, $b, $P);
        $tot->[$_] = ($tot->[$_] + $cf * $pr->[$_]) % $p for 0..$P; } } }
  $tot }
my @lams; for my $n (0..$d) { R2mv::each_partition($n, sub { push @lams, [ @{$_[0]} ]; }); }
my ($n, $bad, $cb, $cc) = (0, 0, 0, 0);
for my $lam (@lams) { my $v = hirota($lam, 0, 1); my @nz = grep { $v->[$_] } 0..$P;
  my $vb = hirota($lam, 0, -1); my @nzb = grep { $vb->[$_] } 0..$P; my $vc = hirota($lam, 1, 0); my @nzc = grep { $vc->[$_] } 0..$P;
  $n++; $bad++ if @nz; $cb++ if @nzb; $cc++ if @nzc;
  printf "lambda=(%s): %s; control tau_0.tau_{-1}: %s; control tau_1.tau_0: %s\n", join(',', @$lam), (@nz ? "NONZERO at x^(@nz[0..($#nz < 4 ? $#nz : 4)])" : "0 through x^$P"),
    (@nzb ? "nonzero (first x^$nzb[0])" : "zero"), (@nzc ? "nonzero (first x^$nzc[0])" : "zero"); }
printf "SUMMARY P_%d p=%d tau=%d d=%d P=%d: %d monomials, %d failures; controls nonzero: %d/%d (tau_{-1}), %d/%d (reversed); %d s\n", $k, $p, $tau, $d, $P, $n, $bad, $cb, $n, $cc, $n, time - $t0;
