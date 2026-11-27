#!/usr/bin/perl
# t_flow_K.pl -- consistency of the fitted flow with the closed form T1:
#   K^flow_t(mu) := log alpha(t) + log u(mu) - log u(phi_t^{-1}(mu))  versus  K_t(mu) = (1+tau) Kappa(t mu) - tau Kappa(t(mu-1)),
# compared on the coefficients of mu^1, mu^2, ... (the only ones that pair with Psi = O(mu^-2)).
# Input: the fitted xi_m (from out/flow_tau*.txt) for orders 1..M.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($file, $taus, $M) = @ARGV; my $tau = R($taus);
my %xi; open my $fh, '<', $file or die; while (<$fh>) { next unless /^xi_(\d+) = (.*?)\s+\(rank (\d+)\/(\d+)/; my ($m, $rest, $r, $n) = ($1, $2, $3, $4); next unless $r == $n; my @c; while ($rest =~ /\(([^)]*)\)l\^(\d+)/g) { $c[$2] = R($1); } $xi{$m} = [ map { $_ // R(0) } @c ]; } close $fh;
$M = (sort { $b <=> $a } grep { my $k = $_; !grep { !exists $xi{$_} } 1..$k } keys %xi)[0] if !$M;
print "using xi_1..xi_$M\n";
# Laurent series in mu: hash power -> coeff; truncate below LOW
my $LOW = -6;
sub lmul { my ($a, $b) = @_; my %c; for my $i (keys %$a) { for my $j (keys %$b) { my $k = $i + $j; next if $k < $LOW; $c{$k} = ($c{$k} // R(0)) + $a->{$i}*$b->{$j}; } } return \%c; }
sub ladd { my ($a, $b, $s) = @_; my %c = %$a; for my $k (keys %$b) { $c{$k} = ($c{$k} // R(0)) + ($s // R(1))*$b->{$k}; } return \%c; }
sub poly2h { my $p = shift; my %h; for my $i (0..$#$p) { $h{$i} = $p->[$i] unless $p->[$i]->is_zero; } return \%h; }
# bivariate: t-order array of hashes
sub bmul { my ($a, $b) = @_; my @c = map { {} } 0..$M; for my $i (0..$M) { for my $j (0..$M-$i) { $c[$i+$j] = ladd($c[$i+$j], lmul($a->[$i], $b->[$j])); } } return \@c; }
sub badd { my ($a, $b, $s) = @_; return [ map { ladd($a->[$_], $b->[$_], $s) } 0..$M ]; }
# delta(mu) := phi^{-1}(mu) - mu : solve phi(mu + delta) = mu  i.e.  delta + sum_m t^m xi_m(mu + delta) = 0  iteratively
my @delta = map { {} } 0..$M;
for my $it (1..$M+1) { my @nd = map { {} } 0..$M;
  for my $m (1..$M) { # xi_m(mu + delta) = sum_p c_p (mu + delta)^p  as bivariate
    my $base = [ {1 => R(1)}, @delta[1..$M] ]; my $pw = [ {0 => R(1)}, map { {} } 1..$M ]; my $acc = [ map { {} } 0..$M ];
    for my $p (0..$#{$xi{$m}}) { $acc = badd($acc, $pw, $xi{$m}[$p]) unless $xi{$m}[$p]->is_zero; $pw = bmul($pw, $base); }
    for my $i (0..$M-$m) { $nd[$i+$m] = ladd($nd[$i+$m], $acc->[$i], R(-1)); } }
  @delta = @nd; }
# log u(mu + delta) - log u(mu) = tau log(1 + delta/(mu-1)) - (1+tau) log(1 + delta/mu);  1/(mu-1) = sum_k mu^{-k-1}
my %inv_m = (-1 => R(1)); my %inv_m1; for my $k (0..10) { $inv_m1{-$k-1} = R(1); }
sub blog1p { my $x = shift; my @res = map { {} } 0..$M; my $pw = [ {0 => R(1)}, map { {} } 1..$M ]; for my $k (1..$M) { $pw = bmul($pw, $x); $res[$_] = ladd($res[$_], $pw->[$_], R((-1)**($k+1))/$k) for 0..$M; } return \@res; }
my $x1 = [ map { lmul($_, \%inv_m1) } @delta ]; my $x2 = [ map { lmul($_, \%inv_m) } @delta ];
my $l1 = blog1p($x1); my $l2 = blog1p($x2);
my $dlogu = badd([ map { lmul($_, {0 => $tau}) } @$l1 ], [ map { lmul($_, {0 => 1+$tau}) } @$l2 ], R(-1));
# alpha(t)
my $Tt = T_of_t_series($M); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $M), $Tt, $M);
my %Pj; for my $j (2..$M+2) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $M), $Tt, $M); my $p = [R(1), map { R(0) } 1..$M]; $p = series_mul($p, $h1, $M) for 1..$j; my $rho = series_mul($hj, series_inv($p, $M), $M);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$M ]; }
my @al = (R(1), map { R(0) } 1..$M);
for my $it (0..$M+1) { my @nw = (R(1), map { R(0) } 1..$M); my $wp = [@al]; for my $j (2..$M+2) { $wp = series_mul($wp, \@al, $M); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $M); $nw[$_] = $nw[$_] + $q->[$_] for 0..$M; } @al = @nw; }
my $logal = series_log(\@al, $M);
# K^flow = log alpha - dlogu  (since log u(phi^{-1}(mu)) - log u(mu) = dlogu)
my @Bn = (R(1)); for my $mm (1..$M+2) { my $s = R(0); for my $k (0..$mm-1) { $s = $s + binom($mm+1, $k) * $Bn[$k]; } $Bn[$mm] = -$s/($mm+1); }
my @kap = (R(0), $tau + R('1/2'), map { $Bn[$_]/($_*fact($_)) } 2..$M+2);
for my $m (1..$M) { my %kf = %{ ladd({0 => $logal->[$m]}, $dlogu->[$m], R(-1)) };
  my @got = map { $kf{$_} // R(0) } 1..$m; my @want = map { my $p = $_; $kap[$m]*((($p == $m) ? 1+$tau : R(0)) - $tau*binom($m,$p)*R((-1)**($m-$p))) } 1..$m;
  my $ok = 1; for my $i (0..$m-1) { $ok = 0 if $got[$i] != $want[$i]; }
  print "t^$m: K^flow coefficients mu^1..mu^$m = (", join(', ', map { "$_" } @got), ")  T1: (", join(', ', map { "$_" } @want), ")  ", ($ok ? "AGREE" : "DIFFER"), "\n"; }
