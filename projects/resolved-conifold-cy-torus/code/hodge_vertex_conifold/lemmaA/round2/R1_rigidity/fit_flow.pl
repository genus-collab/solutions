#!/usr/bin/perl
# fit_flow.pl -- is the state along the leaf curve a flow of the lambda-line?   Uhat_t(lambda) = [ Psi(phi_t(lambda)) ]_-,
#   phi_t(lambda) = lambda + sum_{m>=1} t^m xi_m(lambda),  xi_m polynomial of degree <= m+1 (universal, independent of Psi).
# Fit xi_m order by order from the exact state of the basic solution Psi = u^{N0} (ode_state.pl) and test on other N.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($K, $ts, $N0, @Ntest) = @ARGV; my $tau = R($ts); @Ntest = (4, 5, 6) unless @Ntest;
my $NN = $K + 6;   # y-truncation for series
# u(y) series, Psi = u^N as a series in y (y = 1/lambda)
my $u = series_mul([R(0), R(1)], series_pow([R(1), R(-1)], $tau, $NN), $NN);
sub psiN { my $N = shift; my $p = [R(1), map { R(0) } 1..$NN]; $p = series_mul($p, $u, $NN) for 1..$N; return $p; }
sub state { my $N = shift; my $L0 = $N - 1; my @lines = `perl "$FindBin::Bin/ode_state.pl" $K $ts $L0`; my %st;
  for (@lines) { next unless /^t\^(\d+) Uhat\(y\): (.*)$/; my $m = $1; for my $e (split /\s+/, $2) { next unless $e =~ /^y(\d+):(.*)$/; $st{"$m,$1"} = R($2); } } return \%st; }
# Laurent series in lambda represented as hash power->coeff (powers can be positive and negative).
sub lmul { my ($a, $b, $lo) = @_; my %c; for my $i (keys %$a) { for my $j (keys %$b) { my $k = $i + $j; next if $k < $lo; $c{$k} = ($c{$k} // R(0)) + $a->{$i}*$b->{$j}; } } return \%c; }
sub ladd { my ($a, $b, $s) = @_; my %c = %$a; for my $k (keys %$b) { $c{$k} = ($c{$k} // R(0)) + ($s // R(1))*$b->{$k}; } return \%c; }
# Psi as Laurent series in lambda: y^n -> lambda^{-n}
sub psiL { my $p = shift; my %h; for my $n (0..$#$p) { $h{-$n} = $p->[$n] unless $p->[$n]->is_zero; } return \%h; }
# derivative d/dlambda
sub ld { my $a = shift; my %c; for my $k (keys %$a) { next if $k == 0; $c{$k-1} = ($c{$k-1} // R(0)) + R($k)*$a->{$k}; } return \%c; }
# compose Psi(phi_t) to order t^M given xi_1..xi_M (hashes), via Taylor: sum_k (delta^k/k!) Psi^{(k)},  delta = sum_m t^m xi_m
sub compose { my ($psi, $xis, $M, $lo) = @_;   # returns array over t-order of Laurent hashes
  my @res = map { {} } 0..$M; $res[0] = { %$psi };
  # delta^k as t-series of Laurent hashes
  my @dpow = ([ {0 => R(1)}, map { {} } 1..$M ]); my @der = ($psi); for my $k (1..$M) { $der[$k] = ld($der[$k-1]); }
  for my $k (1..$M) { my @prev = @{$dpow[$k-1]}; my @cur = map { {} } 0..$M;
    for my $i (0..$M) { next unless %{$prev[$i]}; for my $m (1..$M-$i) { next unless $xis->[$m]; $cur[$i+$m] = ladd($cur[$i+$m], lmul($prev[$i], $xis->[$m], $lo - 50)); } }
    $dpow[$k] = \@cur; my $f = R(1)/fact($k);
    for my $i (0..$M) { next unless %{$cur[$i]}; $res[$i] = ladd($res[$i], lmul($cur[$i], $der[$k], $lo), $f); } }
  return \@res; }
my $st0 = state($N0); my $psi0 = psiL(psiN($N0));
my @xis = (undef);
for my $m (1..$K) {
  # unknown xi_m = sum_{p=0}^{m+1} x_p lambda^p ; linear contribution at order t^m: xi_m * Psi'(lambda); rest from lower xis
  my $base = compose($psi0, [ @xis, {} ], $m, -($K+3))->[$m];   # with xi_m = 0
  my $dpsi = ld($psi0);
  my @eq; for my $n (1..$K+2-$m) { my $target = $st0->{"$m,$n"} // R(0); my $rhs = $target - ($base->{-$n} // R(0));
    my @row = map { $dpsi->{-$n - $_} // R(0) } 0..$m+1; push @row, $rhs; push @eq, \@row; }
  my $nc = $m + 2; my @A = map { [ @$_ ] } @eq; my $rk = 0; my @piv;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#A) { unless ($A[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @A[$rk, $pr] = @A[$pr, $rk];
    my $inv = R(1)/$A[$rk][$c]; $A[$rk] = [ map { $_ * $inv } @{$A[$rk]} ]; for my $i (0..$#A) { next if $i == $rk || $A[$i][$c]->is_zero; my $f = $A[$i][$c]; $A[$i] = [ map { $A[$i][$_] - $f*$A[$rk][$_] } 0..$nc ]; } push @piv, $c; $rk++; }
  my $inc = 0; for my $i ($rk..$#A) { $inc++ unless $A[$i][$nc]->is_zero; }
  my %xi; for my $i (0..$rk-1) { $xi{$piv[$i]} = $A[$i][$nc]; } $xis[$m] = \%xi;
  printf "xi_%d = %s   (rank %d/%d, inconsistent %d)\n", $m, join(' + ', map { "(" . ($xi{$_} // 0) . ")l^$_" } 0..$m+1), $rk, $nc, $inc;
  # test on other N
  for my $N (@Ntest) { my $st = state($N); my $pr = compose(psiL(psiN($N)), \@xis, $m, -($K+3))->[$m]; my $bad = 0;
    for my $n (1..$K+2-$m) { $bad++ if ($st->{"$m,$n"} // R(0)) != ($pr->{-$n} // R(0)); } print "    test N=$N at t^$m: ", ($bad ? "$bad mismatches" : "OK"), "\n"; } }
