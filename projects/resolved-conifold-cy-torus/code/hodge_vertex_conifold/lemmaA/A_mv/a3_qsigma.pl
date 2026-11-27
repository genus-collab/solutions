#!/usr/bin/perl
# a3_qsigma.pl -- read lemmaA/data/q_g<h>_r<r>.txt, convert to c = 1 and write sigma^h Q_h(s; j) = Sum_k (-1)^k sigma^(h-k) Q^(k)(s; j)
# (Q^(k) the tau-independent Hodge components, Q^(h) = lambda_g part) as polynomials in s, j1..jr; prints each Q^(k).
# Usage: perl a3_qsigma.pl h r
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($h, $r) = @ARGV;
my $rec = parse_file("$FindBin::Bin/../data/q_g${h}_r${r}.txt"); my $ast = parse_expr($rec->[0][2]);
my @V = ('s', map { "j$_" } 1..$r);
my @taus = map { R($_) } ('1', '2', '3', '1/2', '5', '7/3')[0..$h];
my @pts;
for my $tau (@taus) { my $x = -$tau/(1+$tau); my $c = 1 + $x; my $sig = $tau*(1+$tau);
  my $p = eval_poly($ast, \@V, { x => $x }); my $f = $c ** (3*$h - 3 + 2*$r) * $sig ** $h;
  for my $k (keys %$p) { $p->{$k} *= $f; } push @pts, [$sig, $p]; }
# Lagrange interpolation in sigma, coefficientwise
my %coef;   # coef{j}{key} = coefficient of sigma^j
for my $i (0..$#pts) { my @b = (R(1)); my $den = R(1);
  for my $jj (0..$#pts) { next if $jj == $i; my @nb = map { R(0) } 0..scalar(@b); for my $k (0..$#b) { $nb[$k+1] += $b[$k]; $nb[$k] -= $b[$k] * $pts[$jj][0]; } @b = @nb; $den *= ($pts[$i][0] - $pts[$jj][0]); }
  for my $k (0..$#b) { for my $key (keys %{$pts[$i][1]}) { $coef{$k}{$key} //= R(0); $coef{$k}{$key} += $b[$k] * $pts[$i][1]{$key} / $den; } } }
sub show { my ($p) = @_; my @t;
  for my $key (sort keys %$p) { next if $p->{$key}->is_zero; my @e = split /,/, $key; my $mon = join('*', map { $e[$_] ? ($V[$_] . ($e[$_] > 1 ? "^$e[$_]" : '')) : () } 0..$#e); push @t, $p->{$key} . ($mon ? "*$mon" : ''); }
  return @t ? join(' + ', @t) : '0'; }
for my $j (0..$h) { my $k = $h - $j; print "Q^($k) [coefficient of (-1)^$k sigma^$j]: ", show({ map { $_ => R((-1)**$k) * $coef{$j}{$_} } keys %{$coef{$j}} }), "\n"; }
