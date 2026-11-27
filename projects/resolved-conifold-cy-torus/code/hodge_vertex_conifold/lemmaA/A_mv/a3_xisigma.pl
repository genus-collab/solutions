#!/usr/bin/perl
# a3_xisigma.pl -- components Xi^(k)_h(j1..jr) (tau-independent) of Xi_h = Q_h(s=-1), from lemmaA/data/q_g<h>_r<r>.txt:
#   sigma^h Xi_h = Sum_k (-1)^k sigma^(h-k) Xi^(k)  (k = h: lambda_g part, k = 0: Hurwitz part). Prints each Xi^(k) as a polynomial in j's,
#   and factorises out the lattice-vanishing test on the simplex sum (j_i - 1) <= 2h-3 (j_i >= 1).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($h, $r, $sval) = @ARGV; $sval //= -1; $sval = R($sval);
my $rec = parse_file("$FindBin::Bin/../data/q_g${h}_r${r}.txt"); my $ast = parse_expr($rec->[0][2]);
my @V = map { "j$_" } 1..$r;
my @taus = map { R($_) } ('1', '2', '3', '1/2', '5', '7/3')[0..$h];
my @pts;
for my $tau (@taus) { my $x = -$tau/(1+$tau); my $c = 1 + $x; my $sig = $tau*(1+$tau);
  my $p = eval_poly($ast, \@V, { x => $x, s => $sval }); my $f = $c ** (3*$h - 3 + 2*$r) * $sig ** $h;
  for my $k (keys %$p) { $p->{$k} *= $f; } push @pts, [$sig, $p]; }
my %coef;
for my $i (0..$#pts) { my @b = (R(1)); my $den = R(1);
  for my $jj (0..$#pts) { next if $jj == $i; my @nb = map { R(0) } 0..scalar(@b); for my $k (0..$#b) { $nb[$k+1] += $b[$k]; $nb[$k] -= $b[$k] * $pts[$jj][0]; } @b = @nb; $den *= ($pts[$i][0] - $pts[$jj][0]); }
  for my $k (0..$#b) { for my $key (keys %{$pts[$i][1]}) { $coef{$k}{$key} //= R(0); $coef{$k}{$key} += $b[$k] * $pts[$i][1]{$key} / $den; } } }
sub show { my ($p) = @_; my @t;
  for my $key (sort keys %$p) { next if $p->{$key}->is_zero; my @e = split /,/, $key; my $mon = join('*', map { $e[$_] ? ($V[$_] . ($e[$_] > 1 ? "^$e[$_]" : '')) : () } 0..$#e); push @t, $p->{$key} . ($mon ? "*$mon" : ''); }
  return @t ? join(' + ', @t) : '0'; }
for my $j (0..$h) { my $k = $h - $j; my %p = map { $_ => R((-1)**$k) * $coef{$j}{$_} } keys %{$coef{$j}};
  # lattice test
  my @bad; my $rec2; $rec2 = sub { my ($pre, $left) = @_; if (@$pre == $r) { my $D = 0; $D += $_ - 1 for @$pre; return if $D > 2*$h - 3;
      my $v = poly_eval_at(\%p, [ map { R($_) } @$pre ]); push @bad, "(@$pre)" unless $v->is_zero; return; }
    for my $jv (1..(2*$h)) { $rec2->([@$pre, $jv]); } }; $rec2->([]);
  print "Xi^($k)_$h (s=$sval): ", show(\%p), "\n   lattice zeros on the range: ", (@bad ? "FAIL at @bad" : "all vanish"), "\n"; }
