#!/usr/bin/perl
# Referee: Lemma 5.1 of the report (normalised Toda equation) checked as an EXACT identity of Laurent polynomials in
# w = 1/(1-T), to all orders (not only the leading order at T = oo used for (REC)), for small (g, R).
#   theta D_J Gt_{g,R} = T * ( [beta^{2g}] { E_0(beta) * Sum_{set partitions pi of R} prod_B Ft_{R_B}(beta) } - D_J^2 Gt_{g,R} ),
#   Gt_{g,R} := kappa_R g_{g,R} = e^{-JT} G_{R,g}(T e^{-T}) = Sum_k p_k w^{e+k}   (Lemma 4.2; p_k from the referee's Hurwitz engine),
#   Ft_{R'}(beta) = Sum_{g',k>=1} beta^{2(g'+k-1)} (2/(2k)!) D_{J'}^{2k} Gt_{g',R'},
#   E_0(beta) = exp( Sum_{(g',k) != (0,1)} beta^{2(g'+k-1)} (2/(2k)!) D^{2k} G_{g',empty} ),
#   theta = T d/dT = (w^2 - w) d/dw,  D_J = w (theta + J),  T = 1 - 1/w.
# Unstable inputs: Gt_{0,(j)} = 1/j^2 - T/(j(j+1)), Gt_{0,(j1,j2)} = 1/(j1+j2), D^2 H_0 = T, D H_1 = T^2/(24 (1-T)^2).
# Usage: perl ref_lemma51.pl g j1 j2 ...
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use Math::BigRat; use Math::BigInt;
sub BR { Math::BigRat->new(shift) }
my ($g, @R) = @ARGV; @R = sort { $a <=> $b } @R; my $r = @R;
# Laurent polynomials: hash exponent => BigRat
sub L_add { my ($a, $b, $c) = @_; $c //= BR(1); my %s = map { $_ => $a->{$_}->copy } keys %$a; for (keys %$b) { $s{$_} = ($s{$_} // BR(0)) + $c * $b->{$_} } clean(\%s) }
sub L_mul { my ($a, $b) = @_; my %s; for my $i (keys %$a) { for my $j (keys %$b) { $s{$i + $j} = ($s{$i + $j} // BR(0)) + $a->{$i} * $b->{$j} } } clean(\%s) }
sub L_scal { my ($a, $c) = @_; clean({ map { $_ => $a->{$_} * $c } keys %$a }) }
sub clean { my $h = shift; delete $h->{$_} for grep { $h->{$_}->is_zero } keys %$h; $h }
sub L_theta { my $a = shift; my %s; for my $i (keys %$a) { next if $i == 0; my $c = $a->{$i} * BR($i);   # (w^2 - w) i w^{i-1}
    $s{$i + 1} = ($s{$i + 1} // BR(0)) + $c; $s{$i} = ($s{$i} // BR(0)) - $c } clean(\%s) }
sub L_DJ { my ($a, $J) = @_; my $t = L_add(L_theta($a), $a, BR($J)); my %s = map { ($_ + 1) => $t->{$_} } keys %$t; \%s }
sub L_D { L_DJ($_[0], 0) }
sub L_eq { my ($a, $b) = @_; my $d = L_add($a, $b, BR(-1)); return !%$d }
my $Tw = { 0 => BR(1), -1 => BR(-1) };   # T = 1 - 1/w
my %GT;   # "g'|R'" => Laurent polynomial Gt
sub Gt { my ($gp, @Rp) = @_; my $key = "$gp|" . join(",", @Rp); return $GT{$key} if $GT{$key};
  my $res; my $rp = @Rp; my $Jp = 0; $Jp += $_ for @Rp;
  if ($gp == 0 && $rp == 1) { my $j = $Rp[0]; $res = L_add({ 0 => BR($ENV{PU3} // 1) / BR($j * $j) }, L_scal($Tw, BR(-1) / BR($j * ($j + 1)))) }
  elsif ($gp == 0 && $rp == 2) { $res = { 0 => BR(1) / BR($Jp) } }
  elsif ($rp == 0 && $gp <= 1) { die "Gt: use the explicit D^2k H_0, H_1" }
  else { my ($p0, $ok, $pd, $p) = RefHur::xi0(\@Rp, $gp, 2, 'char'); die "Lambert structure failed for $key" unless $ok;
    my $e = 2 * $gp - 2 + $rp; $res = clean({ map { ($e + $_) => $p->[$_]->copy } 0..$#$p }) }
  $GT{$key} = $res; return $res }
sub Y { my ($gp, $k, @Rp) = @_; my $Jp = 0; $Jp += $_ for @Rp;           # D_{J'}^{2k} Gt_{g',R'} (for R' empty: D^{2k} G_{g',empty})
  if (!@Rp && $gp == 0) { my $x = $Tw; $x = L_D($x) for 1..2 * $k - 2; return $x }                  # D^{2k} H_0 = D^{2k-2} T
  if (!@Rp && $gp == 1) { my $x = { 2 => BR(1) / 24, 1 => BR(-2) / 24, 0 => BR(1) / 24 }; $x = L_scal($x, BR($ENV{PH1} // 1)); $x = L_D($x) for 1..2 * $k - 1; return $x }   # D H_1 = (w-1)^2/24
  my $x = Gt($gp, @Rp); $x = L_DJ($x, $Jp) for 1..2 * $k; return $x }
my $N = $g;   # beta^{2m}, m = 0..g ; series arrays indexed by m
sub S_mul { my ($a, $b) = @_; my @c = map { {} } 0..$N; for my $i (0..$N) { for my $j (0..$N - $i) { $c[$i + $j] = L_add($c[$i + $j], L_mul($a->[$i], $b->[$j])) } } \@c }
# E_0(beta): exponent X_m = Sum_{g'+k-1 = m, (g',k) != (0,1)} (2/(2k)!) D^{2k} G_{g',empty}, m = 1..N (m = 0 term excluded)
my @X0 = ({}); for my $m (1..$N) { my $s = {}; for my $gp (0..$m) { my $k = $m + 1 - $gp; next if $gp == 0 && $k == 1;
    $s = L_add($s, Y($gp, $k), BR(2) / BR(RefHur::fac(2 * $k))) } push @X0, $s }
my @E0 = ({ 0 => BR(1) }, map { {} } 1..$N);   # exp via E0' = X0' E0 in the grading m:  m E0_m = sum_{i} i X0_i E0_{m-i}
for my $m (1..$N) { my $s = {}; for my $i (1..$m) { $s = L_add($s, L_mul($X0[$i], $E0[$m - $i]), BR($i)) } $E0[$m] = L_scal($s, BR(1) / BR($m)) }
# block series Ft_B for every non-empty subset (bitmask) of the labelled parts
my $full = (1 << $r) - 1; my %FT;
for my $mask (1..$full) { my @B = map { $R[$_] } grep { $mask & (1 << $_) } 0..$r - 1; my @s = map { {} } 0..$N;
  for my $m (0..$N) { for my $gp (0..$m) { my $k = $m + 1 - $gp; $s[$m] = L_add($s[$m], Y($gp, $k, @B), BR(2) / BR(RefHur::fac(2 * $k))) } }
  $FT{$mask} = \@s }
my %F = (0 => [{ 0 => BR(1) }, map { {} } 1..$N]);
for my $S (1..$full) { my $low = $S & (-$S); my $rest = $S ^ $low; my $acc = [map { {} } 0..$N]; my $sub = $rest;
  while (1) { my $B = $low | $sub; my $p = S_mul($FT{$B}, $F{$S ^ $B}); $acc->[$_] = L_add($acc->[$_], $p->[$_]) for 0..$N; last if $sub == 0; $sub = ($sub - 1) & $rest }
  $F{$S} = $acc }
my $tot = S_mul(\@E0, $F{$full})->[$g];
my $J = 0; $J += $_ for @R; my $G = Gt($g, @R);
my $lhs = L_theta(L_DJ($G, $J));
my $rhs = L_mul($Tw, L_add($tot, L_DJ(L_DJ($G, $J), $J), BR(-1)));
my $ok = L_eq($lhs, $rhs);
my @ex = sort { $a <=> $b } keys %$lhs;
printf "Lemma 5.1 exact in w: g=%d R=(%s): %s  (LHS has %d terms, w-exponents %s..%s)\n", $g, join(",", @R), ($ok ? "IDENTITY HOLDS" : "FAILS"), scalar(@ex), $ex[0] // '-', $ex[-1] // '-';
unless ($ok) { my $d = L_add($lhs, $rhs, BR(-1)); print "  difference: ", join(" + ", map { "($d->{$_}) w^$_" } sort { $a <=> $b } keys %$d), "\n" }
