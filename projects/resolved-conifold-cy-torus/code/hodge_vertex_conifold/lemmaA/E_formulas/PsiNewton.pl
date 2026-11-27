# PsiNewton.pl nmax : test the "regular Newton" structure of Psi = H/(S1 S2) for two leaves:
#   Psi = Sum_{k>=0} Sum_{n'=k}^{max(2k-1,0)} A_{k,n'}(lambda) (j1)_k (j2)_k / ((J-1)(J-2)...(J-n')),
# (j)_k falling factorial, (j1)_k (j2)_k = prod_{i<k} (P - iJ + i^2); A_{k,n'} universal series (sigma-polynomial coefficients,
# independent of J and P).  Every basis function is regular at the lattice points j1, j2 >= 0 (n' <= 2k-1), so this form makes
# Lemma A for r = 2 evident through the H-form.  Solves at each order lambda^{2n} for the A's and reports consistency.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use SerS; use HFun;
sub R { Math::BigRat->new(shift) }
my $nmax = shift // 4;
my $N = HFun::numerators($FindBin::Bin, $nmax);
my $inv = SerS::sexp(SerS::sscale(SerS::logS1S2(), R(-1)));
sub sp { my $a = shift; my $p = {}; for my $i (0..$#$a) { $p = Poly::add($p, Poly::scale(Poly::pow(Poly::var('s'), $i), $a->[$i])) if $a->[$i] != 0; } return $p; }
my $J = Poly::var('J'); my $P = Poly::var('P');
sub ff2 { my $k = shift; my $p = Poly::const(1); for my $i (0..$k-1) { $p = Poly::mul($p, Poly::add(Poly::sub_($P, Poly::scale($J, $i)), Poly::const($i*$i))); } return $p; }
sub prodJ { my ($a, $b) = @_; my $p = Poly::const(1); for my $l ($a..$b) { $p = Poly::mul($p, Poly::sub_($J, Poly::const($l))); } return $p; }
for my $n (1..$nmax) { my $L = 2*$n - 1;
  # numerator of Psi_n over prod_{l=1}^{L}(J-l)
  my $M = {}; for my $a (0..$n) { my $t = Poly::mul(sp(SerS::coeff($inv, $n - $a)), $N->[$a]); my $La = 2*$a - 1; $La = 0 if $La < 0; $t = Poly::mul($t, prodJ($La + 1, $L)); $M = Poly::add($M, $t); }
  # basis
  my @lab; my @bas; push @lab, "k=0"; push @bas, prodJ(1, $L);
  for my $k (1..$n) { for my $np ($k .. 2*$k - 1) { next if $np > $L; push @lab, "k=$k,n'=$np"; push @bas, Poly::mul(ff2($k), prodJ($np + 1, $L)); } }
  # monomials in J, P
  my %mi; my @mons; for my $b (@bas, $M) { for my $m (keys %$b) { my $e = Poly::mono_parse($m); delete $e->{s}; my $k = Poly::mono_key($e); unless (exists $mi{$k}) { $mi{$k} = scalar @mons; push @mons, $k; } } }
  my $ok = 1; my %A;
  for my $spow (0..$n + 1) { my $tgt = Poly::coeff_of($M, 's', $spow);
    my @rows = map { my $m = $_; [ (map { $_->{$m} // R(0) } @bas), ($tgt->{$m} // R(0)) ] } @mons;
    my $nc = scalar @bas; my $row = 0; my @piv;
    for my $c (0..$nc-1) { my $p; for my $i ($row..$#rows) { if ($rows[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @rows[$row, $p] = @rows[$p, $row];
      my $pv = $rows[$row][$c]; $_ = $_ / $pv for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$c] == 0; my $f = $rows[$i][$c]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $c; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    for my $q (0..$#piv) { my $x = $rows[$q][$nc]; $A{$lab[$piv[$q]]}{$spow} = $x if $x != 0; } }
  print "[lambda^", 2*$n, "] Psi in the regular Newton basis (", scalar(@bas), " functions): ", ($ok ? "CONSISTENT" : "inconsistent"), "\n";
  if ($ok) { for my $l (@lab) { next unless $A{$l}; print "    A($l) = ", join(' + ', map { "($A{$l}{$_})s^$_" } sort { $a <=> $b } keys %{$A{$l}}), "\n"; } } }
