package HodgePt;
# Point evaluation of Xi_h(j_1..j_n) (all components) by forgetting the points one at a time, with leaves lam_j = 1/(1 - j psi):
#   K kappa_beta on Mbar_{h,n} = pi^*(K) (1 - psi_n) prod_{b in beta} (pi^* kappa_b + psi_n^b);
#   pi_*[ g(psi_n) prod_{i<n} lam_{j_i}(psi_i) ] = kappa[g] prod_{i<n} lam_{j_i} + g(0) Sum_{i<n} j_i prod_{i<n} lam_{j_i}   (S lam_j = j lam_j),
#   g(z) = lam_{j_n}(z) (1 - z) z^{|A|},  kappa[g] = Sum_c g_c kappa_{c-1}  (kappa_{-1} = 0, kappa_0 = 2h-2+(n-1)).
# Base: Hodge::T(beta) = Int_{Mbar_h} Omega K kappa_beta.  Values are v-polynomials (v = -1/sigma), i.e. the components Xi^(k).
use strict; use warnings; use Math::BigRat; use Hodge;
sub R { Math::BigRat->new(shift) }
our %memo;
sub WP { my ($beta, $J) = @_; my $h = $Hodge::HCUR;
  my @b = sort { $a <=> $b } @$beta; my @jj = sort { $b <=> $a } @$J; my $n = scalar @jj; my $D = 3*$h - 3 + $n;   # kappa-degree bound on Mbar_{h,n}
  my $sb = 0; $sb += $_ for @b; return [R(0)] if $sb > $D;
  return Hodge::T(@b) if $n == 0;
  if ($h == 1 && $n == 1) { my $res = [R(0)]; for my $a (0..1) { $res = Hodge::vadd($res, Hodge::vscale(Hodge::T1(\@b, $a), R($jj[0]) ** $a)); } return $res; }
  my $key = join(',', @b) . "|" . join(',', @jj); return $memo{$key} if exists $memo{$key};
  my $jn = pop @jj; my $Jp = 0; $Jp += $_ for @jj;
  my $res = [R(0)];
  for my $s (Hodge::subms(@b)) { my ($A, $rest, $mult) = @$s; my $sa = 0; $sa += $_ for @$A; my $srest = 0; $srest += $_ for @$rest;
    for (my $e = 0; ; $e++) { my $idx = $e + $sa - 1; last if $srest + ($idx > 0 ? $idx : 0) > $D - 1;
      my $g = $e == 0 ? R(1) : R($jn) ** $e - R($jn) ** ($e - 1); next if $g == 0;
      my $w = $g * $mult;
      if ($idx < 0) { next }
      elsif ($idx == 0) { $res = Hodge::vadd($res, Hodge::vscale(WP($rest, \@jj), $w * (2*$h - 3 + $n))); }
      else { $res = Hodge::vadd($res, Hodge::vscale(WP([@$rest, $idx], \@jj), $w)); }
      last if $jn == 0 && $e >= 1; last if $jn == 1 && $e >= 0 && $idx >= 0 && 0; }
    if ($sa == 0 && $Jp != 0) { $res = Hodge::vadd($res, Hodge::vscale(WP(\@b, \@jj), R($Jp))); } }
  $memo{$key} = $res; return $res; }
1;
