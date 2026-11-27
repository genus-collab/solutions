#!/usr/bin/perl
# Referee re-implementation of the H' test and a sharpness probe beyond the proved range of b.
#   I(j, b) := Int_{Mbar_{h,r}} c_{N - |b|}(-R pi_* S_D) prod_i psi_i^{b_i},  c(-R pi_* S_D) = c(T(-log)) / prod_i prod_{k=1}^{j_i} (1 - k psi_i),
# computed from the k = 0 component polynomial X_0(j) = Sum_a <c(T(-log)) prod psi^a> prod j^a (comp_g<h>_r<r>.txt, here re-derived by
# ref_ucomp.pl and read from out/ucomp_g<h>_r<r>.txt).  Complete homogeneous symmetric polynomials h_n(1..j) by the generating function.
# Tests, for all j_i >= 0 with J >= 1, D <= 2h-3:
#   (a) H' as stated: 0 <= b_i <= j_i, |b| <= J-1  ->  I = 0 (Theorem H');
#   (b) box form: Xi^(0)_h(m) = 0 for all 0 <= m_i <= j_i, m != 0 (equivalent to (a) by partial fractions);
#   (c) beyond: some b_i = j_i + 1 (and |b| <= J - 1): the proof gives no vanishing; report how many are non-zero.
# Usage: perl ref_Hprime.pl
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rpoly; use Math::BigRat;
sub hcomp { my ($n, $j) = @_; return Math::BigRat->new($n == 0 ? 1 : 0) if $j == 0; # h_n(1, ..., j) = coefficient of t^n in prod_k 1/(1 - k t)
  my @c = (Math::BigRat->new(1), map { Math::BigRat->new(0) } 1 .. $n); for my $k (1 .. $j) { for my $m (1 .. $n) { $c[$m] += $c[$m - 1] * $k; } } return $c[$n]; }
my %tot;
for my $hr ([2, 1], [2, 2], [3, 1], [3, 2], [4, 1], [4, 2]) { my ($h, $r) = @$hr;
  open(my $fh, '<', "$FindBin::Bin/out/ucomp_g${h}_r${r}.txt") or die "run ref_ucomp.pl first"; my $X0; while (<$fh>) { s/\r//; chomp; my @f = split /\s*\|\s*/; $X0 = parse_rp($f[3]) if $f[2] == 0; } close $fh;
  my @mono = map { [ split /,/ ] } keys %$X0; my %coef = map { $_ => $X0->{$_} } keys %$X0;
  my $I = sub { my ($j, $b) = @_; my $S = Math::BigRat->new(0);
    for my $k (keys %coef) { my @a = split /,/, $k; my $w = Math::BigRat->new(1); my $ok = 1;
      for my $i (0 .. $r - 1) { my $nn = $a[$i] - $b->[$i]; if ($nn < 0) { $ok = 0; last; } $w *= hcomp($nn, $j->[$i]); }
      $S += $w * $coef{$k} if $ok; } return $S; };
  my @pts; my $rec; my @cur; $rec = sub { my $i = shift; if ($i == $r) { my $J = 0; $J += $_ for @cur; my $D = 0; $D += $_ - 1 for @cur; push @pts, [@cur] if ($J >= 1 && $D <= 2 * $h - 3); return; }
    for my $jj (($i ? $cur[-1] : 0) .. 2 * $h) { push @cur, $jj; $rec->($i + 1); pop @cur; } }; $rec->(0);
  my ($na, $za, $nb, $zb, $nc, $zc) = (0, 0, 0, 0, 0, 0); my @exc;
  for my $j (@pts) { my $J = 0; $J += $_ for @$j;
    my @bs = ([]); for my $i (0 .. $r - 1) { my @nb; for my $b (@bs) { push @nb, [@$b, $_] for 0 .. $j->[$i] + 1; } @bs = @nb; }
    for my $b (@bs) { my $sb = 0; $sb += $_ for @$b; next if $sb > $J - 1; my $beyond = grep { $b->[$_] > $j->[$_] } 0 .. $r - 1;
      my $v = $I->($j, $b);
      if (!$beyond) { $na++; $za++ unless $v->is_zero; } else { $nc++; unless ($v->is_zero) { $zc++; push @exc, "j=(@$j) b=(@$b): " . $v->bstr if @exc < 4; } } }
    # (b) box form
    my @ms = ([]); for my $i (0 .. $r - 1) { my @nm; for my $m (@ms) { push @nm, [@$m, $_] for 0 .. $j->[$i]; } @ms = @nm; }
    for my $m (@ms) { my $M = 0; $M += $_ for @$m; next if $M == 0; $nb++; my $v = rp_eval($X0, @$m); $zb++ unless $v->is_zero; } }
  print "genus $h r $r: (a) H' pairs $na, non-zero $za;  (b) box points $nb, non-zero $zb;  (c) beyond (some b_i = j_i + 1): $nc pairs, non-zero $zc",
    (@exc ? "  e.g. " . join(' | ', @exc) : ''), "\n";
  $tot{a} += $na; $tot{za} += $za; $tot{c} += $nc; $tot{zc} += $zc; }
print "TOTAL: H' pairs $tot{a} (non-zero $tot{za}); beyond pairs $tot{c} (non-zero $tot{zc})\n";
