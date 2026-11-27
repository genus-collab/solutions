# p1moments.pl h : the first Lambert coefficient of the Hurwitz class, p_1(h,j) = Int_{Mbar_{h,1}} c(T(-log)) kappa_1 lam_j
#   = [eps^1] Int Lambda^vee(1) K_{eps-1} lam_j,  K_{eps-1} <-> dilaton leaves (eps-1) psi^2/(1-psi)  (Route B, Lemma 4.2),
# hence p_1 = Sum_m (m (-1)^{m-1}/m!) Int_{Mbar_{h,1+m}} Lambda^vee(1) lam_j prod_{k=1}^m Tt(psi_{1+k}),  Tt = psi^2/(1-psi).
# Then checks the Mumford reduction of the lambda_1-component (report, Prop. M1):
#   12 Xi^(1)_h(j) = p_1(h,j) - (Xi^(0)_h(j) - Xi^(0)_h(0))/j + (1/2)(j-2h+2)(j-2h+3) Xi^(0)_{h-1}(j)
#                    + Sum_{g1=1}^{h-1} (j-2g1+1) Xi^(0)_{h-g1}(0) Xi^(0)_{g1}(j)
# against the component tables refdata/ref_engine_h<g>_r1.txt, at j = 1..3h+1.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use Pair;
sub Q { Math::BigRat->new(shift) }
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
my $h = shift // 3;
# Lambda^vee(1) as ch-polynomial
my $L1 = Pair::lamvee($h, 1);
# p_1 as a polynomial in j: coefficient of j^a = Sum_mu L1_mu Sum_P weight * m * I(h, mu, [a, D])
my @p1 = map { Q(0) } 0..3*$h;
for my $mu (keys %$L1) { my $c = $L1->{$mu}; next if $c == 0; my @mu = split /,/, $mu; my $dm = RefHodge::cdeg($mu); my $B = 3*$h - 2 - $dm; next if $B < 0;
  for my $P (@{RefHodge::allparts($B)}) { my $m = scalar @$P; next if $m == 0; my $ex = 0; $ex += $_ for @$P; my $t = $B - $ex;
    my @D = map { $_ + 1 } @$P; my %mult; $mult{$_}++ for @D; my $den = Q(1); $den *= fact($_) for values %mult;
    my $w = (($m % 2) ? Q(1) : Q(-1)) * $m / $den;   # m (-1)^{m-1} / (multiset symmetry)
    my $v = RefHodge::I($h, \@mu, [$t, @D]); next if $v == 0; $p1[$t] += $c * $w * $v; } }
# component tables
sub table { my $g = shift; my @T; open my $fh, '<', "$FindBin::Bin/refdata/ref_engine_h${g}_r1.txt" or die "no table $g"; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; $rest //= ''; my @c;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($a, $v) = split /:/, $t; $c[$a] = Q($v); } $T[$k] = \@c; } close $fh; return \@T; }
sub ev { my ($p, $x) = @_; my $s = Q(0); for my $a (0..$#$p) { next unless defined $p->[$a]; $s += $p->[$a] * Q($x) ** $a; } return $s; }
my %X0; for my $g (1..$h) { my $T = table($g); $X0{$g} = $T->[0]; } my $X1 = table($h)->[1];
my ($n, $bad) = (0, 0);
for my $j (1..3*$h + 1) {
  my $lhs = 12 * ev($X1, $j);
  my $rhs = ev(\@p1, $j) - (ev($X0{$h}, $j) - ev($X0{$h}, 0)) / $j;
  $rhs += Q(1)/2 * ($j - 2*$h + 2) * ($j - 2*$h + 3) * ev($X0{$h-1}, $j) if $h >= 2;
  for my $g1 (1..$h-1) { $rhs += ($j - 2*$g1 + 1) * ev($X0{$h - $g1}, 0) * ev($X0{$g1}, $j); }
  $n++; $bad++ if $lhs != $rhs;
  printf "  j=%2d  12 Xi^(1) = %s   p_1 = %s\n", $j, $lhs, ev(\@p1, $j) if $ENV{VERBOSE} || $j <= 2*$h - 2; }
printf "genus %d: Mumford reduction of the lambda_1 component checked at %d values of j: %d mismatches\n", $h, $n, $bad;
print "p_1(h,j) polynomial: ", Pair::pstr(\@p1), "\n";
