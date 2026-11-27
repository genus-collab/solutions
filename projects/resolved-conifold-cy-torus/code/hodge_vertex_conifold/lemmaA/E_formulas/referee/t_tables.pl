# t_tables.pl [gmax] : the referee engine against the Mathematica Hodge-side tables of lemmaA/data (Mumford's formula, hvc_lib.m),
# which the author did NOT use: psi_g<g>.txt (raw <Omega psi^a>), mv_g<g>_l<l>.txt (no dilaton leaves, up to 3 leaves) and
# q_g<g>_r<r>.txt (the full s-polynomial of the interpolated dilaton leaves s Tt, Tt = psi^2/(c - psi)).
# Normalisations: Omega_c = prod Lambda^vee(w)/w^g with w = (c, -x, -1) is used DIRECTLY for the psi tables (no conversion);
# for mv/q the c = 1 engine value is converted by Xi^c = c^{-(3g-3+2r)} Xi^{c=1} (tau = -x/(1+x)), so the psi tables also test
# that conversion indirectly.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use RefParse; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 4; my $t0 = time; $| = 1;
my $DATA = "$FindBin::Bin/../../data";
my @xs = (Q(2), Q('1/3'), Q('-3/5'));
sub readlines { my $f = shift; open my $fh, '<', $f or die "no $f"; local $/; my $t = <$fh>; return grep { /\S/ } split /\n/, $t; }
my ($ntot, $nbad) = (0, 0);
for my $g (1..$gmax) {
  my $maxd = $g == 1 ? 1 : 3*$g - 3;
  # ---- psi tables
  if (-e "$DATA/psi_g$g.txt") { my ($n, $bad) = (0, 0);
    for my $x (@xs) { my $c = 1 + $x; my $om = RefHodge::omega($g, [$c, -$x, Q(-1)], $maxd);
      for my $line (readlines("$DATA/psi_g$g.txt")) { my ($gg, $nn, $a, $e) = split /\s*\|\s*/, $line; my @a = $a =~ /(\d+)/g;
        my $v = Q(0); for my $k (keys %$om) { my @mu = split /,/, $k; $v += $om->{$k} * RefHodge::I($g, \@mu, \@a); }
        my $d = RefParse::evalx($e, { x => $x }); $n++; $bad++ if $v != $d; } }
    printf "genus %d psi table: %d comparisons, %d mismatches [%.1f s]\n", $g, $n, $bad, time - $t0; $ntot += $n; $nbad += $bad; }
  # ---- mv tables (no dilaton leaves)
  for my $l (1..3) { my $f = "$DATA/mv_g${g}_l$l.txt"; next unless -e $f; my ($line) = readlines($f); my ($gg, $ll, $e) = split /\s*\|\s*/, $line, 3;
    my ($n, $bad) = (0, 0);
    for my $x (@xs) { my $c = 1 + $x; my $tau = -$x / $c; my $om = RefHodge::omega($g, [Q(1), $tau, -1 - $tau], $maxd);
      my @pts = $l == 1 ? ([0], [1], [2], [3], [5], [-2], ['1/2']) : $l == 2 ? ([1, 1], [2, 3], [0, 4], [-1, 2], [3, 3]) : ([1, 2, 3], [2, 2, 2], [0, 1, 5], [-1, 2, 2]);
      for my $p (@pts) { my $q = RefHodge::xi($g, $om, $p, 0); my $mine = $q->[0] / ($c ** (3*$g - 3 + 2*$l));
        my %v = (x => $x); $v{"m$_"} = Q($p->[$_ - 1]) for 1..$l; my $d = RefParse::evalx($e, \%v); $n++; $bad++ if $mine != $d; } }
    printf "genus %d mv table l = %d: %d comparisons, %d mismatches [%.1f s]\n", $g, $l, $n, $bad, time - $t0; $ntot += $n; $nbad += $bad; }
  # ---- q tables (s-polynomial)
  for my $r (1..2) { my $f = "$DATA/q_g${g}_r$r.txt"; next unless -e $f; my ($line) = readlines($f); my ($gg, $rr, $e) = split /\s*\|\s*/, $line, 3;
    my ($n, $bad) = (0, 0);
    for my $x (@xs[0..1]) { my $c = 1 + $x; my $tau = -$x / $c; my $om = RefHodge::omega($g, [Q(1), $tau, -1 - $tau], $maxd);
      my @pts = $r == 1 ? ([0], [2], [3], [7], ['-1/2']) : ([2, 2], [0, 3], [1, 4], [3, 5], [-1, 2]);
      for my $p (@pts) { my $q = RefHodge::xi($g, $om, $p);
        for my $s (map { Q($_) } (-1, 0, 1, 2, -2, 3, '1/2', 5, -3, 4, 6, 7, -4)) { my $mine = RefHodge::xi_value($q, $s) / ($c ** (3*$g - 3 + 2*$r));
          my %v = (x => $x, s => $s); $v{"j$_"} = Q($p->[$_ - 1]) for 1..$r; my $d = RefParse::evalx($e, \%v); $n++; $bad++ if $mine != $d; } } }
    printf "genus %d q table r = %d: %d comparisons (13 values of s each), %d mismatches [%.1f s]\n", $g, $r, $n, $bad, time - $t0; $ntot += $n; $nbad += $bad; }
}
printf "TOTAL: %d comparisons, %d mismatches; memo %d; %.1f s\n", $ntot, $nbad, scalar(keys %RefHodge::IM), time - $t0;
