# t_moments.pl h r [outfile] : the components Xi^(k)_h(j_1..j_r) as exact polynomials in the j's, computed with the REFEREE engine
# (dilaton leaves as explicit points, no kappa classes, no K).  Steps:
#  (1) tau-independent moments S_mu(A) = Sum_m (-1)^m/m! Int_{Mbar_{h,r+m}} ch^mu psi^A prod_{k=1}^m Tt(psi_{r+k}), Tt = Sum_{d>=2} psi^d,
#      for every ch-monomial mu and every multiset A of leaf exponents;
#  (2) the Hodge class Omega(tau) = prod_{w in (1,tau,-1-tau)} Lambda^vee(w)/w^h at h+2 values of tau, its ch-coefficients
#      interpolated in v = -1/sigma (degree <= h; the last tau value is a check): Omega = Sum_k v^k Omega^(k);
#  (3) Xi^(k)(A) = Sum_mu Omega^(k)_mu S_mu(A), i.e. the coefficient of j^a (any ordering a of A) in Xi^(k).
# Output in the author's format  "h | r | k | a1,..,ar:coef ; ..."  (coefficient of j1^a1 ... jr^ar).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($h, $r, $out) = @ARGV; $out //= "$FindBin::Bin/ref_engine_h${h}_r${r}.txt"; my $t0 = time;
my @taus = map { Q($_) } qw(1 2 3 1/2 1/3 -1/3 3/2 5 -7/2 7 1/5 -2/7); @taus = @taus[0..$h+1];
my @vs = map { Q(-1) / ($_ * (1 + $_)) } @taus;
my $maxd = $h == 1 ? 1 : 3*$h - 3;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
sub interp { my ($v, $y, $d) = @_; my @C = map { Q(0) } 0..$d;
  for my $i (0..$d) { my @L = (Q(1)); my $den = Q(1);
    for my $m (0..$d) { next if $m == $i; my @N = map { Q(0) } 0..$#L+1; for my $k (0..$#L) { $N[$k+1] += $L[$k]; $N[$k] -= $L[$k] * $v->[$m]; } @L = @N; $den *= ($v->[$i] - $v->[$m]); }
    for my $k (0..$d) { $C[$k] += $y->[$i] * $L[$k] / $den; } }
  return @C; }
# (2) Omega components
my @om = map { RefHodge::omega($h, [Q(1), $_, -1 - $_], $maxd) } @taus;
my %mus; for my $o (@om) { $mus{$_} = 1 for keys %$o; }
my %OMK; my $badO = 0;   # OMK{mu}[k]
for my $mu (keys %mus) { my @y = map { $om[$_]{$mu} // Q(0) } 0..$#taus; my @C = interp(\@vs, [@y[0..$h]], $h);
  my $chk = Q(0); $chk += $C[$_] * $vs[$h+1] ** $_ for 0..$h; $badO++ if $chk != $y[$h+1]; $OMK{$mu} = \@C; }
printf "Omega: %d ch-monomials, v-interpolation check failures %d [%.1f s]\n", scalar(keys %mus), $badO, time - $t0;
# (1) moments and (3) components
my @comp = map { {} } 0..$h; my $nmu = 0;
for my $mu (sort { RefHodge::cdeg($a) <=> RefHodge::cdeg($b) || $a cmp $b } keys %mus) { my @mu = split /,/, $mu; my $dm = RefHodge::cdeg($mu); my $B = 3*$h - 3 + $r - $dm; next if $B < 0;
  my %S;
  for my $P (@{RefHodge::allparts($B)}) { my $m = scalar @$P; my $ex = 0; $ex += $_ for @$P; my $t = $B - $ex;
    my @D = map { $_ + 1 } @$P; my %mult; $mult{$_}++ for @D; my $den = Q(1); $den *= fact($_) for values %mult; my $w = (($m % 2) ? Q(-1) : Q(1)) / $den;
    for my $Pa (RefHodge::parts_le($t, $r)) { my @A = ((0) x ($r - scalar(@$Pa)), reverse @$Pa); my $v = RefHodge::I($h, \@mu, [@A, @D]); next if $v == 0;
      my $key = join(',', @A); my $c = $w * $v; if (exists $S{$key}) { $S{$key} += $c } else { $S{$key} = $c } } }
  for my $A (keys %S) { for my $k (0..$h) { my $c = $OMK{$mu}[$k]; next if $c == 0; my $x = $c * $S{$A}; if (exists $comp[$k]{$A}) { $comp[$k]{$A} += $x } else { $comp[$k]{$A} = $x } } }
  $nmu++; printf "  mu %s (deg %d): %d moments [%.1f s, memo %d]\n", $mu, $dm, scalar(keys %S), time - $t0, scalar(keys %RefHodge::IM) if $ENV{VERBOSE}; }
sub uniqperms { my @x = @_; my %s; my @o; my $rec; $rec = sub { my ($pre, $rest) = @_; if (!@$rest) { my $k = join(',', @$pre); push @o, [@$pre] unless $s{$k}++; return; }
    my %u; for my $i (0..$#$rest) { next if $u{$rest->[$i]}++; my @r2 = @$rest; splice(@r2, $i, 1); $rec->([@$pre, $rest->[$i]], \@r2); } }; $rec->([], \@x); return @o; }
if ($ENV{COMPACT}) {   # multiset format: one line per (k, sorted A), coefficient of j^a for EVERY ordering a of A (file ref_mset_h<h>_r<r>.txt)
  $out = "$FindBin::Bin/ref_mset_h${h}_r${r}.txt"; open my $o, '>', $out or die;
  for my $k (0..$h) { for my $A (sort keys %{$comp[$k]}) { next if $comp[$k]{$A} == 0; print $o "$h | $r | $k | $A:$comp[$k]{$A}\n"; } }
  close $o; }
else {
open my $o, '>', $out or die;
for my $k (0..$h) { my @terms; for my $A (sort keys %{$comp[$k]}) { next if $comp[$k]{$A} == 0; for my $p (uniqperms(split /,/, $A)) { push @terms, join(',', @$p) . ":" . $comp[$k]{$A}; } }
  print $o "$h | $r | $k | ", join(' ; ', @terms), "\n"; }
close $o; }
printf "genus %d r = %d: %d ch-monomials processed, Omega check failures %d; written %s; memo %d [%.1f s]\n", $h, $r, $nmu, $badO, $out, scalar(keys %RefHodge::IM), time - $t0;
