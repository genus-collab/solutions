# data_components.pl : extract the sigma-components Xi^(k)_h(R) from the MATHEMATICA q tables (lemmaA\data\q_g<g>_r<r>.txt,
# Mumford's formula, no Perl engine involved), r = 1 and r = 2, h = 2, 3, 4, and write them as exact polynomials.
#   c-normalisation: Q(s, j, x), w = (c, -x, -1), c = 1 + x;  Xi^(c=1)(tau) = c^(3h-3+2r) Q(-1, j, x),  tau = -x/c,  sigma = -x/c^2,
#   Xi^(c=1) = Sum_k v^k Xi^(k), v = -1/sigma = c^2/x.   Components by exact interpolation in v (h+1 values of x, one more as a check),
#   then as polynomials in j (r = 1: 3h-1 values of j, one more as check) or in (j1, j2) (r = 2: tensor grid).
# Output: referee/out/data_h<h>_r<r>.txt in the format "h | r | k | a[,b]:coef ; ...".
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RParse; use RL;
sub Q { Math::BigRat->new(shift) }
my $D = "$FindBin::Bin/../../../data";
mkdir "$FindBin::Bin/out";
sub interp1 { my ($xs, $ys) = @_; my $n = scalar @$xs; my @C = map { Q(0) } 0..$n-1;   # Lagrange -> monomial coefficients
  for my $i (0..$n-1) { my @L = (Q(1)); my $den = Q(1);
    for my $m (0..$n-1) { next if $m == $i; my @N = map { Q(0) } 0..scalar(@L); for my $k (0..$#L) { $N[$k+1] += $L[$k]; $N[$k] -= $L[$k] * $xs->[$m]; } @L = @N; $den *= ($xs->[$i] - $xs->[$m]); }
    for my $k (0..$n-1) { $C[$k] += $ys->[$i] * ($L[$k] // Q(0)) / $den; } }
  return @C; }
my @xv = map { Q($_) } (1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
for my $h (2, 3, 4) { for my $r (1, 2) {
  my ($g, $rr, $e) = RParse::qexpr("$D/q_g${h}_r${r}.txt"); die "header mismatch" unless $g == $h && $rr == $r;
  my $nx = $h + 2; my @xs = @xv[0..$nx-1]; my @vs = map { (1 + $_) ** 2 / $_ } @xs;
  my $comp = sub { my %pt = @_;   # components at one leaf point
    my @y = map { my $x = $_; my $val = RParse::evalx($e, { %pt, s => -1, x => $x }); $val * (1 + $x) ** (3*$h - 3 + 2*$r) } @xs;
    my @C = interp1([@vs[0..$h]], [@y[0..$h]]); my $chk = Q(0); $chk += $C[$_] * $vs[$h+1] ** $_ for 0..$h; die "v-interpolation check failed" if $chk != $y[$h+1];
    return @C; };
  my $out = "$FindBin::Bin/out/data_h${h}_r${r}.txt"; open my $o, '>', $out or die;
  if ($r == 1) { my $nj = 3*$h - 1; my @js = map { Q($_) } 0..$nj; my @vals = map { [ $comp->(j1 => $_) ] } @js;
    for my $k (0..$h) { my @C = interp1([@js[0..$nj-1]], [ map { $vals[$_][$k] } 0..$nj-1 ]); my $chk = Q(0); $chk += $C[$_] * $js[$nj] ** $_ for 0..$#C; die "j-check failed" if $chk != $vals[$nj][$k];
      print $o "$h | 1 | $k | ", RL::pstr(RL::ptrim([@C])), "\n"; } }
  else { my $nd = 3*$h + 1; my @js = map { Q($_) } 0..$nd; my %val;
    for my $a (@js) { for my $b (@js) { next if $b < $a; my @C = $comp->(j1 => $a, j2 => $b); $val{"$a,$b"} = \@C; $val{"$b,$a"} = \@C; } }
    for my $k (0..$h) { my @rows;   # for each j1, interpolate in j2; then interpolate coefficients in j1
      my @byA; for my $a (@js[0..$nd-1]) { my @C = interp1([@js[0..$nd-1]], [ map { $val{"$a,$_"}[$k] } @js[0..$nd-1] ]); push @byA, \@C; }
      my @terms; for my $bpow (0..$nd-1) { my @C = interp1([@js[0..$nd-1]], [ map { $byA[$_][$bpow] } 0..$nd-1 ]); for my $apow (0..$#C) { push @terms, "$apow,$bpow:$C[$apow]" if $C[$apow] != 0 } }
      print $o "$h | 2 | $k | ", join(' ; ', @terms), "\n"; } }
  close $o; print "genus $h r = $r: components written to $out\n"; } }
