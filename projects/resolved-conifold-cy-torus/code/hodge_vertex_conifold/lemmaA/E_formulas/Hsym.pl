# Hsym.pl nmax : from the slice series C_s(m) (polynomials in m, slice_m*.txt), build
#   H(j1, j2; lambda) = B^{-j1} Sum_s C_s(j1; lambda) / ((J-1)(J-2)...(J-s)),   J = j1 + j2,
# so that sigma^h Xi_h(j1,j2) = prod_{l=1}^{2h-1}(J-l) [lambda^{2h}] B^J H  (equivalent to the slice form);
# prints the numerators N_n(j1,j2) of [lambda^{2n}] H = N_n / prod_{l=1}^{2n-1}(J-l), tests their symmetry j1 <-> j2,
# and rewrites them in J and P = j1 j2 when symmetric.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use SerS;
sub R { Math::BigRat->new(shift) }
my $nmax = shift // 3;
# ---- load C_s(m) values and interpolate in m
my %C; for my $m (1..16) { open my $fh, '<', "$FindBin::Bin/slice_m$m.txt" or next;
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = R($1); } } }
my $S12 = SerS::S1S2(); for my $n (0..5) { my $c = SerS::coeff($S12, $n); for my $p (0..$#$c) { $C{0}{0}{$n}{$p} = $c->[$p] if $c->[$p] != 0; } }
my @ms = (0..16);
sub interp_poly { my ($xs, $ys, $var) = @_; my $N = scalar @$xs; my @dd = @$ys; for my $k (1..$N-1) { for (my $i = $N-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($xs->[$i] - $xs->[$i-$k]); } }
  my $p = {}; my $b = Poly::const(1); for my $k (0..$N-1) { $p = Poly::add($p, Poly::scale($b, $dd[$k])) if $dd[$k] != 0; $b = Poly::mul($b, Poly::sub_(Poly::var($var), Poly::const($xs->[$k]))); } return $p; }
my %Cp;   # Cp{s}{n} = Poly in m, s(sigma)
for my $s (0..2*$nmax) { for my $n (0..$nmax) { my $tot = {};
  for my $p (0..$n+1) { my @y = map { $C{$_}{$s}{$n}{$p} // R(0) } @ms; next unless grep { $_ != 0 } @y;
    my $poly = interp_poly([map { R($_) } @ms], \@y, 'm'); $tot = Poly::add($tot, Poly::mul($poly, Poly::pow(Poly::var('s'), $p))); }
  $Cp{$s}{$n} = $tot; } }
# ---- B^{-m} as series with Poly(m, s) coefficients: exp(-m beta)
my $bet = SerS::beta(); my @Bm = (Poly::const(1)); { my @pw = (SerS::sone()); for my $k (1..$nmax) { push @pw, SerS::smul($pw[-1], $bet); }
  for my $n (1..$nmax) { my $t = {}; for my $k (1..$n) { my $c = SerS::coeff($pw[$k], $n); my $f = R(1); $f *= $_ for 2..$k;
      my $cp = {}; for my $i (0..$#$c) { $cp = Poly::add($cp, Poly::scale(Poly::pow(Poly::var('s'), $i), $c->[$i])) if $c->[$i] != 0; }
      $t = Poly::add($t, Poly::mul($cp, Poly::scale(Poly::pow(Poly::var('m'), $k), (($k % 2) ? R(-1) : R(1)) / $f))); } push @Bm, $t; } }
# ---- numerator of [lambda^{2n}] H over the common denominator prod_{l=1}^{2n-1}(J-l)
my $Jv = Poly::var('J');
for my $n (0..$nmax) { my $num = {}; my $L = 2*$n - 1; $L = 0 if $L < 0;
  for my $s (0..2*$nmax) { for my $a (0..$n) { my $c = $Cp{$s}{$a} // next; next if Poly::is_zero($c); my $b = $Bm[$n - $a];
      my $t = Poly::mul($c, $b); # times prod_{l=s+1}^{L}(J-l)
      next if $s > $L && !Poly::is_zero($t) && do { print "  (warning: s = $s > 2n-1 at n = $n)\n"; 1 };
      for my $l ($s+1 .. $L) { $t = Poly::mul($t, Poly::sub_($Jv, Poly::const($l))); }
      $num = Poly::add($num, $t); } }
  # substitute m -> j1, J -> j1 + j2 ; test symmetry
  my $sub = { m => Poly::var('j1'), J => Poly::add(Poly::var('j1'), Poly::var('j2')) };
  my $N2 = {}; for my $k (keys %$num) { my $e = Poly::mono_parse($k); my $t = Poly::const($num->{$k});
    for my $v (keys %$e) { my $f = exists $sub->{$v} ? $sub->{$v} : Poly::var($v); $t = Poly::mul($t, Poly::pow($f, $e->{$v})); } $N2 = Poly::add($N2, $t); }
  my $sw = {}; for my $k (keys %$N2) { my $e = Poly::mono_parse($k); my %f = %$e; my ($a, $b) = ($f{j1} // 0, $f{j2} // 0); delete $f{j1}; delete $f{j2}; $f{j1} = $b if $b; $f{j2} = $a if $a; my $nk = Poly::mono_key(\%f); $sw->{$nk} = ($sw->{$nk} // R(0)) + $N2->{$k}; }
  my $sym = Poly::equal($N2, Poly::clean($sw));
  print "[lambda^", 2*$n, "] H: numerator over prod_{l=1}^{", $L, "}(J-l) is ", ($sym ? "SYMMETRIC" : "not symmetric"), " in (j1,j2)\n";
  # express in J, P: substitute j2 = J - j1 and then j1^2 = j1 J - P iteratively -> polynomial in J, P (j1-free if symmetric)
  if ($sym) { my $q = CPsub($N2); print "   = ", Poly::to_string($q), "\n"; } }
sub CPsub { my $p = shift; # rewrite symmetric poly in j1,j2 as poly in J,P via elementary symmetric reduction
  my $res = {}; my $cur = { %$p };
  while (!Poly::is_zero($cur)) { # leading monomial j1^a j2^b with a >= b
    my ($best, $ba, $bb); for my $k (keys %$cur) { my $e = Poly::mono_parse($k); my $a = $e->{j1} // 0; my $b = $e->{j2} // 0; next if $a < $b; my $d = $a + $b; if (!defined $best || $d > $ba + $bb || ($d == $ba + $bb && $a > $ba)) { ($best, $ba, $bb) = ($k, $a, $b); } }
    my $e = Poly::mono_parse($best); my $c = Poly::coeff_of(Poly::coeff_of({ $best => $cur->{$best} }, 'j1', $ba), 'j2', $bb);   # coefficient (poly in s)
    my $term = Poly::mul($c, Poly::mul(Poly::pow(Poly::var('P'), $bb), Poly::pow(Poly::var('J'), $ba - $bb)));
    $res = Poly::add($res, $term);
    my $back = Poly::mul($c, Poly::mul(Poly::pow(Poly::mul(Poly::var('j1'), Poly::var('j2')), $bb), Poly::pow(Poly::add(Poly::var('j1'), Poly::var('j2')), $ba - $bb)));
    $cur = Poly::sub_($cur, $back); }
  return $res; }
