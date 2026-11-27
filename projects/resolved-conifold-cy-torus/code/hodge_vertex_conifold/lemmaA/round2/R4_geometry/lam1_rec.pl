# lam1_rec.pl hmax : the lambda_1-component of the one-leaf function for all genera up to hmax, WITHOUT any Hodge engine:
#  (1) Hurwitz one-leaf X_g(j) = Xi^(0)_g(j) = (j-1)_(2g-2) [t^{2g}] S^2 H^j  (Route B, Theorem B);
#  (2) B_g(j) := d/dy Xi^(0)_g(j,y)|_{y=0} from Route B's (REC) for two leaves (j,y), differentiated at y = 0 (report, Prop. M2):
#      (2h+1)(j-2h) B_h = -(2h+1)(j-2h+1) X_h + Sum_{k=2}^{h} c_k [ (x)_(2k)'|_{x=j-2h''} (j-2h''+1) X_{h''} + (j-2h'')_(2k) B_{h''} ]
#                         + c_{h+1} d/dj (j-1)_(2h+1) + Sum_{h1+k1+h2+k2 = h+2} c_{k1} c_{k2} psi_{h1,k1}(j) phi'_{h2,k2}(0),
#      h'' = h+1-k, c_k = 2/(2k)!, psi_{h1,k1}(j) = (j-2h1+1)_(2k1) X_{h1}(j), phi_{h2,k2}(y) = (y-2h2+1)_(2k2) X_{h2}(y)
#      (genus 0: (x+1)_(2k) X_0(x) = (x-1)_(2k-2));
#  (3) p_1 = (2h-1) X_h - B_h  and the Mumford reduction (Prop. M1)
#      12 Xi^(1)_h(j) = p_1 - (X_h(j) - X_h(0))/j + (1/2)(j-2h+2)(j-2h+3) X_{h-1}(j) + Sum_{g1=1}^{h-1} (j-2g1+1) X_{h-g1}(0) X_{g1}(j);
#  (4) comparison with the closed form predicted by (6.1):  Xi^(1)_h(j) = -(1/4) (j-1)_(2h-2) [t^{2h}] ( G Phi + (t/2) G' ),
#      G = S^{2+j} exp(-j t S'/S), Phi = (1 + j/2) t^2 f'' - (j/6) t^3 f''', f = log S;
#  (5) the vanishing Xi^(1)_h(j) = 0 for 2 <= j <= 2h-2;  (6) comparison with the engine tables where available.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use PJ;
sub Q { Math::BigRat->new(shift) }
my $hmax = shift // 8;
$PJ::N = 2*$hmax + 4;
my $f = PJ::logS(); my $tf1 = PJ::stder($f);                                  # t f'
my $logH = PJ::sadd($f, $tf1, Q(-1));                                          # log H = f - t f'
my $GX = PJ::sexp(PJ::sadd(PJ::sscalp($f, PJ::pc(2)), PJ::sscalp($logH, PJ::pj())));   # S^2 H^j
my %X; for my $g (1..$hmax) { $X{$g} = PJ::pmul(PJ::ffj(-1, 2*$g - 2), $GX->[2*$g]); }
my %X0; $X0{$_} = PJ::peval($X{$_}, Q(0)) for 1..$hmax;
sub c { my $k = shift; return Q(2) / PJ::fact(2*$k) }
sub psi { my ($h1, $k1) = @_; return $h1 == 0 ? PJ::ffj(-1, 2*$k1 - 2) : PJ::pmul(PJ::ffj(-2*$h1 + 1, 2*$k1), $X{$h1}) }
sub phid0 { my ($h2, $k2) = @_; my $p = psi($h2, $k2); return $p->[1] // Q(0) }   # d/dy at 0 of the same polynomial in y
my %B;
for my $h (1..$hmax) {
  my $rhs = PJ::pscale(PJ::pmul([Q(-2*$h + 1), Q(1)], $X{$h}), Q(-(2*$h + 1)));
  for my $k (2..$h) { my $hp = $h + 1 - $k;
    $rhs = PJ::padd($rhs, PJ::pscale(PJ::pmul(PJ::dffj(-2*$hp, 2*$k), PJ::pmul([Q(-2*$hp + 1), Q(1)], $X{$hp})), c($k)));
    $rhs = PJ::padd($rhs, PJ::pscale(PJ::pmul(PJ::ffj(-2*$hp, 2*$k), $B{$hp}), c($k))); }
  $rhs = PJ::padd($rhs, PJ::pscale(PJ::dffj(-1, 2*$h + 1), c($h + 1)));
  for my $h1 (0..$h+1) { for my $k1 (1..$h+1) { for my $h2 (0..$h+1) { my $k2 = $h + 2 - $h1 - $k1 - $h2; next if $k2 < 1;
        my $d = phid0($h2, $k2); next if $d == 0; $rhs = PJ::padd($rhs, PJ::pscale(psi($h1, $k1), c($k1) * c($k2) * $d)); } } }
  my ($q, $rem) = PJ::pdiv_lin($rhs, 2*$h); die "genus $h: RHS not divisible by (j - 2h): remainder $rem" if $rem != 0;
  $B{$h} = PJ::pscale($q, Q(1) / (2*$h + 1)); }
# validation of B against the two-leaf engine tables (nr_test/engine_h<h>_r2.txt, k = 0, coefficient of y^1)
for my $h (1..$hmax) { my $fn = "$FindBin::Bin/nr_test/engine_h${h}_r2.txt"; next unless -e $fn; open my $fh, '<', $fn or die; my @b;
  while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 0; for my $t (split /\s*;\s*/, $rest) { my ($m, $v) = split /:/, $t; my ($a, $bb) = split /,/, $m; $b[$a] = Q($v) if $bb == 1; } } close $fh;
  my $bad = 0; my $n = $#b > $#{$B{$h}} ? $#b : $#{$B{$h}}; for my $a (0..$n) { $bad++ if ($b[$a] // Q(0)) != ($B{$h}[$a] // Q(0)); }
  printf "genus %d: B_h from REC vs two-leaf engine table: %s\n", $h, ($bad ? "$bad MISMATCHES" : "equal"); }
# lambda_1 component via Prop. M1
my %X1;
for my $h (1..$hmax) {
  my $p1 = PJ::psub(PJ::pscale($X{$h}, Q(2*$h - 1)), $B{$h});
  my ($qq, $rr) = PJ::pdiv_lin(PJ::psub($X{$h}, PJ::pc($X0{$h})), 0); die "not divisible by j" if $rr != 0;
  my $t = PJ::psub($p1, $qq);
  $t = PJ::padd($t, PJ::pscale(PJ::pmul(PJ::pmul([Q(-2*$h + 2), Q(1)], [Q(-2*$h + 3), Q(1)]), $X{$h-1}), Q(1)/2)) if $h >= 2;
  $t = PJ::padd($t, PJ::pc(Q(1)/2)) if $h == 1;   # genus 1: delta_irr term = (1/2) Int_{Mbar_{0,3}} 1
  for my $g1 (1..$h-1) { $t = PJ::padd($t, PJ::pscale(PJ::pmul([Q(-2*$g1 + 1), Q(1)], $X{$g1}), $X0{$h - $g1})); }
  $X1{$h} = PJ::pscale($t, Q(1)/12); }
# closed form from (6.1)
my $S = PJ::sexp($f); my $f2 = PJ::sderiv(PJ::sderiv($f)); my $f3 = PJ::sderiv($f2);
my $G = PJ::sexp(PJ::sadd(PJ::sscalp($f, [Q(2), Q(1)]), PJ::sscalp($tf1, [Q(0), Q(-1)])));   # S^{2+j} e^{-j t f'}
my $Phi = PJ::sadd(PJ::sscalp(PJ::sshift($f2, 2), [Q(1), Q(1)/2]), PJ::sscalp(PJ::sshift($f3, 3), [Q(0), Q(-1)/6]));
my $Wser = PJ::sadd(PJ::smul($G, $Phi), PJ::sscalp(PJ::stder($G), PJ::pc(Q(1)/2)));
for my $h (1..$hmax) {
  my $cf = PJ::pscale(PJ::pmul(PJ::ffj(-1, 2*$h - 2), $Wser->[2*$h]), Q(-1)/4);
  my $same = PJ::iszero(PJ::psub($cf, $X1{$h}));
  my @vals = map { PJ::peval($X1{$h}, Q($_)) } 2..(2*$h - 2); my $nz = grep { $_ != 0 } @vals;
  printf "genus %2d: Xi^(1) from REC+Mumford %s the (6.1) closed form; in-range values j=2..%d: %s\n", $h, ($same ? "EQUALS" : "DIFFERS FROM"), 2*$h-2, ($nz ? "$nz NON-ZERO" : "all zero"); }
# engine comparison for the lambda_1 component
for my $h (1..$hmax) { my $fn = "$FindBin::Bin/refdata/ref_engine_h${h}_r1.txt"; next unless -e $fn; open my $fh, '<', $fn or die; my @c;
  while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 1; for my $t (split /\s*;\s*/, $rest) { my ($a, $v) = split /:/, $t; $c[$a] = Q($v); } } close $fh;
  my $bad = 0; my $n = $#c > $#{$X1{$h}} ? $#c : $#{$X1{$h}}; for my $a (0..$n) { $bad++ if ($c[$a] // Q(0)) != ($X1{$h}[$a] // Q(0)); }
  printf "genus %d: Xi^(1) from REC+Mumford vs Hodge engine table: %s\n", $h, ($bad ? "$bad MISMATCHES" : "equal"); }
