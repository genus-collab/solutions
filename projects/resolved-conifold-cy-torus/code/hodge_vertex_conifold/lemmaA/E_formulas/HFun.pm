package HFun;
# The symmetric two-leaf correction series H(j1,j2;lambda) of the shortened-Norlund structure:
#   sigma^h Xi_h(j1,j2) = prod_{l=1}^{2h-1}(J-l) [lambda^{2h}] B^J H,    [lambda^{2n}] H = N_n(J,P) / prod_{l=1}^{2n-1}(J-l),
# built from the slice series C_s(m) (slice_m*.txt, m = 0..16, interpolated in m): H = B^{-j1} Sum_s C_s(j1)/((J-1)...(J-s)).
use strict; use warnings; use Math::BigRat; use Poly; use SerS;
sub R { Math::BigRat->new(shift) }
sub interp_poly { my ($xs, $ys, $var) = @_; my $N = scalar @$xs; my @dd = @$ys; for my $k (1..$N-1) { for (my $i = $N-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($xs->[$i] - $xs->[$i-$k]); } }
  my $p = {}; my $b = Poly::const(1); for my $k (0..$N-1) { $p = Poly::add($p, Poly::scale($b, $dd[$k])) if $dd[$k] != 0; $b = Poly::mul($b, Poly::sub_(Poly::var($var), Poly::const($xs->[$k]))); } return $p; }
sub to_JP { my $p = shift; # symmetric poly in j1, j2 (and s) -> poly in J, P
  my $res = {}; my $cur = { %$p };
  while (!Poly::is_zero($cur)) { my ($best, $ba, $bb);
    for my $k (keys %$cur) { my $e = Poly::mono_parse($k); my $a = $e->{j1} // 0; my $b = $e->{j2} // 0; next if $a < $b; my $d = $a + $b; if (!defined $best || $d > $ba + $bb || ($d == $ba + $bb && $a > $ba)) { ($best, $ba, $bb) = ($k, $a, $b); } }
    die "not symmetric" unless defined $best;
    my $c = Poly::coeff_of(Poly::coeff_of({ $best => $cur->{$best} }, 'j1', $ba), 'j2', $bb);
    $res = Poly::add($res, Poly::mul($c, Poly::mul(Poly::pow(Poly::var('P'), $bb), Poly::pow(Poly::var('J'), $ba - $bb))));
    $cur = Poly::sub_($cur, Poly::mul($c, Poly::mul(Poly::pow(Poly::mul(Poly::var('j1'), Poly::var('j2')), $bb), Poly::pow(Poly::add(Poly::var('j1'), Poly::var('j2')), $ba - $bb)))); }
  return $res; }
sub numerators { my ($dir, $nmax) = @_;
  my %C; for my $m (1..16) { open my $fh, '<', "$dir/slice_m$m.txt" or next;
    while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = R($1); } } }
  my $S12 = SerS::S1S2(); for my $n (0..5) { my $c = SerS::coeff($S12, $n); for my $p (0..$#$c) { $C{0}{0}{$n}{$p} = $c->[$p] if $c->[$p] != 0; } }
  my @ms = (0..16); my %Cp;
  for my $s (0..2*$nmax) { for my $n (0..$nmax) { my $tot = {};
    for my $p (0..$n+1) { my @y = map { $C{$_}{$s}{$n}{$p} // R(0) } @ms; next unless grep { $_ != 0 } @y;
      $tot = Poly::add($tot, Poly::mul(interp_poly([map { R($_) } @ms], \@y, 'm'), Poly::pow(Poly::var('s'), $p))); }
    $Cp{$s}{$n} = $tot; } }
  my $bet = SerS::beta(); my @Bm = (Poly::const(1)); { my @pw = (SerS::sone()); for my $k (1..$nmax) { push @pw, SerS::smul($pw[-1], $bet); }
    for my $n (1..$nmax) { my $t = {}; for my $k (1..$n) { my $c = SerS::coeff($pw[$k], $n); my $f = R(1); $f *= $_ for 2..$k;
        my $cp = {}; for my $i (0..$#$c) { $cp = Poly::add($cp, Poly::scale(Poly::pow(Poly::var('s'), $i), $c->[$i])) if $c->[$i] != 0; }
        $t = Poly::add($t, Poly::mul($cp, Poly::scale(Poly::pow(Poly::var('m'), $k), (($k % 2) ? R(-1) : R(1)) / $f))); } push @Bm, $t; } }
  my $Jv = Poly::var('J'); my @N;
  for my $n (0..$nmax) { my $num = {}; my $L = 2*$n - 1; $L = 0 if $L < 0;
    for my $s (0..2*$nmax) { for my $a (0..$n) { my $c = $Cp{$s}{$a} // next; next if Poly::is_zero($c); my $t = Poly::mul($c, $Bm[$n - $a]); next if $s > $L;
        for my $l ($s+1 .. $L) { $t = Poly::mul($t, Poly::sub_($Jv, Poly::const($l))); } $num = Poly::add($num, $t); } }
    my $sub = { m => Poly::var('j1'), J => Poly::add(Poly::var('j1'), Poly::var('j2')) };
    my $N2 = {}; for my $k (keys %$num) { my $e = Poly::mono_parse($k); my $t = Poly::const($num->{$k}); for my $v (keys %$e) { my $f = exists $sub->{$v} ? $sub->{$v} : Poly::var($v); $t = Poly::mul($t, Poly::pow($f, $e->{$v})); } $N2 = Poly::add($N2, $t); }
    push @N, to_JP($N2); }
  return \@N; }
1;
