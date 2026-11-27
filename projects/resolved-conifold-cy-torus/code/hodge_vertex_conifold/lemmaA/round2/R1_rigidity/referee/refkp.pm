package refkp;
# Referee's independent implementation (written from scratch, not copied from the route's kpmod.pm / kp_point.pl).
# Arithmetic: mode 'mod' (prime P, default 1000000007) or 'rat' (Math::BigRat).
#
# KP propagation.  F_G(P) = Sum_T Y(T) P^T/|Aut T| on the slice P_1 = 0 (T multisets of parts >= 2), normalised
# coefficients <f>_R := |Aut R| [P^R] f, so <f g>_R = Sum_{S subset R} binom(R,S) <f>_S <g>_{R-S}.
# Genus 0: g_k(P) := d_1 d_k F_0(P), <g_k>_S = sigma Pi_k Pi_S Xi_0(S+{1,k}), Xi_0(T) = Gamma(J)/Gamma(J-r+3).
# Linearised dispersionless Hirota at every point P (Route A Theorem D for all R), in my own form:
#   (Y_P(w) - Y_P(z)) Sum_{a,b>=1} G_ab z^-a w^-b = V(z) - V(w),   Y_P(z) = z - Sum_k g_k(P) z^-k,
# which is equivalent to   G_{1,n} = G_{n,1} = V_n   and, for m,n >= 1,
#   G_{m+1,n} = G_{m,n+1} + Sum_k g_k (G_{m-k,n} - G_{m,n-k})      (anti-diagonal sweep).
# Here G_ab = d_a d_b F_G, V_n = d_1 d_n F_G; string Y(T+{1}) = J_T Y(T), so <V_n>_R = (J_R+n) Y(R+{n}).
# Unknowns: one-leaf values y_L = Y({L+1}), L = 1..K; every Y(T) is a linear form in them.
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(set_mode set_prime NN FR AD SB ML IV NG IZ EQ ONE ZERO
                 vzero vunit vadd vaxpy vscale vsub viszero veq
                 Pi_of fact_f binom_gen gcd parts_ge2 mkey aut_f subms
                 kp_propagate ser_mul ser_inv ser_comp ser_log ser_exp hyp2f1 Tser);
our $MODE = 'mod'; our $P = 1000000007;
sub set_mode { $MODE = shift; }
sub set_prime { $P = shift; }
sub NN { my $n = shift; return Math::BigRat->new("$n") if $MODE eq 'rat'; $n %= $P; $n += $P if $n < 0; return $n; }
sub AD { return $MODE eq 'rat' ? $_[0] + $_[1] : ($_[0] + $_[1]) % $P; }
sub SB { return $MODE eq 'rat' ? $_[0] - $_[1] : ($_[0] - $_[1]) % $P; }
sub NG { return $MODE eq 'rat' ? -$_[0] : (-$_[0]) % $P; }
sub ML { return $MODE eq 'rat' ? $_[0] * $_[1] : ($_[0] * $_[1]) % $P; }
sub IZ { return $MODE eq 'rat' ? $_[0]->is_zero : ($_[0] % $P == 0); }
sub EQ { return $MODE eq 'rat' ? ($_[0] == $_[1]) : (($_[0] - $_[1]) % $P == 0); }
sub ONE { return NN(1); } sub ZERO { return NN(0); }
sub powmod { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub IV { my $a = shift; if ($MODE eq 'rat') { die "IV: zero" if $a->is_zero; return Math::BigRat->new(1) / $a; } die "IV: zero mod p" if $a % $P == 0; return powmod($a, $P - 2); }
sub FR { my ($a, $b) = @_; $b //= 1; return Math::BigRat->new("$a/$b") if $MODE eq 'rat'; return ML(NN($a), IV(NN($b))); }
# vectors
sub vzero { my $n = shift; return [ map { ZERO() } 1..$n ]; }
sub vunit { my ($n, $i) = @_; my $v = vzero($n); $v->[$i] = ONE(); return $v; }
sub vadd { my ($a, $b) = @_; return [ map { AD($a->[$_], $b->[$_]) } 0..$#$a ]; }
sub vsub { my ($a, $b) = @_; return [ map { SB($a->[$_], $b->[$_]) } 0..$#$a ]; }
sub vscale { my ($c, $a) = @_; return [ map { ML($c, $_) } @$a ]; }
sub vaxpy { my ($y, $c, $x) = @_; return if IZ($c); for my $i (0..$#$y) { $y->[$i] = AD($y->[$i], ML($c, $x->[$i])); } }
sub viszero { my $a = shift; for (@$a) { return 0 unless IZ($_); } return 1; }
sub veq { my ($a, $b) = @_; for my $i (0..$#$a) { return 0 unless EQ($a->[$i], $b->[$i]); } return 1; }
# combinatorics
my %FC; sub fact_f { my $n = shift; return $FC{"$MODE$P,$n"} //= do { my $f = ONE(); $f = ML($f, NN($_)) for 2..$n; $f }; }
sub binom_gen { my ($x, $m) = @_; return ZERO() if $m < 0; my $v = ONE(); for my $i (0..$m-1) { $v = ML($v, SB($x, NN($i))); } return ML($v, IV(fact_f($m))); }
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
sub Pi_of { my ($j, $tau) = @_; return ONE() if $j <= 1; my $v = ONE(); for my $b (1..$j-1) { $v = ML($v, AD(ML(NN($j), $tau), NN($b))); } return ML($v, IV(fact_f($j-1))); }
# multisets of parts >= 2 (descending) with level D = Sum (t-1)
sub parts_ge2 { my ($D, $maxpart) = @_; return ([]) if $D == 0; my @o; my $top = $maxpart < $D + 1 ? $maxpart : $D + 1;
  for my $t (reverse 2..$top) { push @o, [$t, @$_] for parts_ge2($D - ($t - 1), $t); } return @o; }
sub mkey { return join(',', sort { $b <=> $a } @_); }
sub aut_f { my %m; $m{$_}++ for @_; my $a = ONE(); $a = ML($a, fact_f($_)) for values %m; return $a; }
# sub-multisets S of R with weight binom(R,S) (as integers); returns list of [S, R-S, w]
my %SMC;
sub subms { my @R = @_; my $k = join(',', @R); return @{ $SMC{$k} //= do {
  my %m; $m{$_}++ for @R; my @out = ([[], [], 1]);
  for my $v (sort { $b <=> $a } keys %m) { my $n = $m{$v}; my @new;
    for my $o (@out) { for my $i (0..$n) { my $bn = 1; for my $u (1..$i) { $bn = $bn * ($n - $u + 1) / $u; }
        push @new, [ [@{$o->[0]}, ($v) x $i], [@{$o->[1]}, ($v) x ($n - $i)], $o->[2] * $bn ]; } }
    @out = @new; }
  \@out } }; }

# ---------------------------------------------------------------------------------------------------------
# kp_propagate(tau, K): returns { Y => {key => vector}, G => {Rkey => table}, stats }
# vector index L-1 <-> y_L (L = 1..K).
sub kp_propagate {
  my ($tau, $K, %opt) = @_;
  my $sig = ML($tau, AD($tau, ONE()));
  my %PiC; my $Pi = sub { my $j = shift; return $PiC{$j} //= Pi_of($j, $tau); };
  # genus-0 coefficient <g_k>_S
  my %gC; my $g = sub { my ($k, @S) = @_; my $kk = "$k;" . join(',', @S); return $gC{$kk} //= do {
      my $s = 0; $s += $_ for @S; my $m = @S; my $v = ML($sig, $Pi->($k));
      $v = ML($v, $Pi->($_)) for @S;
      # Xi_0(S+{1,k}) = (s+k)!/(s+k+1-m)!
      my $num = fact_f($s + $k); my $den = fact_f($s + $k + 1 - $m);
      ML($v, ML($num, IV($den))) }; };
  my %Y; my %Gt; my $nchk = 0; my $nbad = 0;
  for my $L (1..$K) { $Y{$L+1} = vunit($K, $L - 1); }
  # all R with level <= K, grouped by number of parts
  my @allR; for my $D (0..$K) { push @allR, parts_ge2($D, $D + 1); }
  my %byp; for my $R (@allR) { push @{$byp{scalar @$R}}, $R; }
  my $Ylookup = sub { my @T = @_; my @Tp = sort { $b <=> $a } grep { $_ >= 2 } @T; my $nu = @T - @Tp;
    return vzero($K) unless @Tp;                       # Y(1^n) = 0 (string), Y(empty) only multiplied by 0 below
    my $v = $Y{join(',', @Tp)}; die "Ylookup: missing (@Tp)" unless defined $v;
    my $J = 0; $J += $_ for @Tp; my $f = ONE(); for my $i (0..$nu-1) { $f = ML($f, NN($J + $i)); }
    return vscale($f, $v); };
  for my $p (sort { $a <=> $b } keys %byp) {
    for my $R (sort { join(',', @$a) cmp join(',', @$b) } @{$byp{$p}}) {
      my $DR = 0; $DR += $_ - 1 for @$R; my $JR = 0; $JR += $_ for @$R;
      my $smax = $K + 2 - $DR; my $Rk = join(',', @$R);
      my @subs = grep { @{$_->[0]} > 0 } subms(@$R);   # nonempty S
      my $T = []; # T->[m][n]
      my $Vn = sub { my $n = shift;
        if ($n >= 2) { return vscale(NN($JR + $n), $Ylookup->(@$R, $n)); }
        # n = 1: <V_1>_R = Y(R+{1,1}) = (J_R+1) J_R Y(R)
        return vzero($K) if $JR == 0;
        return vscale(ML(NN($JR + 1), NN($JR)), $Ylookup->(@$R)); };
      my $Gget = sub { my ($tab, $a, $b) = @_; return undef if $a <= 0 || $b <= 0; return $tab->[$a][$b]; };
      for my $s (2..$smax) {
        $T->[1][$s-1] = $Vn->($s - 1);
        for my $m (1..$s-2) { my $n = $s - $m - 1;
          my $val = [ @{$T->[$m][$n+1]} ];
          for my $k (1..$s) {
            # S = empty term: <g_k>_empty (G_R[m-k][n] - G_R[m][n-k])
            my $ge = $g->($k);
            my $x1 = $Gget->($T, $m - $k, $n); my $x2 = $Gget->($T, $m, $n - $k);
            vaxpy($val, $ge, $x1) if $x1; vaxpy($val, NG($ge), $x2) if $x2;
            for my $sb (@subs) { my ($S, $C, $w) = @$sb; my $tabC = $Gt{join(',', @$C)}; die "missing table (@$C)" unless $tabC;
              my $gs = ML(NN($w), $g->($k, @$S));
              my $y1 = $Gget->($tabC, $m - $k, $n); my $y2 = $Gget->($tabC, $m, $n - $k);
              vaxpy($val, $gs, $y1) if $y1; vaxpy($val, NG($gs), $y2) if $y2; } }
          $T->[$m+1][$n] = $val; }
        # consistency: G_{s-1,1} = V_{s-1}
        $nchk++; unless (veq($T->[$s-1][1], $T->[1][$s-1])) { $nbad++; print "INCONSISTENT end of anti-diagonal R=($Rk) s=$s\n" if $nbad <= 5; }
      }
      # symmetry checks
      for my $a (1..$smax-1) { for my $b ($a+1..$smax-$a) { next unless $T->[$a][$b] && $T->[$b][$a]; $nchk++; unless (veq($T->[$a][$b], $T->[$b][$a])) { $nbad++; print "ASYMMETRIC R=($Rk) ($a,$b)\n" if $nbad <= 5; } } }
      $Gt{$Rk} = $T;
      # extract Y(R+{a,b}), a,b >= 2
      for my $a (2..$smax) { for my $b ($a..$smax-$a) { next if $a + $b > $smax; my $Tk = mkey(@$R, $a, $b); my $v = $T->[$a][$b];
          if (exists $Y{$Tk}) { $nchk++; unless (veq($Y{$Tk}, $v)) { $nbad++; print "DECOMPOSITION MISMATCH T=($Tk) from R=($Rk)\n" if $nbad <= 5; } }
          else { $Y{$Tk} = $v; } } }
    }
  }
  return { Y => \%Y, G => \%Gt, nchk => $nchk, nbad => $nbad, sig => $sig, Pi => $Pi };
}
sub k_dummy { return 0; }

# ---------------------------------------------------------------------------------------------------------
# truncated power series in t (arrays, index = power, length N+1)
sub ser_mul { my ($a, $b, $N) = @_; my @c = map { ZERO() } 0..$N;
  for my $i (0..$N) { next if $i > $#$a || IZ($a->[$i]); for my $j (0..$N-$i) { next if $j > $#$b; $c[$i+$j] = AD($c[$i+$j], ML($a->[$i], $b->[$j])); } } return \@c; }
sub ser_inv { my ($a, $N) = @_; my @b = (IV($a->[0]));
  for my $n (1..$N) { my $s = ZERO(); for my $k (1..$n) { next if $k > $#$a; $s = AD($s, ML($a->[$k], $b[$n-$k])); } $b[$n] = NG(ML($s, $b[0])); } return \@b; }
sub ser_comp { my ($f, $g, $N) = @_; my @res = map { ZERO() } 0..$N; my $pw = [ ONE(), map { ZERO() } 1..$N ];
  for my $k (0..$N) { last if $k > $#$f; for my $i (0..$N) { $res[$i] = AD($res[$i], ML($f->[$k], $pw->[$i])); } $pw = ser_mul($pw, $g, $N); } return \@res; }
sub ser_log { my ($a, $N) = @_; # a[0] = 1
  my $ia = ser_inv($a, $N); my @da = map { ($_+1 <= $#$a) ? ML(NN($_+1), $a->[$_+1]) : ZERO() } 0..$N; my $q = ser_mul(\@da, $ia, $N);
  my @b = (ZERO()); for my $n (1..$N) { $b[$n] = ML($q->[$n-1], IV(NN($n))); } return \@b; }
sub ser_exp { my ($a, $N) = @_; my @b = (ONE());
  for my $n (1..$N) { my $s = ZERO(); for my $k (1..$n) { next if $k > $#$a; $s = AD($s, ML(NN($k), ML($a->[$k], $b[$n-$k]))); } $b[$n] = ML($s, IV(NN($n))); } return \@b; }
sub hyp2f1 { my ($a, $b, $c, $N) = @_; my @s = (ONE()); my $t = ONE();
  for my $n (1..$N) { $t = ML($t, ML(ML(AD($a, NN($n-1)), AD($b, NN($n-1))), IV(ML(AD($c, NN($n-1)), NN($n))))); push @s, $t; } return \@s; }
sub Tser { my $N = shift; my @s = (ZERO()); for my $n (1..$N) { push @s, ML(NN(($n % 2) ? 1 : -1), IV(fact_f($n))); } return \@s; }
1;
