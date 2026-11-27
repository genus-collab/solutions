package RefRec;
# Referee's own implementation of the report's recursion (REC), of Theorem B, and of the pure-unit cancellation (Lemma 6.3).
# (REC): (e+1) delta Xi_h(R) = [x^h] Sum_{set partitions pi of the labelled parts} prod_{B in pi} W_B(x)  minus the excluded term,
#   W_B(x) = Sum_{hB >= 0, kB >= 1} x^{hB+kB-1} (2/(2kB)!) (delta_B)_(2kB) Xi_{hB}(R_B),  delta_B = J_B - (2hB - 2 + r_B).
# Implemented by the subset recursion F(S) = Sum_{B subset S, B contains min S} W_B F(S \ B) on bitmasks, with polynomials in x
# truncated at degree h (different code path from the author's explicit enumeration in Rec.pm).
use strict; use warnings; use Math::BigInt; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
sub BI { Math::BigInt->new(shift) }
my @FAC = (BI(1));
sub fac { my $n = shift; while (@FAC <= $n) { push @FAC, $FAC[-1] * BI(scalar @FAC) } return $FAC[$n]->copy }
sub ff { my ($x, $n) = @_; my $p = BR(1); $p *= BR($x - $_) for 0..$n - 1; return $p }   # x(x-1)...(x-n+1), n >= 0

sub xi_genus0 { my @R = @_; my $r = @R; my $J = 0; $J += $_ for @R;
  return BR(1) / BR($J * ($J + 1)) if $r == 1; return BR(1) / BR($J) if $r == 2; return ff($J - 1, $r - 3) }

my %XI;
sub xi { my ($h, @R) = @_; @R = sort { $a <=> $b } @R; my $key = "$h;" . join(",", @R);
  return $XI{$key} if exists $XI{$key};
  my $val;
  if ($h == 0) { $val = xi_genus0(@R) }
  else {
    my $r = @R; my $J = 0; $J += $_ for @R; my $e = 2 * $h - 2 + $r; my $delta = $J - $e;
    if ($delta == 0) { $val = undef }
    else {
      my $full = (1 << $r) - 1;
      # block polynomials W[mask] = [coef of x^0..x^h]
      my %W;
      for my $mask (1..$full) { my @B = map { $R[$_] } grep { $mask & (1 << $_) } 0..$r - 1; my $rB = @B; my $JB = 0; $JB += $_ for @B;
        my @w = map { BR(0) } 0..$h;
        for my $m (0..$h) { for my $hB (0..$m) { my $kB = $m - $hB + 1;
            next if $mask == $full && $hB == $h && $kB == 1;
            my $dB = $JB - (2 * $hB - 2 + $rB);
            next if $dB >= 0 && $dB <= 2 * $kB - 1;
            my $x = xi($hB, @B); die "RefRec: needed undetermined Xi_$hB(@B)" unless defined $x; next if $x->is_zero;
            $w[$m] += BR(2) / BR(fac(2 * $kB)) * ff($dB, 2 * $kB) * $x } }
        $W{$mask} = \@w }
      my %F = (0 => [BR(1), map { BR(0) } 1..$h]);
      for my $S (1..$full) { my $low = $S & (-$S); my $rest = $S ^ $low; my @acc = map { BR(0) } 0..$h;
        # subsets B of S containing the lowest bit: B = low | sub, sub subset of rest
        my $sub = $rest;
        while (1) { my $B = $low | $sub; my $wB = $W{$B}; my $fR = $F{$S ^ $B};
          for my $a (0..$h) { next if $wB->[$a]->is_zero; for my $b (0..$h - $a) { next if $fR->[$b]->is_zero; $acc[$a + $b] += $wB->[$a] * $fR->[$b] } }
          last if $sub == 0; $sub = ($sub - 1) & $rest }
        $F{$S} = \@acc }
      $val = $F{$full}[$h] / BR(($e + 1) * $delta) } }
  $XI{$key} = $val; return $val }

# ---------------------------------------------------------------- power series helpers (arrays of BigRat, truncated at N)
sub s_mul { my ($a, $b, $N) = @_; my @c = map { BR(0) } 0..$N;
  for my $i (0..$N) { next if $a->[$i]->is_zero; for my $j (0..$N - $i) { next if $b->[$j]->is_zero; $c[$i + $j] += $a->[$i] * $b->[$j] } } return \@c }
sub s_inv { my ($a, $N) = @_; die "s_inv" if $a->[0]->is_zero; my @c = (BR(1) / $a->[0]);
  for my $n (1..$N) { my $s = BR(0); for my $k (1..$n) { $s += $a->[$k] * $c[$n - $k] } push @c, -$s / $a->[0] } return \@c }
sub s_log { my ($a, $N) = @_; die "s_log" unless $a->[0] == 1;   # log a, a(0) = 1:  (log a)' = a'/a
  my @da = map { $a->[$_ + 1] * BR($_ + 1) } 0..$N - 1; push @da, BR(0); my $q = s_mul(\@da, s_inv($a, $N), $N);
  return [BR(0), map { $q->[$_ - 1] / BR($_) } 1..$N] }
sub s_exp { my ($a, $N) = @_; die "s_exp" unless $a->[0]->is_zero; my @c = (BR(1));   # c' = a' c
  for my $n (1..$N) { my $s = BR(0); for my $k (1..$n) { $s += BR($k) * $a->[$k] * $c[$n - $k] } push @c, $s / BR($n) } return \@c }

# Theorem B, form 1:  Xi_h(j) = (j-1)_(2h-2) [t^{2h}] S^2 H^j,  S = (t/2)/sinh(t/2),  H = S exp(-t S'/S)  (computed by series algebra)
sub thmB_form1 { my ($h, $j) = @_; my $N = 2 * $h;
  my @sh = map { $_ % 2 == 0 ? BR(1) / (BR(fac($_ + 1)) * BR(BI(2) ** $_)) : BR(0) } 0..$N;   # sinh(t/2)/(t/2) = sum (t/2)^{2n}/(2n+1)!
  my $S = s_inv(\@sh, $N); my $logS = s_log($S, $N);
  my @dS = map { $S->[$_ + 1] * BR($_ + 1) } 0..$N - 1; push @dS, BR(0);
  my $q = s_mul(\@dS, s_inv($S, $N), $N);                  # S'/S
  my @E = (BR(0), map { -$q->[$_ - 1] } 1..$N);             # -t S'/S
  my @ex = map { BR(2 + $j) * $logS->[$_] + BR($j) * $E[$_] } 0..$N;
  my $F = s_exp(\@ex, $N);
  my $pref = $h >= 1 ? ff($j - 1, 2 * $h - 2) : BR(1) / BR($j * ($j + 1));
  return $pref * $F->[$N] }
# Theorem B, form 2 (Lagrange):  Xi_h(j) = (j-1)_(2h-2) [y^j] x(y)^{j-2h},  x = y H(x), computed by fixed-point iteration (h >= 1)
sub thmB_form2 { my ($h, $j) = @_; my $N = 2 * $h + 1;
  my @sh = map { $_ % 2 == 0 ? BR(1) / (BR(fac($_ + 1)) * BR(BI(2) ** $_)) : BR(0) } 0..$N;
  my $S = s_inv(\@sh, $N); my @dS = map { $S->[$_ + 1] * BR($_ + 1) } 0..$N - 1; push @dS, BR(0);
  my $q = s_mul(\@dS, s_inv($S, $N), $N); my @E = (BR(0), map { -$q->[$_ - 1] } 1..$N);
  my $H = s_mul($S, s_exp(\@E, $N), $N);                    # H(t)
  # w(y) := x(y)/y satisfies w = H(y w); iterate
  my @w = (BR(1), map { BR(0) } 1..$N);
  for my $it (0..$N + 1) { my @yw = (BR(0), @w[0..$N - 1]);       # y w(y)
    my @acc = map { BR(0) } 0..$N; my @pw = (BR(1), map { BR(0) } 1..$N);
    for my $k (0..$N) { if (!$H->[$k]->is_zero) { $acc[$_] += $H->[$k] * $pw[$_] for 0..$N } my $np = s_mul(\@pw, \@yw, $N); @pw = @$np }
    @w = @acc }
  my $lw = s_log(\@w, $N); my $pw = s_exp([map { BR($j - 2 * $h) * $_ } @$lw], $N);   # w^{j-2h}
  return ff($j - 1, 2 * $h - 2) * $pw->[2 * $h] }               # [y^j] x^{j-2h} = [y^{2h}] w^{j-2h}

# ---------------------------------------------------------------- Bernoulli numbers and Lemma 6.3
my @BER;
sub bern { my $n = shift; if (!@BER) { push @BER, BR(1) } while (@BER <= $n) { my $m = scalar @BER; my $s = BR(0);
    for my $k (0..$m - 1) { $s += BR(fac($m + 1) / (fac($k) * fac($m + 1 - $k))) * $BER[$k] } push @BER, -$s / BR($m + 1) } return $BER[$n]->copy }
sub L_empty { my ($g, $k) = @_; return BR(2 * $g - 1) * bern(2 * $g) * BR(fac(2 * $g + 2 * $k - 3)) / BR(fac(2 * $g)) }
1;
