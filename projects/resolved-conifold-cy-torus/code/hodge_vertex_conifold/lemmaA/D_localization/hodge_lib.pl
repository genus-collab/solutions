# hodge_lib.pl -- exact Hodge integrals on Mbar_{g,n} (g <= 4) in Perl + Math::BigRat (no Mathematica, no TR).
#   kw($g, @a)            = < tau_{a_1} ... tau_{a_n} >_g            (Witten-Kontsevich, DVV recursion + string/dilaton)
#   hi($g, $k1, $k3, @a)  = < lambda_1^k1 ch_3^k3 prod tau_{a_i} >_g (Mumford's GRR formula, as in hvc_lib.m hiBody;
#                           string and dilaton equations are applied first: lambda_1, ch_3 are pulled back along forgetful maps)
#   omega_class($g, @w)   = prod_{a} Lambda^vee_g(w_a)/w_a^g = prod_a sum_i (-1)^i lambda_i w_a^{-i}, as a hash "k1,k3" => coeff
#                           (lambda_2 = L^2/2, lambda_3 = L^3/6 + 2C, lambda_4 = L^4/24 + 2LC; L = lambda_1, C = ch_3)
#   vertex($g, $cls, @series)        = Int_{Mbar_{g,n}} cls prod_i f_i(psi_i), f_i = [c_0, c_1, ...] power series in psi_i
#   integral_groups($g,$k1,$k3,@grp) = Int L^k1 C^k3 prod over groups [series, mult]: identical series grouped (multisets)
# ASCII only.
use strict; use warnings;
use Math::BigRat;
our (%KW, %HI);
sub Q { Math::BigRat->new(@_) }
sub dfact { my $n = shift; my $r = Q(1); while ($n > 1) { $r *= $n; $n -= 2; } return $r; }   # n!! with (-1)!! = 1
sub fact { my $n = shift; my $r = Q(1); for my $i (2..$n) { $r *= $i } return $r; }
sub binom { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; my $r = Q(1); for my $i (1..$k) { $r *= ($n - $k + $i); $r /= $i; } return $r; }

# ---------------- Witten-Kontsevich ----------------
sub kw {
  my ($g, @a) = @_;
  @a = sort { $a <=> $b } @a;
  my $n = scalar @a;
  return Q(0) if $g < 0;
  return Q(0) if 2*$g - 2 + $n <= 0;
  my $s = 0; $s += $_ for @a;
  return Q(0) if $s != 3*$g - 3 + $n;
  for (@a) { return Q(0) if $_ < 0 }
  my $key = "$g|" . join(",", @a);
  return $KW{$key} if exists $KW{$key};
  my $res;
  if ($g == 0 && $n == 3) { $res = Q(1); }
  elsif ($g == 1 && $n == 1) { $res = Q(1, 24); }
  elsif ($a[0] == 0) {             # string equation
    my @rest = @a[1..$#a]; $res = Q(0);
    for my $j (0..$#rest) { next if $rest[$j] == 0; my @b = @rest; $b[$j]--; $res += kw($g, @b); }
  } elsif ($a[0] == 1) {           # dilaton equation
    my @rest = @a[1..$#a]; $res = Q(2*$g - 2 + scalar(@rest)) * kw($g, @rest);
  } else {                         # DVV on the largest index: a_max = k+1
    my @rest = @a[0..$#a-1]; my $k = $a[-1] - 1; my $m = scalar @rest;
    my $tot = Q(0);
    for my $j (0..$m-1) {
      my @b = @rest; my $dj = $b[$j]; $b[$j] += $k;
      $tot += dfact(2*$k + 2*$dj + 1) / dfact(2*$dj - 1) * kw($g, @b);
    }
    for my $r (0..$k-1) {
      my $s2 = $k - 1 - $r; my $cf = dfact(2*$r + 1) * dfact(2*$s2 + 1) / 2;
      $tot += $cf * kw($g - 1, $r, $s2, @rest);
      for my $mask (0 .. (1 << $m) - 1) {
        my (@I, @J); for my $i (0..$m-1) { if (($mask >> $i) & 1) { push @I, $rest[$i] } else { push @J, $rest[$i] } }
        for my $g1 (0..$g) { $tot += $cf * kw($g1, $r, @I) * kw($g - $g1, $s2, @J); }
      }
    }
    $res = $tot / dfact(2*$k + 3);
  }
  $KW{$key} = $res; return $res;
}

# ---------------- Mumford: lambda_1^k1 ch_3^k3 ----------------
sub dmax { my $g = shift; return 0 if $g <= 0; return 1 if $g == 1; my $x = $g*($g+1)/2; my $y = 3*$g - 3; return $x < $y ? $x : $y; }
sub subset_sum {
  my ($aref, $f) = @_; my @a = @$aref; my $n = scalar @a; my $tot = Q(0);
  for my $mask (0 .. (1 << $n) - 1) {
    my (@S, @T); for my $i (0..$n-1) { if (($mask >> $i) & 1) { push @S, $a[$i] } else { push @T, $a[$i] } }
    $tot += $f->(\@S, \@T);
  }
  return $tot;
}
sub hi {
  my ($g, $k1, $k3, @a) = @_;
  @a = sort { $a <=> $b } @a;
  my $n = scalar @a;
  return Q(0) if $g < 0 || $k1 < 0 || $k3 < 0;
  return Q(0) if 2*$g - 2 + $n <= 0;
  my $s = 0; $s += $_ for @a;
  return Q(0) if $k1 + 3*$k3 + $s != 3*$g - 3 + $n;
  for (@a) { return Q(0) if $_ < 0 }
  return Q(0) if $g == 0 && ($k1 > 0 || $k3 > 0);
  return kw($g, @a) if $k1 == 0 && $k3 == 0;
  return Q(0) if $k1 + 3*$k3 > dmax($g);
  my $key = "$g|$k1|$k3|" . join(",", @a);
  return $HI{$key} if exists $HI{$key};
  my $tot;
  if ($n >= 1 && $a[0] <= 1 && 2*$g - 2 + ($n - 1) > 0) {
    my @rest = @a[1..$#a];
    if ($a[0] == 0) { $tot = Q(0); for my $j (0..$#rest) { next if $rest[$j] == 0; my @b = @rest; $b[$j]--; $tot += hi($g, $k1, $k3, @b); } }
    else { $tot = Q(2*$g - 2 + scalar(@rest)) * hi($g, $k1, $k3, @rest); }
  } elsif ($k1 > 0) {
    $tot = hi($g, $k1 - 1, $k3, @a, 2);
    for my $i (0..$n-1) { my @b = @a; $b[$i]++; $tot -= hi($g, $k1 - 1, $k3, @b); }
    $tot += hi($g - 1, $k1 - 1, $k3, @a, 0, 0) / 2;
    $tot += subset_sum(\@a, sub {
      my ($S, $T) = @_; my $t = Q(0);
      for my $h (0..$g) { for my $j1 (0..$k1-1) { for my $j3 (0..$k3) {
        my $c = binom($k1 - 1, $j1) * binom($k3, $j3); next if $c == 0;
        my $x = hi($h, $j1, $j3, @$S, 0); next if $x == 0;
        $t += $c * $x * hi($g - $h, $k1 - 1 - $j1, $k3 - $j3, @$T, 0);
      } } }
      return $t; }) / 2;
    $tot /= 12;
  } else {
    $tot = hi($g, 0, $k3 - 1, @a, 4);
    for my $i (0..$n-1) { my @b = @a; $b[$i] += 3; $tot -= hi($g, 0, $k3 - 1, @b); }
    for my $j (0..2) { $tot += Q((-1)**$j) * hi($g - 1, 0, $k3 - 1, @a, $j, 2 - $j) / 2; }
    $tot += subset_sum(\@a, sub {
      my ($S, $T) = @_; my $t = Q(0);
      for my $h (0..$g) { for my $j (0..2) { for my $j3 (0..$k3-1) {
        my $c = Q((-1)**$j) * binom($k3 - 1, $j3);
        my $x = hi($h, 0, $j3, @$S, $j); next if $x == 0;
        $t += $c * $x * hi($g - $h, 0, $k3 - 1 - $j3, @$T, 2 - $j);
      } } }
      return $t; }) / 2;
    $tot *= Q(-1, 30) / 24;    # B_4/4!
  }
  $HI{$key} = $tot; return $tot;
}

# ---------------- classes as polynomials in L = lambda_1, C = ch_3 ----------------
sub cls_mul {
  my ($A, $B, $g) = @_; my %R;
  for my $ka (keys %$A) { my ($a1, $a3) = split /,/, $ka;
    for my $kb (keys %$B) { my ($b1, $b3) = split /,/, $kb;
      my ($c1, $c3) = ($a1 + $b1, $a3 + $b3); next if $c1 + 3*$c3 > dmax($g);
      my $k = "$c1,$c3"; $R{$k} = (exists $R{$k} ? $R{$k} : Q(0)) + $A->{$ka} * $B->{$kb}; } }
  return \%R;
}
sub cls_add { my ($A, $B, $cb) = @_; my %R = %$A; $cb = Q(1) unless defined $cb;
  for my $k (keys %$B) { $R{$k} = (exists $R{$k} ? $R{$k} : Q(0)) + $cb * $B->{$k}; } return \%R; }
sub lambda_i {   # lambda_i as class in genus g
  my ($i, $g) = @_;
  return { "0,0" => Q(1) } if $i == 0;
  return {} if $i > $g;
  return { "1,0" => Q(1) } if $i == 1;
  return { "2,0" => Q(1, 2) } if $i == 2;
  return { "3,0" => Q(1, 6), "0,1" => Q(2) } if $i == 3;
  return { "4,0" => Q(1, 24), "1,1" => Q(2) } if $i == 4;
  die "lambda_$i not implemented";
}
sub cvee { my ($u, $g) = @_; my %R;   # c_u(E^vee) = sum_i (-1)^i lambda_i u^i
  for my $i (0..$g) { my $l = lambda_i($i, $g); my $c = Q((-1)**$i) * $u**$i;
    for my $k (keys %$l) { $R{$k} = (exists $R{$k} ? $R{$k} : Q(0)) + $c * $l->{$k}; } }
  return \%R; }
sub omega_class { my ($g, @w) = @_; my $R = { "0,0" => Q(1) };
  for my $w (@w) { $R = cls_mul($R, cvee(Q(1) / $w, $g), $g); } return $R; }
sub cvee_prod { my ($g, @u) = @_; my $R = { "0,0" => Q(1) };
  for my $u (@u) { $R = cls_mul($R, cvee($u, $g), $g); } return $R; }

# ---------------- integrals ----------------
sub tuples_sum {
  my ($lens, $target) = @_; my @out; my $n = scalar @$lens;
  my $rec; $rec = sub { my ($i, $left, $cur) = @_;
    if ($i == $n) { push @out, [@$cur] if $left == 0; return; }
    my $max = $lens->[$i] - 1; $max = $left if $left < $max;
    for my $a (0..$max) { $rec->($i + 1, $left - $a, [@$cur, $a]); } };
  $rec->(0, $target, []); return @out;
}
sub vertex {
  my ($g, $cls, @f) = @_;
  my $n = scalar @f; my $dim = 3*$g - 3 + $n; my $tot = Q(0);
  return $tot if 2*$g - 2 + $n <= 0;
  my @lens = map { scalar @$_ } @f;
  for my $k (keys %$cls) { my ($k1, $k3) = split /,/, $k; my $c = $cls->{$k}; next if $c == 0;
    my $t = $dim - $k1 - 3*$k3; next if $t < 0;
    for my $tp (tuples_sum(\@lens, $t)) {
      my $p = $c; for my $i (0..$n-1) { $p *= $f[$i][$tp->[$i]]; last if $p == 0; } next if $p == 0;
      $tot += $p * hi($g, $k1, $k3, @$tp);
    } }
  return $tot;
}
# groups: list of [series_ref, multiplicity]; identical series in a group are summed over multisets of exponents
sub integral_groups {
  my ($g, $k1, $k3, @grp) = @_;
  my $n = 0; $n += $_->[1] for @grp;
  return Q(0) if 2*$g - 2 + $n <= 0;
  my $t = 3*$g - 3 + $n - $k1 - 3*$k3; return Q(0) if $t < 0;
  my $tot = Q(0);
  my $rec; $rec = sub {
    my ($gi, $left, $exps, $wt) = @_;
    if ($gi == scalar @grp) { return unless $left == 0; $tot += $wt * hi($g, $k1, $k3, @$exps); return; }
    my ($s, $m) = @{$grp[$gi]}; my $maxa = scalar(@$s) - 1;
    # nonincreasing sequences of length m with entries <= maxa, sum <= left
    my $inner; $inner = sub {
      my ($k, $cap, $lft, $cur, $w) = @_;
      if ($k == $m) {
        my %mult; $mult{$_}++ for @$cur; my $perm = fact($m); $perm /= fact($_) for values %mult;
        $rec->($gi + 1, $lft, [@$exps, @$cur], $wt * $w * $perm); return; }
      my $top = $cap < $lft ? $cap : $lft;
      for my $a (reverse 0..$top) { my $c = $s->[$a]; next if $c == 0; $inner->($k + 1, $a, $lft - $a, [@$cur, $a], $w * $c); }
    };
    $inner->(0, $maxa, $left, [], Q(1));
  };
  $rec->(0, $t, [], Q(1));
  return $tot;
}
# power series helpers (truncated at degree N)
sub ser_geom { my ($scale, $N) = @_; my @s; my $p = Q(1); for my $i (0..$N) { push @s, $p; $p *= $scale; } return \@s; }   # 1/(1 - scale psi)
sub ser_T { my ($c, $N) = @_; my @s = (Q(0), Q(0)); for my $a (2..$N) { push @s, -(Q(1) / $c**($a - 1)); } return \@s; }  # -psi^2/(c - psi)
sub ser_const { my ($x, $N) = @_; my @s = (Q($x)); push @s, Q(0) for 1..$N; return \@s; }
1;
