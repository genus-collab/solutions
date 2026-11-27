package kpmod;
# kpmod.pm -- modular (p = 2^31 - 1) engine for the MV-normalised KP propagation at one framing value,
# the closed-form (polar-order) relation rows, and the regularity conditions at special framings.
# Same mathematics as kp_point.pl (exact) -- see the header there; cross-checked against it by t_kpmod.pl.
use strict; use warnings;
use Exporter 'import';
our @EXPORT = qw($P md mm mi mfrac parts_list kp_setup kp_rows set_prime);
our $P = 2147483647;
sub md { my $x = shift; $x %= $P; $x += $P if $x < 0; return $x; }
sub mm { return ($_[0] * $_[1]) % $P; }
sub mi { my $a = md(shift); die "mi: zero" unless $a; my ($t, $nt, $r, $nr) = (0, 1, $P, $a); while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q*$nt); ($r, $nr) = ($nr, $r - $q*$nr); } return md($t); }
sub mfrac { my ($n, $d) = @_; return mm(md($n), mi($d)); }
my %PC;
sub parts_list { my ($n, $max) = @_; return ([]) if $n == 0; my $k = "$n,$max"; return @{ $PC{$k} //= do { my @o; for my $k1 (reverse 1..($max < $n ? $max : $n)) { push @o, [$k1, @$_] for parts_list($n - $k1, $k1); } \@o } }; }
our %FACT; sub set_prime { $P = shift; %FACT = (); } sub fct { my $n = shift; return $FACT{$n} //= do { my $f = 1; $f = mm($f, $_) for 2..$n; $f }; }
sub binomm { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; return mm(fct($n), mi(mm(fct($k), fct($n-$k)))); }
# kp_setup(tau, D, Lmax): returns hashref with Y{key} = arrayref (linear form in y_D..y_Lmax, mod P), Pi, etc.
sub kp_setup {
  my ($tau, $D, $Lmax) = @_;   # tau mod P
  my $sig = mm($tau, md($tau + 1)); my $NU = $Lmax - $D + 1;
  my %Pi; my $Pi = sub { my $j = shift; return 1 if $j <= 1; return $Pi{$j} //= do { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } mm($v, mi(fct($j-1))) }; };
  my $X0 = sub { my @T = @_; my $J = 0; $J += $_ for @T; my $r = @T; if ($r >= 3) { my $v = 1; $v = mm($v, md($J - $_)) for 1..$r-3; return $v; } elsif ($r == 2) { return mi($J); } else { return mi(mm($J, $J+1)); } };
  my %SubC; my $subms = sub { my $S = shift; my $k = join(',', @$S); return @{ $SubC{$k} //= do {
      my %m; $m{$_}++ for @$S; my @v = sort { $b <=> $a } keys %m; my @out = ([[], [], 1]);
      for my $val (@v) { my $n = $m{$val}; my @new; for my $o (@out) { for my $i (0..$n) { push @new, [ [@{$o->[0]}, ($val) x $i], [@{$o->[1]}, ($val) x ($n-$i)], mm($o->[2], binomm($n,$i)) ]; } } @out = @new; }
      \@out } }; };
  my %E0C; my $E0; $E0 = sub { my ($m, $n, $S) = @_; my $k = "$m;$n;" . join(',', @$S); return $E0C{$k} if exists $E0C{$k};
    my $res;
    if ($m == 0 && $n == 0) { $res = @$S ? 0 : 1; } elsif ($m <= 0 || $n <= 0) { $res = 0; }
    else { $res = 0;
      for my $sb ($subms->($S)) { my ($K, $C, $w) = @$sb;
        for my $a1 (1..$m) { for my $b1 (1..$n) { my $e = $E0->($m-$a1, $n-$b1, $C); next unless $e;
            $res = ($res + mm(mm($w, mm($a1, mm($sig, mm(mm($Pi->($a1), $Pi->($b1)), $X0->(@$K, $a1, $b1))))), $e)) % $P; } } }
      $res = mm($res, mi($m)); }
    $E0C{$k} = $res; return $res; };
  my %Y; for my $L ($D..$Lmax) { my @v = (0) x $NU; $v[$L-$D] = 1; $Y{$L+1} = \@v; }
  my $lvl = sub { my $s = 0; $s += $_ - 1 for @_; return $s; };
  my $YG = sub { my @L = @_; my @Tp = sort { $b <=> $a } grep { $_ >= 2 } @L; my $nb = @L - @Tp; return undef unless @Tp;
    my $lv = $lvl->(@Tp); return undef if $lv < $D; die "level $lv > Lmax" if $lv > $Lmax;
    my $J = 0; $J += $_ for @Tp; my $f = 1; $f = mm($f, $J + $_) for 0..$nb-1;
    my $v = $Y{join(',', @Tp)}; die "undefined Y(@Tp)" unless defined $v; return [ map { mm($_, $f) } @$v ]; };
  my @allT;
  for my $L ($D..$Lmax) {
    my @Ts = grep { @$_ >= 2 } map { [ map { $_ + 1 } @$_ ] } parts_list($L, $L); @Ts = sort { @$a <=> @$b } @Ts;
    for my $T (@Ts) { my $kT = join(',', @$T);
      # canonical relation: R = T minus its two largest parts (al <= be)
      my @t = @$T; my $be = shift @t; my $al = shift @t; my @Rm = @t; my $cc = $al + $be - 1; my $JR = 0; $JR += $_ for @Rm;
      my @val = (0) x $NU; my $add = sub { my ($c, $v) = @_; return unless defined $v && $c; for my $i (0..$NU-1) { $val[$i] = ($val[$i] + mm($c, $v->[$i])) % $P if $v->[$i]; } };
      $add->(md($JR + $cc), $YG->(@Rm, $cc));
      for my $sb ($subms->(\@Rm)) { my ($K, $C, $w) = @$sb; my $PiC = 1; $PiC = mm($PiC, $Pi->($_)) for @$C;
        for my $a (1..$al-1) { for my $b (1..$be-1) { my $e = $E0->($al-$a, $be-$b, $C); next unless $e;
            my @arg = (@$K, $a, $b); my @Tp = grep { $_ >= 2 } @arg; next unless @Tp; next if $lvl->(@Tp) < $D;
            $add->(md(-mm($w, mm($PiC, $e))), $YG->(@arg)); } } }
      $Y{$kT} = \@val; push @allT, $T; } }
  return { tau => $tau, sig => $sig, D => $D, Lmax => $Lmax, NU => $NU, Y => \%Y, allT => \@allT, Pi => $Pi };
}
# closed-form relation rows (orders k = D..Lmax) for a setup
sub kp_rows {
  my $S = shift; my ($tau, $D, $Lmax, $NU) = @$S{qw(tau D Lmax NU)}; my $Pi = $S->{Pi};
  my $M = $Lmax - $D + 1;
  # T = 1 - e^-t series; 2F1 series; rho_j = H_j/H_1^j
  my @Tt = (0); for my $n (1..$M+1) { push @Tt, mm(md(($n % 2) ? 1 : -1), mi(fct($n))); }
  my $mul = sub { my ($a, $b, $N) = @_; my @c = (0) x ($N+1); for my $i (0..$N) { next unless $a->[$i]; for my $j (0..$N-$i) { $c[$i+$j] = ($c[$i+$j] + mm($a->[$i], $b->[$j])) % $P if $b->[$j]; } } return \@c; };
  my $inv = sub { my ($a, $N) = @_; my @b = (mi($a->[0])); for my $n (1..$N) { my $s = 0; for my $k (1..$n) { $s = ($s + mm($a->[$k], $b[$n-$k])) % $P; } $b[$n] = md(-mm($s, $b[0])); } return \@b; };
  my $compose = sub { my ($f, $g, $N) = @_; my @res = (0) x ($N+1); my $pw = [1, (0) x $N]; for my $k (0..$N) { for my $i (0..$N) { $res[$i] = ($res[$i] + mm($f->[$k], $pw->[$i])) % $P; } $pw = $mul->($pw, $g, $N); } return \@res; };
  my $hyp = sub { my ($a, $b, $c, $N) = @_; my @s = (1); my $t = 1; for my $n (1..$N) { $t = mm($t, mm(mm(md($a + $n - 1), md($b + $n - 1)), mi(mm(md($c + $n - 1), $n)))); push @s, $t; } return \@s; };
  my $N = $M + 1;
  my $h1 = $compose->($hyp->(md(-$tau), md(-$tau), 1, $N), \@Tt, $N);
  my %rho;
  for my $j (2..$Lmax+1) { my $hj = $compose->($hyp->(md(-$j*$tau), md(-$j*$tau), $j, $N), \@Tt, $N);
    my $p = [1, (0) x $N]; $p = $mul->($p, $h1, $N) for 1..$j; $rho{$j} = $mul->($hj, $inv->($p, $N), $N); }
  my @rows;
  for my $k ($D..$Lmax) { my @eq = (0) x $NU;
    for my $LT ($D..$k) { for my $q (parts_list($LT, $LT)) { my @T = map { $_ + 1 } @$q; my $top = 1; my $pr = [1, (0) x $M];
        for my $t (@T) { $top = mm($top, mm(md(($t % 2) ? -1 : 1), mm($Pi->($t), mi($t)))); $pr = $mul->($pr, $rho{$t}, $M); }
        my %m; $m{$_}++ for @T; my $aut = 1; $aut = mm($aut, fct($_)) for values %m;
        my $co = mm(mm($top, $pr->[$k - $LT]), mi($aut)); next unless $co; my $v = $S->{Y}{join(',', @T)};
        for my $i (0..$NU-1) { $eq[$i] = ($eq[$i] + mm($co, $v->[$i])) % $P if $v->[$i]; } } }
    push @rows, [$k, \@eq]; }
  return \@rows;
}
1;
