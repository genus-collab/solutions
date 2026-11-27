package rmv;
# rmv.pm -- referee's own implementation (independent of mvlib.pm) of the Schur side of the Marino-Vafa formula modulo a prime p,
# with unit legs p_1 = zeta*x, cumulants over labelled non-unit parts, and the Laguerre transform to the dilaton point.
#   R = sum_nu exp((tau+1/2) kappa_nu x/2) prod_box 1/(2 sinh(h x/2)) s_nu(p),   W_nu = x^{-|nu|} w_nu(x),
#   w_nu = (1/H_nu) exp( (tau+1/2)(sum of contents) x + sum_k c_k P_2k(nu) x^2k ),  c_k = -B_2k/(2k (2k)!),  P_2k = sum_box h^2k.
#   [x^(2g-2+l)][p_mu] log R = sigma^(g+l-1) Pi_mu M_g(mu)/|Aut mu|   (the route's (1.1); re-checked against the Hodge tables).
#   Xi_h(B) = sum_{k=0}^{N} (-1)^k/k! L_{N-k}^{(2h-2+r+k)}(-J-k) M_h(B+1^k),  N = 3h-3+r.
use strict; use warnings;
use Math::BigRat;
our $p = 1000000007;
sub setp { $p = shift; }
sub md { my $x = $_[0] % $p; return $x < 0 ? $x + $p : $x; }
sub mi { my $a = md($_[0]); die "rmv: zero not invertible" if $a == 0; my ($t, $nt, $r, $nr) = (0, 1, $p, $a);
  while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q * $nt); ($r, $nr) = ($nr, $r - $q * $nr); } return md($t); }
sub rat2mod { my $q = Math::BigRat->new($_[0]); my $n = $q->numerator->copy->bmod($p)->numify; my $d = $q->denominator->copy->bmod($p)->numify; return ($n * mi($d)) % $p; }
my @fac = (1); sub fac { my $n = shift; for my $i (scalar(@fac)..$n) { $fac[$i] = ($fac[$i-1] * $i) % $p; } return $fac[$n]; }
sub reset_cache { @fac = (1); }
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; return (fac($n) * mi((fac($k) * fac($n - $k)) % $p)) % $p; }
# exact Bernoulli numbers (B_1 = -1/2) as BigRat
my @B = (Math::BigRat->new(1));
sub bern { my $n = shift; for my $m (scalar(@B)..$n) { my $s = Math::BigRat->new(0); my $bin = Math::BigInt->new(1);
    for my $k (0..$m-1) { $s += Math::BigRat->new($bin) * $B[$k]; $bin = $bin * ($m + 1 - $k) / ($k + 1); } $B[$m] = -$s / ($m + 1); } return $B[$n]; }
# partitions of n (descending arrays), cached
my %PC;
sub partitions { my ($n, $mx) = @_; $mx //= $n; $mx = $n if $mx > $n; my $key = "$n,$mx"; return $PC{$key} if $PC{$key};
  my @o; if ($n == 0) { @o = ([]); } else { for my $k (reverse 1..$mx) { for my $q (@{ partitions($n - $k, $k) }) { push @o, [$k, @$q]; } } }
  $PC{$key} = \@o; return \@o; }
sub conj { my ($nu) = @_; my @c; for my $i (0..$#$nu) { $c[$_]++ for 0..$nu->[$i]-1; } return \@c; }
sub hooks { my ($nu) = @_; my $c = conj($nu); my @h; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { push @h, $nu->[$i] - $j + $c->[$j] - $i - 1; } } return @h; }
sub content_sum { my ($nu) = @_; my $s = 0; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { $s += $j - $i; } } return $s; }
my %FC;
sub fdim { my ($nu) = @_; my $key = join(',', @$nu); return $FC{$key} if exists $FC{$key}; my $n = 0; $n += $_ for @$nu; my $H = 1; $H = ($H * $_) % $p for hooks($nu);
  return $FC{$key} = (fac($n) * mi($H)) % $p; }
# character chi_nu(B + 1^rest) by Murnaghan-Nakayama (beta sets), B = list of non-unit parts
my %CC;
sub chi { my ($nu, $Bp) = @_; my @Bs = @$Bp; return fdim($nu) unless @Bs; my $key = join(',', @$nu) . '|' . join(',', @Bs); return $CC{$key} if exists $CC{$key};
  my $bsz = shift @Bs; my $L = scalar @$nu; my @beta = map { $nu->[$_] + $L - 1 - $_ } 0..$L-1; my %occ = map { $_ => 1 } @beta; my $res = 0;
  for my $i (0..$L-1) { my $nb = $beta[$i] - $bsz; next if $nb < 0 || $occ{$nb}; my $ht = 0; for my $bb (@beta) { $ht++ if $bb > $nb && $bb < $beta[$i]; }
    my @n2 = sort { $b <=> $a } map { $_ == $i ? $nb : $beta[$_] } 0..$L-1; my @lam; for my $q (0..$L-1) { my $v = $n2[$q] - ($L - 1 - $q); push @lam, $v if $v > 0; }
    my $v = chi(\@lam, \@Bs); $res = ($ht % 2) ? md($res - $v) : md($res + $v); }
  return $CC{$key} = $res; }
sub clear_chi { %CC = (); }
# w_nu(x) mod p to order P, given taup = tau mod p
my @ck;   # c_k mod p
sub prep_ck { my $P = shift; @ck = (); for (my $k = 1; 2*$k <= $P; $k++) { my $c = -bern(2*$k) / (Math::BigRat->new(2*$k) * Math::BigRat->new(Math::BigInt->new(2*$k)->bfac)); $ck[$k] = rat2mod($c); } }
sub wser { my ($nu, $taup, $P) = @_; my @h = hooks($nu); my %hm; $hm{$_}++ for @h; my $H = 1; $H = ($H * $_) % $p for @h;
  my @e = (0) x ($P + 1); my $half = mi(2);
  $e[1] = (((2*$taup + 1) % $p) * md(content_sum($nu)) % $p) * $half % $p if $P >= 1;
  my @h2 = map { [ ($_ * $_) % $p, $hm{$_} ] } keys %hm; my @pw = map { 1 } @h2;
  for (my $k = 1; 2*$k <= $P; $k++) { my $s = 0; for my $q (0..$#h2) { $pw[$q] = ($pw[$q] * $h2[$q][0]) % $p; $s = ($s + $pw[$q] * $h2[$q][1]) % $p; } $e[2*$k] = ($ck[$k] * $s) % $p; }
  my @nz = grep { $e[$_] } 1..$P; my @me = map { ($_ * $e[$_]) % $p } 0..$P;
  my @f = (mi($H)); for my $n (1..$P) { my $s = 0; for my $m (@nz) { last if $m > $n; $s = ($s + $me[$m] * $f[$n - $m]) % $p; } $f[$n] = ($s * mi($n)) % $p; }
  return \@f; }
# bivariate series helpers: A[k][m], k <= K, m <= P
sub bzero { my ($K, $P) = @_; return [ map { [ (0) x ($P + 1) ] } 0..$K ]; }
sub bmul { my ($A, $Bb, $K, $P) = @_; my $C = bzero($K, $P);
  for my $k1 (0..$K) { my $a = $A->[$k1]; for my $m1 (0..$P) { my $x = $a->[$m1]; next unless $x;
      for my $k2 (0..$K - $k1) { my $b = $Bb->[$k2]; my $c = $C->[$k1 + $k2]; for my $m2 (0..$P - $m1) { next unless $b->[$m2]; $c->[$m1 + $m2] = ($c->[$m1 + $m2] + $x * $b->[$m2]) % $p; } } } }
  return $C; }
sub badd { my ($A, $Bb, $s, $K, $P) = @_; my $C = bzero($K, $P); for my $k (0..$K) { for my $m (0..$P) { $C->[$k][$m] = md($A->[$k][$m] + $s * $Bb->[$k][$m]); } } return $C; }
sub binv0 { my ($A, $K, $P) = @_;   # inverse of A with A[0] = 1 + O(x)? here A[0][0] = 1 and A[0][m>0] possibly non-zero: invert as series in zeta with x-series coefficients
  my $a0 = $A->[0]; die "binv0: constant term" unless $a0->[0] == 1;
  my @i0 = (1); for my $m (1..$P) { my $s = 0; for my $j (1..$m) { $s = ($s + $a0->[$j] * $i0[$m - $j]) % $p; } $i0[$m] = md(-$s); }   # 1/a0(x)
  my $G = bzero($K, $P); $G->[0] = [ @i0 ];
  my $smul = sub { my ($u, $v) = @_; my @w = (0) x ($P + 1); for my $i (0..$P) { next unless $u->[$i]; for my $j (0..$P - $i) { $w[$i + $j] = ($w[$i + $j] + $u->[$i] * $v->[$j]) % $p if $v->[$j]; } } return \@w; };
  for my $k (1..$K) { my @s = (0) x ($P + 1); for my $j (1..$k) { my $t = $smul->($A->[$j], $G->[$k - $j]); $s[$_] = ($s[$_] + $t->[$_]) % $p for 0..$P; }
    my $t = $smul->(\@s, \@i0); $G->[$k] = [ map { md(-$_) } @$t ]; }
  return $G; }
# set partitions of 0..n-1
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @o; for my $sp (setparts($n - 1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n - 1; push @o, \@c; } my @c = map { [@$_] } @$sp; push @c, [$n - 1]; push @o, \@c; } return @o; }
sub autf { my @b = @_; my %m; $m{$_}++ for @b; my $a = 1; $a = ($a * fac($m{$_})) % $p for keys %m; return $a; }
sub zB { my @b = @_; my %m; $m{$_}++ for @b; my $z = 1; for my $v (keys %m) { for (1..$m{$v}) { $z = ($z * $v) % $p; } $z = ($z * fac($m{$v})) % $p; } return $z; }
sub Pi { my ($j, $taup) = @_; return 1 if $j <= 1; my $v = 1; for my $b (1..$j-1) { $v = ($v * md($j * $taup + $b)) % $p; } return ($v * mi(fac($j - 1))) % $p; }
sub lag { my ($n, $al, $x) = @_; my $s = 0; my $xp = 1; for my $i (0..$n) { my $t = (binom($n + $al, $n - $i) * $xp) % $p; $t = ($t * mi(fac($i))) % $p; $s = ($i % 2) ? md($s - $t) : ($s + $t) % $p; $xp = ($xp * md($x)) % $p; } return $s; }
# main driver: for genus h and a list of multisets Bs (non-unit parts >= 2), at tau (mod p): returns { "B" => { M => [M_h(B+1^k)]_k, Xi => value } }
# also returns M_g for g = 0..h if requested (for table checks)
sub compute { my ($h, $Bs, $taup, %o) = @_; my $gmax = $o{gmax} // $h;
  my $sig = ($taup * md($taup + 1)) % $p; die "sigma = 0 mod p" unless $sig;
  # needed sub-multisets (as sorted lists) for the cumulants
  my %need; my ($Kmax, $Pmax, $Smax) = (0, 0, 0);
  for my $B (@$Bs) { my @b = @$B; my $r = scalar @b; my $J = 0; $J += $_ for @b; my $N = 3*$h - 3 + $r;
    $Kmax = $N if $N > $Kmax; my $P = 2*$h - 2 + $r + 2*$N + $J; $Pmax = $P if $P > $Pmax;
    for my $mask (0..(1 << $r) - 1) { my @s = sort { $b <=> $a } map { $b[$_] } grep { ($mask >> $_) & 1 } 0..$r-1; $need{join(',', @s)} = \@s; } }
  for my $k (keys %need) { my $J = 0; $J += $_ for @{$need{$k}}; $Smax = $J + $Kmax if $J + $Kmax > $Smax; }
  prep_ck($Pmax); fac($Smax + $Pmax + 10);
  my %Rs; for my $k (keys %need) { $Rs{$k} = bzero($Kmax, $Pmax); }
  my %zinv; for my $k (keys %need) { $zinv{$k} = mi(zB(@{$need{$k}})); }
  for my $n (0..$Smax) { for my $nu (@{ partitions($n) }) { my $w;
      for my $k (keys %need) { my $Bn = $need{$k}; my $J = 0; $J += $_ for @$Bn; my $kk = $n - $J; next if $kk < 0 || $kk > $Kmax;
        my $ch = chi($nu, $Bn); next unless $ch; $w //= wser($nu, $taup, $Pmax);
        my $cf = ((($ch * $zinv{$k}) % $p) * mi(fac($kk))) % $p; my $row = $Rs{$k}[$kk];
        for my $m (0..$Pmax) { $row->[$m] = ($row->[$m] + $cf * $w->[$m]) % $p if $w->[$m]; } } } clear_chi(); }
  my $Rinv = binv0($Rs{''}, $Kmax, $Pmax);
  my %mt; for my $k (keys %need) { next if $k eq ''; my @b = @{$need{$k}}; my $au = autf(@b); my $T = bmul($Rs{$k}, $Rinv, $Kmax, $Pmax);
    for my $kk (0..$Kmax) { $_ = ($_ * $au) % $p for @{$T->[$kk]}; } $mt{$k} = $T; }
  # cumulants kappa(S) = d^S log R (derivative-normalised), by the multiset recursion
  #   kappa(S) = m(S) - sum_{T subset S, s0 in T, T != S} #(labelled T of that type containing s0) kappa(T) m(S\T)
  my %kap; my $kapf; $kapf = sub { my @s = @_; my $key = join(',', @s); return $kap{$key} if exists $kap{$key};
    my %m; $m{$_}++ for @s; my @vv = sort { $b <=> $a } keys %m; my $v0 = $vv[0]; my $acc = [ map { [ @$_ ] } @{$mt{$key}} ];
    my @types = ([]); for my $v (@vv) { my @nt; my $lo = ($v == $v0) ? 1 : 0; for my $t (@types) { push @nt, [ @$t, $_ ] for $lo..$m{$v}; } @types = @nt; }
    for my $t (@types) { my $tot = 0; $tot += $_ for @$t; next if $tot == scalar @s;
      my $cnt = 1; for my $q (0..$#vv) { my $v = $vv[$q]; $cnt = ($cnt * ($v == $v0 ? binom($m{$v} - 1, $t->[$q] - 1) : binom($m{$v}, $t->[$q]))) % $p; }
      my @T = map { ($vv[$_]) x $t->[$_] } 0..$#vv; my @C = map { ($vv[$_]) x ($m{$vv[$_]} - $t->[$_]) } 0..$#vv;
      my $prod = bmul($kapf->(@T), $mt{join(',', @C)}, $Kmax, $Pmax); $acc = badd($acc, $prod, md(-$cnt), $Kmax, $Pmax); }
    return $kap{$key} = $acc; };
  my %res;
  for my $B (@$Bs) { my @b = @$B; my $r = scalar @b; my $J = 0; $J += $_ for @b; my $N = 3*$h - 3 + $r;
    my $C;
    if ($o{setparts}) { $C = bzero($Kmax, $Pmax);
      for my $sp (setparts($r)) { my $nb = scalar @$sp; my $prod;
        for my $blk (@$sp) { my $key = join(',', sort { $b <=> $a } map { $b[$_] } @$blk); $prod = defined $prod ? bmul($prod, $mt{$key}, $Kmax, $Pmax) : $mt{$key}; }
        my $cf = fac($nb - 1); $cf = md(-$cf) if ($nb - 1) % 2; $C = badd($C, $prod, $cf, $Kmax, $Pmax); } }
    else { $C = $kapf->(@b); }
    my $ai = mi(autf(@b)); my $PiB = 1; $PiB = ($PiB * Pi($_, $taup)) % $p for @b;
    my %Mg;
    for my $g (0..$gmax) { my @M; for my $kk (0..$N) { my $idx = 2*$g - 2 + $r + 2*$kk + $J; if ($idx < 0 || $idx > $Pmax) { push @M, undef; next; }
        my $v = ($C->[$kk][$idx] * $ai) % $p; $v = ($v * autf(@b)) % $p; $v = ($v * fac($kk)) % $p;   # Aut(B) k! [..]C_B with C_B = (1/Aut) sum ... : net factor k!
        my $den = 1; for (1..($g + $r + $kk - 1)) { $den = ($den * $sig) % $p; } $den = ($den * $PiB) % $p;
        push @M, ($v * mi($den)) % $p; } $Mg{$g} = \@M; }
    my $xi = 0; my $Mh = $Mg{$h};
    for my $kk (0..$N) { my $t = (lag($N - $kk, 2*$h - 2 + $r + $kk, md(-$J - $kk)) * $Mh->[$kk]) % $p; $t = ($t * mi(fac($kk))) % $p; $xi = ($kk % 2) ? md($xi - $t) : ($xi + $t) % $p; }
    # s-coefficients of Q_h(s;B): [s^m]Q = sum_{a<=m} M_a/a! (-1)^(m-a) L_{m-a}^{(2h-3+r+a)}(-(J+a))
    my @qs; for my $m (0..$N + 2) { my $s = 0; for my $a (0..($m < $N ? $m : $N)) { my $t = (lag($m - $a, 2*$h - 3 + $r + $a, md(-$J - $a)) * $Mh->[$a]) % $p; $t = ($t * mi(fac($a))) % $p; $s = (($m - $a) % 2) ? md($s - $t) : ($s + $t) % $p; } push @qs, $s; }
    $res{join(',', @b)} = { Xi => $xi, M => \%Mg, Q => \@qs, N => $N }; }
  return \%res; }
1;
