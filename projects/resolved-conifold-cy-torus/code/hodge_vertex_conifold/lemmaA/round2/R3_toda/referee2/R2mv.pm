package R2mv;
# Second referee's own implementation (written from scratch, shares no code with mvlib.pm or rmv.pm) of the Schur side of the
# Marino-Vafa formula modulo a prime p, with unit legs, and the Laguerre transform to the dilaton point.
#   Z = sum_nu W_nu(x) s_nu(p),  W_nu = exp((tau+1/2) cs(nu) x) prod_box 1/(2 sinh(h x/2)),  cs = sum of contents = kappa/2,
#   w_nu := x^{|nu|} W_nu = (1/H_nu) exp((tau+1/2) cs x + sum_k c_k P_2k(nu) x^2k),  c_k = -B_2k/(2k (2k)!),  P_2k = sum_box h^2k.
#   Zt_S[k][m] = [x^m] x^{|S|+k} [p_S p_1^k] Z = [x^m] sum_{nu |- |S|+k} chi_nu(S u 1^k) w_nu / (z_S k!).
#   Unit legs zeta = p_1/x; connected part over labelled non-unit parts by the set-partition cumulant formula.
#   (1.1): [x^(2g-2+l)][p_mu] log Z = sigma^(g+l-1) Pi_mu M_g(mu)/|Aut mu|,   M_g(mu) = <Omega prod lam_{mu_i}>_g.
#   Xi_g(B) = sum_{a=0}^{N} (-1)^a/a! L_{N-a}^{(2g-2+r+a)}(-J-a) M_g(B+1^a),  N = 3g-3+r  (Route A Prop. 2.1),  g >= 1.
# Characters: Murnaghan-Nakayama on beta-sets (rim hooks of the non-unit parts), then f^lambda by the hook length formula.
use strict; use warnings;
use Math::BigRat; use Math::BigInt;

our $p = 1000000007;
our %fcache; our @fac = (1); our @ifac = (1);
sub setp { $p = shift; %fcache = (); @fac = (1); @ifac = (1); }
sub md { my $x = $_[0] % $p; $x += $p if $x < 0; $x }
sub inv { my $a = md($_[0]); die "R2mv::inv(0)" unless $a; my ($t, $nt, $r, $nr) = (0, 1, $p, $a);
  while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q * $nt); ($r, $nr) = ($nr, $r - $q * $nr); } md($t) }
sub fac { my $n = shift; for my $i (scalar(@fac)..$n) { $fac[$i] = ($fac[$i-1] * $i) % $p; } $fac[$n] }
sub ifac { my $n = shift; fac($n); for my $i (scalar(@ifac)..$n) { $ifac[$i] = inv($fac[$i]); } $ifac[$n] }
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n || $n < 0; fac($n) * ifac($k) % $p * ifac($n - $k) % $p }
sub rat { my $q = Math::BigRat->new("$_[0]"); my $n = $q->numerator->copy->bmod($p)->numify; my $d = $q->denominator->copy->bmod($p)->numify; md($n) * inv($d) % $p }
sub pw { my ($b, $e) = @_; my $r = 1; $b = md($b); while ($e > 0) { $r = $r * $b % $p if $e & 1; $b = $b * $b % $p; $e >>= 1; } $r }
sub Pi { my ($j, $tau) = @_; return 1 if $j <= 1; my $v = 1; for my $a (1..$j-1) { $v = $v * md($j * $tau + $a) % $p; } $v * ifac($j - 1) % $p }

# Bernoulli numbers (exact), B_1 = -1/2
my @bern = (Math::BigRat->new(1));
sub bern { my $n = shift; for my $m (scalar(@bern)..$n) { my $s = Math::BigRat->new(0); my $c = Math::BigInt->new(1);
    for my $k (0..$m-1) { $s += Math::BigRat->new($c) * $bern[$k]; $c = $c * ($m + 1 - $k) / ($k + 1); } $bern[$m] = -$s / ($m + 1); } $bern[$n] }
sub ck { my $k = shift; my $f = Math::BigInt->new(2*$k)->bfac; rat(-bern(2*$k) / (Math::BigRat->new(2*$k) * Math::BigRat->new($f))) }

# partition generator (descending), callback on each partition (array ref, do not keep it)
sub each_partition { my ($n, $cb) = @_; my @a; my $rec; $rec = sub { my ($rem, $mx) = @_;
    if ($rem == 0) { $cb->(\@a); return; } for my $k (reverse 1..($mx < $rem ? $mx : $rem)) { push @a, $k; $rec->($rem - $k, $k); pop @a; } };
  $rec->($n, $n); }
sub conjugate { my ($nu) = @_; my @c; for my $i (0..$#$nu) { $c[$_]++ for 0..$nu->[$i]-1; } \@c }
sub cmp_part { my ($a, $b) = @_; my $n = @$a > @$b ? @$a : @$b; for my $i (0..$n-1) { my $x = $a->[$i] // 0; my $y = $b->[$i] // 0; return $x <=> $y if $x != $y; } 0 }
sub hooks { my ($nu, $c) = @_; my @h; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { push @h, $nu->[$i] - $j + $c->[$j] - $i - 1; } } @h }
sub contsum { my ($nu) = @_; my $s = 0; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { $s += $j - $i; } } $s }

# f^lambda mod p, lambda given as descending list; cached
sub flam { my @l = grep { $_ > 0 } @_; my $key = join(',', @l); return $fcache{$key} if exists $fcache{$key};
  my $n = 0; $n += $_ for @l; my $c = conjugate(\@l); my $H = 1; $H = $H * $_ % $p for hooks(\@l, $c);
  $fcache{$key} = fac($n) * inv($H) % $p }
# chi_nu(parts u 1^rest) via beta sets; beta = descending distinct non-negative integers; parts descending
my %cmemo;
sub chi_beta { my ($beta, $parts, $pi) = @_; $pi //= 0;
  if ($pi > $#$parts) { my $L = scalar @$beta; my @l = map { $beta->[$_] - ($L - 1 - $_) } 0..$L-1; return flam(@l); }
  my $key = join(',', @$beta) . '|' . $pi . '|' . join(',', @$parts); if ($pi > 0) { return $cmemo{$key} if exists $cmemo{$key}; }
  my $s = $parts->[$pi]; my %occ; $occ{$_} = 1 for @$beta; my $sum = 0;
  for my $i (0..$#$beta) { my $bd = $beta->[$i]; my $t = $bd - $s; next if $t < 0 || $occ{$t};
    my $ht = 0; for my $bb (@$beta) { $ht++ if $bb > $t && $bb < $bd; }
    my @nb = sort { $b <=> $a } map { $_ == $bd ? $t : $_ } @$beta;
    my $v = chi_beta(\@nb, $parts, $pi + 1); $sum = ($ht & 1) ? $sum - $v : $sum + $v; }
  my $res = md($sum); $cmemo{$key} = $res if $pi > 0; $res }
sub chi { my ($nu, $parts) = @_; my $L = scalar @$nu; my @beta = map { $nu->[$_] + $L - 1 - $_ } 0..$L-1; chi_beta(\@beta, [sort { $b <=> $a } @$parts]) }
sub clear_chi { %cmemo = (); }

# w_nu series mod p to order P (from hooks multiset and content sum)
my @ckc;
sub prep { my $P = shift; @ckc = (); for (my $k = 1; 2*$k <= $P; $k++) { $ckc[$k] = ck($k); } fac(4*$P + 200); ifac(4*$P + 200); }
sub wser { my ($hm, $H, $cs, $tau, $P) = @_;   # hm: hash hook -> multiplicity
  my @e = (0) x ($P + 1);
  $e[1] = md((2*$tau + 1) % $p * inv(2) % $p * md($cs)) if $P >= 1;
  my @hv = keys %$hm; my @h2 = map { $_ * $_ % $p } @hv; my @mu = map { $hm->{$_} } @hv; my @pw = (1) x scalar(@hv);
  for (my $k = 1; 2*$k <= $P; $k++) { my $s = 0; for my $q (0..$#hv) { $pw[$q] = $pw[$q] * $h2[$q] % $p; $s = ($s + $pw[$q] * $mu[$q]) % $p; } $e[2*$k] = $ckc[$k] * $s % $p; }
  my @nz = grep { $e[$_] } 1..$P; my @me = map { $_ * $e[$_] % $p } 0..$P;
  my @f = (inv($H)); my @invn = map { $_ ? inv($_) : 0 } 0..$P;
  for my $n (1..$P) { my $s = 0; for my $m (@nz) { last if $m > $n; $s = ($s + $me[$m] * $f[$n - $m]) % $p; } $f[$n] = $s * $invn[$n] % $p; }
  \@f }

# bivariate series [k][m]
sub bz { my ($K, $P) = @_; [ map { [ (0) x ($P + 1) ] } 0..$K ] }
sub bmul { my ($A, $B, $K, $P) = @_; my $C = bz($K, $P);
  for my $k1 (0..$K) { my $a = $A->[$k1]; my @ia = grep { $a->[$_] } 0..$P; next unless @ia;
    for my $k2 (0..$K - $k1) { my $b = $B->[$k2]; my @ib = grep { $b->[$_] } 0..$P; next unless @ib; my $c = $C->[$k1 + $k2];
      for my $i (@ia) { my $x = $a->[$i]; for my $j (@ib) { last if $i + $j > $P; $c->[$i + $j] = ($c->[$i + $j] + $x * $b->[$j]) % $p; } } } }
  $C }
sub bscal { my ($A, $s) = @_; [ map { [ map { $_ * $s % $p } @$_ ] } @$A ] }
sub badd { my ($A, $B) = @_; [ map { my $k = $_; [ map { ($A->[$k][$_] + $B->[$k][$_]) % $p } 0..$#{$A->[$k]} ] } 0..$#$A ] }
sub sinv { my ($a, $P) = @_; die "sinv: a0 != 1" unless $a->[0] == 1; my @g = (1); for my $m (1..$P) { my $s = 0; for my $j (1..$m) { $s = ($s + $a->[$j] * $g[$m - $j]) % $p if $a->[$j]; } $g[$m] = md(-$s); } \@g }
sub smul { my ($a, $b, $P) = @_; my @c = (0) x ($P + 1); for my $i (0..$P) { next unless $a->[$i]; for my $j (0..$P - $i) { $c[$i + $j] = ($c[$i + $j] + $a->[$i] * $b->[$j]) % $p if $b->[$j]; } } \@c }
sub binv { my ($A, $K, $P) = @_;   # A[0] = series with constant 1
  my $i0 = sinv($A->[0], $P); my $G = bz($K, $P); $G->[0] = $i0;
  for my $k (1..$K) { my @s = (0) x ($P + 1); for my $j (1..$k) { my $t = smul($A->[$j], $G->[$k - $j], $P); $s[$_] = ($s[$_] + $t->[$_]) % $p for 0..$P; }
    my $t = smul(\@s, $i0, $P); $G->[$k] = [ map { md(-$_) } @$t ]; }
  $G }
sub blog { my ($A, $K, $P) = @_;   # log of A with A[0] = 1 (as series in x), A[k] for k >= 1; returns L with L[0] = log A[0] = 0 (A[0] must be exactly 1)
  for my $m (1..$P) { die "blog: A[0] not 1" if $A->[0][$m]; }
  my $L = bz($K, $P);
  for my $k (1..$K) { my @s = map { $_ * $k % $p } @{$A->[$k]}; for my $j (1..$k-1) { my $t = smul($L->[$j], $A->[$k - $j], $P); $s[$_] = md($s[$_] - $j * $t->[$_]) for 0..$P; }
    my $ik = inv($k); $L->[$k] = [ map { $_ * $ik % $p } @s ]; }
  $L }
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @o; for my $sp (setparts($n - 1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n - 1; push @o, \@c; } my @c = map { [@$_] } @$sp; push @c, [$n - 1]; push @o, \@c; } @o }
sub autm { my %m; $m{$_}++ for @_; my $a = 1; $a = $a * fac($m{$_}) % $p for keys %m; $a }
sub zS { my %m; $m{$_}++ for @_; my $z = 1; for my $v (keys %m) { $z = $z * pw($v, $m{$v}) % $p * fac($m{$v}) % $p; } $z }
sub lag { my ($n, $al, $x) = @_; my $s = 0; my $xp = 1; for my $i (0..$n) { my $t = binom($n + $al, $n - $i) * $xp % $p * ifac($i) % $p; $s = ($i & 1) ? md($s - $t) : ($s + $t) % $p; $xp = $xp * md($x) % $p; } $s }

# main: compute(gmax, [B lists], tau mod p, %opt) -> { "b1,b2,.." => { M => {g => [M_g(B+1^a)]}, Xi => {g => value} } }
# option empty => 1 also returns key '' with M_g(1^a) (from log Zt_empty).
sub compute { my ($gmax, $Bs, $tau, %o) = @_;
  my $sig = $tau * md($tau + 1) % $p; die "sigma = 0 mod p" unless $sig;
  my (%need, $K, $P); $K = 0; $P = 0;
  for my $B (@$Bs) { my @b = sort { $b <=> $a } @$B; my $r = @b; my $J = 0; $J += $_ for @b; my $N = 3*$gmax - 3 + $r; $N = 0 if $N < 0;
    $K = $N if $N > $K; my $PP = 2*$gmax - 2 + $r + 2*$N + $J; $P = $PP if $PP > $P;
    for my $mask (0..(1 << $r) - 1) { my @s = sort { $b <=> $a } map { $b[$_] } grep { ($mask >> $_) & 1 } 0..$r-1; $need{join(',', @s)} = \@s; } }
  $need{''} = [];
  prep($P); my $t0 = time;
  my %Zt; my %Jk; my %izk;
  for my $k (keys %need) { $Zt{$k} = bz($K, $P); my $J = 0; $J += $_ for @{$need{$k}}; $Jk{$k} = $J; $izk{$k} = inv(zS(@{$need{$k}})); }
  my $Smax = 0; for my $k (keys %need) { $Smax = $Jk{$k} + $K if $Jk{$k} + $K > $Smax; }
  my $npart = 0;
  for my $n (0..$Smax) {
    my @keys = grep { my $kk = $n - $Jk{$_}; $kk >= 0 && $kk <= $K } keys %need; next unless @keys;
    each_partition($n, sub { my ($nu) = @_; my $c = conjugate($nu); my $cmp = cmp_part($nu, $c); return if $cmp < 0;
      my @h = hooks($nu, $c); my %hm; $hm{$_}++ for @h; my $H = 1; $H = $H * $_ % $p for @h; my $cs = contsum($nu);
      my $w = wser(\%hm, $H, $cs, $tau, $P); $npart++;
      for my $pass (0, 1) { next if $pass == 1 && $cmp == 0; my $nn = $pass ? $c : $nu;
        for my $k (@keys) { my $kk = $n - $Jk{$k}; my $ch = chi($nn, $need{$k}); next unless $ch;
          my $cf = $ch * $izk{$k} % $p * ifac($kk) % $p; my $row = $Zt{$k}[$kk];
          if ($pass) { for my $m (0..$P) { next unless $w->[$m]; my $v = $cf * $w->[$m] % $p; $row->[$m] = ($m & 1) ? md($row->[$m] - $v) : ($row->[$m] + $v) % $p; } }
          else { for my $m (0..$P) { $row->[$m] = ($row->[$m] + $cf * $w->[$m]) % $p if $w->[$m]; } } } } });
    clear_chi();
    printf STDERR "#   n=%d done (%d partitions so far, %d s)\n", $n, $npart, time - $t0 if $o{verbose}; }
  my $Inv = binv($Zt{''}, $K, $P);
  my %Y; for my $k (keys %need) { next if $k eq ''; $Y{$k} = bscal(bmul($Zt{$k}, $Inv, $K, $P), autm(@{$need{$k}})); }
  my %res;
  my @todo = map { [sort { $b <=> $a } @$_] } @$Bs; push @todo, [] if $o{empty};
  for my $bb (@todo) { my @b = @$bb; my $r = @b; my $J = 0; $J += $_ for @b; my $LB;
    if ($r == 0) { $LB = blog($Zt{''}, $K, $P); }
    else { $LB = bz($K, $P);
      for my $sp (setparts($r)) { my $nb = scalar @$sp; my $prod;
        for my $blk (@$sp) { my $key = join(',', sort { $b <=> $a } map { $b[$_] } @$blk); $prod = defined $prod ? bmul($prod, $Y{$key}, $K, $P) : $Y{$key}; }
        my $cf = fac($nb - 1); $cf = md(-$cf) if ($nb - 1) & 1; $LB = badd($LB, bscal($prod, $cf)); }
      $LB = bscal($LB, inv(autm(@b))); }
    my $PiB = 1; $PiB = $PiB * Pi($_, $tau) % $p for @b; my $au = autm(@b);
    my (%M, %Xi);
    for my $g (0..$gmax) { my $N = 3*$g - 3 + $r; my @Mg;
      my $Na = $N + ($o{extra} // 0); for my $a (0..($Na > $K ? $K : ($Na < 0 ? 0 : $Na))) { my $idx = 2*$g - 2 + $r + 2*$a + $J; if ($idx < 0 || $idx > $P) { push @Mg, undef; next; }
        my $den = pw($sig, $g + $r + $a - 1) * $PiB % $p; push @Mg, $LB->[$a][$idx] * $au % $p * fac($a) % $p * inv($den) % $p; }
      $M{$g} = \@Mg;
      if ($g >= 1 && ($g >= 2 || $r >= 1) && $N >= 0 && $N <= $K) { my $xi = 0;
        for my $a (0..$N) { my $t = lag($N - $a, 2*$g - 2 + $r + $a, md(-$J - $a)) * $Mg[$a] % $p * ifac($a) % $p; $xi = ($a & 1) ? md($xi - $t) : ($xi + $t) % $p; }
        $Xi{$g} = $xi; } }
    $res{join(',', @b)} = { M => \%M, Xi => \%Xi }; }
  printf STDERR "# R2mv::compute gmax=%d K=%d P=%d Smax=%d: %d partitions (conjugate pairs) in %d s\n", $gmax, $K, $P, $Smax, $npart, time - $t0 if $o{verbose} || $o{report};
  \%res }
1;
