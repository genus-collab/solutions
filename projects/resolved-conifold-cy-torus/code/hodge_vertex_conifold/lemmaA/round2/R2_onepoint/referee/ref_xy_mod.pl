#!/usr/bin/perl
# ref_xy_mod.pl TAU G N JMAX [ENGINE_TABLE|-] [PRIME_INDEX]
# Referee's independent implementation of Theorem 1 of the R2 report (ABDKS x-y swap closed formula on C_tau), in the FORMAL form
# that is valid for every rational tau (no finite products):
#     Q_mu(v) = exp(-mu Sc(mu hb d) xhat),  xhat = tau S(hb d/tau) log v - e2 S(hb d/e2) log(1+v),
# expanded at v = oo by exp-log:  Q_mu / v^mu = exp( e2 mu log(1+w)
#     + Sum_{k>=1} (2k-1)! hb^{2k} [ phi_k(tau mu) tau^{-2k} w^{2k} - phi_k(e2 mu) e2^{-2k} w^{2k} (1+w)^{-2k} ] ),  w = 1/v,
# with sinh(N t/2)/sinh(t/2) = Sum_k phi_k(N) t^{2k}.  Iterated residue: [v_1^-1] ... [v_n^-1] (v_1 first) of
#     prod Q_{mu_i}(v_i) Sum_{n-cycles} prod_i 1/(v_i - v_s(i) + k_i hb),  k_i = (mu_i + mu_s(i))/2,
# using  [v^-1] G(v) prod_{t=1}^m 1/(v - c_t) = Sum_a g_a h_{a-m+1}(c)  (complete homogeneous polynomials of the nodes).
# Then [x^mu]F_{g,n} = (-1)^{n-1}/prod mu_i [hb^{2g-2+n}] (...), and Xi = sigma^{g-1} [x^mu]F / prod Pi_{mu_i} (dictionary P2).
# All arithmetic modulo a prime p ~ 2^31 (exact native integers); engine rationals are reduced mod p.  A "match" mod p is a
# necessary condition for equality; a mismatch mod p is a certain mismatch.
use strict; use warnings;
my ($TAUS, $G, $N, $JMAX, $TABLE, $PIDX) = @ARGV;
die "usage: TAU G N JMAX [TABLE|-] [PIDX]\n" unless defined $JMAX;
$PIDX //= 0;
# ---------- primes ----------
sub isprime { my $q = shift; return 0 if $q < 2; for (my $d = 2; $d * $d <= $q; $d++) { return 0 if $q % $d == 0; } return 1; }
my @PR; { my $q = 2147483647; while (@PR < 4) { push @PR, $q if isprime($q); $q -= 2; } }
my $P = $PR[$PIDX];
sub mm { return ($_[0] * $_[1]) % $P; }
sub ad { my $s = $_[0] + $_[1]; return $s >= $P ? $s - $P : $s; }
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub iv { my $a = shift; die "division by zero mod p" if $a % $P == 0; return pw($a, $P - 2); }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; $s =~ s/\s+//g; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), iv(decmod($2))); } die "bad rational '$s'" unless $s =~ /^-?\d+$/; return decmod($s); }
sub im { my $k = shift; return $k >= 0 ? $k % $P : ng((-$k) % $P); }
# ---------- parameters ----------
my $tau = ratmod($TAUS); my $e2 = ad($tau, 1); my $sig = mm($tau, $e2);
die "sigma = 0" if $sig == 0;
my $H = 2*$G - 2 + $N; die "unstable" if $H <= 0;
my $HT = $H;   # hb truncation for the polynomial arithmetic
my $HB = $N;   # variable index of hb (v_1..v_n have indices 0..n-1)
my $BITS = 6;
sub vkey { my ($var, $e) = @_; return $e << ($BITS * $var); }
sub hbdeg { return ($_[0] >> ($BITS * $HB)) & 63; }
sub vdeg { my ($k, $var) = @_; return ($k >> ($BITS * $var)) & 63; }
my @fact = (1); push @fact, mm($fact[-1], $_) for 1..300;
my @ifact = map { iv($_) } @fact;
sub binm { my ($a, $b) = @_; return 0 if $b < 0 || $b > $a; return mm($fact[$a], mm($ifact[$b], $ifact[$a-$b])); }
# phi_k(Nm): coefficients of sinh(N t/2)/sinh(t/2) in t^{2k}, N given mod p
sub phis { my ($Nm, $K) = @_; my $half = iv(2); my (@num, @den);
  for my $i (0..$K) { $num[$i] = mm(pw(mm($Nm, $half), 2*$i+1), $ifact[2*$i+1]); $den[$i] = mm(pw($half, 2*$i+1), $ifact[2*$i+1]); }
  my @q; my $id0 = iv($den[0]);
  for my $i (0..$K) { my $s = $num[$i]; for my $j (1..$i) { $s = ad($s, ng(mm($den[$j], $q[$i-$j]))); } $q[$i] = mm($s, $id0); }
  return \@q; }
# Laurent coefficients q_a(hb) of Q_mu/v^mu, a = 0..A, hb-degree <= HTq ; returns ref to array of arrays [a][p]
my %QC;
sub qser { my ($mu, $A, $HTq) = @_; my $key = "$mu|$A|$HTq"; return $QC{$key} if $QC{$key};
  my $K = int($HTq / 2);
  my $Na = mm($tau, im($mu)); my $Nb = mm($e2, im($mu));
  my $pa = phis($Na, $K); my $pb = phis($Nb, $K);
  my @L = map { [ (0) x ($HTq + 1) ] } 0..$A;   # L[a][p]
  for my $a (1..$A) { my $c = mm($Nb, iv($a)); $c = ng($c) if $a % 2 == 0; $L[$a][0] = ad($L[$a][0], $c); }
  for my $k (1..$K) { my $f = $fact[2*$k-1];
    my $ca = mm($f, mm($pa->[$k], iv(pw($tau, 2*$k))));
    my $cb = mm($f, mm($pb->[$k], iv(pw($e2, 2*$k))));
    $L[2*$k][2*$k] = ad($L[2*$k][2*$k], $ca) if 2*$k <= $A;
    for my $b (0..$A - 2*$k) { my $cf = binm(2*$k + $b - 1, $b); $cf = ng($cf) if $b % 2 == 1;   # (1+w)^{-2k}
      $L[2*$k + $b][2*$k] = ad($L[2*$k + $b][2*$k], ng(mm($cb, $cf))); } }
  # exp: E_0 = 1, E_a = (1/a) Sum_{j=1}^a j L_j E_{a-j}
  my @E = ([ 1, (0) x $HTq ]);
  for my $a (1..$A) { my @acc = (0) x ($HTq + 1);
    for my $j (1..$a) { my $Lj = $L[$j]; my $Ea = $E[$a-$j];
      for my $p (0..$HTq) { next unless $Lj->[$p]; my $c = mm($Lj->[$p], $j);
        for my $q (0..$HTq - $p) { next unless $Ea->[$q]; $acc[$p+$q] = ($acc[$p+$q] + $c * $Ea->[$q]) % $P; } } }
    my $ia = iv($a); $E[$a] = [ map { mm($_, $ia) } @acc ]; }
  $QC{$key} = \@E; return \@E; }
# ---------- polynomials: hash key (packed exponents) -> residue ----------
sub pmul { my ($A, $B) = @_; my %C;
  for my $ka (keys %$A) { my $ca = $A->{$ka}; my $ha = hbdeg($ka);
    for my $kb (keys %$B) { next if $ha + hbdeg($kb) > $HT; my $k = $ka + $kb; $C{$k} = (($C{$k} // 0) + $ca * $B->{$kb}) % $P; } }
  delete $C{$_} for grep { !$C{$_} } keys %C; return \%C; }
sub padd_into { my ($A, $B, $s) = @_; $s //= 1; for my $k (keys %$B) { $A->{$k} = (($A->{$k} // 0) + $s * $B->{$k}) % $P; } }
sub pconst { return { 0 => $_[0] % $P }; }
# ---------- engine table ----------
my %ENG; my $have = 0;
if (defined $TABLE && $TABLE ne '-') { open my $fh, '<', $TABLE or die "no table $TABLE";
  while (my $line = <$fh>) { chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; my ($h, $r, $k, $rest) = @f;
    next unless $h == $G && $r == $N; $have = 1;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /^([\d,]+):(\S+)$/; push @{$ENG{$k}}, [ [ split /,/, $1 ], ratmod($2) ]; } } close $fh; }
sub xi_eng { my @mu = @_; return undef unless $have; my $mis = ng(iv($sig)); my $v = 0;
  for my $k (keys %ENG) { my $s = 0; for my $t (@{$ENG{$k}}) { my ($ex, $c) = @$t; my $m = $c; for my $i (0..$#mu) { $m = mm($m, pw(im($mu[$i]), $ex->[$i])); } $s = ad($s, $m); }
    $v = ad($v, mm($s, pw($mis, $k))); } return $v; }
# ---------- Pi_mu ----------
sub Pi { my $mu = shift; my $r = 1; for my $a (1..$mu-1) { $r = mm($r, ad(mm(im($mu), $tau), im($a))); } return mm($r, $ifact[$mu-1]); }
# ---------- n-cycles ----------
sub ncycles { my $n = shift; return ([0]) if $n == 1; my @out;
  my $rec; $rec = sub { my ($path, $rest) = @_; if (!@$rest) { my @s; for my $i (0..$#$path) { $s[$path->[$i]] = $path->[($i+1) % @$path]; } push @out, \@s; return; }
    for my $i (0..$#$rest) { my @r2 = @$rest; my ($x) = splice(@r2, $i, 1); $rec->([@$path, $x], \@r2); } };
  $rec->([0], [1..$n-1]); return @out; }
my @CYC = ncycles($N);
my $half = iv(2);
my $AG = $N * $JMAX + $N + $JMAX + 3;   # global bound for the Laurent order
# ---------- the iterated residue for one tuple (n >= 2) ----------
sub resid { my @mu = @_; die "n = 1 is handled separately" if $N == 1;
  my @q = map { qser($_, $AG, $HT) } @mu;
  my $tot = 0;
  for my $sg (@CYC) {
    my @kk = map { mm(im($mu[$_] + $mu[$sg->[$_]]), $half) } 0..$N-1;
    my $Pp = pconst(1);
    for my $i (0..$N-1) {
      my @nodes; my $sgn = 1;
      for my $l (0..$N-1) { my $s = $sg->[$l];
        if ($l == $i && $s > $i) { push @nodes, [ $s, ng($kk[$l]) ]; }               # 1/(v_i - (v_s - k hb))
        elsif ($s == $i && $l > $i) { push @nodes, [ $l, $kk[$l] ]; $sgn = -$sgn; }  # -1/(v_i - (v_l + k hb))
      }
      my $m = scalar @nodes;
      my %R; my $emax = 0;
      for my $k (keys %$Pp) { my $e = vdeg($k, $i); my $k2 = $k - vkey($i, $e); $R{$e}{$k2} = $Pp->{$k}; $emax = $e if $e > $emax; }
      my $Smax = $emax + $mu[$i] + 1 - $m; $Smax = 0 if $Smax < 0;
      my @c = map { my ($j, $kh) = @$_; my %p = (vkey($j, 1) => 1); $p{vkey($HB, 1)} = $kh if $kh; \%p } @nodes;
      my @h;
      if ($m == 0) { @h = (pconst(1)); push @h, {} for 1..$Smax; }
      elsif ($m == 1) { @h = (pconst(1)); push @h, pmul($h[-1], $c[0]) for 1..$Smax; }
      else { my @p2 = (pconst(1)); push @p2, pmul($p2[-1], $c[1]) for 1..$Smax;
        @h = (pconst(1)); for my $s (1..$Smax) { my $t = pmul($h[-1], $c[0]); padd_into($t, $p2[$s]); delete $t->{$_} for grep { !$t->{$_} } keys %$t; push @h, $t; } }
      my %NP;
      for my $e (keys %R) {
        my %Ge; my $top = $e + $mu[$i] + 1 - $m;
        for my $a (0..$top) { my $s = $top - $a; die "A too small ($a > $AG)" if $a > $AG; my $qa = $q[$i][$a];
          my %qp; for my $p (0..$HT) { $qp{vkey($HB, $p)} = $qa->[$p] if $qa->[$p]; } next unless %qp; next unless %{$h[$s]};
          padd_into(\%Ge, pmul(\%qp, $h[$s])); }
        next unless %Ge;
        padd_into(\%NP, pmul($R{$e}, \%Ge), $sgn > 0 ? 1 : $P - 1); }
      delete $NP{$_} for grep { !$NP{$_} } keys %NP;
      $Pp = \%NP; }
    for my $k (keys %$Pp) { die "leftover v" if $k & ((1 << ($BITS * $HB)) - 1); $tot = ad($tot, $Pp->{$k}) if hbdeg($k) == $H; } }
  return $tot; }
# ---------- main loop ----------
my @tuples; { my $rec; $rec = sub { my ($pre, $max) = @_; if (@$pre == $N) { push @tuples, [@$pre]; return; } for my $x (reverse 1..$max) { $rec->([@$pre, $x], $x); } }; $rec->([], $JMAX); }
my ($nt, $nok, $nbad, $nnd, $nrange, $nrangez) = (0, 0, 0, 0, 0, 0);
for my $t (@tuples) { my @mu = @$t; my $sum = 0; $sum += $_ for @mu;
  my $F;
  if ($N == 1) { my $mu = $mu[0]; my $q = qser($mu, $mu + 2, $H + 1);   # [x^mu]F_{g,1} = mu^-2 [hb^{2g}] q_{mu+1}
    $F = mm($q->[$mu+1][$H+1], iv(mm($mu, $mu))); }
  else { my $r = resid(@mu); my $pm = 1; $pm = mm($pm, $_) for @mu; $F = mm($r, iv($pm)); $F = ng($F) if $N % 2 == 0; }
  my $pp = 1; $pp = mm($pp, Pi($_)) for @mu;
  if ($pp == 0) { # Pi_mu = 0 (negative tau): P2 then forces [x^mu]F = 0; check that and skip the engine comparison
    printf "tau=%s g=%d mu=(%s): Pi = 0, [x^mu]F %s\n", $TAUS, $G, join(',', @mu), ($F == 0 ? "= 0 (consistent)" : "!= 0 INCONSISTENT"); $nbad++ if $F != 0; next; }
  my $sgp = ($G >= 1) ? pw($sig, $G - 1) : iv($sig);   # sigma^{g-1} (g = 0: sigma^{-1}; pw() only takes exponents >= 0)
  my $xi = mm($sgp, mm($F, iv($pp)));
  my $xe = xi_eng(@mu); my $range = ($sum <= 2*$G - 3 + $N);
  $nt++; $nrange++ if $range; $nrangez++ if $range && $xi == 0;
  my $tag; if (!defined $xe) { $tag = "(no data)"; $nnd++; } elsif ($xe == $xi) { $tag = "OK"; $nok++; } else { $tag = "MISMATCH"; $nbad++; }
  printf "tau=%s g=%d mu=(%s): ref %d  engine %s  %s%s\n", $TAUS, $G, join(',', @mu), $xi, (defined $xe ? $xe : '-'), $tag, ($range ? "  [range]" : ""); }
printf "SUMMARY tau=%s g=%d n=%d jmax=%d p=%d: tuples %d, compared %d, OK %d, MISMATCH %d, no-data %d; range points %d (ref zero at %d)\n",
  $TAUS, $G, $N, $JMAX, $P, $nt, $nok + $nbad, $nok, $nbad, $nnd, $nrange, $nrangez;
