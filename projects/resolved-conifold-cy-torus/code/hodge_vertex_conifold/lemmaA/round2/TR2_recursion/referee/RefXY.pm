package RefXY;
# Referee (TR2, round 2): the ABDKS formula (17) for C_tau, integrated against x^(-mu) (integration by parts in every
# variable, Res x^(-mu) (d 1/dX)^r eta = mu^r Res x^(-mu) eta, hence Sum_r mu^r [u^r] f(u) = f(mu)):
#   n = 1:   [x^j] F_{g,1} = -(1/j^2) [hbar^(2g)] Res_{v=oo} Psi_j(v) dv
#   n >= 2:  [x^mu] F_{g,n} = -(prod 1/mu_i) [hbar^(2g-2+n)] Res_{v_1=oo}..Res_{v_n=oo} prod_i Psi_{mu_i}(v_i)
#                              Sum_{sigma n-cycle} prod_i dv_i / (v_i - v_sigma(i) + (mu_i + mu_sigma(i)) hbar/2)
#   with the iterated residues computed in C((t_1))((t_2))...((t_n)), t_i = 1/v_i (|t_n| << ... << |t_1|),
#   Psi_mu = x^(-mu) E(u = mu),  E = exp(-u (S(u hbar d_v) xhat - X)),  xhat = Sum_a (1/S(alpha_a hbar d_v)) log(v - a)/alpha_a,
#   (a, 1/alpha_a) = (0, tau), (-1, -(tau+1)),  S(u) = u^(-1)(e^(u/2) - e^(-u/2)).
# Expansion: Psi_mu(t) = t^(-mu) (1+t)^(mu(tau+1)) exp( Sum_k (hbar t)^(2k) [A_k + B_k (1+t)^(-2k)] ),
#   A_k = (2k-1)! mu tau c_k(mu, 1/tau), B_k = -(2k-1)! mu (tau+1) c_k(mu, -1/(tau+1)), c_k(u,al) = [y^(2k)] S(u y)/S(al y).
# Written by the referee; it shares no code with the author's XY.pm / two_point.pl / xy_point.pl.
use strict; use warnings; use integer;
use RefMod;

sub new { my ($tn,$td,$G,$KT) = @_;
  my $tau = ratm($tn,$td); my $tau1 = addm($tau,1);
  my $X = { tau=>$tau, tau1=>$tau1, sig=>mulm($tau,$tau1), G=>$G, KT=>$KT, psi=>{} };
  my $K = $G + 3; my @fac = (1); push @fac, mulm($fac[-1],$_) for 1..2*$K+2;
  my @s = map { invm(mulm(powm(4,$_), $fac[2*$_+1])) } 0..$K;       # S(y) = Sum s_k y^(2k)
  my @si = (1); for my $n (1..$K) { my $a = 0; $a = addm($a, mulm($s[$_],$si[$n-$_])) for 1..$n; push @si, negm($a); }  # 1/S
  $X->{s} = \@s; $X->{si} = \@si; $X->{fac} = \@fac;
  $X->{al} = [ invm($tau), negm(invm($tau1)) ];
  $X }
sub ck { my ($X,$u,$al,$k) = @_; my $c = 0;
  for my $p (0..$k) { $c = addm($c, mulm(mulm($X->{s}[$p], powm($u,2*$p)), mulm($X->{si}[$k-$p], powm($al,2*($k-$p))))); } $c }
# P_m(t), m = 0..G+1: [hbar^(2m)] Psi_mu = t^(2m-mu) P_m(t)
sub psi { my ($X,$mu) = @_; return $X->{psi}{$mu} if $X->{psi}{$mu};
  my $KT = $X->{KT}; my $MM = $X->{G} + 2; my $tau = $X->{tau}; my $tau1 = $X->{tau1};
  my @Q = (undef);
  for my $k (1..$MM) { my $f = $X->{fac}[2*$k-1];
    my $A = mulm(mulm($f, mulm($mu,$tau)), ck($X,$mu,$X->{al}[0],$k));
    my $B = negm(mulm(mulm($f, mulm($mu,$tau1)), ck($X,$mu,$X->{al}[1],$k)));
    my $pw = pspow([1,1,(0) x ($KT-1)], -2*$k, $KT); my $q = psscal($pw,$B); $q->[0] = addm($q->[0],$A); $Q[$k] = $q; }
  my @E = ([1,(0) x $KT]);
  for my $m (1..$MM) { my $acc = [(0) x ($KT+1)];
    for my $k (1..$m) { my $t = psmul($Q[$k],$E[$m-$k],$KT); $acc = psadd($acc, psscal($t, md($k)), $KT); }
    $E[$m] = psscal($acc, invm($m)); }
  my $M = mulm($mu,$tau1); my $onep = [ map { binm($M,$_) } 0..$KT ];
  my @Pm = map { psmul($onep,$E[$_],$KT) } 0..$MM;
  $X->{psi}{$mu} = \@Pm; \@Pm }

sub xF1 { my ($X,$g,$j) = @_; my $P = psi($X,$j); my $o = $j + 1 - 2*$g; return 0 if $o < 0; die "KT too small\n" if $o > $X->{KT};
  mulm($P->[$g][$o], invm(mulm($j,$j))) }

sub cycles { my ($n) = @_; my @perms = ([]); for my $e (1..$n-1) { my @np; for my $pp (@perms) { for my $pos (0..@$pp) { my @q = @$pp; splice(@q,$pos,0,$e); push @np, \@q; } } @perms = @np; }
  my @cyc; for my $pp (@perms) { my @ord = (0,@$pp); my @img; for my $i (0..$n-1) { $img[$ord[$i]] = $ord[($i+1) % $n]; } push @cyc, \@img; } @cyc }

# [x^mu] F_{g,n}, n >= 1
sub xF { my ($X,$g,@mu) = @_; my $n = @mu; return xF1($X,$g,$mu[0]) if $n == 1;
  my $H = 2*$g - 2 + $n; my $KT = $X->{KT}; my @P = map { psi($X,$_) } @mu; my $tot = 0;
  for my $sg (cycles($n)) {
    my @edges = map { my ($i,$j) = ($_, $sg->[$_]); { i=>$i, j=>$j, s=>($i>$j?$i:$j), b=>($i>$j?$j:$i), sgn=>($i>$j?1:0), c=>mulm(md($mu[$i]+$mu[$j]), invm(2)) } } 0..$n-1;
    my @sE; my @bE; for my $e (0..$#edges) { push @{$sE[$edges[$e]{s}]}, $e; push @{$bE[$edges[$e]{b}]}, $e; }
    # hbar distributions: m_v (Psi, weight 2) and me_e (edges, weight 1), total H
    my @dist; my $rec; $rec = sub { my ($pre,$left,$idx) = @_;
      if ($idx == 2*$n) { push @dist, [@$pre] if $left == 0; return; }
      my $w = $idx < $n ? 2 : 1; for (my $m = 0; $m*$w <= $left; $m++) { $rec->([@$pre,$m], $left - $m*$w, $idx+1); } };
    $rec->([], $H, 0);
    for my $d (@dist) { my @m = @$d[0..$n-1]; my @me = @$d[$n..2*$n-1];
      my $ecoef = 1; for my $e (0..$#edges) { $ecoef = mulm($ecoef, powm(negm($edges[$e]{c}), $me[$e])); }
      next unless $ecoef;
      my @l = (0) x $n; my $sum = 0;
      my $go; $go = sub { my ($v,$acc) = @_;
        if ($v < 0) { $sum = addm($sum,$acc); return; }
        my $base = 1 - 2*$m[$v] + $mu[$v]; for my $e (@{$bE[$v] // []}) { $base += $l[$e]; }
        my @se = @{$sE[$v] // []}; for my $e (@se) { $base -= $me[$e] + 1; }
        return if $base < 0;
        my $fin = sub { my ($kv,$acc2) = @_; die "KT too small ($kv)\n" if $kv > $KT; my $c = $P[$v][$m[$v]][$kv] or return; $go->($v-1, mulm($acc2,$c)); };
        my $efac = sub { my ($e) = @_; my $c = binm(md($me[$e]+$l[$e]), $l[$e]); $c = negm($c) if !$edges[$e]{sgn} && ($me[$e] % 2 == 0); $c };
        if (@se == 0) { $fin->($base,$acc); }
        elsif (@se == 1) { my $e = $se[0]; for my $le (0..$base) { $l[$e] = $le; $fin->($base-$le, mulm($acc,$efac->($e))); } }
        else { my ($e1,$e2) = @se; for my $l1 (0..$base) { $l[$e1] = $l1; my $f1 = mulm($acc,$efac->($e1));
            for my $l2 (0..$base-$l1) { $l[$e2] = $l2; $fin->($base-$l1-$l2, mulm($f1,$efac->($e2))); } } } };
      $go->($n-1, $ecoef);
      $tot = addm($tot,$sum); } }
  # -(prod 1/mu)(-1)^n * tot
  my $pre = 1; $pre = mulm($pre, invm($_)) for @mu; $pre = negm($pre) if $n % 2 == 0;   # -(-1)^n = +1 for odd n, -1 for even n
  mulm($pre,$tot) }
1;
