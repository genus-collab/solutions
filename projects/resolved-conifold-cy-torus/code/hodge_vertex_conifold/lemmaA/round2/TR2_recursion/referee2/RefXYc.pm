package RefXYc;
# Referee 2: closed formulas of the TR2 report, implemented independently, modulo a prime p.
#  Psi_j(v) = exp(-j S(j hbar d) xhat) = x^(-j) exp( sum_{k>=1} hbar^(2k) sum_a (j/alpha_a) (2k-1)! c_k(j,alpha_a) (v-a)^(-2k) ),
#  a in {0,-1}, 1/alpha_0 = tau, 1/alpha_{-1} = -(tau+1), c_k(u,al) = [t^(2k)] S(ut)/S(al t), S(t) = 2 sinh(t/2)/t.
#  In w = 1/v:  [hbar^(2k)] Psi_j = w^(2k-j) P_{j,k}(w),  P_{j,k}(w) = sum_{m+q=k} AA_m BB_q (1+w)^(M-2q),  M = j(tau+1),
#  AA = exp-coefficients of sum a_k t^(2k), BB = exp-coefficients of sum b_k t^(2k) (computed from first principles here).
#  General n (report s.5):  [x^mu] F_{g,n} = -(prod 1/mu_i) [hbar^(2g-2+n)] IterRes prod_i Psi_{mu_i}(v_i)
#        * sum_{sigma n-cycle} prod_i dv_i/(v_i - v_sigma(i) + (mu_i+mu_sigma(i)) hbar/2),   domain |v_n| > ... > |v_1|.
use strict; use warnings; use integer;
use RefTR;
*md = \&RefTR::md; *inv = \&RefTR::inv; *pw = \&RefTR::pw; *gbin = \&RefTR::gbin;

sub new { my ($class,$tau,$p,$G)=@_; my $S = bless { p=>$p, tau=>$tau, t1=>($tau+1)%$p, G=>$G }, $class;
  # s_m = [t^(2m)] S(t) = 1/(4^m (2m+1)!),  sp_m = [t^(2m)] 1/S(t)
  my @s; my $f = 1; for my $m (0..$G+2) { $f = 1; for my $i (1..2*$m+1) { $f = ($f*$i) % $p; } push @s, inv(pw(4,$m,$p)*$f % $p, $p); }
  my @sp = (1); for my $n (1..$G+2) { my $acc = 0; for my $k (1..$n) { $acc = ($acc + $s[$k]*$sp[$n-$k]) % $p; } push @sp, md(-$acc,$p); }
  $S->{s} = \@s; $S->{sp} = \@sp; return $S; }
sub ck { my ($S,$u,$al,$k)=@_; my $p = $S->{p}; my $r = 0;
  for my $m (0..$k) { my $q = $k - $m; $r = ($r + pw($u,2*$m,$p) * $S->{s}[$m] % $p * $S->{sp}[$q] % $p * pw($al,2*$q,$p)) % $p; } $r }
sub expser { my ($S,$c,$K)=@_; my $p = $S->{p}; my @F = (1);   # F(T) = exp(sum_{k>=1} c[k] T^k)
  for my $n (1..$K) { my $s = 0; for my $k (1..$n) { $s = ($s + $k * $c->[$k] % $p * $F[$n-$k]) % $p; } push @F, $s * inv($n,$p) % $p; } \@F }
sub AB { my ($S,$j)=@_; return @{$S->{AB}{$j}} if $S->{AB}{$j}; my $p = $S->{p}; my $tau = $S->{tau}; my $t1 = $S->{t1}; my $G = $S->{G} + 1;
  my $al0 = inv($tau,$p); my $al1 = md(-inv($t1,$p),$p); my (@a,@b) = ((0),(0)); @a = (0); @b = (0);
  for my $k (1..$G) { my $f = 1; for my $i (1..2*$k-1) { $f = ($f*$i) % $p; }
    push @a, $j % $p * $tau % $p * $f % $p * $S->ck($j % $p,$al0,$k) % $p;
    push @b, md(-($j % $p * $t1 % $p * $f % $p * $S->ck($j % $p,$al1,$k) % $p), $p); }
  my $A = $S->expser(\@a,$G); my $B = $S->expser(\@b,$G); $S->{AB}{$j} = [$A,$B]; ($A,$B) }
# Norlund closed forms (report (5.2)): AA_m = (-N)_{2m} [y^{2m}] S1^(1-N) tau^(-2m), BB_q = (M)_{2q} [y^{2q}] S1^(1+M) (tau+1)^(-2q)
sub s1pow { my ($S,$e,$K)=@_; my $p = $S->{p}; # S1(y)^e = (1/S(y))^e = S(y)^(-e); coefficients in T = y^2
  my $c = -$e; my @f = @{$S->{s}}; my @g = (1); # g = f^c: n g_n = sum_k (c k - (n-k)) f_k g_{n-k}
  for my $n (1..$K) { my $acc = 0; for my $k (1..$n) { $acc = ($acc + md($c*$k - ($n-$k), $p) * $f[$k] % $p * $g[$n-$k]) % $p; } push @g, $acc * inv($n,$p) % $p; } \@g }
sub falling { my ($a,$n,$p)=@_; my $r = 1; for my $i (0..$n-1) { $r = $r * md($a - $i,$p) % $p; } $r }
sub AB_norlund { my ($S,$j)=@_; my $p = $S->{p}; my $tau = $S->{tau}; my $t1 = $S->{t1}; my $G = $S->{G} + 1;
  my $N = $j % $p * $tau % $p; my $M = $j % $p * $t1 % $p;
  my $sa = $S->s1pow(md(1-$N,$p),$G); my $sb = $S->s1pow(md(1+$M,$p),$G); my (@A,@B);
  for my $m (0..$G) { push @A, falling(md(-$N,$p),2*$m,$p) * $sa->[$m] % $p * pw($tau,-2*$m,$p) % $p;
                      push @B, falling($M,2*$m,$p) * $sb->[$m] % $p * pw($t1,-2*$m,$p) % $p; } (\@A,\@B) }
# P_{j,k}(w) up to w^L
sub Pjk { my ($S,$j,$k,$L)=@_; my $key = "$j,$k,$L"; return $S->{P}{$key} if $S->{P}{$key}; my $p = $S->{p};
  my ($A,$B) = $S->AB($j); my $M = $j % $p * $S->{t1} % $p; my @P = (0) x ($L+1);
  for my $q (0..$k) { my $m = $k - $q; my $cf = $A->[$m] * $B->[$q] % $p; next unless $cf; my $e = md($M - 2*$q,$p);
    for my $i (0..$L) { $P[$i] = ($P[$i] + $cf * gbin($e,$i,$p)) % $p; } }
  $S->{P}{$key} = \@P; \@P }
# one-point (5.1)-(5.2): [x^j] F_{g,1} = (1/j^2) P_{j,g}[j+1-2g]
sub onept { my ($S,$g,$j)=@_; my $p = $S->{p}; my $i = $j + 1 - 2*$g; return 0 if $i < 0; my $P = $S->Pjk($j,$g,$i); $P->[$i] * inv($j*$j % $p,$p) % $p }
# two-point (5.4) with the domain |v2| > |v1|
sub twopt { my ($S,$g,$j1,$j2)=@_; my $p = $S->{p}; my $c = ($j1+$j2) % $p * inv(2,$p) % $p; my $tot = 0;
  for my $k1 (0..$g) { for my $k2 (0..$g-$k1) { my $m = $g - $k1 - $k2; my $c2m = pw($c,2*$m,$p);
      my $Lmax = $j2 - 1 - 2*$k2 - 2*$m; next if $Lmax < 0;
      my $P1 = $S->Pjk($j1,$k1,$j1+1+$Lmax+2*$g+2); my $P2 = $S->Pjk($j2,$k2,$j2+2*$g+2);
      for my $l (0..$Lmax) { my $i1 = $j1 + 1 - 2*$k1 + $l; next if $i1 < 0; my $i2 = $j2 - 1 - 2*$k2 - 2*$m - $l;
        $tot = ($tot + $c2m * gbin(2*$m+1+$l,$l,$p) % $p * $P1->[$i1] % $p * $P2->[$i2]) % $p; } } }
  $tot * inv($j1*$j2 % $p,$p) % $p }
# general n: integrated formula over n-cycles, domain |v_n| > ... > |v_1| (indices 0..n-1 here)
sub npt { my ($S,$g,@mu)=@_; my $p = $S->{p}; my $n = @mu; my $H = 2*$g - 2 + $n; my $tot = 0;
  my @perms = ([]); for my $e (1..$n-1) { my @np; for my $q (@perms) { for my $pos (0..@$q) { my @r = @$q; splice(@r,$pos,0,$e); push @np, \@r; } } @perms = @np; }
  my @cycles; for my $q (@perms) { my @ord = (0,@$q); my @img; for my $i (0..$n-1) { $img[$ord[$i]] = $ord[($i+1) % $n]; } push @cycles, \@img; }
  my $Lbig = 0; $Lbig += $_ + 2 for @mu; $Lbig += 2*$H + 4;
  for my $sg (@cycles) {
    my @F; for my $i (0..$n-1) { my $j = $sg->[$i]; my ($o,$in) = $i > $j ? ($i,$j) : ($j,$i); my $sgn = ($i == $o) ? 1 : -1;
      push @F, { o=>$o, in=>$in, s=>$sgn, c=>(($mu[$i]+$mu[$j]) % $p) * inv(2,$p) % $p }; }
    # factors sorted by outer index descending, for the l-enumeration
    my @ordF = sort { $F[$b]{o} <=> $F[$a]{o} } 0..$#F;
    # enumerate k (per variable) and m (per factor)
    my @kv; my $reck; $reck = sub { my ($i,$left,$acc)=@_; if ($i == $n) { push @kv, [@$acc]; return; } for my $k (0..$left) { $reck->($i+1,$left-$k,[@$acc,$k]); } };
    $reck->(0, $H/2, []);
    for my $k (@kv) { my $sk = 0; $sk += $_ for @$k; my $Mt = $H - 2*$sk; next if $Mt < 0;
      my @Pv = map { $S->Pjk($mu[$_],$k->[$_],$Lbig) } 0..$n-1;
      my @mv; my $recm; $recm = sub { my ($f,$left,$acc)=@_; if ($f == $n-1) { push @mv, [@$acc,$left]; return; } for my $m (0..$left) { $recm->($f+1,$left-$m,[@$acc,$m]); } };
      $recm->(0,$Mt,[]);
      for my $m (@mv) { my $pref = 1;
        for my $f (0..$#F) { $pref = $pref * pw(md(-$F[$f]{c},$p),$m->[$f],$p) % $p; $pref = md(-$pref,$p) if $F[$f]{s} < 0 && ($m->[$f]+1) % 2; }
        next unless $pref;
        # E_i base: sum over factors with outer i of (m_f + 1)
        my @E0 = (0) x $n; for my $f (0..$#F) { $E0[$F[$f]{o}] += $m->[$f] + 1; }
        # recursive l-enumeration in the order @ordF; budget for variable i: E_i <= mu_i + 1 - 2 k_i
        my $recl; $recl = sub { my ($t,$E,$cf)=@_;
          if ($t == @ordF) { my $v = $cf; for my $i (0..$n-1) { my $idx = $mu[$i] + 1 - 2*$k->[$i] - $E->[$i]; return if $idx < 0 || $idx > $Lbig; $v = $v * $Pv[$i][$idx] % $p; return unless $v; }
            $tot = md($tot + ($n % 2 ? -$v : $v), $p); return; }   # (-1)^n from the n residues  -Res
          my $f = $ordF[$t]; my ($o,$in) = ($F[$f]{o},$F[$f]{in}); my $mf = $m->[$f];
          # outer budget: once all factors with outer o are chosen the constraint at o is final; use a safe bound
          my $bud = $mu[$o] + 1 - 2*$k->[$o] - $E->[$o]; # remaining room at o (E already includes m+1 of all its outer factors)
          # inner-l of o (factors where o is inner, partner outer > o) have already been chosen (descending order), so E[o] is final apart from its own outer l's
          for (my $l = 0; $l <= $bud; $l++) { my @E2 = @$E; $E2[$o] += $l; $E2[$in] -= $l;
            $recl->($t+1, \@E2, $cf * gbin($mf+$l,$l,$p) % $p); } };
        $recl->(0,\@E0,$pref); } } }
  my $den = 1; $den = $den * $_ % $p for @mu; md(-$tot * inv($den,$p), $p) }
1;
