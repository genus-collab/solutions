package RefTR;
# Referee (TR2, round 2): an independent exact (mod P) Eynard-Orantin recursion on the Hodge base curve C_tau
#   X = log x = tau log v - (tau+1) log(1+v),  omega_{0,1} = y dX with y = -v,  B = dv dv'/(v-v')^2,  branch point v = tau.
# Differences from the author's TRtau.pm/TRmod.pm: no local coordinate z with X = X(tau) - z^2/(2 sigma) and no basis e_k.
# We use e = v - tau and the global involution sigma(v) = tau + s(e), s(e) = -e + O(e^2), solved order by order from
# X(tau + s) = X(tau + e); correlators are stored in the pole basis prod_i dv_i/(v_i - tau)^{k_i}, k_i >= 2.
# Kernel: K(v0, v) = (1/2) [1/(v0 - v) - 1/(v0 - sigma v)] dv0 / ((y(v) - y(sigma v)) dX(v)),
#   omega_{g,n+1}(v0, J) = Res_{v=tau} K(v0,v) [omega_{g-1,n+2}(v, sigma v, J) + Sum' omega_{h}(v, I) omega_{g-h}(sigma v, J\I)].
# x-expansion at v = oo (t = 1/v, x = t (1+t)^(-tau-1)):
#   [x^mu] Int_oo^v dv/(v - tau)^k = -(1/mu) [t^(mu+1-k)] (1 - tau t)^(-k) (1+t)^(mu(tau+1))   (Lagrange inversion).
use strict; use warnings; use integer;
use RefMod;

sub new_curve { my ($tn,$td,$AMAX) = @_;
  my $tau = ratm($tn,$td); my $tau1 = addm($tau,1); my $sig = mulm($tau,$tau1);
  die "tau = 0 or -1\n" unless $tau && $tau1;
  my $NR = 2*$AMAX + 6; my $NN = $NR + 2;
  my ($it,$it1) = (invm($tau), invm($tau1));
  # X(tau+e) - X(tau) = Sum_{k>=2} pk[k] e^k, pk[k] = (-1)^(k+1)/k (tau^(1-k) - (tau+1)^(1-k))
  my @pk = (0,0);
  for my $k (2..$NN+2) { my $x = mulm(subm(powm($it,$k-1), powm($it1,$k-1)), invm($k)); $x = negm($x) if $k % 2 == 0; push @pk, $x; }
  die "p2 check failed\n" unless $pk[2] == negm(invm(mulm(2,$sig)));
  # s(e): [e^(j+1)] (P(s) - P(e)) = 0 determines s_j (the s_j-dependence is -2 pk[2] s_j)
  my @s = (0, negm(1), (0) x ($NN-1)); my $i2p2 = invm(mulm(2,$pk[2]));
  for my $j (2..$NN) { my $T = $j+1; my $sj = [@s[0..$j]];
    my $pw = psmul($sj,$sj,$T); my $c = mulm($pk[2], $pw->[$T] // 0);
    for my $k (3..$T) { $pw = psmul($pw,$sj,$T); $c = addm($c, mulm($pk[$k], $pw->[$T] // 0)); }
    $s[$j] = mulm(subm($c,$pk[$T]), $i2p2); }
  # verify P(s) - P(e) = O(e^(NN+2)) on the first orders
  { my $T = 12; my $sj=[@s[0..$T]]; my $pw = psmul($sj,$sj,$T); my @Ps = map { mulm($pk[2], $pw->[$_]) } 0..$T;
    for my $k (3..$T) { $pw = psmul($pw,$sj,$T); $Ps[$_] = addm($Ps[$_], mulm($pk[$k],$pw->[$_])) for 0..$T; }
    for my $i (2..$T) { die "involution check failed at order $i\n" unless $Ps[$i] == $pk[$i]; } }
  my @r  = map { $s[$_+1] } 0..$NR;                      # r = s/e
  my @sp = map { mulm($_+1, $s[$_+1]) } 0..$NR;           # s'(e)
  my $te = psmul([$tau,1],[$tau1,1],$NR);                 # (tau+e)(tau+1+e)
  my $half = invm(2);
  my $C = { tau=>$tau, tau1=>$tau1, sig=>$sig, AMAX=>$AMAX, NR=>$NR, r=>\@r, sp=>\@sp, W=>{} };
  # g_l = (1/2)(1 + r + ... + r^(l-1)) (tau+e)(tau+1+e),  l = 1..2 AMAX
  my @g; my $G = [1,(0) x $NR]; my $rp = [1,(0) x $NR];
  for my $l (1..2*$AMAX+2) { $g[$l] = psscal(psmul($G,$te,$NR),$half); $rp = psmul($rp,\@r,$NR); $G = psadd($G,$rp,$NR); }
  $C->{g} = \@g;
  # h_q = r^q s'
  my %h; my $rinv = psinv(\@r,$NR); my $cur = [@sp];
  $h{0} = $cur; for my $q (1..2*$AMAX+2) { $cur = psmul($cur,\@r,$NR); $h{$q} = $cur; }
  $cur = [@sp]; for my $q (1..2*$AMAX+2) { $cur = psmul($cur,$rinv,$NR); $h{-$q} = $cur; }
  $C->{h} = \%h;
  # B(v, sigma v) = s'/(e^2 (1-r)^2): k-series
  my $omr = [ map { negm($_) } @r ]; $omr->[0] = addm($omr->[0],1);
  $C->{kser} = psmul(\@sp, pspow($omr,-2,$NR), $NR);
  $C->{prod} = {};
  return $C; }

# RES(l,p,q) = Res_{e=0} G_l(e) e^p s(e)^q s'(e) de,  G_l = (1/2)(e^l - s^l)/((s - e) X'(tau+e)) = e^(l-2) g_l
sub RES { my ($C,$l,$p,$q) = @_; my $o = 1 - $l - $p - $q; return 0 if $o < 0;
  die "RES order $o > NR\n" if $o > $C->{NR};
  my $key = "$l,$q"; my $pr = $C->{prod}{$key};
  unless ($pr) { my $gl = $C->{g}[$l] or die "g_$l missing\n"; my $hq = $C->{h}{$q} or die "h_$q missing\n";
    $pr = psmul($gl,$hq,$C->{NR}); $C->{prod}{$key} = $pr; }
  $pr->[$o] }

sub stable { my ($g,$n) = @_; 2*$g - 2 + $n > 0 }

sub compute { my ($C,$g,$n1) = @_; my $n = $n1 - 1; my $W = $C->{W}; my %acc; my $AM = $C->{AMAX};
  my $add = sub { my ($k,$v) = @_; return unless $v; my $o = $acc{$k} // 0; $acc{$k} = addm($o,$v); };
  if ($g >= 1) {
    if ($g == 1 && $n == 0) { for my $l (1..3) { my $o = 3 - $l; $add->("".($l+1), psmul($C->{g}[$l],$C->{kser},$C->{NR})->[$o]); } }
    else { my $Wa = $W->{($g-1).",".($n+2)} or die "missing ".($g-1).",".($n+2)."\n";
      for my $key (keys %$Wa) { my $c = $Wa->{$key}; my @k = split /,/, $key; my $a = shift @k; my $b = shift @k; my $rest = join(",",@k);
        for my $l (1..$a+$b+1) { my $r = RES($C,$l,-$a,-$b) or next; $add->(($rest eq '' ? ($l+1) : ($l+1).",".$rest), mulm($c,$r)); } } } }
  for my $h (0..$g) { for my $mask (0..(1<<$n)-1) {
      my @I = grep { ($mask>>$_) & 1 } 0..$n-1; my @Ic = grep { !(($mask>>$_) & 1) } 0..$n-1;
      my ($na,$nb) = (1+@I, 1+@Ic); my $hb = $g - $h;
      next if ($h == 0 && $na == 1) || ($hb == 0 && $nb == 1);
      next unless ($h == 0 && $na == 2) || stable($h,$na);
      next unless ($hb == 0 && $nb == 2) || stable($hb,$nb);
      my $mk = sub { my ($hh,$nn) = @_; my @L;
        if ($hh == 0 && $nn == 2) { for my $m (0..$AM) { push @L, [$m, md($m+1), [$m+2]]; } }
        else { my $Wx = $W->{"$hh,$nn"} or die "missing $hh,$nn\n"; for my $key (keys %$Wx) { my @k = split /,/, $key; my $a = shift @k; push @L, [-$a, $Wx->{$key}, \@k]; } }
        \@L };
      my $LA = $mk->($h,$na); my $LB = $mk->($hb,$nb);
      for my $A (@$LA) { my ($p,$ca,$kA) = @$A;
        for my $B (@$LB) { my ($q,$cb,$kB) = @$B; my $lmax = 1 - $p - $q; next if $lmax < 1;
          my @kJ; my ($ia,$ib) = (0,0); for my $j (0..$n-1) { if (($mask>>$j) & 1) { push @kJ, $kA->[$ia++]; } else { push @kJ, $kB->[$ib++]; } }
          my $rest = join(",",@kJ); my $cab = mulm($ca,$cb);
          for my $l (1..$lmax) { my $r = RES($C,$l,$p,$q) or next; $add->(($rest eq '' ? ($l+1) : ($l+1).",".$rest), mulm($cab,$r)); } } } } }
  for my $k (keys %acc) { delete $acc{$k} unless $acc{$k}; }
  $W->{"$g,$n1"} = \%acc; }

# run(C, G, N): all (g,n) needed for (G, N) and below: n <= N + G - g
sub run { my ($C,$G,$N,%o) = @_; my @ord;
  for my $chi (1..2*$G-2+$N) { for my $g (0..$G) { my $n = $chi - 2*$g + 2; next if $n < 1 || $n > $N + $G - $g; push @ord, [$g,$n]; } }
  for my $gn (@ord) { my ($g,$n) = @$gn; my $t0 = time;
    if ($g == 0 && $n == 3) { compute($C,0,3); } elsif ($g == 1 && $n == 1) { compute($C,1,1); } else { compute($C,$g,$n); }
    print STDERR "  ($g,$n): ".scalar(keys %{$C->{W}{"$g,$n"}})." terms, ".(time-$t0)."s\n" if $o{verbose}; } }

sub xc { my ($C,$k,$mu) = @_; my $key = "$k,$mu"; return $C->{xc}{$key} if exists $C->{xc}{$key};
  my $tau = $C->{tau}; my $M = mulm($mu, $C->{tau1}); my $top = $mu + 1 - $k; my $s = 0;
  for my $i (0..$top) { $s = addm($s, mulm(mulm(binm(md($k+$i-1),$i), powm($tau,$i)), binm($M,$top-$i))) if $top >= 0; }
  my $r = negm(mulm($s, invm($mu))); $C->{xc}{$key} = $r; $r }
sub xcoef { my ($C,$g,@mu) = @_; my $Wg = $C->{W}{"$g,".scalar(@mu)} or die "no W $g,".scalar(@mu)."\n"; my $tot = 0;
  for my $key (keys %$Wg) { my @k = split /,/, $key; my $v = $Wg->{$key}; for my $i (0..$#mu) { $v = mulm($v, xc($C,$k[$i],$mu[$i])); last unless $v; } $tot = addm($tot,$v); } $tot }
# P2: Xi_h(mu) = sigma^(h-1) [x^mu] F_{h,r} / prod C(mu_i (tau+1) - 1, mu_i - 1)
sub xi_p2 { my ($C,$g,@mu) = @_; my $v = mulm(xcoef($C,$g,@mu), powm($C->{sig},$g-1));
  for my $m (@mu) { my $b = binm(subm(mulm($m,$C->{tau1}),1), $m-1); return undef unless $b; $v = mulm($v, invm($b)); } $v }
1;
