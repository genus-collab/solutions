use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RegLib; use Math::BigRat;
# Test 7, part B: the stored identity row solves the T_F Picard-Fuchs equation in WKB form, exactly, as an identity in Q(s).
# BS eq:pf on T_F (alpha_j = alpha = E/3): [(z th + alpha)^3 + z^3 y (3th)(3th+1)(3th+2)] I = 0, th = y d/dy.
# WKB: I = e^{u/z} S, th u = Y = L.  On T_F, with s = 1/(L+alpha):  Y = 1/s - alpha,  y = -1/(27 q^3), q = 1 - alpha s,
# th = v d/ds, v = q/(3 alpha).  Conjugated operator: Op = (a + z th)^3 - q^-3 [N^3 + z N^2 + (2/9) z^2 N], a = 1/s, N = Y + z th.
# We apply q^3 Op (Laurent-polynomial coefficients) to S(z) = Sum_k S_k z^k (S_k from form.pl) and check every z^m, m <= K.
# Controls: (i) perturb S_c by delta s^j; (ii) structural: Op1 (s sigma) = -sigma' (transport), q^3 Op2, q^3 Op3 on the basis s^n.
my ($E1,$E2,$PERT) = @ARGV; $E1 = QQ($E1); $E2 = QQ($E2); my $E = $E1+$E2; my $al = $E/3;
my @S; open my $in, '<', "$FindBin::Bin/Sk_${E1}_${E2}.txt" or die; while (<$in>) { chomp; my ($k,@t) = split ' '; my %p; for (@t) { my ($e,$c) = split /:/; $p{$e} = QQ($c) } $S[$k] = \%p } close $in;
my $K = $#S;
if (defined $PERT) { my ($c,$j,$d) = split /,/, $PERT; $S[$c] = lp_add($S[$c], {$j => QQ($d)}); print "CONTROL: S_$c += $d s^$j\n"; }
my $a = {-1 => QQ(1)}; my $Y = {-1 => QQ(1), 0 => -$al}; my $q = {0 => QQ(1), 1 => -$al};
my $q3 = lp_mul($q, lp_mul($q,$q));
sub opX { my ($c,$F) = @_; # (c + z th) on a z-series
  my @R; for my $m (0..$K) { $R[$m] = lp_mul($c, $F->[$m]); $R[$m] = lp_add($R[$m], lp_theta($F->[$m-1],$al)) if $m >= 1; } return \@R; }
my $A1 = opX($a,\@S); my $A2 = opX($a,$A1); my $A3 = opX($a,$A2);
my $N1 = opX($Y,\@S); my $N2 = opX($Y,$N1); my $N3 = opX($Y,$N2);
my $allzero = 1;
for my $m (0..$K) {
  my $r = lp_add(lp_mul($q3,$A3->[$m]), lp_scale($N3->[$m],QQ(-1)));
  $r = lp_add($r, lp_scale($N2->[$m-1],QQ(-1))) if $m >= 1;
  $r = lp_add($r, lp_scale($N1->[$m-2],QQ(-2)/9)) if $m >= 2;
  my $z = (keys %$r) ? 0 : 1; $allzero &&= $z;
  print "[z^$m] q^3 Op S = ", ($z ? "0" : lp_str($r)), "\n";
}
print "PF/WKB identity through z^$K: ", ($allzero ? "EXACT ZERO" : "FAILS"), "\n";
exit if defined $PERT;
# structural operator data: Op1 = coefficient of z^1, Op2, Op3 as operators on a single function f
sub th { lp_theta($_[0],$al) }
sub Op1 { my $f=shift; my $fp = th($f);
  my $A = lp_add(lp_scale(lp_mul(lp_mul($a,$a),$fp),QQ(3)), lp_scale(lp_mul(lp_mul($a,th($a)),$f),QQ(3)));
  my $B = lp_add(lp_scale(lp_mul(lp_mul($Y,$Y),$fp),QQ(3)), lp_scale(lp_mul(lp_mul($Y,th($Y)),$f),QQ(3)), lp_mul(lp_mul($Y,$Y),$f));
  return (lp_mul($q3,$A), $B); } # q^3 Op1 f = q3*A - B
sub Op2 { my $f=shift; my $fp = th($f); my $fpp = th($fp);
  my $A = lp_add(lp_scale(lp_mul($a,$fpp),QQ(3)), lp_scale(lp_mul(th($a),$fp),QQ(3)), lp_mul(th(th($a)),$f));
  my $B = lp_add(lp_scale(lp_mul($Y,$fpp),QQ(3)), lp_scale(lp_mul(th($Y),$fp),QQ(3)), lp_mul(th(th($Y)),$f),
                 lp_scale(lp_mul($Y,$fp),QQ(2)), lp_mul(th($Y),$f), lp_scale(lp_mul($Y,$f),QQ(2)/9));
  return lp_add(lp_mul($q3,$A), lp_scale($B,QQ(-1))); }
sub Op3 { my $f=shift; my $f1 = th($f); my $f2 = th($f1); my $f3 = th($f2);
  return lp_add(lp_mul($q3,$f3), lp_scale(lp_add($f3,$f2,lp_scale($f1,QQ(2)/9)),QQ(-1))); }
# divide a Laurent poly by q^3 = (1 - alpha s)^3 if possible: returns (quotient, remainder-is-zero)
sub divq3 { my $p = shift; return ({},1) unless keys %$p; my $lo = lp_minexp($p); my $hi = lp_maxexp($p);
  my %r = %$p; my %quo; my $c3 = -$al**3; # leading coeff of q^3 in s^3
  for (my $e = $hi; $e >= $lo+3; $e--) { next unless exists $r{$e} && $r{$e} != 0; my $c = $r{$e}/$c3; $quo{$e-3} = $c;
    my $sub = lp_scale(lp_shift($q3,$e-3),$c); for my $kk (keys %$sub) { $r{$kk} = ($r{$kk}//QQ(0)) - $sub->{$kk} } }
  lp_clean(\%r); return (\%quo, (keys %r) ? 0 : 1); }
print "\nStructure of the transport recursion (exact, on the basis s^n):\n";
for my $n (1..6) {
  my $sigma = {$n-1 => QQ(1)}; my ($P,$Bq) = Op1(lp_shift($sigma,1));
  my $op1 = lp_add($P, lp_scale($Bq,QQ(-1)));  # = q^3 Op1(s sigma)
  my $want = lp_scale(lp_mul($q3, {$n-2 => QQ($n-1)}),QQ(-1)); # q^3 * (-sigma')
  my $chk = lp_add($op1, lp_scale($want,QQ(-1)));
  my ($o2q,$o2ok) = divq3(Op2({$n => QQ(1)})); my ($o3q,$o3ok) = divq3(Op3({$n => QQ(1)}));
  print " n=$n: Op1(s * s^", $n-1, ") + (s^", $n-1, ")' = ", ((keys %$chk)?"NONZERO":"0"),
        "; q^3 | q^3 Op2(s^$n): ", ($o2ok?'yes':'no'), "; q^3 | q^3 Op3(s^$n): ", ($o3ok?'yes':'no'), "\n";
  print "     Op3(s^$n) = ", ($o3ok ? lp_str($o3q) : 'non-polynomial'), "\n";
  print "     Op2(s^$n) = ", ($o2ok ? lp_str($o2q) : 'has poles at q=0 (orbifold)'), "\n";
}
# q-pole part of Op2 on a general element s*sigma: is Op2 S_{k-1} + Op3 S_{k-2} polynomial for the stored rows?
print "\nSources src_k = -(Op2 S_{k-1} + Op3 S_{k-2}) for the stored rows (transport: sigma_k' = -src_k... sign as computed):\n";
for my $k (1..$K) {
  my $src = Op2($S[$k-1]); $src = lp_add($src, Op3($S[$k-2])) if $k >= 2;
  my ($qq,$ok) = divq3($src);
  my $sig = {}; for my $e (keys %{$S[$k]}) { $sig->{$e-1} = $S[$k]{$e} }
  my $sigp = {}; for my $e (keys %$sig) { next if $e == 0; $sigp->{$e-1} = $sig->{$e}*$e } # d sigma/ds
  my $chk = $ok ? lp_add($qq, lp_scale($sigp,QQ(-1))) : {x=>1};
  printf " k=%d: q^3 divides q^3*(Op2 S_{k-1}+Op3 S_{k-2}): %s; quotient s-exponents %s..%s; (Op2 S_{k-1}+Op3 S_{k-2}) - d sigma_k/ds = %s\n",
     $k, ($ok?'yes':'no'), ($ok? (lp_minexp($qq)//'-') : '-'), ($ok ? (lp_maxexp($qq)//'-') : '-'), ((keys %$chk)?'NONZERO':'0');
}
