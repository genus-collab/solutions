#!/usr/bin/perl
# Referee check for S2: the interface (I1)-(I4) built from the KEY FORMULA, against the exact tables of leaf_structure.md
# (fitted there from the stored z^9 factors, i.e. from the code's EdN/DilN directly).
#   dilaton:  [t^k] (-psi g(psi)) = delta_k(psi), k = 0..3  (k = 0: q_0 = c - c^2 lam_1)
#   loops:    [t^k] rho^2 E_11 from (I4)  [-x c^3/t^2 - <sigma(psi), sigma(psi')>]/(psi + psi')   (uses K_1 = 1)
#             and from the flow (V8) with Theta = -x^3 tau s*^2, against the e^(k)_{il} tables, k = 0, 1, 2
#   E11 leg:  [t^1] rho f_{1,0} = (1-x)^2 x^3/c^3 [-(1/12) lam_1 + (2/3) lam_2]
#   symplectic constant from the closed form: t^2 <sigma(psi), sigma(-psi)> = -x c^3
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $N = @ARGV ? $ARGV[0] : 8;           # t-order
my ($ntest, $nfail) = (0, 0);
sub ok { my ($c, $msg) = @_; $ntest++; unless ($c) { $nfail++; print "FAIL: $msg\n"; } }
sub z0 { [ map { R(0) } 0..$N ] }
sub c0 { my $a = z0(); $a->[0] = R($_[0]); $a }
sub add { my ($a,$b)=@_; [ map { $a->[$_] + $b->[$_] } 0..$N ] }
sub sub_ { my ($a,$b)=@_; [ map { $a->[$_] - $b->[$_] } 0..$N ] }
sub sc { my ($s,$a)=@_; [ map { $a->[$_]*$s } 0..$N ] }
sub mul { my ($a,$b)=@_; my $c = z0(); for my $i (0..$N) { next if $a->[$i]->is_zero; for my $k (0..$N-$i) { $c->[$i+$k] += $a->[$i]*$b->[$k]; } } $c }
sub inv { my $a = shift; my $b = z0(); $b->[0] = 1/$a->[0]; for my $n (1..$N) { my $s = R(0); $s += $a->[$_]*$b->[$n-$_] for 1..$n; $b->[$n] = -$s/$a->[0]; } $b }
sub sqrt1 { my $a = shift; die unless $a->[0] == 1; my $b = z0(); $b->[0] = R(1); for my $n (1..$N) { my $s = R(0); $s += $b->[$_]*$b->[$n-$_] for 1..$n-1; $b->[$n] = ($a->[$n]-$s)/2; } $b }
sub expo { my $a = shift; die unless $a->[0]->is_zero; my $e = z0(); $e->[0] = R(1); for my $n (1..$N) { my $s = R(0); $s += $a->[$_]*$e->[$n-$_]*$_ for 1..$n; $e->[$n] = $s/$n; } $e }
sub der { my $a = shift; my $b = z0(); $b->[$_-1] = $a->[$_]*$_ for 1..$N; $b }
sub integ { my $a = shift; my $b = z0(); $b->[$_] = $a->[$_-1]/$_ for 1..$N; $b }
sub up { my $a = shift; my $b = z0(); $b->[$_] = $a->[$_-1] for 1..$N; $b }        # times t
sub down { my $a = shift; die "down: nonzero const" unless $a->[0]->is_zero; my $b = z0(); $b->[$_-1] = $a->[$_] for 1..$N; $b }   # divided by t (top lost)
sub comp { my ($f, $g) = @_; my $r = z0(); my $p = c0(1); for my $k (0..$N) { $r = add($r, sc($f->[$k], $p)); $p = mul($p, $g); } $r }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }

for my $x (R(2), R("3/5"), R(-3), R("5/2")) { my $c = 1 + $x;
  my $tt = z0(); $tt->[$_] = R((($_ % 2) ? -1 : 1))/fact($_ + 1) for 0..$N;    # tau/t
  my $tau = up($tt);
  my $dc = sqrt1(sub_(c0(1), sc(4*$x/$c**2, $tau)));                         # delta/c
  my $cs = sc(R(2), inv(add(c0(1), $dc)));                                   # c s*
  my $ss = sc(1/$c, $cs);
  my $dphi = integ(sub_(c0(1), sc($x, $ss)));                                # d(dphi)/dt = 1 - x s*  (S2.2(e))
  my $alg = mul(sqrt1($dc), inv($cs));
  my $rho = mul($tt, inv($alg));
  my $tL1 = sub_(sc($c, inv($tt)), sc($x, up($ss)));                          # t L_1
  my $Theta = sc(-$x**3, mul($tau, mul($ss, $ss)));
  # G(z) = c/(c+z) (tau/t) F(z, tau) Exp[-dphi/z] at a rational z
  my $G = sub { my $z = shift; my @F; for my $n (0..$N) { my $num = R(1); my $den = fact($n)*$z**$n;
        $num *= (1 + $_*$z)**2 for 1..$n; $den *= ($c + $_*$z) for 2..$n+1; push @F, $num/$den; }
      sc($c/($c + $z), mul($tt, mul(comp(\@F, $tau), expo(sc(-1/$z, $dphi))))) };
  my $lam = sub { my ($j, $p) = @_; 1/($c - $j*$p) };
  # tables of leaf_structure.md
  my @delta = (
    sub { my $p = shift; $c - $c**2*$lam->(1,$p) },
    sub { my $p = shift; ((1-$x)**2/2)*($lam->(2,$p) - $lam->(1,$p)) },
    sub { my $p = shift; (-(2 - 4*$x + 3*$x**2 - 4*$x**3 + 2*$x**4)*$lam->(1,$p) + 6*(1-$x)**4*$lam->(2,$p) - (2-$x)**2*(1-2*$x)**2*$lam->(3,$p))/(12*$c**2) },
    sub { my $p = shift; (1-$x)**2*(-(3 + 5*$x**2 + 3*$x**4)/72*$lam->(1,$p) + (7 - 20*$x + 22*$x**2 - 20*$x**3 + 7*$x**4)/24*$lam->(2,$p)
          - (2-$x)**2*(1-2*$x)**2/8*$lam->(3,$p) + (3-$x)**2*(1-3*$x)**2/36*$lam->(4,$p))/$c**4 } );
  my %e; $e{0} = { "11" => -$x**3/2 };
  $e{1} = { "11" => -(1-$x)**2*$x**3/(2*$c**2), "12" => 2*(1-$x)**2*$x**3/(3*$c**2), "21" => 2*(1-$x)**2*$x**3/(3*$c**2) };
  $e{2} = { "11" => -$x**3*(7 - 20*$x + 21*$x**2 - 20*$x**3 + 7*$x**4)/(24*$c**4), "12" => (1-$x)**4*$x**3/$c**4, "21" => (1-$x)**4*$x**3/$c**4,
            "13" => -3*(2-$x)**2*(1-2*$x)**2*$x**3/(16*$c**4), "31" => -3*(2-$x)**2*(1-2*$x)**2*$x**3/(16*$c**4), "22" => -(1-$x)**4*$x**3/$c**4 };
  my $loopT = sub { my ($k, $p, $q) = @_; my $s = R(0); for my $il (keys %{$e{$k}}) { my ($i, $l) = split //, $il; $s += $e{$k}{$il}*$lam->($i,$p)*$lam->($l,$q); } $s };
  my $pair = sub { my ($a, $b) = @_; add(sc($x**2, add(mul($a->[0],$b->[1]), mul($a->[1],$b->[0]))), sc(-$c*$x, mul($a->[1],$b->[1]))) };
  my @psis = grep { my $r = $c/$_; !($r->is_int && abs($r) <= $N + 4) } (R("2/7"), R("-5/3"), R("4/9"), R("3/7"), R("-1/5"));
  my %g; my %Y;
  my $mkY = sub { my $p = shift; my $gg = $p->is_zero ? $rho : $G->(-$p);   # g(psi) = G(-psi); G(0; t) = rho r_{0,1}[0] = rho
    # Y := t sigma(psi) = (t g, (t L_1 - psi) g - psi t g')   (I1)
    ($gg, [ up($gg), sub_(mul(sub_($tL1, c0($p)), $gg), sc($p, up(der($gg)))) ]) };
  for my $p (@psis, R(0)) { ($g{$p}, $Y{$p}) = $mkY->($p); }
  for my $p (@psis) { my $q = sc(-$p, $g{$p});
    for my $k (0..3) { ok($q->[$k] == $delta[$k]->($p), "dilaton t^$k (x=$x, psi=$p)"); } }
  for my $p (@psis) { my (undef, $Ym) = $mkY->(-$p);
    my $K = $pair->($Y{$p}, $Ym);
    ok($K->[$_] == ($_ == 0 ? -$x*$c**3 : 0), "K_1 = 1: t^2 <sigma(psi), sigma(-psi)> = -x c^3 at t^$_ (x=$x, psi=$p)") for 0..$N-2; }
  # loops
  for my $p (@psis) { for my $q (@psis) { next if ($p + $q)->is_zero;
      my $num = sub_(c0(-$x*$c**3), $pair->($Y{$p}, $Y{$q}));               # t^2 rho^2 E_11 (psi + psi')
      my $E = sc(1/($p + $q), $num);
      ok($E->[0]->is_zero && $E->[1]->is_zero, "t^2 rho^2 E_11 = O(t^2) (x=$x, $p, $q)");
      for my $k (0..2) { ok($E->[$k+2] == $loopT->($k, $p, $q), "loop t^$k via (I4) (x=$x, psi=$p, psi'=$q)"); }
      # flow (V8): (k+2) E^(k) = [t^{k+1}] Theta (t g(psi))' (t g(psi'))'
      my $fl = mul($Theta, mul(der(up($g{$p})), der(up($g{$q}))));
      for my $k (0..$N-3) { ok($fl->[$k+1] == ($k+2)*$E->[$k+2], "flow (V8) at t^$k (x=$x, $p, $q)"); } } }
  # E11 leg b = 0 at t^1: rho f_{1,0} = rho^-1 [psi'^0] rho^2 E_11 = rho^-1 E(psi, 0)/t^2
  for my $p (@psis) {
    my $E = sc(1/$p, sub_(c0(-$x*$c**3), $pair->($Y{$p}, $Y{R(0)})));
    my $leg = mul(inv($rho), down(down($E)));
    ok($leg->[0] == -($x**3/2)*$lam->(1,$p)/$c, "E11 leg b=0 t^0 (x=$x)");
    ok($leg->[1] == (1-$x)**2*$x**3/$c**3*(-$lam->(1,$p)/12 + 2*$lam->(2,$p)/3), "E11 leg b=0 t^1 (x=$x, psi=$p)"); }
  print "x=$x: $ntest checks so far, $nfail failures\n";
}
print "ref_S2_leaves.pl N=$N: $ntest checks, $nfail failures\n";
