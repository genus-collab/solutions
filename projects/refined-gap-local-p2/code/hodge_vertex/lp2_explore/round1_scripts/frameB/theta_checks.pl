use strict; use warnings; use Math::BigRat;
# theta_checks.pl -- exact checks of the normalisation identities used in the frame-B plan (z^0..z^N):
#  (1) Theta(z,-1)(1-z) = Theta(z)                                   [BS exmpl:stirling, used in thm:cfasym]
#  (2) Upsilon(z/e1,z/e2) * hinv_w(z) = 1, w=(E,-e1,-e2), hinv as in edges_sing.pl (frame B)   [Upsilon = h_w]
#  (3) exp(-c_0(z)) (BS Normalisation, 5 weights at P_0 of a generic torus) = prod_{j=0..2} Theta(-z/w_j) * N^-1-part,
#      i.e. exp(-c_0) / prod_{base} Theta(-z/w_j) = Theta(z/e1) Theta(z/e2) = sqrt(e1e2)/N_0  (BS eq:normNLR, z-part)
#  (4) leading conifold row: E*Upsilon(z/e1,z/e2)/(E+z) times hinv_w = E/(E+z)  (F2 from eq:Rcf, frame B)
my $N = $ARGV[0] // 9;
sub Q { Math::BigRat->new("$_[0]") }
my @B = (Q(1)); for my $n (1..$N+2) { my $s = Q(0); my $c = Q(1);
  for my $k (0..$n-1) { $s += binom($n+1,$k)*$B[$k]; } push @B, -$s/Q($n+1); }   # B_1 = -1/2
sub binom { my ($n,$k)=@_; my $r=Q(1); for my $i (1..$k){ $r = $r*Q($n-$i+1)/Q($i);} $r }
sub Bpoly { my ($n,$x)=@_; my $s=Q(0); for my $k (0..$n){ $s += binom($n,$k)*$B[$k]*($x**($n-$k)); } $s }
sub sexp { my @a=@_; my @e=(Q(1),(Q(0))x$N); for my $n (1..$N){ my $v=Q(0); for my $k (1..$n){ $v += Q($k)*$a[$k]*$e[$n-$k]; } $e[$n]=$v/Q($n);} @e }
sub smul { my ($a,$b)=@_; my @c=(Q(0))x($N+1); for my $i (0..$N){ for my $j (0..$N-$i){ $c[$i+$j]+=$a->[$i]*$b->[$j]; } } @c }
# log Theta(c z, x) = -Sum_k B_{k+1}(x) (-c z)^k/(k(k+1))
sub logTheta { my ($c,$x)=@_; my @l=(Q(0)); for my $k (1..$N){ push @l, -Bpoly($k+1,$x)*((-$c)**$k)/Q($k*($k+1)); } @l }
sub Theta { sexp(logTheta(@_)) }
sub iszero1 { my @a=@_; my $bad=0; $bad++ unless $a[0]==1; for my $i (1..$N){ $bad++ unless $a[$i]==0; } $bad }
# (1)
{ my @t1 = Theta(Q(1),Q(-1)); my @t0 = Theta(Q(1),Q(0)); my @omz=(Q(1),Q(-1),(Q(0))x($N-1));
  my @lhs = smul(\@t1,\@omz); my $bad=0; for my $i (0..$N){ $bad++ unless $lhs[$i]==$t0[$i]; }
  print "(1) Theta(z,-1)(1-z) = Theta(z) to z^$N: ", ($bad?"FAIL":"OK"), "\n"; }
for my $pr ([7,3],[5,2],[4,9],['1/3','1']) {
  my ($e1,$e2)=(Q($pr->[0]),Q($pr->[1])); my $E=$e1+$e2;
  # hinv_w exactly as edges_sing.pl
  my @Bb=(Q('1/6'),Q('-1/30'),Q('1/42'),Q('-1/30'),Q('5/66')); my @w=($E,-$e1,-$e2); my @lg=(Q(0))x($N+1);
  for my $k (1..5){ my $m=2*$k-1; next if $m>$N; my $S=Q(0); $S+=Q(1)/($_**$m) for @w; $lg[$m]=$Bb[$k-1]/Q(2*$k*(2*$k-1))*$S; }
  my @hinv = sexp(@lg);
  my @U = smul([smul([Theta(Q(1)/$e1,Q(0))],[Theta(Q(1)/$e2,Q(0))])],[Theta(-Q(1)/$E,Q(0))]);   # Theta(-z/E)=1/Theta(z/E)
  my @c2 = smul(\@U,\@hinv);
  print "(2) [e1,e2]=[$e1,$e2] Upsilon(z/e1,z/e2)*h_w^-1 = 1 to z^$N: ", (iszero1(@c2)?"FAIL":"OK"), "\n";
  # (3) generic torus at P_0 with weights alpha=(a0,a1,a2): tangent weights 3a0, a1-a0, a2-a0, -e1, -e2
  my ($a0,$a1,$a2)=(Q('2/7'),Q('5/3'),Q('-4/11'));
  my @W5=(3*$a0,$a1-$a0,$a2-$a0,-$e1,-$e2);
  my @mc=(Q(0)); for my $m (1..$N){ if($m%2){ my $n=($m+1)/2; my $S=Q(0); $S+=Q(1)/($_**$m) for @W5; push @mc, -$B[2*$n]/Q(2*$n*(2*$n-1))*$S; } else { push @mc, Q(0);} }
  my @emc = sexp(@mc);                                     # exp(-c_0(z))
  my @prodbase=(Q(1),(Q(0))x$N); for my $wj (@W5[0..2]) { @prodbase = smul(\@prodbase,[Theta(-Q(1)/$wj,Q(0))]); }
  my @rhs = smul(\@prodbase, [smul([Theta(Q(1)/$e1,Q(0))],[Theta(Q(1)/$e2,Q(0))])]);
  my $bad=0; for my $i (0..$N){ $bad++ unless $emc[$i]==$rhs[$i]; }
  print "(3) exp(-c_0) = prod_base Theta(-z/w) * Theta(z/e1)Theta(z/e2) to z^$N: ", ($bad?"FAIL":"OK"), "\n";
  # (4) E*Upsilon/(E+z) * hinv vs E/(E+z): equivalent to (2); print first coefficients of E/(E+z)
}
