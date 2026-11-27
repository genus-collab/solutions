package XY;
# Closed formulas for the TR correlators of C_tau from the x-y swap theorem for log topological recursion:
#   A. Alexandrov, B. Bychkov, P. Dunin-Barkowski, M. Kazarian, S. Shadrin, "Log topological recursion through the prism
#   of x-y swap", IMRN 2024 (arXiv:2312.16950), Theorem 2.8 and Proposition 3.1 (Sigma = CP^1, y an affine coordinate).
# Their convention: omega^(0)_1 = -y dx.  Our TR (P2 normalisation): omega_{0,1} = y_ours dX with y_ours = -v, so in their
# convention  x = X = tau log v - (tau+1) log(1+v),  y = v  (z = v the affine coordinate).
# Log singularities of x at a_0 = 0 (residue 1/alpha_0 = tau) and a_1 = -1 (residue 1/alpha_1 = -(tau+1)) are the
# LogTR-vital points of the dual curve; xhat = Sum_i (1/S(alpha_i hbar d_y)) log(v - a_i)/alpha_i,  S(t) = 2 sinh(t/2)/t.
# n = 1:  omega_{g,1} = -[hbar^{2g}] Sum_{r>=0} (d (1/dX))^r [u^{r+1}] E dv,   E = exp(-u (S(u hbar d_y) xhat - X)).
# All series are expanded at v = oo in t = 1/v (the expansion point of P2).
use strict; use warnings; use Math::BigRat;
use Exporter 'import';
our @EXPORT_OK = qw(R xi_formula_n1 Elist);
sub R { Math::BigRat->new(@_) }
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $k (0..$i-1){ $r *= ($a-$k); } $r/fact($i) }

# returns hashref: E->[m]{udeg} = t-series (arrayref, length T+1) of [hbar^(2m) u^udeg] E
sub Elist { my ($tau,$G,$T)=@_;
  my $zero = sub { [ map { R(0) } 0..$T ] };
  my $smul = sub { my ($a,$b)=@_; my @c=map { R(0) } 0..$T; for my $i (0..$T){ next if $a->[$i]->is_zero; for my $j (0..$T-$i){ next if $b->[$j]->is_zero; $c[$i+$j] += $a->[$i]*$b->[$j]; } } \@c };
  my $sadd = sub { my ($a,$b,$s)=@_; $s //= R(1); [ map { $a->[$_] + $s*$b->[$_] } 0..$T ] };
  my $K = $G + 1; my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$K;
  my @sp = (R(1)); for my $n (1..$K){ my $acc=R(0); for my $k (1..$n){ $acc += $s[$k]*$sp[$n-$k]; } push @sp, -$acc; }
  my @a = (R(0), R(-1)); my @al = (R(1)/$tau, R(-1)/($tau+1));
  my $vpow_inv = sub { my ($aa,$m)=@_; my $r=$zero->(); for my $i (0..$T-$m){ $r->[$m+$i] = binomr(R(-$m),$i)*((-$aa)**$i); } $r };
  my @L; for my $k (1..$G){ my %P;
    for my $i (0..1){ my $vi = $vpow_inv->($a[$i],2*$k); my $f = fact(2*$k-1)/$al[$i];
      for my $p (0..$k){ my $q=$k-$p; my $cf = $f*$s[$p]*$sp[$q]*($al[$i]**(2*$q)); my $deg = 2*$p+1;
        $P{$deg} = exists $P{$deg} ? $sadd->($P{$deg},$vi,$cf) : [ map { $cf*$_ } @$vi ]; } }
    $L[$k]=\%P; }
  my $pmul = sub { my ($A,$B)=@_; my %C; for my $i (keys %$A){ for my $j (keys %$B){ my $t=$smul->($A->{$i},$B->{$j}); $C{$i+$j} = exists $C{$i+$j} ? $sadd->($C{$i+$j},$t) : $t; } } \%C };
  my $one=$zero->(); $one->[0]=R(1);
  my @E = ({ 0 => $one });
  for my $m (1..$G){ my %acc; for my $k (1..$m){ my $t = $pmul->($L[$k],$E[$m-$k]); for my $d (keys %$t){ $acc{$d} = exists $acc{$d} ? $sadd->($acc{$d},$t->{$d},R($k)/$m) : [ map { (R($k)/$m)*$_ } @{$t->{$d}} ]; } } $E[$m]=\%acc; }
  return \@E; }

# xi_formula_n1(tau, G, M): returns hash "g,mu" => Xi_g(mu) (P2 normalisation) for 1 <= g <= G, 1 <= mu <= M
sub xi_formula_n1 { my ($tau,$G,$M)=@_; $tau = R($tau) unless ref $tau; my $sig=$tau*($tau+1); my $T=$M+4;
  my $E = Elist($tau,$G,$T);
  my $zero = sub { [ map { R(0) } 0..$T ] };
  my $smul = sub { my ($a,$b)=@_; my @c=map { R(0) } 0..$T; for my $i (0..$T){ next if $a->[$i]->is_zero; for my $j (0..$T-$i){ next if $b->[$j]->is_zero; $c[$i+$j] += $a->[$i]*$b->[$j]; } } \@c };
  my $pre = do { my $n=$zero->(); $n->[0]=R(1); $n->[1]=R(1); my $di=$zero->(); for my $i (0..$T){ $di->[$i]=$tau**$i; } $smul->($n,$di) }; # (1+t)/(1-tau t)
  my $Dop = sub { my ($f)=@_; my $g=$zero->(); for my $i (1..$T){ $g->[$i] = $i*$f->[$i]; } $smul->($g,$pre) };     # t f'(t) (1+t)/(1-tau t)
  my %out;
  for my $g (1..$G){ my $Eg = $E->[$g]; my $F = $zero->();
    if (exists $Eg->{1}) { my $e=$Eg->{1}; for my $i (0..$T){ next if $e->[$i]->is_zero; my $pw=$i-2; die "log term" if $pw==-1; $F->[$pw+1] += -$e->[$i]/($pw+1) if $pw+1<=$T; } }
    for my $r (1..3*$g){ my $e=$Eg->{$r+1} or next; my $h=$smul->($e,$pre); die "pole" unless $h->[0]->is_zero; my $f=$zero->(); for my $i (1..$T){ $f->[$i-1] = -$h->[$i]; }
      for (1..$r-1){ $f = $Dop->($f); } $F->[$_] += $f->[$_] for 0..$T; }
    # sign: their y = v  ->  overall factor -1 relative to the y = -v evaluation
    $F->[$_] = -$F->[$_] for 0..$T;
    for my $mu (1..$M){ my $c=R(0); for my $k (1..$mu){ next if $F->[$k]->is_zero; $c += $k*$F->[$k]*binomr(($tau+1)*$mu,$mu-$k); } $c /= $mu;
      $out{"$g,$mu"} = $c*($sig**($g-1))/binomr($mu*($tau+1)-1,$mu-1); } }
  return \%out; }
1;
