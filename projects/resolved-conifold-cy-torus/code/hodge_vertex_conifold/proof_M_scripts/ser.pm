package ser;
use strict; use warnings;
use Math::BigRat;
our $N = 12;
sub R { Math::BigRat->new(@_) }
sub zero { [ map { R(0) } 0..$N ] }
sub mul { my ($a,$b)=@_; my $c=zero(); for my $i (0..$N){ next if $a->[$i]->is_zero; for my $j (0..$N-$i){ $c->[$i+$j] += $a->[$i]*$b->[$j]; } } $c }
sub add { my ($a,$b)=@_; [ map { $a->[$_]+$b->[$_] } 0..$N ] }
sub sub_ { my ($a,$b)=@_; [ map { $a->[$_]-$b->[$_] } 0..$N ] }
sub scal { my ($c,$a)=@_; [ map { $c*$a->[$_] } 0..$N ] }
sub inv { my ($a)=@_; die "inv: a0=0" if $a->[0]->is_zero; my $b=zero(); $b->[0]=R(1)/$a->[0];
  for my $n (1..$N){ my $s=R(0); for my $k (1..$n){ $s += $a->[$k]*$b->[$n-$k]; } $b->[$n] = -$s/$a->[0]; } $b }
sub shiftdown { my ($a,$v)=@_; my $c=zero(); for my $i ($v..$N){ $c->[$i-$v]=$a->[$i]; } $c }  # divide by t^v (assumes low coeffs zero)
sub shiftup { my ($a,$v)=@_; my $c=zero(); for my $i (0..$N-$v){ $c->[$i+$v]=$a->[$i]; } $c }
sub deriv { my ($a)=@_; my $c=zero(); for my $i (1..$N){ $c->[$i-1]=$a->[$i]*$i; } $c }
sub sqrt1 { my ($a)=@_; die "sqrt: a0 != 1" unless $a->[0] == 1; my $b=zero(); $b->[0]=R(1);
  for my $n (1..$N){ my $s=R(0); for my $k (1..$n-1){ $s += $b->[$k]*$b->[$n-$k]; } $b->[$n] = ($a->[$n]-$s)/2; } $b }
sub compose { my ($f,$g)=@_; # f(g(t)), g0 = 0
  die "compose: g0 != 0" unless $g->[0]->is_zero; my $r=zero(); my $p=zero(); $p->[0]=R(1);
  for my $k (0..$N){ $r = add($r, scal($f->[$k],$p)); $p = mul($p,$g); } $r }
sub revert { my ($w)=@_; # w = t + ..., return t(w)
  my $t=zero(); $t->[1]=R(1);
  for my $it (1..$N+1){ my $wt = compose($w,$t); my $err = sub_($wt, do{my $x=zero(); $x->[1]=R(1); $x}); $t = sub_($t,$err); } $t }
sub coef { my ($a,$i)=@_; $a->[$i] }
sub show { my ($a,$m)=@_; $m//=$N; join(", ", map { $a->[$_]->bstr } 0..$m) }
1;
