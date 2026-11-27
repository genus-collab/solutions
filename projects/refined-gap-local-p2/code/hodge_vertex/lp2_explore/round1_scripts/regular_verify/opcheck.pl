# check the claimed closed forms of Op1, Op2, Op3 on s^n: compute q^3 Op_j(s^n) directly from
# q^3 Op = q^3 (1/s + z th)^3 - (N^3 + z N^2 + 2/9 z^2 N) and compare with q^3 * (closed form).  usage: perl opcheck.pl E1 E2 NMAX
use strict; use warnings; use Math::BigRat; sub R { Math::BigRat->new(@_) }
my ($E1,$E2,$NM)=@ARGV; my $al=(R($E1)+R($E2))/3;
sub lp_add { my %r=%{$_[0]}; for my $n (keys %{$_[1]}){ $r{$n}=($r{$n}//R(0))+$_[1]{$n} } for(keys %r){delete $r{$_} if $r{$_}->is_zero} \%r }
sub lp_sc { my($p,$c)=@_; my %r; for(keys %$p){ my $v=$p->{$_}*$c; $r{$_}=$v unless $v->is_zero } \%r }
sub lp_mul { my($p,$q)=@_; my %r; for my $i (keys %$p){ for my $j (keys %$q){ $r{$i+$j}=($r{$i+$j}//R(0))+$p->{$i}*$q->{$j} } } for(keys %r){delete $r{$_} if $r{$_}->is_zero} \%r }
sub lp_der { my($p)=@_; my %r; for(keys %$p){ next if $_==0; $r{$_-1}=$p->{$_}*$_ } \%r }
my $Lp={-1=>R(1),0=>-$al}; my $qp={0=>R(1),1=>-$al}; my $q3=lp_mul($qp,lp_mul($qp,$qp));
sub th { lp_sc(lp_mul($qp,lp_der($_[0])), R(1)/(3*$al)) }
my $K=3;
sub sadd { my($x,$y,$c)=@_; [map{ lp_add($x->[$_], lp_sc($y->[$_],$c)) }0..$K] }
sub zth { my($s)=@_; [ {}, map{ th($s->[$_]) }0..$K-1 ] }
sub smulp { my($x,$p)=@_; [map{ lp_mul($x->[$_],$p) }0..$K] }
sub zsh { my($x,$m)=@_; [ (map{{}}1..$m), @{$x}[0..$K-$m] ] }
sub Aop { sadd(smulp($_[0],{-1=>R(1)}), zth($_[0]), R(1)) } sub Nop { sadd(smulp($_[0],$Lp), zth($_[0]), R(1)) }
my $all=1;
for my $n (1..$NM){ my $S=[{$n=>R(1)},{},{},{}];
  my $A3=Aop(Aop(Aop($S))); my $N1=Nop($S); my $N2=Nop($N1); my $N3=Nop($N2);
  my $r=smulp($A3,$q3); $r=sadd($r,$N3,R(-1)); $r=sadd($r,zsh($N2,1),R(-1)); $r=sadd($r,zsh($N1,2),R(-2)/9);
  my $c0=$r->[0];
  my $op1={ $n-2 => R(-$n)+1 };             # Op1(s^n) = -n s^(n-2) + s^(n-2)
  my $op2={ $n-2 => (-1)*R(($n-1)*(2*$n-1))/(3*$al), $n-1 => R(3*$n*$n-3*$n+1)/9 };
  my $op3={ $n-2 => (-1)*R($n*($n-1)**2)/(9*$al**2), $n-1 => R($n*(3*$n*$n-3*$n+1))/(27*$al), $n => (-1)*R($n**3)/27 };
  my @ok=(!%$c0); my $j=1; for my $op ($op1,$op2,$op3){ my $d=lp_add($r->[$j], lp_sc(lp_mul($q3,$op),R(-1))); push @ok, !%$d; $j++ }
  $all&&=!grep{!$_}@ok; print "n=$n: z^0 vanishes: ".($ok[0]?1:0)."; Op1,Op2,Op3 closed forms: @ok[1..3]\n" }
print "ALL: ".($all?"OK":"FAIL")."\n";
