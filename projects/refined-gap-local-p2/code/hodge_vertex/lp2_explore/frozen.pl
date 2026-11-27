use strict; use warnings;
# Frame-B identity row of K_P2 x C^2 (T_F) at the singular idempotent, conifold limit L = LLL[0][y] -> infinity.
# Reads cPP_z<K>.m (list of identity-column coefficients, rational in eps1, eps2, L, sq), evaluates at large L,
# normalises by the z^0 entry, multiplies by h_w^(+-1) for w = (E, -e1, -e2), prints coefficient k times E^k.
my $f = shift; local $/; open my $fh,'<',$f or die; my $s=<$fh>; close $fh;
$s =~ s/\s+//g; $s =~ s/^\{//; $s =~ s/\}$//;
my @parts; my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
print "entries: ",scalar(@parts),"\n";
our ($e1,$e2)=(@ARGV>=2?@ARGV[0,1]:(2,1));
my $Lval = (@ARGV>=3?$ARGV[2]:1e7);
our $E=$e1+$e2; our $sq=sqrt($e1*$e2/$E**3);
my $bs = chr(92);
sub ev { my ($x,$L)=@_; my $t=$x; $t =~ s/\Q$bs\E//g;
  $t =~ s/\Q[Epsilon][1]\E/\$e1/g; $t =~ s/\Q[Epsilon][2]\E/\$e2/g; $t =~ s/\QLLL[0][y]\E/\$L/g;
  $t =~ s/\QSqrt[3]\E/sqrt(3)/g; $t =~ s/\bsq\b/\$sq/g; $t =~ s/\^/**/g; $t =~ tr/[]/()/;
  my $v = eval $t; die "$@ in $t" if $@; return $v;}
my @rho = map { ev($_,$Lval) } @parts; my @r = map { $_/$rho[0] } @rho;
printf "L*cPP0 = %.10g  (1/(sqrt3 sq E^2) = %.10g)\n", $rho[0]*$Lval, 1/(sqrt(3)*$sq*$E**2);
my @B=(1/6,-1/30,1/42,-1/30,5/66);
my @w=($E,-$e1,-$e2); my $N=scalar(@parts)-1;
my @logh=(0)x($N+1);
for my $k(1..5){ my $m=2*$k-1; next if $m>$N; my $S=0; $S+=$_**(-$m) for @w; $logh[$m]=-$B[$k-1]/(2*$k*(2*$k-1))*$S; }
sub expser { my @a=@_; my @e=(1,(0)x$N); for my $n(1..$N){ my $v=0; for my $k(1..$n){ $v+=$k*$a[$k]*$e[$n-$k]; } $e[$n]=$v/$n;} return @e;}
my @h=expser(@logh); my @hinv=expser(map{-$_}@logh);
sub mul { my ($a,$b)=@_; my @c=(0)x($N+1); for my $i(0..$N){for my $j(0..$N-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} return @c;}
for my $sg (1,-1){ my @rs = map { $r[$_]*($sg**$_) } 0..$N;
  for my $hp (1,-1){ my @hh = $hp==1?@h:@hinv;
   my @p=mul(\@rs,\@hh); printf "row sign %+d, h^%+d: %s\n",$sg,$hp, join(' ', map { sprintf "%.8g", $p[$_]*$E**$_ } 0..$N); } }
