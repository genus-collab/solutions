# trace_qle.pl: check of the global quadratic loop equation in trace form for C_tau (n = 0 case):
#   Tr_x Phi_g = 0,  Phi_g(v) = Brk_g(v,v)/(2 dX^2) - v omega_{g,1}(v)/dX,  Brk_g(v,v') = omega_{g-1,2}(v,v') + Sum_h omega_{h,1}(v) omega_{g-h,1}(v').
# Tr Phi_g = 0 is equivalent to: Phi_g has poles only at v = tau, vanishes at v = -1 (automatic), and the sigma-even part of its
# Laurent expansion in the local coordinate z at v = tau is holomorphic.  We test the last condition (mod p).
# usage: perl trace_qle.pl tn td G
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin; use TRmod;
my ($tn,$td,$G)=@ARGV; my $p=2147483629;
my $C = TRmod::new_curve($tn,$td,3*$G-1,$p); TRmod::run_tr($C,$G,2);
my $E=$C->{E}; my $sig=$C->{sig};
# Laurent series helpers on {lo, c}
sub lmul { my ($a,$b,$hi)=@_; my $lo=$a->{lo}+$b->{lo}; my $L=$hi-$lo; my @c=(0) x ($L+1); for my $i (0..$#{$a->{c}}){ last if $i>$L; my $ai=$a->{c}[$i] or next; for my $j (0..$#{$b->{c}}){ last if $i+$j>$L; my $bj=$b->{c}[$j] or next; $c[$i+$j]=($c[$i+$j]+$ai*$bj)%$p; } } { lo=>$lo, c=>\@c } }
sub lacc { my ($acc,$a,$f)=@_; for my $i (0..$#{$a->{c}}){ my $q=$a->{lo}+$i; $acc->{$q} = (($acc->{$q}//0) + $f*$a->{c}[$i])%$p; } }
my $HI = 2;
for my $g (2..$G){ my %phi;   # power -> coefficient of Phi(z)
  # omega_{g,1}(z)/dz = Sum C_k E_k(z)
  my %om; for my $key (keys %{$C->{W}{"$g,1"}}){ lacc(\%om, $E->[$key], $C->{W}{"$g,1"}{$key}); }
  # Brk(z,z)/dz^2
  my %br; for my $key (keys %{$C->{W}{($g-1).",2"}}){ my ($a,$b)=split /,/,$key; lacc(\%br, lmul($E->[$a],$E->[$b],$HI), $C->{W}{($g-1).",2"}{$key}); }
  for my $h (1..$g-1){ for my $k1 (keys %{$C->{W}{"$h,1"}}){ for my $k2 (keys %{$C->{W}{($g-$h).",1"}}){ lacc(\%br, lmul($E->[$k1],$E->[$k2],$HI), $C->{W}{"$h,1"}{$k1}*$C->{W}{($g-$h).",1"}{$k2}%$p); } } }
  # Phi = sig^2/(2 z^2) Brk + sig v(z)/z * om ,  v(z) = tau + u(z)
  my $i2 = TRmod::inv(2,$p); my $s2 = $sig*$sig%$p*$i2%$p;
  for my $q (keys %br){ $phi{$q-2} = (($phi{$q-2}//0) + $s2*$br{$q})%$p; }
  my @v = ($C->{tau}, @{$C->{u}}[1..40]);
  for my $q (keys %om){ for my $m (0..$#v){ my $pw=$q+$m-1; last if $pw > $HI; $phi{$pw} = (($phi{$pw}//0) + $sig*$v[$m]%$p*$om{$q})%$p; } }
  my @even_polar = grep { $_ < 0 && $_ % 2 == 0 && $phi{$_} } keys %phi; my @odd_polar = grep { $_ < 0 && $_ % 2 != 0 && $phi{$_} } keys %phi;
  print "tau=$tn/$td g=$g: Phi_g polar part: even powers non-zero: ".(@even_polar ? join(",",sort { $a<=>$b } @even_polar) : "none")."; odd polar powers present: ".scalar(@odd_polar)."\n"; }
