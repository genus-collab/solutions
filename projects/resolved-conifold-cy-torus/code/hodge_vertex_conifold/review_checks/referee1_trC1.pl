# Independent exact EO recursion on C_1 in the odd coordinate zeta = (v-1)/(v+1):
#   x = (1-zeta^2)/4, y = -(1+zeta)/(1-zeta), B = dz dz'/(z-z')^2, K = (1-z^2)^2 dz0 / (8 z (z0^2-z^2) dz).
# omega_{g,n} = Sum C[k1..kn] prod e_{k_i}(z_i), e_k(z) = dz/z^(2k+2).
use strict; use warnings; use Math::BigRat;
my $GN = shift // 6; my $GMAX = shift // 3; my $SMAX = shift // 99;
sub R { Math::BigRat->new(@_) }
my @c = (R(1), R(-2), R(1));
my %W;  # W{"g,n"}{"k1,...,kn"} = coeff
sub add { my ($h,$k,$v)=@_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }
$W{"0,3"} = { "0,0,0" => R("-1/4") };
$W{"1,1"} = { "1" => R("-1/32"), "0" => R("1/16") };
sub stable { my ($g,$n)=@_; return 2*$g-2+$n > 0 }
sub compute { my ($g,$np1)=@_; my $n=$np1-1; my %res; my %D;
  # term A
  if ($g>=1 && stable($g-1,$n+2)) { my $Wa=$W{($g-1).",".($n+2)};
    for my $key (keys %$Wa){ my @k=split /,/,$key; my $a=shift @k; my $b=shift @k; add(\%D, join(",",$a,$b,@k), $Wa->{$key}); } }
  # term B: stable products
  for my $h (0..$g){ for my $mask (0..(2**$n)-1){ my @I=grep { ($mask>>$_)&1 } 0..$n-1; my @Ic=grep { !(($mask>>$_)&1) } 0..$n-1;
      next unless stable($h,1+@I) && stable($g-$h,1+@Ic);
      my $W1=$W{"$h,".(1+@I)} or die "missing $h,".(1+@I); my $W2=$W{($g-$h).",".(1+@Ic)} or die "missing";
      for my $k1 (keys %$W1){ my @a=split /,/,$k1; my $a0=shift @a;
        for my $k2 (keys %$W2){ my @b=split /,/,$k2; my $b0=shift @b; my @kJ;
          my ($ia,$ib)=(0,0); for my $j (0..$n-1){ if (($mask>>$j)&1) { push @kJ,$a[$ia++]; } else { push @kJ,$b[$ib++]; } }
          add(\%D, join(",",$a0,$b0,@kJ), $W1->{$k1}*$W2->{$k2}); } } } }
  for my $key (keys %D){ my @k=split /,/,$key; my $a=shift @k; my $b=shift @k;
    for my $i (0..2){ my $j=$a+$b+2-$i; next if $j<0; add(\%res, join(",",$j,@k), R("-1/8")*$c[$i]*$D{$key}); } }
  # term C: omega_{0,2}(z,z_q) x omega_{g,n}
  if ($n>=1 && stable($g,$n)) { my $Wc=$W{"$g,$n"} or die "missing $g,$n";
    for my $q (0..$n-1){ for my $key (keys %$Wc){ my @k=split /,/,$key; my $b=shift @k;   # k = indices of J\q in order
        for my $j (0..$b+1){ for my $p (0..$b+1-$j){ my $i=$b+1-$j-$p; next if $i>2;
            my @kJ=@k; splice(@kJ,$q,0,$p);
            add(\%res, join(",",$j,@kJ), R(-(2*$p+1))/4*$c[$i]*$Wc->{$key}); } } } } }
  $W{"$g,$np1"}=\%res; }
my @order; for my $chi (2..$GN){ for my $g (0..$GMAX){ my $n=$chi-2*$g+2; next if $n<1; next if $g+$n > $SMAX; next if ($g==0 && $n<=3) || ($g==1 && $n==1); push @order,[$g,$n]; } }
for my $o (@order){ compute(@$o); }
# symmetry check
for my $gn (sort keys %W){ my $h=$W{$gn}; my $bad=0; for my $key (keys %$h){ my $s=join(",",sort { $a <=> $b } split /,/,$key);
    my $acc=R(0); } ; }
sub symcheck { my $h=shift; my %s; for my $key (keys %$h){ my @k=split /,/,$key; my $sk=join(",",sort { $a<=>$b } @k); push @{$s{$sk}}, $h->{$key}; }
  my $bad=0; for my $sk (keys %s){ my @v=@{$s{$sk}}; for my $v (@v){ $bad++ if $v != $v[0]; } } $bad }
# x-expansion: [x^mu] F = Sum_k C[k] prod (-2/mu) 4^(mu-1) (k+3/2)_(mu-1)/(mu-1)!
sub poch { my ($a,$m)=@_; my $r=R(1); for my $i (0..$m-1){ $r *= $a+$i; } $r }
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binom { my ($n,$k)=@_; return R(0) if $k<0 || $k>$n; fact($n)/(fact($k)*fact($n-$k)) }
sub leafcoef { my ($k,$mu)=@_; R(-2)/$mu * R(4)**($mu-1) * poch(R($k)+R("3/2"),$mu-1)/fact($mu-1) }
sub xcoef { my ($gn,@mu)=@_; my $h=$W{$gn}; my $tot=R(0);
  for my $key (keys %$h){ my @k=split /,/,$key; my $p=$h->{$key}->copy; for my $i (0..$#mu){ $p *= leafcoef($k[$i],$mu[$i]); } $tot += $p; } $tot }
sub Xi { my ($gn,@mu)=@_; my $v=xcoef($gn,@mu); for my $m (@mu){ $v /= binom(2*$m-1,$m-1); } $v }
binmode STDOUT;
for my $gn (sort keys %W){ my ($g,$n)=split /,/,$gn; next unless $g>=2 || ($g==1);
  my $nt=scalar keys %{$W{$gn}}; my $bad=symcheck($W{$gn});
  print "($g,$n): ", $nt, " terms, symmetry violations: $bad\n"; }
# degree test: all mu with mu_i>=1 and |mu| <= 2g-2+n+1 (one beyond the predicted first nonzero layer)
for my $gn (sort keys %W){ my ($g,$n)=split /,/,$gn; next unless $g>=2;
  my @out; my $lim=2*$g-2+$n+($n==1?4:1);
  my @stack=([]); my @tuples;
  my $rec; $rec = sub { my ($pre,$left,$min)=@_; if (@$pre==$n){ push @tuples,[@$pre]; return; } for my $m ($min..$left-($n-@$pre-1)){ $rec->([@$pre,$m],$left-$m,$m); } };
  $rec->([], $lim, 1);
  for my $t (@tuples){ my $s=0; $s+=$_ for @$t; my $v=Xi($gn,@$t); my $pred = ($s <= 2*$g-3+$n) ? "zero" : "nonzero";
    my $ok = ($pred eq "zero") ? $v->is_zero : !$v->is_zero;
    push @out, "(".join(",",@$t).")=".$v->bstr.($ok?"":" <<VIOLATION>>"); }
  print "g=$g n=$n: ", join("  ",@out), "\n"; }
