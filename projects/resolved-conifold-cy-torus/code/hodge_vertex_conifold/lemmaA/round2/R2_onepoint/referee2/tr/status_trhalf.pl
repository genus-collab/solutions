# status_trhalf.pl: copy of referee1_trhalf.pl with the odd-j range of the n = 1 test set by env JMAX.
# Exact EO recursion on C_{-1/2}: zeta = v + 1/2, x^2 = 1/(zeta^2-1/4), y = 1/2 - zeta, K = ((zeta^2-1/4)/2) dz0/(zeta (z0^2-zeta^2) dzeta).
use strict; use warnings; use Math::BigRat;
sub Qr { Math::BigRat->new(@_) }
my ($GN,$GMAX,$SMAX) = (shift//6, shift//3, shift//99);
my @p = (Qr("-1/8"), Qr("1/2"));
my %W;
sub add { my ($h,$k,$v)=@_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }
$W{"0,3"} = { "0,0,0" => Qr(-2)*$p[0] };
$W{"1,1"} = { "1" => Qr("-1/4")*$p[0], "0" => Qr("-1/4")*$p[1] };
sub stable { 2*$_[0]-2+$_[1] > 0 }
sub compute { my ($g,$np1)=@_; my $n=$np1-1; my %res; my %D;
  if ($g>=1 && stable($g-1,$n+2)) { my $Wa=$W{($g-1).",".($n+2)}; for my $key (keys %$Wa){ add(\%D,$key,$Wa->{$key}); } }
  for my $h (0..$g){ for my $mask (0..(2**$n)-1){ my @I=grep { ($mask>>$_)&1 } 0..$n-1; my @Ic=grep { !(($mask>>$_)&1) } 0..$n-1;
      next unless stable($h,1+@I) && stable($g-$h,1+@Ic);
      my $W1=$W{"$h,".(1+@I)}; my $W2=$W{($g-$h).",".(1+@Ic)};
      for my $k1 (keys %$W1){ my @a=split /,/,$k1; my $a0=shift @a; for my $k2 (keys %$W2){ my @b=split /,/,$k2; my $b0=shift @b; my @kJ; my ($ia,$ib)=(0,0);
          for my $j (0..$n-1){ push @kJ, (($mask>>$j)&1) ? $a[$ia++] : $b[$ib++]; } add(\%D, join(",",$a0,$b0,@kJ), $W1->{$k1}*$W2->{$k2}); } } } }
  for my $key (keys %D){ my @k=split /,/,$key; my $a=shift @k; my $b=shift @k;
    for my $i (0..$#p){ my $j=$a+$b+2-$i; next if $j<0; add(\%res, join(",",$j,@k), -$p[$i]*$D{$key}); } }
  if ($n>=1 && stable($g,$n)) { my $Wc=$W{"$g,$n"};
    for my $q (0..$n-1){ for my $key (keys %$Wc){ my @k=split /,/,$key; my $b=shift @k;
        for my $j (0..$b+1){ for my $pp (0..$b+1-$j){ my $i=$b+1-$j-$pp; next if $i>$#p; my @kJ=@k; splice(@kJ,$q,0,$pp);
            add(\%res, join(",",$j,@kJ), Qr(-2*(2*$pp+1))*$p[$i]*$Wc->{$key}); } } } } }
  $W{"$g,$np1"}=\%res; }
for my $chi (2..$GN){ for my $g (0..$GMAX){ my $n=$chi-2*$g+2; next if $n<1 || $g+$n>$SMAX; next if ($g==0 && $n<=3) || ($g==1 && $n==1); compute($g,$n); } }
sub fact { my $n=shift; my $r=Qr(1); $r*=$_ for 1..$n; $r }
sub gbinom { my ($a,$k)=@_; my $r=Qr(1); for my $i (0..$k-1){ $r *= ($a-$i); } $r/fact($k) }
sub leaf { my ($k,$mu)=@_; my $d=$mu-2*$k-1; return Qr(0) if $d<0 || $d%2; my $m=$d/2; my $a=Qr(0)-Qr($k)-Qr("1/2"); my $r=Qr(-1)/Qr(2*$k+1); $r=$r*gbinom($a,$m); $r=$r/(Qr(4)**$m); return $r; }
sub Xi { my ($gn,@mu)=@_; my $h=$W{$gn}; my $tot=Qr(0); for my $key (keys %$h){ my @k=split /,/,$key; my $t=$h->{$key}->copy; $t*=leaf($k[$_],$mu[$_]) for 0..$#mu; $tot+=$t; }
  for my $m (@mu){ $tot /= gbinom(Qr($m)/Qr(2)-Qr(1),$m-1); } $tot }
use File::Basename (); use File::Spec ();
my $dir = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../../../..";
sub getexpr { my ($file,$tag)=@_; open my $f,'<',"$dir/$file" or die; while(<$f>){ if(index($_,$tag)>=0){ my $e=substr($_,index($_,$tag)+length($tag)); chomp $e; $e=~s/\r//g; return $e; } } die }
sub ev { my ($e,$x,@j)=@_; my $s=$e; $s=~s/jj\[(\d)\]/\$J[$1]/g; $s=~s/\bx\b/\$X/g; $s=~s/\^/**/g; my @J=(undef,map { Qr($_) } @j); my $X=Qr($x); my $v=eval $s; die $@ if $@; $v }
my %F = ("2,2"=>getexpr("hvc_multiprobe_out.txt","genus 2 symbolic r=2: "), "3,2"=>getexpr("hvc_multiprobe_out.txt","genus 3 symbolic r=2: "),
         "4,1"=>getexpr("hvc_multiprobe_g4_out.txt","genus 4 r=1 symbolic: "));
for my $gn (sort keys %W){ my ($g,$n)=split /,/,$gn; next unless $g>=2 && $n<=2; my @out;
  my @tups = $n==1 ? (map {[$_]} grep { $_%2 } 1..($ENV{JMAX}//(2*$g+5))) : (map { my $a=$_; map { [$a,$_] } grep { $_%2 && $_>=$a } 1..2*$g+3 } grep { $_%2 } 1..2*$g+3);
  for my $t (@tups){ my $s=0; $s+=$_ for @$t; my $v=Xi($gn,@$t); my $str="(".join(",",@$t).")=".$v->bstr;
    if ($s<=2*$g-3+$n && !$v->is_zero){ $str.=" <<VIOLATION>>"; }
    if (exists $F{$gn}) { my $m=ev($F{$gn},1,@$t); $str .= " [M=".$m->bstr.(!$m->is_zero ? ", ratio ".($v/$m)->bstr : "")."]"; }
    push @out,$str; }
  print "tau=-1/2 g=$g n=$n: ", join("  ",@out), "\n"; }
