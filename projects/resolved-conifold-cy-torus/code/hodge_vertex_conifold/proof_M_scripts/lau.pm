package lau;   # truncated Laurent series in u as hashes power => BigRat, powers <= $HI kept
use strict; use warnings; use Math::BigRat;
our $HI = 24;
sub R { Math::BigRat->new(@_) }
sub mono { my ($p,$c)=@_; return { $p => (defined $c ? $c : R(1)) }; }
sub add { my ($a,$b)=@_; my %c=%$a; for my $k (keys %$b){ $c{$k} = exists $c{$k} ? $c{$k}+$b->{$k} : $b->{$k}->copy; } clean(\%c) }
sub scal { my ($s,$a)=@_; my %c; for my $k (keys %$a){ $c{$k}=$a->{$k}*$s; } clean(\%c) }
sub clean { my $a=shift; for my $k (keys %$a){ delete $a->{$k} if $a->{$k}->is_zero || $k > $HI; } $a }
sub mul { my ($a,$b)=@_; my %c; for my $i (keys %$a){ for my $j (keys %$b){ my $p=$i+$j; next if $p > $HI; my $v=$a->{$i}*$b->{$j}; $c{$p} = exists $c{$p} ? $c{$p}+$v : $v; } } clean(\%c) }
sub val { my $a=shift; my @k = sort { $a <=> $b } keys %$a; return @k ? $k[0] : undef }
sub inv { my ($a)=@_; my $v = val($a); my $a0 = $a->{$v};
  # a = a0 u^v (1 + r), r has positive powers
  my %r; for my $k (keys %$a){ next if $k==$v; $r{$k-$v} = $a->{$k}/$a0; }
  my $lim = $HI - (-$v); my $out = { 0 => R(1) }; my $pw = { 0 => R(1) }; my $rr = \%r;
  local $HI = $HI + 2*abs($v) + 2;
  for my $n (1..$HI){ $pw = mul($pw, $rr); last unless %$pw; $out = add($out, scal(R((-1)**$n), $pw)); }
  my %o; for my $k (keys %$out){ $o{$k-$v} = $out->{$k}/$a0; } clean(\%o) }
sub res { my $a=shift; exists $a->{-1} ? $a->{-1} : R(0) }
sub fromarr { my ($arr,$sh)=@_; $sh//=0; my %c; for my $i (0..$#$arr){ $c{$i+$sh}=$arr->[$i] unless $arr->[$i]->is_zero; } \%c }
sub show { my $a=shift; join(" + ", map { $a->{$_}->bstr."*u^$_" } sort { $a <=> $b } keys %$a) }
1;
