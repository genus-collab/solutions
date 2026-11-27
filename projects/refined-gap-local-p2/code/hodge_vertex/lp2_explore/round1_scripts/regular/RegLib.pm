package RegLib;
# Shared helpers for the regular-colour tests (exact, Math::BigRat).
# Laurent polynomials in s are hashrefs {exponent => BigRat}; z-series are arrayrefs of those (index = z-power).
# Q(omega) numbers are [a,b] = a + b*omega, omega^2 = -1 - omega.
use strict; use warnings; use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(QQ parse_cpp ev_part lp_add lp_scale lp_mul lp_theta lp_shift lp_clean lp_str lp_eval
                 lp_minexp lp_maxexp w_add w_sub w_mul w_inv w_scale w_str w_isint w_eq w_fromq w_pow);
sub QQ { Math::BigRat->new(shift) }
our ($e1,$e2,$s3,$sqv,$L);
my $bs = chr(92);
sub parse_cpp {
  my $f = shift; local $/; open my $fh,'<',$f or die "$f: $!"; my $s=<$fh>; close $fh;
  $s =~ s/\s+//g; $s =~ s/^\{//; $s =~ s/\}$//;
  my @parts; my $d=0; my $cur='';
  for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
    if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
  return @parts;
}
# exact evaluation of one stored coefficient at eps1 = $E1, eps2 = $E2, LLL = $Lv, Sqrt[3] -> $s3v, sq -> 1
sub ev_part { my ($x,$E1,$E2,$Lv,$s3v)=@_; my $t=$x; $t =~ s/\Q$bs\E//g;
  $t =~ s/\Q[Epsilon][1]\E/\$e1/g; $t =~ s/\Q[Epsilon][2]\E/\$e2/g; $t =~ s/\QLLL[0][y]\E/\$L/g;
  $t =~ s/\QSqrt[3]\E/\$s3/g; $t =~ s/\bsq\b/\$sqv/g; $t =~ s/\^/**/g; $t =~ tr/[]/()/;
  $t =~ s/(?<![\$\w.])(\d+)/QQ('$1')/g;
  local $e1 = $E1; local $e2 = $E2; local $L = $Lv; local $s3 = $s3v; local $sqv = QQ(1);
  my $v = eval $t; die "$@ in $t" if $@; return $v; }
# ---- Laurent polynomials in s ----
sub lp_clean { my $p=shift; for my $k (keys %$p){ delete $p->{$k} if $p->{$k}==0 } return $p }
sub lp_add { my %r; for my $p (@_){ for my $k (keys %$p){ $r{$k} = ($r{$k}//QQ(0)) + $p->{$k} } } return lp_clean(\%r) }
sub lp_scale { my ($p,$c)=@_; my %r; $r{$_} = $p->{$_}*$c for keys %$p; return lp_clean(\%r) }
sub lp_mul { my ($p,$q)=@_; my %r; for my $i (keys %$p){ for my $j (keys %$q){ $r{$i+$j} = ($r{$i+$j}//QQ(0)) + $p->{$i}*$q->{$j} } } return lp_clean(\%r) }
sub lp_shift { my ($p,$n)=@_; my %r; $r{$_+$n} = $p->{$_} for keys %$p; return \%r }
# theta_y = v d/ds, v = (1 - alpha s)/(3 alpha)
sub lp_theta { my ($p,$al)=@_; my %r;
  for my $n (keys %$p){ next if $n==0; my $c = $p->{$n}*$n/(3*$al);
    $r{$n-1} = ($r{$n-1}//QQ(0)) + $c; $r{$n} = ($r{$n}//QQ(0)) - $c*$al; } return lp_clean(\%r) }
sub lp_minexp { my $p=shift; my @k = sort {$a<=>$b} keys %$p; return @k ? $k[0] : undef }
sub lp_maxexp { my $p=shift; my @k = sort {$a<=>$b} keys %$p; return @k ? $k[-1] : undef }
sub lp_str { my $p=shift; my @k = sort {$a<=>$b} keys %$p; return '0' unless @k; join(' + ', map { "($p->{$_}) s^$_" } @k) }
# ---- Q(omega) ----
sub w_fromq { [QQ($_[0]), QQ(0)] }
sub w_add { [ $_[0][0]+$_[1][0], $_[0][1]+$_[1][1] ] }
sub w_sub { [ $_[0][0]-$_[1][0], $_[0][1]-$_[1][1] ] }
sub w_scale { [ $_[0][0]*$_[1], $_[0][1]*$_[1] ] }
sub w_mul { my ($x,$y)=@_; my ($a,$b)=@$x; my ($c,$d)=@$y;   # (a+bw)(c+dw) = ac + (ad+bc)w + bd w^2, w^2 = -1-w
  [ $a*$c - $b*$d, $a*$d + $b*$c - $b*$d ] }
sub w_inv { my $x=shift; my ($a,$b)=@$x; my $n = $a*$a - $a*$b + $b*$b; die "zero" if $n==0;
  [ ($a-$b)/$n, -$b/$n ] }   # conj(a+bw) = a + b w^2 = (a-b) - b w
sub w_pow { my ($x,$n)=@_; my $r=[QQ(1),QQ(0)]; if($n<0){ $x = w_inv($x); $n=-$n } $r = w_mul($r,$x) for 1..$n; return $r }
sub w_eq { $_[0][0]==$_[1][0] && $_[0][1]==$_[1][1] }
sub w_isint { 1 }
sub w_str { my $x=shift; "($x->[0]) + ($x->[1]) w" }
sub lp_eval { my ($p,$sw)=@_; my $r=[QQ(0),QQ(0)]; for my $n (keys %$p){ $r = w_add($r, w_scale(w_pow($sw,$n), $p->{$n})) } return $r }
1;
