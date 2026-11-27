package RefMod;
# Referee (TR2, round 2): modular arithmetic and truncated power series mod a prime P < 2^31.
# Written independently of the author's TRmod.pm.  ASCII only.
use strict; use warnings; use integer;
use Exporter 'import';
our @EXPORT = qw(setp md mulm addm subm negm invm powm ratm binm psmul psinv pspow psadd psscal prime_list);
our $P;
sub setp { $P = shift; }
sub md   { my $a = $_[0] % $P; $a += $P if $a < 0; $a }
sub mulm { ($_[0] * $_[1]) % $P }
sub addm { my $s = $_[0] + $_[1]; $s -= $P if $s >= $P; $s }
sub subm { my $s = $_[0] - $_[1]; $s += $P if $s < 0; $s }
sub negm { $_[0] ? $P - $_[0] : 0 }
sub invm { my $a = md($_[0]); die "inverse of 0 mod $P\n" unless $a; my ($t,$nt,$r,$nr) = (0,1,$P,$a);
  while ($nr) { my $q = $r / $nr; ($t,$nt) = ($nt, $t - $q*$nt); ($r,$nr) = ($nr, $r - $q*$nr); } $t += $P if $t < 0; $t }
sub powm { my ($a,$e) = @_; $a = md($a); if ($e < 0) { $a = invm($a); $e = -$e; } my $r = 1;
  while ($e) { $r = ($r*$a) % $P if $e & 1; $e >>= 1; $a = ($a*$a) % $P if $e; } $r }
sub ratm { my ($n,$d) = @_; mulm(md($n), invm($d)) }
# generalized binomial C(a,k), a a residue, k >= 0 an integer
sub binm { my ($a,$k) = @_; return 0 if $k < 0; my ($num,$den) = (1,1);
  for my $i (0..$k-1) { $num = mulm($num, subm(md($a), md($i))); $den = mulm($den, md($i+1)); } mulm($num, invm($den)) }
sub psmul { my ($a,$b,$N) = @_; my @c = (0) x ($N+1);
  for my $i (0..$#$a) { last if $i > $N; my $ai = $a->[$i] or next; my $lim = $N - $i; $lim = $#$b if $#$b < $lim;
    for my $j (0..$lim) { my $bj = $b->[$j] or next; $c[$i+$j] = ($c[$i+$j] + $ai*$bj) % $P; } } \@c }
sub psinv { my ($a,$N) = @_; my $i0 = invm($a->[0]); my @b = ($i0);
  for my $n (1..$N) { my $s = 0; for my $k (1..$n) { last if $k > $#$a; my $ak = $a->[$k] or next; $s = ($s + $ak*$b[$n-$k]) % $P; }
    push @b, mulm(negm($s), $i0); } \@b }
sub pspow { my ($a,$m,$N) = @_; my $r = [1, (0) x $N]; my $b = $a; if ($m < 0) { $b = psinv($a,$N); $m = -$m; }
  while ($m) { $r = psmul($r,$b,$N) if $m & 1; $m >>= 1; $b = psmul($b,$b,$N) if $m; } $r }
sub psadd { my ($a,$b,$N) = @_; [ map { addm(($_ <= $#$a ? $a->[$_] : 0), ($_ <= $#$b ? $b->[$_] : 0)) } 0..$N ] }
sub psscal { my ($a,$c) = @_; [ map { mulm($_, $c) } @$a ] }
sub prime_list { my ($k) = @_; my @p; my $n = 2147483647;
  while (@p < $k) { my $ok = 1; for (my $d = 3; $d*$d <= $n; $d += 2) { if ($n % $d == 0) { $ok = 0; last; } } push @p, $n if $ok && $n % 2; $n -= 2; } @p }
1;
