# Hsym3.pl : symmetry test of H for three leaves. From the slice series of prefixes (m1,m2) (slice_fitR.pl output files
# sliceR_<m1>_<m2>.txt), H_{(m1,m2)}(J; lambda) = B^{-(m1+m2)} Sum_s C_s^{(m1,m2)}(lambda)/((J-1)...(J-s)) evaluated at J = m1+m2+j.
# For a triple (a,b,c) the two prefixes (a,b) [free c] and (a,c) [free b] must give the same series if H is symmetric.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use SerS;
sub R { Math::BigRat->new(shift) }
$SerS::K = 4;
sub loadC { my $f = shift; my %C; open my $fh, '<', $f or die "no $f";
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$s}{$n}{$2} = R($1); } } return \%C; }
sub Hser { my ($C, $M, $J) = @_; my $res = SerS::szero();
  for my $s (keys %$C) { my $den = R(1); $den *= ($J - $_) for 1..$s; my $ser = SerS::szero(); for my $n (keys %{$C->{$s}}) { next if $n > $SerS::K; my @p; for my $p (keys %{$C->{$s}{$n}}) { $p[$p] = $C->{$s}{$n}{$p}; } $ser->[$n] = [ map { $_ // R(0) } @p ]; }
    $res = SerS::sadd($res, SerS::sscale($ser, R(1)/$den)); }
  my $Bm = SerS::sexp(SerS::sscale(SerS::beta(), R(-$M))); return SerS::smul($Bm, $res); }
sub show { my $x = shift; return join(' | ', map { my $n = $_; "l^" . 2*$n . ": " . join(' ', map { "($x->[$n][$_])s^$_" } grep { ($x->[$n][$_] // 0) != 0 } 0..$#{$x->[$n]}) } 0..$SerS::K); }
for my $t (@ARGV) { my ($a, $b, $c) = split /,/, $t; my $J = $a + $b + $c;
  my $H1 = Hser(loadC("$FindBin::Bin/sliceR_${a}_${b}.txt"), $a + $b, $J);
  my $H2 = Hser(loadC("$FindBin::Bin/sliceR_${a}_${c}.txt"), $a + $c, $J);
  my $same = 1; for my $n (0..$SerS::K) { for my $p (0..6) { $same = 0 if (($H1->[$n][$p] // R(0)) != ($H2->[$n][$p] // R(0))); } }
  print "triple ($a,$b,$c): H via prefix ($a,$b) ", ($same ? "EQUALS" : "differs from"), " H via prefix ($a,$c) through lambda^", 2*$SerS::K, "\n";
  print "   H = ", show($H1), "\n" unless $same; print "   H' = ", show($H2), "\n" unless $same; }
