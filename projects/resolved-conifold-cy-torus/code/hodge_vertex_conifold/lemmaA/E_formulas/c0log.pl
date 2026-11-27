# c0log.pl : log(C_0(m)/S1S2) and (C_s(m)/C(m,s+1)) structure, as polynomials in m (from the slice files m = 0..11)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use SerS;
sub R { Math::BigRat->new(shift) }
my %C; for my $m (1..16) { open my $fh, '<', "$FindBin::Bin/slice_m$m.txt" or next;
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = R($1); } } }
my $S12 = SerS::S1S2(); for my $n (0..5) { my $c = SerS::coeff($S12, $n); for my $p (0..$#$c) { $C{0}{0}{$n}{$p} = $c->[$p] if $c->[$p] != 0; } }
$SerS::K = 5;
sub ser { my ($m, $s) = @_; my $x = SerS::szero(); for my $n (0..5) { my @p = map { $C{$m}{$s}{$n}{$_} // R(0) } 0..5; $x->[$n] = \@p; } return $x; }
sub slog { my $a = shift; my $X = SerS::sadd($a, SerS::sscale(SerS::sone(), R(-1))); my $L = SerS::szero(); my $pw = SerS::sone();
  for my $n (1..5) { $pw = SerS::smul($pw, $X); $L = SerS::sadd($L, SerS::sscale($pw, R((($n % 2) ? 1 : -1)) / $n)); } return $L; }
my $lS = slog($S12);
for my $m (0..8) { my $L = SerS::sadd(slog(ser($m, 0)), SerS::sscale($lS, R(-1)));
  print "m=$m log(C_0/S1S2): ", join(' | ', map { my $n = $_; "l^" . 2*$n . ": " . join(' ', map { "($L->[$n][$_])s^$_" } grep { ($L->[$n][$_] // 0) != 0 } 0..5) } 1..5), "\n"; }
