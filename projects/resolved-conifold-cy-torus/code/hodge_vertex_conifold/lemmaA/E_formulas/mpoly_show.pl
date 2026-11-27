# mpoly_show.pl : print [lambda^{2n} sigma^p] C_s(m) as polynomials in m, divided by C(m, s+1) for s >= 1
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use SerS;
sub R { Math::BigRat->new(shift) }
my %C; for my $m (1..16) { open my $fh, '<', "$FindBin::Bin/slice_m$m.txt" or next;
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = R($1); } } }
my $S12 = SerS::S1S2(); for my $n (0..5) { my $c = SerS::coeff($S12, $n); for my $p (0..$#$c) { $C{0}{0}{$n}{$p} = $c->[$p] if $c->[$p] != 0; } }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
my ($smax, $nmax) = @ARGV; $smax //= 3; $nmax //= 3;
for my $s (0..$smax) { for my $n (0..$nmax) { for my $p (0..$n) {
  my @x; my @y; for my $m (0..11) { next if $s >= 1 && $m <= $s; my $v = $C{$m}{$s}{$n}{$p} // R(0); my $b = $s >= 1 ? binom($m, $s+1) : R(1); push @x, R($m); push @y, $v / $b; }
  # Lagrange -> monomial coefficients
  my $N = scalar @x; my @dd = @y; for my $k (1..$N-1) { for (my $i = $N-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($x[$i] - $x[$i-$k]); } }
  my @p = (R(0)) x $N; my @basis = (R(1)); for my $k (0..$N-1) { for my $i (0..$#basis) { $p[$i] += $dd[$k] * $basis[$i]; } my @nb = (R(0)) x (scalar(@basis) + 1); for my $i (0..$#basis) { $nb[$i+1] += $basis[$i]; $nb[$i] -= $basis[$i] * $x[$k]; } @basis = @nb; }
  my @t = map { "($p[$_])m^$_" } grep { $p[$_] != 0 } 0..$#p; next unless @t;
  printf "[l^%d s^%d] C_%d(m)%s = %s\n", 2*$n, $p, $s, ($s >= 1 ? "/C(m,".($s+1).")" : ""), join(' + ', @t); } } }
