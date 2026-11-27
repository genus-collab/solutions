# t_lines.pl h : report 0(f)(1) on the referee data: for r = 2, which lines J = l divide Xi^(k)_h(j1,j2)?  Unlike the author's
# explore_lines.pl (integer l in [-3, 3h+2] only) we find ALL complex l: P(t, l - t) == 0 in t means every t-coefficient of
# P(t, l - t), a polynomial in l, vanishes; we compute the gcd of these polynomials in l (over Q) and report its degree and rational roots.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $h = shift;
my @C; open my $fh, '<', "$FindBin::Bin/ref_engine_h${h}_r2.txt" or die; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
  for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; my ($a, $b) = split /,/, $m; $C[$k]{"$a,$b"} = Q($cf); } } close $fh;
sub binom { my ($n, $k) = @_; my $x = Q(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } $x }
# polynomials in l as arrays (index = power)
sub ptrim { my $p = shift; pop @$p while @$p && $p->[-1] == 0; return $p; }
sub pmod { my ($a, $b) = @_; my @r = @$a; my @d = @$b; ptrim(\@r); ptrim(\@d); while (@r >= @d && @r) { my $f = $r[-1] / $d[-1]; my $sh = @r - @d; for my $i (0..$#d) { $r[$sh + $i] -= $f * $d[$i]; } ptrim(\@r); } return \@r; }
sub pgcd { my ($a, $b) = @_; $a = ptrim([@$a]); $b = ptrim([@$b]); while (@$b) { my $r = pmod($a, $b); ($a, $b) = ($b, $r); } return $a; }
for my $k (0..$h) { my $P = $C[$k] or next;
  # P(t, l - t) = Sum_{a,b} c_ab t^a (l - t)^b ; coefficient of t^n is a polynomial in l
  my %coef;   # n -> array in l
  for my $key (keys %$P) { my ($a, $b) = split /,/, $key; my $c = $P->{$key};
    for my $i (0..$b) { # (l - t)^b = Sum_i C(b,i) l^{b-i} (-t)^i
      my $n = $a + $i; my $x = $c * binom($b, $i) * (($i % 2) ? -1 : 1); $coef{$n}[$b - $i] = ($coef{$n}[$b - $i] // Q(0)) + $x; } }
  my $g; for my $n (sort { $a <=> $b } keys %coef) { my $p = [ map { $_ // Q(0) } @{$coef{$n}} ]; ptrim($p); next unless @$p; $g = defined $g ? pgcd($g, $p) : $p; }
  my $deg = @$g - 1; my @roots;
  for my $num (-3*$h .. 6*$h) { for my $den (1..4) { my $l = Q($num) / $den; next if $den > 1 && $l->is_int(); my $v = Q(0); my $pw = Q(1); for my $c (@$g) { $v += $c * $pw; $pw *= $l; } push @roots, "$l" if $v == 0; } }
  printf "genus %d, k = %d: the lines J = l dividing Xi^(k) are the roots of a polynomial of degree %d; rational roots found: {%s}\n", $h, $k, $deg, join(' ', @roots); }
