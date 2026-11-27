package RefR2;
# Referee 2 (round 2, route R2) -- small exact library, written from scratch (no code shared with the author's scripts).
#  * Bernoulli numbers, factorials, binomials (Math::BigRat)
#  * formula (6.1) as an exact polynomial in (sigma, j):  sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lam^{2h}] S(e1 lam)^{1+j e2} S(e2 lam)^{1+j e1},
#    e1 + e2 = 1, e1 e2 = -sigma, S(y) = (y/2)/sinh(y/2), log S(y) = -Sum_k B_{2k} y^{2k}/(2k (2k)!)
#  * engine-table parser (format  h | r | k | a1,..,ar:coef ; ...  = coefficient of j1^a1..jr^ar in Xi^(k), Xi = Sum_k (-1/sigma)^k Xi^(k))
use strict; use warnings; use Math::BigRat; use Math::BigInt;
our @EXPORT_OK;
sub Q { Math::BigRat->new(shift) }
my @BN = (Q(1));
sub bern { my $n = shift; while ($#BN < $n) { my $m = scalar @BN; my $s = Q(0); for my $k (0..$m-1) { $s += Q(binomI($m + 1, $k)) * $BN[$k]; } push @BN, -$s / ($m + 1); } return $BN[$n]; }
sub binomI { my ($n, $k) = @_; return Math::BigInt->new(0) if $k < 0 || $k > $n; my $r = Math::BigInt->new(1); for my $i (1..$k) { $r->bmul($n - $k + $i); $r->bdiv($i); } return $r; }
sub fact { my $n = shift; my $r = Math::BigInt->new(1); $r->bmul($_) for 2..$n; return Q($r); }
# bivariate polynomials in (s = sigma, j): hash "a,b" => coef  (s^a j^b)
sub padd { my ($p, $q, $c) = @_; $c //= Q(1); my %r = map { $_ => $p->{$_}->copy } keys %$p;
  for my $k (keys %$q) { my $x = $q->{$k} * $c; if (exists $r{$k}) { $r{$k} += $x } else { $r{$k} = $x } } for (keys %r) { delete $r{$_} if $r{$_}->is_zero } return \%r; }
sub pmul { my ($p, $q) = @_; my %r;
  for my $a (keys %$p) { my ($a1, $a2) = split /,/, $a; for my $b (keys %$q) { my ($b1, $b2) = split /,/, $b; my $k = ($a1 + $b1) . "," . ($a2 + $b2);
      my $x = $p->{$a} * $q->{$b}; if (exists $r{$k}) { $r{$k} += $x } else { $r{$k} = $x } } }
  for (keys %r) { delete $r{$_} if $r{$_}->is_zero } return \%r; }
sub pconst { my $c = shift; return { "0,0" => Q($c) }; }
sub pscal { my ($p, $c) = @_; my %r = map { $_ => $p->{$_} * $c } keys %$p; for (keys %r) { delete $r{$_} if $r{$_}->is_zero } return \%r; }
# power sums p_m(sigma) = e1^m + e2^m as polynomials in s (j-degree 0)
sub powersums { my $M = shift; my @p = ({ "0,0" => Q(2) }, { "0,0" => Q(1) });
  for my $m (2..$M) { push @p, padd($p[$m-1], pmul({ "1,0" => Q(1) }, $p[$m-2])); } return @p; }
# f_n(s,j) = [lam^{2n}] S(e1 lam)^{1+j e2} S(e2 lam)^{1+j e1}, n = 0..H
sub fseries { my $H = shift; my @p = powersums(2*$H + 1);
  # exponent coefficients e_k = c_k (p_{2k} - j s p_{2k-1}), c_k = -B_{2k}/(2k (2k)!)
  my @e = ({});
  for my $k (1..$H) { my $c = -bern(2*$k) / (Q(2*$k) * fact(2*$k)); my $t = padd($p[2*$k], pmul({ "1,1" => Q(-1) }, $p[2*$k-1])); push @e, pscal($t, $c); }
  my @f = (pconst(1));
  for my $n (1..$H) { my $acc = {}; for my $k (1..$n) { $acc = padd($acc, pmul($e[$k], $f[$n-$k]), Q($k)); } push @f, pscal($acc, Q(1)/$n); }
  return @f; }
# sigma^h Xi_h(j) by (6.1), as a polynomial in (s, j)
sub rhs61 { my $h = shift; my @f = fseries($h); my $P = $f[$h];
  for my $l (1..2*$h-2) { $P = pmul($P, { "0,1" => Q(1), "0,0" => Q(-$l) }); } return $P; }
# components Xi^(k)(j) from (6.1): Xi^(k) = (-1)^k [s^{h-k}] rhs  ->  hash k => { b => coef of j^b }
sub comps61 { my $h = shift; my $P = rhs61($h); my %C;
  for my $key (keys %$P) { my ($a, $b) = split /,/, $key; my $k = $h - $a; die "negative component index" if $k < 0; $C{$k}{$b} = $P->{$key} * (($k % 2) ? Q(-1) : Q(1)); }
  return \%C; }
sub eval_poly { my ($P, $s, $j) = @_; my $v = Q(0); for my $key (keys %$P) { my ($a, $b) = split /,/, $key; $v += $P->{$key} * $s**$a * $j**$b; } return $v; }
# formula (6.1) value Xi_h(j) at numeric tau
sub xi61 { my ($h, $tau, $j, $cache) = @_; my $P = ($cache->{$h} //= rhs61($h)); my $s = $tau * (1 + $tau); return eval_poly($P, $s, Q($j)) / $s**$h; }
# engine table: returns hash k => { "a1,..,ar" => coef }
sub load_table { my ($file, $h, $r) = @_; my %C; open my $fh, '<', $file or die "cannot open $file";
  while (my $line = <$fh>) { $line =~ s/\r//g; chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; my ($hh, $rr, $k, $rest) = @f;
    next unless $hh == $h && $rr == $r; $C{$k} //= {};
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /^([\d,]+):(\S+)$/; $C{$k}{$1} = Q($2); } }
  close $fh; return \%C; }
sub table_value { my ($C, $tau, @j) = @_; my $sig = $tau * (1 + $tau); my $v = Q(0);
  for my $k (keys %$C) { my $s = Q(0); for my $key (keys %{$C->{$k}}) { my @a = split /,/, $key; my $t = $C->{$k}{$key}->copy; $t *= Q($j[$_])**$a[$_] for 0..$#a; $s += $t; }
    $v += $s * (Q(-1) / $sig)**$k; }
  return $v; }
sub Pi { my ($j, $tau) = @_; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
1;
