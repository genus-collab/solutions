# ref_c61.pl table h... : referee's own check of formula (6.1) against a one-leaf component table, as an EXACT polynomial identity in (j, sigma):
#   sigma^h Xi_h(j) = (j-1)_(2h-2) [l^{2h}] exp( Sum_m c_m l^{2m} (p_{2m} - j sigma p_{2m-1}) ),   c_m = -B_{2m}/(2m (2m)!),
#   p_n = e1^n + e2^n (e1 = -tau, e2 = 1+tau): p_0 = 2, p_1 = 1, p_n = p_{n-1} + sigma p_{n-2};
#   left side: Sum_k (-1)^k sigma^{h-k} Xi^(k)_h(j) from the table.
# Also: component k = 0 vs Theorem B, k = h vs the lambda_g formula b_h (j-1)...(j-2h+2), b_h = [z^2h] (z/2)/sin(z/2).
# Usage: perl ref_c61.pl <dir with ref_engine_h<h>_r1.txt> h1 h2 ...
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL;
sub Q { Math::BigRat->new(shift) }
my ($dir, @hs) = @ARGV;
for my $h (@hs) {
  my $f = "$dir/ref_engine_h${h}_r1.txt"; open my $fh, '<', $f or die "no $f"; my @C;
  while (my $l = <$fh>) { chomp $l; next unless $l =~ /\|/; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; die "wrong genus" unless $hh == $h; my @p = (Q(0)); $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($a, $v) = split /:/, $t; $p[$a] = Q($v); } $_ //= Q(0) for @p; $C[$k] = RL::ptrim(\@p); } close $fh;
  # left side as bivariate in (j, sigma): index = power of sigma
  my $lhs = [ map { [Q(0)] } 0..$h ];
  for my $k (0..$h) { next unless $C[$k]; $lhs->[$h - $k] = RL::padd($lhs->[$h - $k], $C[$k], Q((-1) ** $k)) }
  # right side
  my @pn = ([Q(2)], [Q(1)]); for my $n (2..2*$h) { $pn[$n] = RL::padd($pn[$n-1], [Q(0), @{$pn[$n-2]}]) }   # polynomials in sigma
  # exponent E(l) = Sum_m l^{2m} e_m(j, sigma), e_m bivariate (index sigma power, entries j-polys)
  my @e; for my $m (1..$h) { my $cm = -RL::bern(2*$m) / (Q(2*$m) * RL::fact(2*$m));
    my $b = [ map { [$_ * $cm] } @{$pn[2*$m]} ];                              # c_m p_{2m}
    my $c2 = [ [Q(0)], map { [Q(0), -$_ * $cm] } @{$pn[2*$m - 1]} ];           # - c_m j sigma p_{2m-1}
    $e[$m] = RL::badd($b, $c2); }
  # exp via the recursion n a_n = Sum_m m e_m a_{n-m} (in powers of l^2)
  my @a = ([[Q(1)]]); for my $n (1..$h) { my $s = [[Q(0)]]; for my $m (1..$n) { $s = RL::badd($s, RL::bsc(RL::bmul($e[$m], $a[$n - $m]), Q($m))) } $a[$n] = RL::bsc($s, Q(1) / $n) }
  my $rhs = RL::bmulj($a[$h], RL::ffj(-1, 2*$h - 2));
  my $d = RL::badd($lhs, $rhs, Q(-1)); my $bad = 0; for my $p (@$d) { $bad++ unless RL::pzero($p) }
  # components k = 0 and k = h
  $RL::NT = 2*$h + 4; my $S2 = RL::sexp(RL::ssc(RL::logS(), Q(2))); my $Xh = RL::pmul(RL::ffj(-1, 2*$h - 2), RL::coeff_expj($S2, RL::logH(), 2*$h));
  my $k0 = RL::peq($Xh, $C[0] // [Q(0)]);
  # b_h = [z^{2h}] (z/2)/sin(z/2) = (-1)^h [z^2h] S(z)... use S(iz): (z/2)/sin(z/2) = exp(Sum_m c_m (-1)^m z^{2m})
  my $ls = RL::logS(); my $lsin = RL::szero(); for my $m (1..$h) { $lsin->[2*$m] = $ls->[2*$m] * Q((-1) ** $m) } my $bh = RL::sexp($lsin)->[2*$h];
  my $kh = RL::peq(RL::psc(RL::ffj(-1, 2*$h - 2), $bh), $C[$h] // [Q(0)]);
  printf "genus %d (%s): (6.1) as a polynomial identity in (j, sigma): %s;  k=0 vs Theorem B: %s;  k=h vs lambda_g formula: %s\n", $h, $f, ($bad ? "$bad sigma-coefficients DIFFER" : "holds"), ($k0 ? "equal" : "DIFFER"), ($kh ? "equal" : "DIFFER"); }
