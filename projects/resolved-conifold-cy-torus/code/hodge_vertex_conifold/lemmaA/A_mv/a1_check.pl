#!/usr/bin/perl
# a1_check.pl -- (A1): match the Schur side of the Marino-Vafa formula EXACTLY against lemmaA/data/mv_g<g>_l<l>.txt.
# Usage: perl a1_check.pl g l maxpart [maxsize]
# Convention tested (x = sqrt(-1) lambda, q = e^x, c = 1, sigma = tau(1+tau)):
#   [x^{2g-2+l}] [p_mu] log sum_nu e^{(tau+1/2) kappa_nu x/2} prod_box 1/(2 sinh(h x/2)) s_nu(p)
#        = (-1)^g sigma^{l-1} Pi_mu(tau)/|Aut mu| * I_g(mu),   I_g(mu) = Int Lambda(1)Lambda(tau)Lambda(-1-tau)/prod(1 - mu_i psi_i),
#   Pi_mu = prod_i prod_{a=1}^{mu_i-1} (mu_i tau + a)/(mu_i - 1)!.
# Data side (c-normalisation, w = (c,-x,-1), c = 1+x, tau = -x/(1+x)):  P_g(mu; x) = c^{-(3g-3+2l)} (-sigma)^{-g} I_g(mu).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($g, $l, $maxpart, $maxsize) = @ARGV; $maxsize //= 1000;
my $file = "$FindBin::Bin/../data/mv_g${g}_l${l}.txt";
my $recs = parse_file($file); my ($gg, $ll, $expr) = @{$recs->[0]}; die "header" unless $gg == $g && $ll == $l;
my $t0 = time; my $ast = parse_expr($expr); printf "parsed %s (%d chars) in %d s\n", $file, length($expr), time - $t0;
my @V = map { "m$_" } 1..$l;
my @taus = map { R($_) } ('2/7', '3/2', '5', '1/4', '-7/3', '-5/2', '7/11');   # 7 values outside [-1,0]: distinct non-zero sigma, Pi_mu != 0
# all sorted l-tuples with parts <= maxpart, size <= maxsize
my @mus; my $rec; $rec = sub { my ($pre, $max, $k) = @_; if ($k == 0) { my $s = 0; $s += $_ for @$pre; push @mus, [@$pre] if $s <= $maxsize; return; }
  for my $p (reverse 1..$max) { $rec->([@$pre, $p], $p, $k - 1); } }; $rec->([], $maxpart, $l);
my $maxn = 0; for my $mu (@mus) { my $n = 0; $n += $_ for @$mu; $maxn = $n if $n > $maxn; }
my ($nok, $nbad) = (0, 0);
for my $tau (@taus) {
  my $x = -$tau / (1 + $tau); my $c = 1 + $x; my $sig = $tau * (1 + $tau);
  my $t1 = time; my $poly = eval_poly($ast, \@V, { x => $x }); my $te = time - $t1;
  my @line;
  for my $mu (@mus) { my $n = 0; $n += $_ for @$mu;
    my $P = 2*$g - 2 + $l + $maxn; my $s = mvlib::mv_conn($mu, $tau, $P);
    my $Pi = R(1); for my $j (@$mu) { for my $b (1..$j-1) { $Pi *= ($j*$tau + $b); } $Pi /= fact($j-1); }
    my $Ig = R((-1)**$g) * $s->[2*$g - 2 + $l + $n] * mvlib::autf(@$mu) / ($sig ** ($l - 1) * $Pi);
    my $pred = $Ig / ((-$sig) ** $g) / ($c ** (3*$g - 3 + 2*$l));
    my $dat = poly_eval_at($poly, [ map { R($_) } @$mu ]);
    if ($pred == $dat) { $nok++; } else { $nbad++; push @line, "MISMATCH mu=(@$mu) pred=$pred data=$dat ratio=" . ($dat->is_zero ? 'inf' : $pred/$dat); } }
  printf "tau=%s: %d partitions, data eval %d s%s\n", $tau, scalar @mus, $te, (@line ? "\n  " . join("\n  ", @line[0..($#line < 5 ? $#line : 5)]) : " all match");
}
print "SUMMARY g=$g l=$l maxpart=$maxpart maxsize=$maxsize: $nok matches, $nbad mismatches (", scalar(@taus), " tau values x ", scalar(@mus), " partitions)\n";
