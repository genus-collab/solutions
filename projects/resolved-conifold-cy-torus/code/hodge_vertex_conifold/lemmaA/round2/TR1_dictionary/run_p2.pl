#!/usr/bin/perl
# run_p2.pl tau G N [ycase] [e]
#   Exact EO recursion on C_tau (TRtau.pm) for all (g,n) with g + n <= G + N, then for every (g,n) with g <= G, n <= N:
#   - symmetry of the computed correlators (the recursion treats the first argument specially);
#   - conversion to the Laplace basis d xi_k (xi_k = D^k xi_0, [x^mu] xi_k = mu^k Pi_mu) and comparison, coefficient by
#     coefficient (all ordered k with |k| <= 3g-3+n, and no coefficient outside), with the Hodge side:
#       ycase base : C^xi[k] = sigma^(1-g) <Omega prod psi^k>^{q_0}_{g,n}                     (P2)
#       ycase bkmp : C^xi[k] = sigma^(g-1+n) <Omega prod psi^k>_{g,n}  (no dilaton leaves)     (BKMP, Zhou)
#       ycase zhou : C^xi[k] = (-1)^n sigma^(g-1+n) <Omega prod psi^k>_{g,n}                 (Zhou's eq. (6) literally)
#       ycase mix  : (y = log(1+w) - e/w) C^xi[k] = sigma^(g-1+n) (1+e sigma)^(2-2g-n) Q_g(-e sigma/(1+e sigma); k)
#   - for ycase base and g >= 1: spot checks of the x-expansion  sigma^(g-1) [x^mu] F_{g,n} / prod Pi_{mu_i} = Xi_g(mu)
#     computed through the global form of e_k and Lagrange inversion (independent of the basis change).
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use TRtau; use HodgeSide;
sub R { Math::BigRat->new(@_) }
my ($tau, $G, $N, $ycase, $e) = @ARGV; $ycase //= 'base'; $e //= 0;
die "usage: run_p2.pl tau G N [base|bkmp|zhou|mix] [e]\n" unless defined $N;
my $t0 = time;
my $kmax = 3*$G - 3 + $N; $kmax = 1 if $kmax < 1;
TRtau::setup(tau => $tau, kmax => $kmax, y => $ycase, e => $e);
my $sigma = $TRtau::sigma; my $tauR = R($tau); my $eR = R($e);
printf "tau = %s, sigma = %s, y = %s%s; G = %d, N = %d, kmax = %d, series length %d; setup %.1f s\n",
  $tauR->bstr, $sigma->bstr, $ycase, ($ycase eq 'mix' ? " (e = " . $eR->bstr . ")" : ""), $G, $N, $kmax, $TRtau::L, time - $t0;
printf "omega_{0,3} = %s e0 e0 e0;  omega_{1,1} = %s\n", $TRtau::W{"0,3"}{"0,0,0"} // 'undef', 'computed below' if 0;
my @order = TRtau::run_all($G + $N, $G);
printf "recursion: %d correlators, %.1f s\n", scalar(@order), time - $t0;
sub tuples { my ($n, $maxsum) = @_; my @out; my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x (0..$left) { $rec->([@$pre, $x], $left - $x); } }; $rec->([], $maxsum); return @out; }
my $totbad = 0; my $totcmp = 0;
for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %TRtau::W) {
  my ($g, $n) = split /,/, $gn; next if $g > $G || $n > $N;
  my $We = $TRtau::W{$gn};
  # symmetry
  my %s; my $symbad = 0;
  for my $key (keys %$We) { my $sk = join(",", sort { $a <=> $b } split /,/, $key); push @{$s{$sk}}, $We->{$key}; }
  for my $sk (keys %s) { my @v = @{$s{$sk}}; my $nperm = 1; for (@v) { $symbad++ if $_ != $v[0]; } }
  # also: every permutation present
  for my $key (keys %$We) { my @k = split /,/, $key; my @perms = perms(@k); for my $p (@perms) { $symbad++ unless exists $We->{join(",", @$p)}; } }
  my $Wx = TRtau::to_xi($We);
  my $dim = 3*$g - 3 + $n; my $bad = 0; my $cmp = 0; my @ex;
  for my $key (keys %$Wx) { my $sum = 0; $sum += $_ for split /,/, $key; if ($sum > $dim) { $bad++; push @ex, "outside:$key"; } }
  my $avail = 1;
  for my $k (tuples($n, $dim)) {
    my $got = $Wx->{join(",", @$k)} // R(0); my $exp;
    if ($ycase eq 'base') { unless (HodgeSide::have_q0($g, $n)) { $avail = 0; last; }
      $exp = $sigma ** (1 - $g) * HodgeSide::hodge_q0($g, $n, $k, $sigma); }
    else { my $mvval;
      if ($g == 0) { # genus 0 MV: (n-3)!/prod k_i! if |k| = n-3
        my $K = 0; $K += $_ for @$k; $mvval = R(0); if ($K == $n - 3) { $mvval = HodgeSide::fact($n - 3); $mvval /= HodgeSide::fact($_) for @$k; } }
      if ($ycase eq 'mix') { my $s = -$eR * $sigma / (1 + $eR * $sigma);
        my $Q;
        if ($g == 0) { # Q_0(s;k) = Sum_m s^m/m! <prod psi^k Tt^m>_0, Tt = psi^2/(1-psi); use H(s') = Sum_m s'^m/m! <prod psi^k T^m>
          $Q = genus0_family($n, $k, $s); }
        else { my $P = HodgeSide::datapoly('q', $g, $n, $tauR->bstr); unless ($P) { $avail = 0; last; }
          $Q = R(0); for my $m (0..3*$g-3+$n) { my $cf = HodgeSide::polycoef($P, 'j', $k, $m); $Q += $cf * $s ** $m unless $cf->is_zero; } }
        $exp = $sigma ** ($g - 1 + $n) * (1 + $eR * $sigma) ** (2 - 2*$g - $n) * $Q; }
      else {
        if ($g >= 1) { my $P = HodgeSide::datapoly('mv', $g, $n, $tauR->bstr); unless ($P) { $avail = 0; last; } $mvval = HodgeSide::polycoef($P, 'm', $k, 0); }
        $exp = $sigma ** ($g - 1 + $n) * $mvval; $exp = -$exp if $ycase eq 'zhou' && $n % 2; } }
    $cmp++; if ($got != $exp) { $bad++; push @ex, "k=(" . join(",", @$k) . ") TR=$got Hodge=$exp" if @ex < 3; } }
  my $line = sprintf "(%d,%d): %4d e-coefficients, symmetry violations %d", $g, $n, scalar(keys %$We), $symbad;
  if ($avail) { $line .= sprintf "; xi-basis vs Hodge: %d coefficients compared, mismatches %d", $cmp, $bad; $totbad += $bad + $symbad; $totcmp += $cmp; }
  else { $line .= "; (no Hodge table for this (g,n))"; $totbad += $symbad; }
  $line .= "  [" . join(" | ", @ex) . "]" if @ex;
  # sample: leading coefficients
  my @samp; for my $kk (sort keys %$Wx) { push @samp, "$kk:" . $Wx->{$kk}->bstr; last if @samp >= 2; }
  print "$line\n";
}
# x-expansion spot checks (ycase base): sigma^(g-1) [x^mu] F / prod Pi = Xi_g(mu) = Sum_k <..>^{q_0} prod mu^k
if ($ycase eq 'base') { my $nx = 0; my $xbad = 0; my @lines;
  for my $gn (sort keys %TRtau::W) { my ($g, $n) = split /,/, $gn; next if $g > $G || $n > $N || $g < 1; next unless HodgeSide::have_q0($g, $n);
    my $dim = 3*$g - 3 + $n; my @mus = mutuples($n, 2*$g - 2 + $n + 2);
    my @vals;
    for my $mu (@mus) { my $pi = R(1); $pi *= TRtau::Pi_mu($_) for @$mu; next if $pi->is_zero;
      my $tr = $sigma ** ($g - 1) * TRtau::xcoef($TRtau::W{$gn}, @$mu) / $pi;
      my $xi = R(0); for my $k (tuples($n, $dim)) { my $h = HodgeSide::hodge_q0($g, $n, $k, $sigma); next if $h->is_zero; my $p = R(1); $p *= R($mu->[$_]) ** $k->[$_] for 0..$n-1; $xi += $h * $p; }
      $nx++; if ($tr != $xi) { $xbad++; }
      push @vals, "(" . join(",", @$mu) . ")=" . $tr->bstr if @vals < 6; }
    push @lines, "  ($g,$n) Xi^TR(mu): " . join("  ", @vals); }
  printf "x-expansion check (P2 at lattice points mu, |mu| <= 2g+n): %d points, mismatches %d\n", $nx, $xbad; $totbad += $xbad;
  print "$_\n" for @lines; }
printf "TOTAL: %d coefficient comparisons, %d failures; time %.1f s\n", $totcmp, $totbad, time - $t0;

sub perms { my @x = @_; return ([]) unless @x; my @out; my %seen;
  for my $i (0..$#x) { next if $seen{$x[$i]}++; my @rest = @x; my ($h) = splice(@rest, $i, 1); push @out, [$h, @$_] for perms(@rest); } return @out; }
sub mutuples { my ($n, $maxsum) = @_; my @out; my $rec; $rec = sub { my ($pre, $left, $min) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x ($min..$left - ($n - @$pre - 1)) { $rec->([@$pre, $x], $left - $x, $x); } }; $rec->([], $maxsum, 1); return @out; }
# genus-0 family: Q_0(s;k) = Sum_m (s^m/m!) <prod psi^{k_i} Tt^m>_{0,n+m},  Tt = psi^2/(1-psi) = Sum_{d>=2} psi^d
sub genus0_family { my ($n, $k, $s) = @_; my $K = 0; $K += $_ for @$k; my $tot = R(0);
  for my $m (0..$n) { my $deg = $n + $m - 3 - $K; next if $deg < 2*$m;   # degrees of the m leaves sum to deg, each >= 2
    # Sum over compositions d_1..d_m >= 2 of deg of <prod tau_k prod tau_d>_0 = (n+m-3)!/(prod k! prod d!)
    my $sum = R(0); my $rec; $rec = sub { my ($left, $cnt, $acc) = @_; if ($cnt == 0) { $sum += $acc if $left == 0; return; }
      for my $d (2..$left) { $rec->($left - $d, $cnt - 1, $acc / HodgeSide::fact($d)); } };
    $rec->($deg, $m, R(1));
    my $v = HodgeSide::fact($n + $m - 3) * $sum; $v /= HodgeSide::fact($_) for @$k;
    $tot += $s ** $m / HodgeSide::fact($m) * $v; }
  return $tot; }
