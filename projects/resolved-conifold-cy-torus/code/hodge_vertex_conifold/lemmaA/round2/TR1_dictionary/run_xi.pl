#!/usr/bin/perl
# run_xi.pl tau G N [ycase] [e] [nmin]
#   Exact EO recursion on C_tau in the Laplace basis (TRxi.pm) for all (g,n) with g + n <= G + N, then for every (g,n)
#   with g <= G, nmin <= n <= N: symmetry check and coefficient-by-coefficient comparison (all ordered k, |k| <= 3g-3+n, and
#   nothing outside) with the Hodge side (see run_p2.pl for the formulas: base = P2, bkmp/zhou = Zhou's theorem, mix = family).
#   Env SAVE=file writes the computed xi-coefficients (g | n | k:coef ; ...) for later use.
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use TRxi; use HodgeSide;
sub R { Math::BigRat->new(@_) }
my ($tau, $G, $N, $ycase, $e, $nmin) = @ARGV; $ycase //= 'base'; $e //= 0; $nmin //= 1;
die "usage: run_xi.pl tau G N [base|bkmp|zhou|mix] [e] [nmin]\n" unless defined $N;
my $t0 = time;
my $kmax = 3*$G - 3 + $N; $kmax = 1 if $kmax < 1;
TRxi::setup(tau => $tau, kmax => $kmax, y => $ycase, e => $e);
my $sigma = $TRxi::sigma; my $tauR = R($tau); my $eR = R($e);
printf "tau = %s, sigma = %s, y = %s%s; G = %d, N = %d, kmax = %d, series length %d; setup %.1f s\n",
  $tauR->bstr, $sigma->bstr, $ycase, ($ycase eq 'mix' ? " (e = " . $eR->bstr . ")" : ""), $G, $N, $kmax, $TRxi::L, time - $t0;
my @order = TRxi::run_all($G + $N, $G);
printf "recursion: %d correlators, %.1f s\n", scalar(@order), time - $t0;
sub tuples { my ($n, $maxsum) = @_; my @out; my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x (0..$left) { $rec->([@$pre, $x], $left - $x); } }; $rec->([], $maxsum); return @out; }
sub perms { my @x = @_; return ([]) unless @x; my @out; my %seen;
  for my $i (0..$#x) { next if $seen{$x[$i]}++; my @rest = @x; my ($h) = splice(@rest, $i, 1); push @out, [$h, @$_] for perms(@rest); } return @out; }
my $totbad = 0; my $totcmp = 0; my $sv;
if ($ENV{SAVE}) { open $sv, '>', $ENV{SAVE} or die; print $sv "# tau = $tau, y = $ycase, e = $e; coefficients of prod d xi_k in omega_{g,n} (TRxi.pm)\n"; }
for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %TRxi::W) {
  my ($g, $n) = split /,/, $gn; next if $g > $G || $n > $N || $n < $nmin;
  my $Wx = $TRxi::W{$gn};
  if ($sv) { print $sv "$g | $n | ", join(" ; ", map { "$_:" . $Wx->{$_}->bstr } sort keys %$Wx), "\n"; }
  my %s; my $symbad = 0;
  for my $key (keys %$Wx) { my $sk = join(",", sort { $a <=> $b } split /,/, $key); push @{$s{$sk}}, $Wx->{$key}; }
  for my $sk (keys %s) { my @v = @{$s{$sk}}; for (@v) { $symbad++ if $_ != $v[0]; } my @pp = perms(split /,/, $sk); $symbad++ if scalar(@pp) != scalar(@v); }
  my $dim = 3*$g - 3 + $n; my $bad = 0; my $cmp = 0; my @ex;
  for my $key (keys %$Wx) { my $sum = 0; $sum += $_ for split /,/, $key; if ($sum > $dim) { $bad++; push @ex, "outside:$key"; } }
  my $avail = 1;
  for my $k (tuples($n, $dim)) {
    my $got = $Wx->{join(",", @$k)} // R(0); my $exp;
    if ($ycase eq 'base') { unless (HodgeSide::have_q0($g, $n)) { $avail = 0; last; }
      $exp = $sigma ** (1 - $g) * HodgeSide::hodge_q0($g, $n, $k, $sigma); }
    elsif ($ycase eq 'mix') { my $s = -$eR * $sigma / (1 + $eR * $sigma); my $Q;
      if ($g == 0) { $Q = genus0_family($n, $k, $s); }
      else { my $P = HodgeSide::datapoly('q', $g, $n, $tauR->bstr); unless ($P) { $avail = 0; last; }
        $Q = R(0); for my $m (0..3*$g-3+$n) { my $cf = HodgeSide::polycoef($P, 'j', $k, $m); $Q += $cf * $s ** $m unless $cf->is_zero; } }
      $exp = $sigma ** ($g - 1 + $n) * (1 + $eR * $sigma) ** (2 - 2*$g - $n) * $Q; }
    else { my $mvval;
      if ($g == 0) { my $K = 0; $K += $_ for @$k; $mvval = R(0); if ($K == $n - 3) { $mvval = HodgeSide::fact($n - 3); $mvval /= HodgeSide::fact($_) for @$k; } }
      else { my $P = HodgeSide::datapoly('mv', $g, $n, $tauR->bstr); unless ($P) { $avail = 0; last; } $mvval = HodgeSide::polycoef($P, 'm', $k, 0); }
      $exp = $sigma ** ($g - 1 + $n) * $mvval; $exp = -$exp if $ycase eq 'zhou' && $n % 2; }
    $cmp++; if ($got != $exp) { $bad++; push @ex, "k=(" . join(",", @$k) . ") TR=$got Hodge=$exp" if @ex < 3; } }
  my $line = sprintf "(%d,%d): %5d xi-coefficients, symmetry violations %d", $g, $n, scalar(keys %$Wx), $symbad;
  if ($avail) { $line .= sprintf "; vs Hodge: %d compared, mismatches %d", $cmp, $bad; $totbad += $bad + $symbad; $totcmp += $cmp; }
  else { $line .= "; (no Hodge table)"; $totbad += $symbad; }
  $line .= "  [" . join(" | ", @ex) . "]" if @ex;
  print "$line\n"; }
close $sv if $sv;
printf "TOTAL: %d coefficient comparisons, %d failures; time %.1f s\n", $totcmp, $totbad, time - $t0;
sub genus0_family { my ($n, $k, $s) = @_; my $K = 0; $K += $_ for @$k; my $tot = R(0);
  for my $m (0..$n) { my $deg = $n + $m - 3 - $K; next if $deg < 2*$m;
    my $sum = R(0); my $rec; $rec = sub { my ($left, $cnt, $acc) = @_; if ($cnt == 0) { $sum += $acc if $left == 0; return; }
      for my $d (2..$left) { $rec->($left - $d, $cnt - 1, $acc / HodgeSide::fact($d)); } };
    $rec->($deg, $m, R(1));
    my $v = HodgeSide::fact($n + $m - 3) * $sum; $v /= HodgeSide::fact($_) for @$k;
    $tot += $s ** $m / HodgeSide::fact($m) * $v; }
  return $tot; }
