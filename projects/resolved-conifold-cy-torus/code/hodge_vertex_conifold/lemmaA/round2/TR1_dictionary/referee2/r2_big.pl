#!/usr/bin/perl
# r2_big.pl TAU S GMAX [SAVEFILE]      (referee 2, TR1 round 2)
# Runs the author's exact TR on C_tau, y = -1/w (TRxi.pm, copied UNCHANGED from ..\TRxi.pm, md5 7de3a1e9...), for every
# (g,n) with g + n <= S, g <= GMAX, n >= 1, and compares EVERY xi-coefficient with Hodge data that do not use TR:
#   genus 0 : X_{0,n}[k] = [j^k] prod_{i=1}^{n-3} (J - i)                    (Route A Lemma 3.1)
#   g >= 1  : E_formulas/engine_h<g>_r<n>.txt (Route E engine) and E_formulas/referee/ref_engine_h<g>_r<n>.txt
#             (Route E referee's independent engine): X[k] = Sum_c (-1/sigma)^c [j^k] Xi^(c)
#   n = 1, no table: formula (6.1) (Hodge-side identity for h <= 7 by Route E's referee; conjectural for h >= 8)
# P2 prediction:  C^TR_{g,n}[k] = sigma^(1-g) X_{g,n}[k]  for all ordered k (|k| <= 3g-3+n) and nothing outside.
# Also reported: the TR degree bound at the zero w = 0 of x, i.e. [x^mu]F_{g,n} = Sum_k C[k] prod mu_i^k_i Pi_mu_i for all
# mu_i >= 1 with |mu| <= 2g-3+n (Theorem 3.1 at tau = 1; Lemma D / Lemma A at other tau), and the first layer |mu| = 2g-2+n.
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use TRxi;
sub R { Math::BigRat->new(@_) }
$| = 1;
my ($tau, $S, $GMAX, $save) = @ARGV; die "usage: r2_big.pl tau S GMAX [savefile]\n" unless defined $GMAX;
my $t0 = time;
my $kmax = 0; for my $g (0..$GMAX) { my $n = $S - $g; next if $n < 1; my $d = 3*$g - 3 + $n; $kmax = $d if $d > $kmax; }
$kmax = 1 if $kmax < 1;
TRxi::setup(tau => $tau, kmax => $kmax, y => 'base');
my $sig = $TRxi::sigma; my $tauR = R($tau);
printf "tau = %s, sigma = %s, S = %d, GMAX = %d, kmax = %d, series length %d; setup %.1f s\n", $tauR->bstr, $sig->bstr, $S, $GMAX, $kmax, $TRxi::L, time - $t0;
my @order = TRxi::run_all($S, $GMAX);
printf "recursion: %d correlators, %.1f s\n", scalar(@order), time - $t0;
if ($save) { open my $sv, '>', $save or die; print $sv "# tau = $tau; xi-coefficients of omega_{g,n} on C_tau, y = -1/w (TRxi.pm via referee2/r2_big.pl)\n";
  for my $gn (sort keys %TRxi::W) { my ($g, $n) = split /,/, $gn; my $Wx = $TRxi::W{$gn};
    print $sv "$g | $n | ", join(" ; ", map { "$_:" . $Wx->{$_}->bstr } sort keys %$Wx), "\n"; } close $sv; }

my $EDIR = "$FindBin::Bin/../../../E_formulas";
sub load_table { my ($f) = @_; return undef unless -e $f; open my $fh, '<', $f or die "open $f"; my %comp;
  while (my $line = <$fh>) { chomp $line; $line =~ s/\r//g; next unless $line =~ /\S/; my ($hh, $rr, $c, $body) = split /\s*\|\s*/, $line, 4;
    $comp{$c} //= {}; next unless defined $body && $body =~ /\S/;
    for my $t (split /\s*;\s*/, $body) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; $comp{$c}{$m} = R($v); } }
  close $fh; return \%comp; }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; return $r; }
sub tuples { my ($n, $maxsum) = @_; my @out; my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x (0..$left) { $rec->([@$pre, $x], $left - $x); } }; $rec->([], $maxsum); return @out; }
sub g0X { my ($n, $k) = @_; my $K = 0; $K += $_ for @$k; return R(0) if $n < 3 || $K > $n - 3;
  my @c = (R(1)); for my $i (1..$n-3) { my @d = map { R(0) } 0..@c; for my $e (0..$#c) { $d[$e+1] += $c[$e]; $d[$e] -= $i * $c[$e]; } @c = @d; }
  my $m = fact($K); $m /= fact($_) for @$k; return ($c[$K] // R(0)) * $m; }
# (6.1): sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lambda^{2h}] S(e1 l)^(1+j e2) S(e2 l)^(1+j e1)
my @BN = (R(1)); sub binomR { my ($n, $k) = @_; my $r = R(1); for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i; } return $r; }
sub bern { my $n = shift; while ($#BN < $n) { my $m = scalar @BN; my $s = R(0); $s += binomR($m + 1, $_) * $BN[$_] for 0..$m-1; push @BN, -$s / ($m + 1); } return $BN[$n]; }
sub pred61 { my ($h, $j) = @_; my $e1 = -$tauR; my $e2 = 1 + $tauR;
  my @c = (R(0)); for my $m (1..$h) { my $ls = -bern(2*$m) / (2*$m * fact(2*$m)); push @c, $ls * ((1 + $j*$e2) * $e1**(2*$m) + (1 + $j*$e1) * $e2**(2*$m)); }
  my @r = (R(1)); for my $nn (1..$h) { my $s = R(0); $s += $_ * $c[$_] * $r[$nn - $_] for 1..$nn; push @r, $s / $nn; }
  my $p = R(1); $p *= ($j - $_) for 1..2*$h-2; return $p * $r[$h]; }   # = sigma^h Xi_h(j)
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tauR + $t) / $t; } return $r; }
sub comps { my ($n, $tot) = @_; return () if $n == 0; return ([$tot]) if $n == 1; my @o; for my $a (1..$tot-$n+1) { push @o, [$a, @$_] for comps($n-1, $tot-$a); } return @o; }

my ($TOT, $BAD) = (0, 0);
for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %TRxi::W) {
  my ($g, $n) = split /,/, $gn; next if $n < 1; my $Wx = $TRxi::W{$gn}; my $dim = 3*$g - 3 + $n;
  # symmetry and support
  my ($symbad, $outside) = (0, 0); my %s;
  for my $key (keys %$Wx) { my @kk = split /,/, $key; my $sum = 0; $sum += $_ for @kk; $outside++ if $sum > $dim;
    my $sk = join(",", sort { $a <=> $b } @kk); push @{$s{$sk}}, $Wx->{$key}; }
  for my $sk (keys %s) { my @v = @{$s{$sk}}; for (@v) { $symbad++ if $_ != $v[0]; } }
  my @lines;
  my @sources;
  if ($g == 0) { push @sources, ['genus0', undef]; }
  else { for my $f (["engine", "$EDIR/engine_h${g}_r${n}.txt"], ["refengine", "$EDIR/referee/ref_engine_h${g}_r${n}.txt"]) {
      my $T = load_table($f->[1]); push @sources, [$f->[0], $T] if $T; } }
  my @K = tuples($n, $dim);
  for my $src (@sources) { my ($name, $T) = @$src; my ($cmp, $bad) = (0, 0); my @ex;
    for my $k (@K) { my $key = join(",", @$k); my $got = $Wx->{$key} // R(0); my $X;
      if ($name eq 'genus0') { $X = g0X($n, $k); }
      else { $X = R(0); my $ms = R(-1) / $sig; for my $c (keys %$T) { my $v = $T->{$c}{$key}; $X += $v * ($ms ** $c) if defined $v; } }
      my $exp = $sig ** (1 - $g) * $X; $cmp++; if ($got != $exp) { $bad++; push @ex, "k=($key) TR=$got pred=$exp" if @ex < 2; } }
    push @lines, sprintf("vs %s: %d compared, %d mismatches%s", $name, $cmp, $bad, @ex ? " [" . join(" | ", @ex) . "]" : "");
    $TOT += $cmp; $BAD += $bad; }
  if ($n == 1 && $g >= 2) { my ($cmp, $bad) = (0, 0); my @ex;
    for my $jj (0..3*$g+1) { my $j = R($jj); my $tr = R(0); for my $key (keys %$Wx) { $tr += $Wx->{$key} * $j ** $key; }
      my $lhs = $sig ** (2*$g - 1) * $tr;   # = sigma^g Xi_g(j) by P2
      my $p = pred61($g, $j); $cmp++; if ($lhs != $p) { $bad++; push @ex, "j=$jj" if @ex < 3; } }
    push @lines, sprintf("vs (6.1)%s: %d j-values, %d mismatches%s", ($g >= 8 ? " [conjectural at h>=8]" : ""), $cmp, $bad, @ex ? " [" . join(",", @ex) . "]" : "");
    if ($g <= 7) { $TOT += $cmp; $BAD += $bad; } }
  # degree bound at w = 0 (x-expansion), g >= 2
  if ($g >= 2) { my $zbad = 0; my $zc = 0; my $first_nz = 0; my $first_c = 0;
    my %PI; my $maxm = 2*$g - 2 + $n; $PI{$_} = Pi_mu($_) for 1..$maxm;
    for my $tot ($n..$maxm) { for my $mu (comps($n, $tot)) { my $v = R(0);
        for my $key (keys %$Wx) { my @kk = split /,/, $key; my $t = $Wx->{$key}; for my $i (0..$n-1) { $t *= R($mu->[$i]) ** $kk[$i]; } $v += $t; }
        for my $i (0..$n-1) { $v *= $PI{$mu->[$i]}; }
        if ($tot <= 2*$g - 3 + $n) { $zc++; $zbad++ unless $v->is_zero; } else { $first_c++; $first_nz++ unless $v->is_zero; } } }
    push @lines, sprintf("x-degree: %d coefficients with |mu| <= 2g-3+n, %d non-zero; layer |mu| = 2g-2+n: %d of %d non-zero", $zc, $zbad, $first_nz, $first_c); }
  printf "(%d,%d): %d xi-coefficients, symmetry violations %d, outside |k|<=dim %d; %s\n", $g, $n, scalar(keys %$Wx), $symbad, $outside, join("; ", @lines);
  $BAD += $symbad + $outside; }
printf "TOTAL: %d comparisons with Hodge-side data (h<=7), %d failures; time %.1f s\n", $TOT, $BAD, time - $t0;
