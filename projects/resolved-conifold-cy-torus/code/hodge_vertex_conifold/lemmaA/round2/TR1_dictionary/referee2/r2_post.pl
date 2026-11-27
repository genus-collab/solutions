#!/usr/bin/perl
# r2_post.pl TAU SAVEFILE [MAXKEYS]   (referee 2) : fast post-processing of a save file written by r2_big.pl.
# For every saved (g,n), n >= 1: symmetry; comparison of every xi-coefficient with sigma^(1-g) X_{g,n}[k] from the Route E engine
# tables and the Route E referee tables (genus 0: [j^k] prod (J-i)); n = 1: formula (6.1); and, for correlators with at most
# MAXKEYS (default 3000) stored coefficients, the x-degree statement [x^mu]F = 0 for |mu| <= 2g-3+n (evaluated by contracting
# one variable at a time) and the number of non-zero coefficients in the first layer |mu| = 2g-2+n.
use strict; use warnings; use FindBin; use Math::BigRat; use Time::HiRes qw(time);
sub R { Math::BigRat->new(@_) }
$| = 1; my ($tau, $file, $MAXKEYS) = @ARGV; $MAXKEYS //= 100000; my $t0 = time;
my $tauR = R($tau); my $sig = $tauR * (1 + $tauR);
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
my @BN = (R(1)); sub binomR { my ($n, $k) = @_; my $r = R(1); for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i; } return $r; }
sub bern { my $n = shift; while ($#BN < $n) { my $m = scalar @BN; my $s = R(0); $s += binomR($m + 1, $_) * $BN[$_] for 0..$m-1; push @BN, -$s / ($m + 1); } return $BN[$n]; }
sub pred61 { my ($h, $j) = @_; my $e1 = -$tauR; my $e2 = 1 + $tauR;
  my @c = (R(0)); for my $m (1..$h) { my $ls = -bern(2*$m) / (2*$m * fact(2*$m)); push @c, $ls * ((1 + $j*$e2) * $e1**(2*$m) + (1 + $j*$e1) * $e2**(2*$m)); }
  my @r = (R(1)); for my $nn (1..$h) { my $s = R(0); $s += $_ * $c[$_] * $r[$nn - $_] for 1..$nn; push @r, $s / $nn; }
  my $p = R(1); $p *= ($j - $_) for 1..2*$h-2; return $p * $r[$h]; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tauR + $t) / $t; } return $r; }
open my $fh, '<', $file or die "no $file"; my %W;
while (my $line = <$fh>) { next if $line =~ /^#/; chomp $line; $line =~ s/\r//g; my ($g, $n, $body) = split /\s*\|\s*/, $line, 3; next unless defined $body;
  my %c; for my $t (split /\s*;\s*/, $body) { next unless $t =~ /:/; my ($k, $v) = split /:/, $t; $c{$k} = R($v); } $W{"$g,$n"} = \%c; }
close $fh;
my ($TOT, $BAD) = (0, 0);
for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %W) {
  my ($g, $n) = split /,/, $gn; next if $n < 1; my $Wx = $W{$gn}; my $dim = 3*$g - 3 + $n; my @lines;
  my ($symbad, $outside) = (0, 0); my %s;
  for my $key (keys %$Wx) { my @kk = split /,/, $key; my $sum = 0; $sum += $_ for @kk; $outside++ if $sum > $dim; push @{$s{join(",", sort { $a <=> $b } @kk)}}, $Wx->{$key}; }
  for my $sk (keys %s) { my @v = @{$s{$sk}}; for (@v) { $symbad++ if $_ != $v[0]; } }
  my @src; if ($g == 0) { push @src, ['genus0', undef]; }
  else { for my $f (["engine", "$EDIR/engine_h${g}_r${n}.txt"], ["refengine", "$EDIR/referee/ref_engine_h${g}_r${n}.txt"]) { my $T = load_table($f->[1]); push @src, [$f->[0], $T] if $T; } }
  for my $src (@src) { my ($name, $T) = @$src; my ($cmp, $bad) = (0, 0);
    for my $k (tuples($n, $dim)) { my $key = join(",", @$k); my $X;
      if ($name eq 'genus0') { $X = g0X($n, $k); } else { $X = R(0); my $ms = R(-1) / $sig; for my $c (keys %$T) { my $v = $T->{$c}{$key}; $X += $v * $ms ** $c if defined $v; } }
      $cmp++; my $got = $Wx->{$key} // R(0); $bad++ if $got != $sig ** (1 - $g) * $X; }
    push @lines, "vs $name: $cmp compared, $bad mismatches"; $TOT += $cmp; $BAD += $bad; }
  if ($n == 1 && $g >= 2) { my ($cmp, $bad) = (0, 0); my $zeros = '';
    for my $jj (0..3*$g+1) { my $j = R($jj); my $tr = R(0); for my $key (keys %$Wx) { $tr += $Wx->{$key} * $j ** $key; }
      $zeros .= " $jj" if $tr->is_zero && $jj >= 1; my $lhs = $sig ** (2*$g - 1) * $tr; $cmp++; $bad++ if $lhs != pred61($g, $j); }
    push @lines, sprintf("vs (6.1)%s: %d j-values, %d mismatches; TR zeros at j =%s", ($g >= 8 ? " [h>=8: (6.1) conjectural]" : ""), $cmp, $bad, $zeros);
    if ($g <= 7) { $TOT += $cmp; $BAD += $bad; } }
  if ($g >= 2 && scalar(keys %$Wx) <= $MAXKEYS) {
    # contract variables: F(mu_1..mu_n) = Sum_k C[k] prod mu_i^k_i Pi_mu_i, for mu_i in 1..M, M = 2g-2+n - (n-1)
    my $M = 2*$g - 2 + 1; my %cur = %$Wx; my %PI = map { $_ => Pi_mu($_) } 1..$M;
    for my $var (reverse 0..$n-1) { my %nxt;
      for my $key (keys %cur) { my @kk = split /,/, $key; my $e = $kk[$var]; my @pre = @kk[0..$var-1]; my @post = @kk[$var+1..$#kk];
        for my $m (1..$M) { my $v = $cur{$key} * R($m) ** $e * $PI{$m}; next if $v->is_zero; my $nk = join(",", @pre, "m$m", @post);
          if (exists $nxt{$nk}) { $nxt{$nk} += $v } else { $nxt{$nk} = $v } } }
      %cur = %nxt; }
    my ($zc, $zbad, $lc, $lnz) = (0, 0, 0, 0);
    my @mus = tuples($n, 2*$g - 2);   # mu_i - 1 >= 0 with Sum (mu_i - 1) <= 2g-2
    for my $t (@mus) { my @mu = map { $_ + 1 } @$t; my $tot = 0; $tot += $_ for @mu; next if $tot > 2*$g - 2 + $n;
      my $v = $cur{join(",", map { "m$_" } @mu)} // R(0);
      if ($tot <= 2*$g - 3 + $n) { $zc++; $zbad++ unless $v->is_zero; } else { $lc++; $lnz++ unless $v->is_zero; } }
    push @lines, sprintf("x-degree: %d coefficients with |mu| <= 2g-3+n, %d non-zero; layer |mu| = 2g-2+n: %d of %d non-zero", $zc, $zbad, $lnz, $lc); }
  printf "(%d,%d): %d coefficients, symmetry violations %d, outside %d; %s\n", $g, $n, scalar(keys %$Wx), $symbad, $outside, join("; ", @lines);
  $BAD += $symbad + $outside; }
printf "TOTAL: %d comparisons with Hodge-side data (h <= 7), %d failures; time %.1f s\n", $TOT, $BAD, time - $t0;
