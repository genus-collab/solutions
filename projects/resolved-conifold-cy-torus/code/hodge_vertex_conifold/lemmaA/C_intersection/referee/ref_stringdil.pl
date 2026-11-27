#!/usr/bin/perl
# Referee: component-wise string and dilaton relations (used in Prop. L and in the C3(f) basis), checked as polynomial identities
#   X_k^{r=2}(1, j) = j * X_k^{r=1}(j),     X_k^{r=2}(0, j) = (j - 1 - 2h + 2) * X_k^{r=1}(j),   k = 0..h, h = 2, 3, 4,
# and the lambda_g component X_h^{r}(j) = b_h (J-1)...(J-N), N = 2h-3+r, b_h = [z^2h] (z/2)/sin(z/2).
# Components from the referee's own extraction out/ucomp_g<h>_r<r>.txt.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rpoly; use Math::BigRat; use Math::BigInt;
sub readc { my ($h, $r) = @_; open(my $fh, '<', "$FindBin::Bin/out/ucomp_g${h}_r${r}.txt") or die; my @X; while (<$fh>) { s/\r//; chomp; my @f = split /\s*\|\s*/; $X[$f[2]] = parse_rp($f[3]); } close $fh; return \@X; }
# substitute j1 = value in a 2-variable rp and rename j2 -> j1
sub sub1 { my ($p, $v) = @_; my %r; for my $k (keys %$p) { my ($a, $b, $c) = split /,/, $k; my $t = $p->{$k} * (Math::BigRat->new($v) ** $a); my $nk = join(',', $b, 0, 0); $r{$nk} = exists $r{$nk} ? $r{$nk} + $t : $t; } for (keys %r) { delete $r{$_} if $r{$_}->is_zero; } return \%r; }
my $j = { '1,0,0' => Math::BigRat->new(1) };
# b_h from (z/2)/sin(z/2) = 1/(sum_{n} (-1)^n (z/2)^{2n}/(2n+1)!)
sub bh { my $h = shift; my @s = map { Math::BigRat->new((-1) ** $_) / (Math::BigRat->new(4) ** $_ * Math::BigRat->new(Math::BigInt->new(2 * $_ + 1)->bfac->bstr)) } 0 .. $h;
  my @inv = (Math::BigRat->new(1)); for my $n (1 .. $h) { my $t = Math::BigRat->new(0); $t -= $s[$_] * $inv[$n - $_] for 1 .. $n; $inv[$n] = $t; } return $inv[$h]; }
my $bad = 0;
for my $h (2, 3, 4) { my $X1 = readc($h, 1); my $X2 = readc($h, 2);
  for my $k (0 .. $h) {
    my $s = rp_add(sub1($X2->[$k], 1), rp_scale(rp_mul($j, $X1->[$k]), -1));
    my $d = rp_add(sub1($X2->[$k], 0), rp_scale(rp_mul(rp_add($j, { '0,0,0' => Math::BigRat->new(1 - 2 * $h) }), $X1->[$k]), -1));
    printf "h=%d k=%d: string %s, dilaton %s\n", $h, $k, (rp_iszero($s) ? 'OK' : 'FAIL'), (rp_iszero($d) ? 'OK' : 'FAIL'); $bad++ unless (rp_iszero($s) && rp_iszero($d)); }
  for my $r (1, 2) { my $X = $r == 1 ? $X1 : $X2; my $N = 2 * $h - 3 + $r; my $ok = 1;
    for my $pt ([1, 1], [2, 5], [3, 4], [7, 2], [0, 9]) { my @v = $r == 1 ? ($pt->[0] + $pt->[1]) : @$pt; my $J = 0; $J += $_ for @v; my $pred = bh($h); $pred *= ($J - $_) for 1 .. $N; $ok = 0 unless rp_eval($X->[$h], @v) == $pred; }
    printf "h=%d r=%d: lambda_g component = b_h (J-1)...(J-N), b_h = %s: %s\n", $h, $r, bh($h)->bstr, ($ok ? 'OK' : 'FAIL'); $bad++ unless $ok; } }
print "SUMMARY: ", ($bad ? "$bad failures" : "all identities hold"), "\n";
