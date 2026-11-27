#!/usr/bin/perl
# q_extract.pl h : read lemmaA/data/q_g<h>_r1.txt (c-normalisation, symbolic x, j1, s) and print
#   sigma^h Q_h(s;(j)) (c = 1 normalisation) as a polynomial in (sig, j, s).
# Conversion: Q^{c=1} = c^{3h-3+2r} Q^{c-norm}, c = 1+x, sigma = -x/c^2.  Interpolation in sigma (deg <= h).
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use mvlib; use MP;
my $h = shift // 2; my $r = 1;
my $file = "$FindBin::Bin/../../data/q_g${h}_r1.txt";
my $recs = parse_file($file); my $expr = $recs->[0][2];
my $node = parse_expr($expr);
my @xs = map { Q($_) } (1..($h+1));
my (@sig, @vals);
for my $x (@xs) {
  my $p = eval_poly($node, ['j1', 's'], { x => $x });          # hash "e_j,e_s" => BigRat
  my $c = 1 + $x; my $sg = -$x / ($c ** 2);
  my $f = ($c ** (3*$h - 3 + 2*$r)) * ($sg ** $h);
  my %q; for my $k (keys %$p) { $q{$k} = $p->{$k} * $f; }
  push @sig, $sg; push @vals, \%q; }
# Lagrange interpolation in sigma: result keys "e_sig,e_j,e_s"
my %res; my $n = scalar @sig;
for my $i (0..$n-1) {
  # basis polynomial l_i(sig) coefficients
  my @l = (Q(1)); my $den = Q(1);
  for my $m (0..$n-1) { next if $m == $i; my @nl = map { Q(0) } 0..$#l+1;
    for my $d (0..$#l) { $nl[$d+1] += $l[$d]; $nl[$d] -= $l[$d] * $sig[$m]; } @l = @nl; $den *= ($sig[$i] - $sig[$m]); }
  for my $k (keys %{$vals[$i]}) { for my $d (0..$#l) { next if $l[$d]->is_zero; my $kk = "$d,$k"; my $v = $vals[$i]{$k} * $l[$d] / $den;
      if (exists $res{$kk}) { $res{$kk} += $v; } else { $res{$kk} = $v; } } } }
for my $k (keys %res) { delete $res{$k} if $res{$k}->is_zero; }
print "h=$h  sigma^h Q_h(s;(j)) =\n", mp_str(\%res, ['sig','j','s']), "\n";
# value at s = -1 and s = 0, by sigma power
for my $sv (-1, 0) { print "s=$sv:\n"; for my $d (0..$h) { my %t; for my $k (keys %res) { my ($es, $ej, $ess) = split /,/, $k; next unless $es == $d;
      my $kk = "$ej"; my $v = $res{$k} * (Q($sv) ** $ess); $t{$kk} = (exists $t{$kk} ? $t{$kk} : Q(0)) + $v; }
    my @ks = sort { $a <=> $b } keys %t; print "  sig^$d: ", join(' + ', map { "($t{$_})*j^$_" } grep { !$t{$_}->is_zero } @ks), "\n"; } }
