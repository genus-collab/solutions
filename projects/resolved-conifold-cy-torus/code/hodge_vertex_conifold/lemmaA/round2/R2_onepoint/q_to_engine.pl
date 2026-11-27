#!/usr/bin/perl
# q_to_engine.pl h r : convert lemmaA/data/q_g<h>_r<r>.txt (c-normalisation, symbolic x, leaves j1..jr, s) into Route E's
# engine format (h | r | k | a1,..,ar:coef ; ...) for Xi_h(R) = Q_h(s = -1; R) (c = 1), Xi = Sum_k (-1/sigma)^k Xi^(k).
# Conversion as in q_extract.pl: Q^{c=1} = c^{3h-3+2r} Q^{c-norm}, c = 1+x, sigma = -x/c^2; sigma^h Q^{c=1} is interpolated in sigma
# (degree <= h, h+1 sample values of x), then evaluated at s = -1; Xi^(k) = (-1)^k [sigma^(h-k)] (sigma^h Xi).
# Output: qeng_h<h>_r<r>.txt in this folder; use it as the engine file of xy_npoint.pl (direct check of Theorem 1 against the data tables).
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use mvlib; use MP;
my ($h, $r) = @ARGV; $h //= 2; $r //= 1;
my $file = "$FindBin::Bin/../../data/q_g${h}_r${r}.txt";
my $recs = parse_file($file); my $node = parse_expr($recs->[0][2]);
my @V = ((map { "j$_" } 1..$r), 's');
my @xs = map { Q($_) } (1..($h+1));
my (@sig, @vals);
for my $x (@xs) { my $p = eval_poly($node, \@V, { x => $x }); my $c = 1 + $x; my $sg = -$x / ($c ** 2);
  my $f = ($c ** (3*$h - 3 + 2*$r)) * ($sg ** $h);
  my %q; for my $k (keys %$p) { my @e = split /,/, $k; my $es = pop @e; my $v = $p->{$k} * $f * (Q(-1) ** $es); my $kk = join(',', @e);
    $q{$kk} = exists $q{$kk} ? $q{$kk} + $v : $v; }
  push @sig, $sg; push @vals, \%q; }
my %res; my $n = scalar @sig;   # keys "d|a1,..,ar" -> coefficient of sigma^d j^a in sigma^h Xi
for my $i (0..$n-1) { my @l = (Q(1)); my $den = Q(1);
  for my $m (0..$n-1) { next if $m == $i; my @nl = map { Q(0) } 0..$#l+1;
    for my $d (0..$#l) { $nl[$d+1] += $l[$d]; $nl[$d] -= $l[$d] * $sig[$m]; } @l = @nl; $den *= ($sig[$i] - $sig[$m]); }
  for my $k (keys %{$vals[$i]}) { for my $d (0..$#l) { next if $l[$d]->is_zero; my $kk = "$d|$k"; my $v = $vals[$i]{$k} * $l[$d] / $den;
      $res{$kk} = exists $res{$kk} ? $res{$kk} + $v : $v; } } }
my $out = "$FindBin::Bin/qeng_h${h}_r${r}.txt"; open my $o, '>', $out or die;
for my $k (0..$h) { my $d = $h - $k; my @terms;
  for my $kk (sort keys %res) { my ($dd, $mono) = split /\|/, $kk; next unless $dd == $d; next if $res{$kk}->is_zero; push @terms, "$mono:" . ((Q(-1) ** $k) * $res{$kk}); }
  print $o "$h | $r | $k | ", join(' ; ', @terms), "\n"; }
close $o; print "wrote $out\n";
