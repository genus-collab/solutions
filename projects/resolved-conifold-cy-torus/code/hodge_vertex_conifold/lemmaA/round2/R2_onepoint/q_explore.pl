#!/usr/bin/perl
# q_explore.pl h : sigma^h Q_h(s;(j)) from the data tables, rewritten in t = 1+s, grouped by sigma^k and t^m
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use mvlib; use MP;
my $h = shift // 2; my $r = 1;
my $file = "$FindBin::Bin/../../data/q_g${h}_r1.txt";
my $recs = parse_file($file); my $node = parse_expr($recs->[0][2]);
my @xs = map { Q($_) } (1..($h+1)); my (@sig, @vals);
for my $x (@xs) { my $p = eval_poly($node, ['j1', 's'], { x => $x }); my $c = 1 + $x; my $sg = -$x / ($c ** 2);
  my $f = ($c ** (3*$h - 3 + 2*$r)) * ($sg ** $h); my %q; for my $k (keys %$p) { $q{$k} = $p->{$k} * $f; } push @sig, $sg; push @vals, \%q; }
my %res; my $n = scalar @sig;
for my $i (0..$n-1) { my @l = (Q(1)); my $den = Q(1);
  for my $m (0..$n-1) { next if $m == $i; my @nl = map { Q(0) } 0..$#l+1; for my $d (0..$#l) { $nl[$d+1] += $l[$d]; $nl[$d] -= $l[$d] * $sig[$m]; } @l = @nl; $den *= ($sig[$i] - $sig[$m]); }
  for my $k (keys %{$vals[$i]}) { for my $d (0..$#l) { next if $l[$d]->is_zero; my $kk = "$d,$k"; my $v = $vals[$i]{$k} * $l[$d] / $den; $res{$kk} = (exists $res{$kk} ? $res{$kk} : Q(0)) + $v; } } }
for my $k (keys %res) { delete $res{$k} if $res{$k}->is_zero; }
# substitute s = t - 1 (variable index 2)
my $tm1 = mp_add(mp_var(2, 3), mp_const(-1, 3));
my $P = mp_subst1(\%res, 2, $tm1, 3);
for my $d (0..$h) { print "sigma^$d:\n";
  for my $m (0..3*$h) { my %t; for my $k (keys %$P) { my ($es, $ej, $et) = split /,/, $k; next unless $es == $d && $et == $m; $t{$ej} = $P->{$k}; }
    next unless %t; my @ks = sort { $a <=> $b } keys %t; print "   t^$m: ", join(' + ', map { "($t{$_})j^$_" } @ks), "\n"; } }
