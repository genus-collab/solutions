#!/usr/bin/perl
# Genus-7 test of the leading-order mKP relation M_{7,eps}(R), R = (13), (14), modulo a prime p at random tau:
#   genus-7 values Xi_7(R), Xi_7(R+2) from the Marino-Vafa Schur side + Laguerre transform (referee's rmv.pm, mod p; no Hodge data),
#   genus <= 6 values from the exact Route E tables (XiData, reduced mod p), genus 0 closed form.
# Usage: perl g7_mkp_modp.pl prime seed
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use lib "$FindBin::Bin/../../A_mv/referee"; use Math::BigRat; use rmv; use XiData;
$| = 1;
my ($prime, $seed, $GG, $tlist) = @ARGV; $prime //= 1000000007; $seed //= 7; $GG //= 7; $tlist //= "13;14"; rmv::setp($prime); my $p = $prime; srand($seed);
my $taup = 2 + int(rand($p - 5)); my $sig = ($taup * (($taup + 1) % $p)) % $p;
my $t0 = time;
my @targets = map { [split /,/, $_] } split /;/, $tlist;
my @Bs = map { ([sort { $b <=> $a } @$_], [sort { $b <=> $a } (@$_, 2)]) } @targets;
my %XH;
for my $gg (7..$GG) { my $res = rmv::compute($gg, \@Bs, $taup, gmax => $gg);
  printf "# rmv genus %d done (%d s), tau = %d mod %d\n", $gg, time - $t0, $taup, $p;
  for my $B (@Bs) { my $k = join(',', sort { $a <=> $b } @$B); $XH{"$gg|$k"} = $res->{join(',', @$B)}{Xi}; print "# Xi_$gg($k) = $XH{\"$gg|$k\"} mod p\n"; } }
sub md { rmv::md(@_) } sub mi { rmv::mi(@_) }
sub r2m { rmv::rat2mod("$_[0]") }
my $vv = md(-mi($sig));
sub X { my ($g, $T) = @_; my @T = sort { $a <=> $b } @$T; my $k = join(',', @T);
  if ($g >= 7) { die "no genus-$g value for $k" unless exists $XH{"$g|$k"}; return $XH{"$g|$k"}; }
  if ($g == 0) { return r2m(XiData::xi0(\@T)); }
  my $v = XiData::xi_v($g, \@T); my $s = 0; my $pw = 1; for my $c (@$v) { $s = ($s + r2m($c) * $pw) % $p; $pw = ($pw * $vv) % $p; } return $s; }
sub fct { my $n = shift; my $f = 1; $f = ($f * $_) % $p for 2..$n; $f }
sub ff { my ($x, $k) = @_; my $q = 1; for my $i (0..$k-1) { $q = ($q * md($x - $i)) % $p; } $q }
sub Pi { rmv::Pi(@_) }
sub pw { my ($b, $e) = @_; my $r = 1; $r = ($r * $b) % $p for 1..$e; $r }
sub Phi { my ($g, $k, $T) = @_;
  if (!@$T) { return 0 if $g > 0; return $taup if $k == 1; return ($taup * $taup % $p) * mi(2) % $p if $k == 2; return 0; }
  my $J = 0; $J += $_ for @$T; my $D = $J - scalar(@$T); my $f = ff($D - 2*$g + 2, $k); return 0 unless $f;
  return pw($sig, $g) * pw($taup, $k) % $p * $J % $p * $f % $p * X($g, $T) % $p * mi(fct($k)) % $p; }
sub rel { my ($G, $Rr, $par) = @_; my @Rr = @$Rr; my $J = 0; $J += $_ for @Rr; my $D = $J - scalar(@Rr); my $n = 2*$G + $par; my @R2 = (@Rr, 2); my $D2 = $D + 1;
  my $P2 = Pi(2, $taup); my $t = 0;
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1; my $f = ff($D2 - 2*$g + 2, $k); next unless $f;
    $t = md($t - 2 * $P2 % $p * pw($sig, $g) % $p * pw($taup, $k) % $p * $f % $p * X($g, \@R2) % $p * mi(fct($k))); }
  $t = ($t + 2 * pw($sig, $G+1) % $p * $J % $p * ($J+1) % $p * X($G, \@Rr)) % $p if $par == 0;
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1; my $f = ff($D - 2*$g + 2, $k); next unless $f;
    $t = ($t + pw($sig, $g+1) * pw($taup, $k) % $p * $J % $p * ($J+1) % $p * $f % $p * X($g, \@Rr) % $p * mi(fct($k))) % $p; }
  for my $s ([[@Rr], []], [[], [@Rr]]) { my ($T1, $T2) = @$s;
    for my $g1 (0..int(($n+1)/2)) { for my $k1 (1..$n+1-2*$g1) { my $rest = $n + 2 - 2*$g1 - $k1; my $f1 = Phi($g1, $k1, $T1); next unless $f1;
        for my $g2 (0..int(($rest-1)/2)) { my $k2 = $rest - 2*$g2; next if $k2 < 1; my $f2 = Phi($g2, $k2, $T2); next unless $f2; $t = ($t + $f1 * $f2) % $p; } } } }
  return $t; }
for my $R (@targets) { for my $par (0, 1) { my $v = rel($GG, $R, $par);
    printf "G=$GG R=(%s) %s: residual %d mod p %s\n", join(',', @$R), ($par ? "ODD " : "EVEN"), $v, ($v == 0 ? "OK" : "FAIL"); } }
# control: the same with a perturbed genus-7 value must fail
my $ck = join(",", sort { $a <=> $b } (@{$targets[0]}, 2)); $XH{"$GG|$ck"} = ($XH{"$GG|$ck"} + 1) % $p; my $vc = rel($GG, $targets[0], 0); printf "control (Xi_G(R+2)+1): residual %d %s\n", $vc, ($vc ? "(nonzero, as it must be)" : "ZERO?!");
printf "# total %d s\n", time - $t0;
