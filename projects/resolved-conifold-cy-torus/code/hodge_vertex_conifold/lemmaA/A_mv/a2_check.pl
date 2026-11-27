#!/usr/bin/perl
# a2_check.pl -- (A2): Q_h(s;R) and Xi_h(R) from Marino-Vafa integrals with extra unit legs.
#   M_a := <Omega lam_R lam_1^a>_{h, r+a}  (c = 1),   computed from the Schur side (mv_unit),
#   [s^m] Q_h(s;R) = Sum_{a=0}^m (M_a/a!) (-1)^{m-a} L_{m-a}^{(2h-3+r+a)}(-(J+a)),
#   Xi_h(R) = Q_h(-1;R) = Sum_{a=0}^{N} (-1)^a/a! L_{N-a}^{(2h-2+r+a)}(-J-a) M_a,   N = 3h-3+r.
# Checks: [s^m]Q = 0 for m = N+1, N+2 (degree bound); [s^m]Q against lemmaA/data/q_g<h>_r<r>.txt (c-normalisation:
# Q^(c) = c^-(3h-3+2r) Q^(c=1)); Xi against formula (6.1) (r = 1) and hvc_multiprobe_out.txt (r = 2, genus 2, 3).
# Usage: perl a2_check.pl h r "j1,j2;j1,j2;..."  [taus]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($h, $r, $Rlist, $taulist) = @ARGV;
my @Rs = map { [ split /,/, $_ ] } split /;/, $Rlist;
my @taus = map { R($_) } split /,/, ($taulist // '2/7,3/2,5,-7/3,7/11');
my $N = 3*$h - 3 + $r;
# q-data
my $qfile = "$FindBin::Bin/../data/q_g${h}_r${r}.txt"; my $qast;
if (-e $qfile) { my $rec = parse_file($qfile); $qast = parse_expr($rec->[0][2]); print "q-data: $qfile\n"; } else { print "q-data: $qfile not present (skipped)\n"; }
# multiprobe symbolic r = 2
my $mpast;
if ($r == 2) { open my $fh, '<', "$FindBin::Bin/../../hvc_multiprobe_out.txt" or die; while (<$fh>) { if (/^genus $h symbolic r=2: (.*)$/) { $mpast = parse_expr($1); } } close $fh;
  print "multiprobe genus $h r = 2: ", ($mpast ? "found" : "absent"), "\n"; }
sub f61 { my ($tau, $j, $h) = @_; my $e1 = -$tau; my $e2 = 1 + $tau; my $K = $h;
  my @c = (R(0)); for my $k (1..$K) { my $lS = -bern(2*$k)/(2*$k*fact(2*$k)); push @c, (1 + $j*$e2)*$lS*($e1**(2*$k)) + (1 + $j*$e1)*$lS*($e2**(2*$k)); }
  my $f = ser_exp(\@c, $K); my $p = R(1); $p *= ($j - $_) for 1..(2*$h-2); my $sig = $tau*(1+$tau); return $p * $f->[$h] / $sig**$h; }
my ($nchk, $nbad) = (0, 0);
for my $tau (@taus) { my $x = -$tau/(1+$tau); my $c = 1 + $x; my $sig = $tau*(1+$tau); my $cf = $c ** (3*$h - 3 + 2*$r);
  my $qpoly = $qast ? eval_poly($qast, ['s', map { "j$_" } 1..$r], { x => $x }) : undef;
  my $mppoly = $mpast ? eval_poly($mpast, ['jj1', 'jj2'], { x => $x }) : undef;
  for my $Rp (@Rs) { my $J = 0; $J += $_ for @$Rp; my $D = $J - $r;
    my $t0 = time; my $M = mv_unit($Rp, $tau, $h, $N + 2); my $dt = time - $t0;
    my @Ma = map { $M->[$h][$_] } 0..$N+2;
    my @q; for my $m (0..$N+2) { my $s = R(0); for my $a (0..$m) { $s += $Ma[$a] / fact($a) * R((-1)**($m-$a)) * lag($m-$a, 2*$h-3+$r+$a, -($J+$a)); } $q[$m] = $s; }
    my $xi = R(0); $xi += R((-1)**$_) * $q[$_] for 0..$N;
    my $xi2 = R(0); for my $a (0..$N) { $xi2 += R((-1)**$a) / fact($a) * lag($N-$a, 2*$h-2+$r+$a, -$J-$a) * $Ma[$a]; }
    my @msg;
    push @msg, "deg-bound " . (($q[$N+1]->is_zero && $q[$N+2]->is_zero) ? "ok" : "FAIL($q[$N+1],$q[$N+2])");
    $nchk++; $nbad++ unless $q[$N+1]->is_zero && $q[$N+2]->is_zero;
    push @msg, "laguerre " . ($xi == $xi2 ? "ok" : "FAIL"); $nchk++; $nbad++ unless $xi == $xi2;
    if ($qpoly) { my @dq = map { R(0) } 0..$N+3;
      for my $k (keys %$qpoly) { my @e = split /,/, $k; my $t = $qpoly->{$k}->copy; $t *= R($Rp->[$_-1]) ** $e[$_] for 1..$r; $dq[$e[0]] += $t; }
      my $ok = 1; for my $m (0..$N+3) { $ok = 0 unless $dq[$m] * $cf == ($m <= $N+2 ? $q[$m] : 0); }
      push @msg, "q-data " . ($ok ? "ok" : "FAIL"); $nchk++; $nbad++ unless $ok;
      if (!$ok) { for my $m (0..$N) { push @msg, " m=$m mv=$q[$m] data=" . ($dq[$m]*$cf); } } }
    if ($r == 1) { my $f = f61($tau, R($Rp->[0]), $h); push @msg, "(6.1) " . ($f == $xi ? "ok" : "FAIL($f)"); $nchk++; $nbad++ unless $f == $xi; }
    if ($mppoly) { my $v = poly_eval_at($mppoly, [ R($Rp->[0]), R($Rp->[1]) ]) * $cf; push @msg, "multiprobe " . ($v == $xi ? "ok" : "FAIL($v)"); $nchk++; $nbad++ unless $v == $xi; }
    printf "h=%d tau=%s R=(%s) D=%d: Xi=%s  [%s] (%d s)\n", $h, $tau, join(',', @$Rp), $D, $xi, join('; ', @msg), $dt;
    print "   q_m (c=1): ", join(', ', map { "$_:" . $q[$_] } 0..$N), "\n" if $ENV{SHOWQ}; } }
print "SUMMARY h=$h r=$r: $nchk checks, $nbad failures\n";
