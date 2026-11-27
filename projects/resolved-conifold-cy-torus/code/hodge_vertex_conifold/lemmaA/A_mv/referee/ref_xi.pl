#!/usr/bin/perl
# ref_xi.pl -- referee: Xi_h(B) directly from the Marino-Vafa formula + Laguerre transform (rmv.pm, mod p), no KP, no TR.
# Usage:
#   perl ref_xi.pl check h prime seed ntau "B1;B2;.."   -> compare with lemmaA/data (mv_g*_l*, q_g*_r*) evaluated mod p
#   perl ref_xi.pl run   h prime seed ntau "B1;B2;.."   -> print Xi_h(B) mod p at ntau random tau, and sigma^h Xi_h(B) interpolated
# B = comma-separated parts, e.g. "2;3;4;2,2;3,2".
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rmv; use rparse;
my ($mode, $h, $prime, $seed, $ntau, $Blist) = @ARGV;
rmv::setp($prime); my $p = $prime; srand($seed);
my @Bs = map { [ sort { $b <=> $a } split /,/, $_ ] } split /;/, $Blist;
my $data = "$FindBin::Bin/../../data";
my %tab;
sub tabexpr { my ($kind, $g, $l) = @_; my $f = "$data/${kind}_g${g}_" . ($kind eq 'mv' ? 'l' : 'r') . "$l.txt"; return undef unless -e $f;
  $tab{$f} //= do { my ($row) = read_table($f); $row->[2]; }; return $tab{$f}; }
my $t0 = time; my ($nok, $nbad) = (0, 0); my %vals;
for my $it (1..$ntau) {
  my $taup = 2 + int(rand($p - 5)); my $x = rmv::md(-$taup * rmv::mi($taup + 1)); my $c = rmv::mi($taup + 1); my $sig = ($taup * ($taup + 1)) % $p;
  my $res = rmv::compute($h, \@Bs, $taup, gmax => $h);
  for my $B (@Bs) { my $key = join(',', @$B); my $rr = $res->{$key}; my $r = scalar @$B;
    push @{$vals{$key}}, [ $sig, $rr->{Xi} ];
    if ($mode eq 'check') {
      # (a) M_g(B + 1^k) against mv tables (l = r + k <= 3), g = 1..h (data exist for g <= 4)
      for my $g (1..$h) { for my $k (0..3 - $r) { my $l = $r + $k; my $e = tabexpr('mv', $g, $l); next unless defined $e; my $Mv = $rr->{M}{$g}[$k]; next unless defined $Mv;
          my %v = (x => $x); my @m = (@$B, (1) x $k); $v{"m$_"} = $m[$_-1] for 1..$l;
          my $d = rparse::eval_expr_mod($e, \%v, $p); my $cf = 1; $cf = ($cf * $c) % $p for 1..(3*$g - 3 + 2*$l); my $pred = ($Mv * rmv::mi($cf)) % $p;
          if ($pred == $d) { $nok++; } else { $nbad++; print "MISMATCH mv g=$g mu=(@m) tau=$taup: $pred vs $d\n"; } } }
      # (b) Xi_h(B) and Q_h(s;B) (s = 2, 5) against q tables
      my $e = tabexpr('q', $h, $r);
      if (defined $e) { my $cf = 1; $cf = ($cf * $c) % $p for 1..(3*$h - 3 + 2*$r);
        for my $s0 (-1, 2, 5) { my %v = (x => $x, s => rmv::md($s0)); $v{"j$_"} = $B->[$_-1] for 1..$r; my $d = rparse::eval_expr_mod($e, \%v, $p);
          my $q = 0; my $sp = 1; for my $m (0..$rr->{N}) { $q = ($q + $rr->{Q}[$m] * $sp) % $p; $sp = ($sp * rmv::md($s0)) % $p; }
          my $pred = ($q * rmv::mi($cf)) % $p; my $pxi = ($rr->{Xi} * rmv::mi($cf)) % $p;
          if ($pred == $d && ($s0 != -1 || $pxi == $d)) { $nok++; } else { $nbad++; print "MISMATCH q h=$h B=($key) s=$s0 tau=$taup: $pred (Xi $pxi) vs $d\n"; } } }
    } } }
if ($mode eq 'check') { printf "CHECK h=%d p=%d: %d matches, %d mismatches (%d s)\n", $h, $p, $nok, $nbad, time - $t0; }
else {
  for my $B (@Bs) { my $key = join(',', @$B); my $L = $vals{$key}; my $nz = grep { $_->[1] } @$L;
    my $J = 0; $J += $_ for @$B; my $D = $J - scalar(@$B);
    printf "h=%d B=(%s) D=%d %s: Xi mod p non-zero at %d of %d random tau%s\n", $h, $key, $D, ($D <= 2*$h - 3 ? "[in range]" : "[outside]"), $nz, scalar @$L,
      ($nz ? "  e.g. Xi=" . (grep { $_->[1] } @$L)[0][1] : ""); }
  printf "(p = %d, %d s)\n", $p, time - $t0; }
