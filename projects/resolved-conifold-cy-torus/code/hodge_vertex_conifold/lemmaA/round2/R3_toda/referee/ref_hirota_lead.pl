#!/usr/bin/perl
# Referee re-implementation (from scratch) of the dilaton-point leading order of a Hirota equation P(D) tau_0.tau_1 = 0,
# tau_0 = Z_MV, tau_1 = Z_MV(e^{tau x} o p), as an exact POLYNOMIAL IN TAU (so every checked case is a proof of that case).
#   P(D) tau0.tau1/(tau0 tau1) = P(d_s) e^{Phi}|_0 = Sum_{set partitions} prod_B d^B Phi,  d^B Phi = d^B S (|B| even), -d^B A (|B| odd),
#   S = 2F + A, A = Sum_{k>=1} (tau x)^k E^k F/k!,  d_{t_i} = i d_{p_i} = i x sigma Pi_i d_{q_i},  F = Sum_g sigma^{g-1} x^{2g-2} F_g.
# Block leading coefficient (U = T_B + {i in B, i >= 2}, v = #{i in B: i = 1}):
#   sign * prod_{i in B} i * [2 if even block, k = 0] * prod_{i in B, i>=2} Pi_i * sigma^{|B|-1} tau^k/k! * Lead(g,k,U,v),
#   Lead = (J_U)^(v) (delta_g(U))_k sigma^g X_g(U)  (U nonempty),   Lead(empty) computed HERE from the closed forms of
#   G^(0)_0 = u - 3u^2/4 + u^3/6 and D G^(1)_0 = (w-1)(w-1-1/sigma)/24 by applying D and d_y (no table is assumed); 0 for g >= 2.
# Data: X_g(T) from Route E's exact tables E_formulas/engine_h<g>_r<r>.txt (own parser) or Route E's point engine; genus 0 closed form.
# Modes:  data  -> print the leading-order relation for every (G, R, parity) and test it is the zero polynomial;
#         random-> replace every sigma^g X_g(T) by a random rational and compare the general expression for P = D_1^2 + D_2 with the
#                  report's Theorem M formula (transcribed separately below): they must agree identically.
# Usage: perl ref_hirota_lead.pl MODE EQ Gmin Gmax rmax Jmax      EQ in {M (=D1^2+D2), P2, P3, P4, P5, JM2}
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use lib "$FindBin::Bin/../../../E_formulas"; use Math::BigRat; use ref_lib;
use Hodge; use HodgePt;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($MODE, $EQ, $Gmin, $Gmax, $rmax, $Jmax) = @ARGV; $MODE //= 'data'; $EQ //= 'M';
my $EDIR = "$FindBin::Bin/../../../E_formulas";
my %RTAB = (1 => 3, 2 => 4, 3 => 4, 4 => 5, 5 => 5, 6 => 3);   # use tables only up to these r (bigger files skipped)
$RTAB{$_} = $ENV{"RTAB$_"} for grep { defined $ENV{"RTAB$_"} } 1..9;
my %EXTRA;   # extra exact values "g|parts" -> [components], e.g. from g7_points_A.out
if (defined $ENV{EXTRA}) { for my $f (split /;/, $ENV{EXTRA}) { open(my $fh, '<', $f) or die "$f"; while (my $l = <$fh>) { $l =~ s/\r//g; chomp $l; next if $l =~ /^#/;
      my @p = split /\s*\|\s*/, $l; next unless @p == 3; my @c = map { R($_) } split /\s*;\s*/, $p[2]; $EXTRA{"$p[0]|" . join(',', sort { $a <=> $b } split /,/, $p[1])} = \@c; } close $fh; } }
# ---------- polynomials in tau ----------
sub ptrim { my $p = shift; pop @$p while @$p && $p->[-1] == 0; $p }
sub padd { my ($a, $b, $c) = @_; $c //= R(1); my @r; my $n = @$a > @$b ? @$a : @$b; for my $i (0..$n-1) { $r[$i] = ($a->[$i] // R(0)) + $c * ($b->[$i] // R(0)); } ptrim(\@r) }
sub pmul { my ($a, $b) = @_; return [] unless @$a && @$b; my @r = map { R(0) } 0..(@$a + @$b - 2);
  for my $i (0..$#$a) { next if $a->[$i] == 0; for my $j (0..$#$b) { $r[$i+$j] += $a->[$i] * $b->[$j]; } } ptrim(\@r) }
sub pscal { my ($a, $c) = @_; ptrim([ map { $_ * $c } @$a ]) }
sub pconst { my $c = R(shift); $c == 0 ? [] : [$c] }
sub pstr { my $p = shift; return "0" unless @$p; join(' + ', map { "($p->[$_])t^$_" } grep { $p->[$_] != 0 } 0..$#$p) }
my $TAU = [R(0), R(1)]; my $SIG = [R(0), R(1), R(1)];
sub ppow { my ($a, $e) = @_; my $r = [R(1)]; $r = pmul($r, $a) for 1..$e; $r }
sub fact { ref_lib::fact(@_) }
sub ff { my ($x, $k) = @_; my $p = R(1); $p *= ($x - $_) for 0..$k-1; $p }
sub rf { my ($x, $k) = @_; my $p = R(1); $p *= ($x + $_) for 0..$k-1; $p }
my %PIc; sub Pi { my $k = shift; return $PIc{$k} //= do { my $p = [R(1)]; for my $a (1..$k-1) { $p = pmul($p, [R($a), R($k)]); } pscal($p, R(1) / fact($k-1)) } }
# ---------- data ----------
my %tab;
sub load { my ($g, $r) = @_; my $key = "$g|$r"; return $tab{$key} if exists $tab{$key};
  my $f = "$EDIR/engine_h${g}_r${r}.txt"; if (!-e $f || $r > ($RTAB{$g} // 0)) { return $tab{$key} = 0; }
  open(my $fh, '<', $f) or die $f; my @rows;
  while (my $l = <$fh>) { $l =~ s/\r//g; chomp $l; my @p = split /\s*\|\s*/, $l, 4; next unless @p == 4;
    for my $t (split /\s*;\s*/, $p[3]) { next unless $t =~ /:/; my ($e, $c) = split /:/, $t; push @rows, [$p[2], [split /,/, $e], R($c)]; } }
  close $fh; return $tab{$key} = \@rows; }
my %xc; my $eng_g = -1; my ($nt, $ne) = (0, 0);
sub comps { my ($g, $T) = @_; my @T = sort { $a <=> $b } @$T; my $key = "$g|" . join(',', @T); return $xc{$key} if exists $xc{$key};
  return $xc{$key} = $EXTRA{$key} if exists $EXTRA{$key};
  my $r = scalar @T; my $tb = load($g, $r); my @c;
  if ($tb) { for my $row (@$tb) { my ($k, $e, $c) = @$row; my $m = $c->copy; for my $i (0..$r-1) { $m *= R($T[$i]) ** $e->[$i] if $e->[$i]; } $c[$k] = ($c[$k] // R(0)) + $m; }
    @c = map { $_ // R(0) } @c; $nt++; }
  else { die "engine needed at genus $g for (@T): refused (NCH/time)" if $g >= 7;
    if ($eng_g != $g) { Hodge::setup($g); %HodgePt::memo = (); $eng_g = $g; } my $v = HodgePt::WP([], [@T]); @c = map { R("$_") } @$v; $ne++; }
  return $xc{$key} = \@c; }
my %rnd; my $RSEED = 12345;
sub Sg {   # sigma^g X_g(T) as a polynomial in tau (T nonempty, parts >= 2)
  my ($g, $T) = @_; my @T = sort { $a <=> $b } @$T; my $key = "$g|" . join(',', @T);
  if ($MODE eq 'random') { return $rnd{$key} //= do { my $h = 0; $h = ($h * 131 + ord($_)) % 1000003 for split //, "$key|$RSEED"; pconst(R(($h % 997) - 498, 1 + ($h % 13))) } }
  if ($g == 0) { my $r = @T; my $J = 0; $J += $_ for @T; return pconst($r == 1 ? R(1, $J*($J+1)) : $r == 2 ? R(1, $J) : do { my $p = R(1); $p *= ($J - $_) for 1..$r-3; $p }); }
  my $c = comps($g, \@T); my $res = [];
  for my $k (0..$#$c) { next if $c->[$k] == 0; $res = padd($res, pscal(ppow($SIG, $g - $k), (($k % 2) ? R(-1) : R(1)) * $c->[$k])); }
  return $res; }
# ---------- empty-profile leading coefficients from closed forms ----------
# function = [J, {m => poly-in-tau coefficient of w^m}]  meaning e^{J u} Sum_m c_m w^m
sub opD { my ($f) = @_; my ($J, $c) = @$f; my %n;
  for my $m (keys %$c) { my $v = $c->{$m}; $n{$m} = padd($n{$m} // [], $v, R(-$J)); $n{$m+1} = padd($n{$m+1} // [], $v, R($J - $m)); $n{$m+2} = padd($n{$m+2} // [], $v, R($m)); }
  return [$J, \%n]; }
sub opDy { my ($f) = @_; my ($J, $c) = @$f; my %n;
  for my $m (keys %$c) { my $v = $c->{$m}; $n{$m+1} = padd($n{$m+1} // [], $v, R($J)); $n{$m+2} = padd($n{$m+2} // [], $v, R($m)); }
  return [$J + 1, \%n]; }
sub lead_at { my ($f, $ord) = @_; my ($J, $c) = @$f; my @nz = sort { $a <=> $b } grep { @{$c->{$_}} } keys %$c;
  die "empty-profile function more singular than expected: lowest w^$nz[0] < w^$ord" if @nz && $nz[0] < $ord; return $c->{$ord} // []; }
# G^(0)_0 = u - 3u^2/4 + u^3/6, u = 1 - 1/w
my $G00 = do { my %c; my @u = ([0, R(1)], [-1, R(-1)]);   # u as {power => coef}
  my %u1 = (0 => R(1), -1 => R(-1)); my %u2 = (0 => R(1), -1 => R(-2), -2 => R(1)); my %u3 = (0 => R(1), -1 => R(-3), -2 => R(3), -3 => R(-1));
  for my $m (-3..0) { my $v = ($u1{$m} // R(0)) - R(3,4) * ($u2{$m} // R(0)) + R(1,6) * ($u3{$m} // R(0)); $c{$m} = pconst($v); } [0, \%c] };
# sigma * D G^(1)_0 = (w-1)(sigma (w-1) - 1)/24 = [sigma w^2 - (2 sigma + 1) w + (sigma + 1)]/24 ;  sigma * d_y G^(1)_0 = e^u w (sigma(w - 1) - 1)/24
my $sDG10 = [0, { 2 => pscal($SIG, R(1,24)), 1 => pscal(padd(pscal($SIG, R(2)), [R(1)]), R(-1,24)), 0 => pscal(padd($SIG, [R(1)]), R(1,24)) }];
my $sDyG10 = [1, { 2 => pscal($SIG, R(1,24)), 1 => pscal(padd($SIG, [R(1)]), R(-1,24)) }];
my %L0c;
sub Lead0 {   # sigma^g * leading coefficient of e^{-vu} d_y^v D^k G^(g)_0 at w^{2g-2+k+v}  (v >= 1 or k >= 1)
  my ($g, $k, $v) = @_; my $key = "$g|$k|$v"; return $L0c{$key} if exists $L0c{$key};
  return $L0c{$key} = [] if $g >= 2;
  die "undifferentiated empty-profile function" if $k == 0 && $v == 0;
  my $f;
  if ($g == 0) { $f = $G00; $f = opD($f) for 1..$k; $f = opDy($f) for 1..$v; }
  else { if ($k >= 1) { $f = $sDG10; $f = opD($f) for 2..$k; $f = opDy($f) for 1..$v; } else { $f = $sDyG10; $f = opDy($f) for 2..$v; } }
  return $L0c{$key} = lead_at($f, 2*$g - 2 + $k + $v); }
# ---------- equations ----------
my %EQS = (M => [[R(1), [1,1]], [R(1), [2]]], JM2 => [[R(1), [1,1,1]], [R(3), [1,2]], [R(-4), [3]]]);
for my $k (2..6) { $EQS{"P$k"} = [[R(1), [$k]], ref_lib::pk_minus_Dtilde($k)]; }
my $P = $EQS{$EQ} or die "unknown equation $EQ";
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @out;
  for my $sp (setparts($n-1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n-1; push @out, \@c; } my @c = map { [@$_] } @$sp; push @c, [$n-1]; push @out, \@c; }
  return @out; }
sub blockval { my ($B, $g, $k, $T) = @_;   # $B = index values
  my $sz = scalar @$B; my $odd = $sz % 2; return [] if $odd && $k < 1;
  my @U = (@$T, grep { $_ >= 2 } @$B); my $v = grep { $_ == 1 } @$B;
  my $c = R($odd ? -1 : 1); $c *= $_ for @$B; $c *= 2 if !$odd && $k == 0; $c /= fact($k);
  my $p = pmul(ppow($SIG, $sz - 1), ppow($TAU, $k)); $p = pmul($p, Pi($_)) for grep { $_ >= 2 } @$B;
  my $L;
  if (@U) { my $J = 0; $J += $_ for @U; my $D = $J - scalar(@U); my $f = ff(R($D - 2*$g + 2), $k) * rf(R($J), $v); return [] if $f == 0; $L = pscal(Sg($g, \@U), $f); }
  else { $L = Lead0($g, $k, $v); }
  return pscal(pmul($p, $L), $c); }
sub general { my ($n, $R) = @_; my @R = @$R; my $r = @R; my $tot = [];
  for my $t (@$P) { my ($c, $idx) = @$t; my $m = scalar @$idx;
    for my $sp (setparts($m)) { my @blocks = map { [ map { $idx->[$_] } @$_ ] } @$sp; my $q = @blocks;
      for my $code (0..$q**$r - 1) { my @T = map { [] } 1..$q; my $cc = $code; for my $i (0..$r-1) { push @{$T[$cc % $q]}, $R[$i]; $cc = int($cc / $q); }
        my $rec; $rec = sub { my ($b, $left, $acc) = @_;
          if ($b == $q) { $tot = padd($tot, pscal($acc, $c)) if $left == 0; return; }
          my $sz = scalar @{$blocks[$b]}; my $kmin = ($sz % 2) ? 1 : 0;
          for my $g (0..int(($left + 2)/2)) { for my $k ($kmin..($left + 2 - $sz - 2*$g)) { my $cont = 2*$g - 2 + $sz + $k; next if $cont < 0 || $cont > $left;
              my $bv = blockval($blocks[$b], $g, $k, $T[$b]); next unless @$bv; $rec->($b + 1, $left - $cont, pmul($acc, $bv)); } } };
        $rec->(0, $n, [R(1)]); } } }
  return $tot; }
# the report's Theorem M, transcribed separately (report s.4), for P = D_1^2 + D_2
sub PhiM { my ($g, $k, $T) = @_;
  if (!@$T) { return [] if $g > 0; return $TAU if $k == 1; return pscal(ppow($TAU, 2), R(1,2)) if $k == 2; return []; }
  my $J = 0; $J += $_ for @$T; my $D = $J - @$T; my $f = ff(R($D - 2*$g + 2), $k); return [] if $f == 0;
  return pscal(pmul(ppow($TAU, $k), Sg($g, $T)), $J * $f / fact($k)); }
sub thmM { my ($G, $eps, $R) = @_; my @R = @$R; my $n = 2*$G + $eps; my $J = 0; $J += $_ for @R; my $D = $J - @R; my @R2 = (@R, 2); my $tot = [];
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1; my $f = ff(R($D + 1 - 2*$g + 2), $k); next if $f == 0;
    $tot = padd($tot, pmul(Pi(2), pmul(ppow($TAU, $k), Sg($g, \@R2))), R(-2) * $f / fact($k)); }
  if ($eps == 0) { if (@R) { $tot = padd($tot, pmul($SIG, Sg($G, \@R)), R(2) * $J * ($J+1)); }
    elsif ($G == 0) { $tot = padd($tot, $SIG); } elsif ($G == 1) { $tot = padd($tot, pmul($SIG, padd($SIG, [R(1)])), R(-1,12)); } }
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1;
    if (@R) { my $f = ff(R($D - 2*$g + 2), $k); next if $f == 0; $tot = padd($tot, pmul($SIG, pmul(ppow($TAU, $k), Sg($g, \@R))), $J * ($J+1) * $f / fact($k)); }
    elsif ($g == 0 && ($k == 1 || $k == 2)) { $tot = padd($tot, pmul($SIG, ppow($TAU, $k)), R(1) / fact($k)); } }
  my $r = @R;
  for my $mask (0..(1<<$r)-1) { my (@A, @B); for my $i (0..$r-1) { if ($mask & (1<<$i)) { push @A, $R[$i] } else { push @B, $R[$i] } }
    for my $g1 (0..$n+2) { for my $k1 (1..$n+2) { my $rest = $n + 2 - 2*$g1 - $k1; next if $rest < 1; my $f1 = PhiM($g1, $k1, \@A); next unless @$f1;
        for my $g2 (0..$rest) { my $k2 = $rest - 2*$g2; next if $k2 < 1; my $f2 = PhiM($g2, $k2, \@B); next unless @$f2; $tot = padd($tot, pmul($f1, $f2)); } } } }
  return $tot; }
# ---------- driver ----------
my $wt = 0; for my $t (@$P) { my $w = 0; $w += $_ for @{$t->[1]}; $wt = $w; }
my @Rs = ([]); my $gen; $gen = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @p = (@$pre, $j); push @Rs, [@p]; $gen->(\@p, $j, $J + $j) if @p < $rmax; } }; $gen->([], 2, 0);
if (defined $ENV{RLIST}) { @Rs = map { [ $_ eq '' ? () : split /,/, $_ ] } split /;/, $ENV{RLIST}, -1; }
my ($N, $bad) = (0, 0);
for my $G ($Gmin..$Gmax) { for my $R (@Rs) { for my $eps (0, 1) { my $n = 2*$G + $eps;
  my $gv = general($n, $R); my $J = 0; $J += $_ for @$R; my $DRk = $J - @$R + $wt - 1;
  my $line = sprintf "%s G=%d R=(%s) %s [target level D(R+%d)=%d, excess %d]", $EQ, $G, join(',', @$R), ($eps ? 'ODD ' : 'EVEN'), $wt, $DRk, 2*$G - 2 - $DRk;
  my $ok;
  if ($MODE eq 'random') { die "random mode needs EQ = M" unless $EQ eq 'M'; my $tm = thmM($G, $eps, $R); my $d = padd($gv, $tm, R(-1)); $ok = !@$d; $line .= $ok ? " general == ThmM (nonzero: " . (@$gv ? 'yes' : 'no') . ")" : " MISMATCH: diff = " . pstr($d); }
  else { $ok = !@$gv; $line .= $ok ? " zero polynomial" : " NONZERO: " . pstr($gv);
    if ($EQ eq 'M') { my $tm = thmM($G, $eps, $R); my $d = padd($gv, $tm, R(-1)); $line .= @$d ? " [ThmM transcription differs!]" : " [= ThmM]"; } }
  $N++; $bad++ unless $ok; print "$line\n"; } } }
printf "SUMMARY mode=%s eq=%s: %d relations, %d failures (table values used %d, engine values %d)\n", $MODE, $EQ, $N, $bad, $nt, $ne;
