#!/usr/bin/perl
# Second referee: mechanical check of the LOGIC of the induction (no numerical data involved).
# (L1) Corollary I_P for P = P_k (k = 2 is the author's Corollary I_2): for every G, every level D (1 <= D <= 2G-3) and every R at level
#      D-k+1, enumerate ALL terms of the dilaton-point leading order of P_k(D) tau_0.tau_1 at (x^{2G}, q_R) and classify each:
#        target  : one block of genus G with U = R+k, (|B|,k_B) = (1,1), single monomial D_{t_k};
#        FF      : a non-empty block with 0 <= delta < k_B (falling factorial (delta)_k = 0);
#        L0      : an empty block whose leading coefficient vanishes (computed exactly, R2rel::L0);
#        H1      : a non-empty block of genus 2 <= g < G with delta <= -1 (in the Lemma-A range at genus g);
#        H2      : a non-empty block of genus G with D(U) < D;
#        SURVIVOR: anything else (would break the corollary).  Also flags any genus-G block with D(U) > D (circularity).
# (L2) Theorem B step of proof (i): after I_2 at level D, every T at level D without a part 2 is reached by Theorem B at
#      R = T - m + (m-1) (m = min T) with all other genus-G terms already known (part 2 / fewer parts / smaller minimum), and every
#      cut / join term in the range by excess counting (string reduction of unit parts).
# Usage: perl r2_logic.pl k Gmin Gmax [rmax]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2rel;
$| = 1;
my ($kP, $Gmin, $Gmax, $rmax) = @ARGV; $rmax //= 999;
my $P = R2rel::hirota_Pk($kP); my $Q = R2rel::ring_Q(); my $sigtest = Math::BigRat->new('2/7') * Math::BigRat->new('9/7');
print "# P_$kP = ", R2rel::show_P($P), "\n";
# multisets of parts >= 2 at level L (sum (j-1) = L) with at most rmax parts: partitions of L, parts +1
sub parts_of { my ($n, $mx) = @_; $mx //= $n; return ([]) if $n == 0; my @o; for my $k (reverse 1..($mx < $n ? $mx : $n)) { push @o, [$k, @$_] for parts_of($n - $k, $k); } @o }
my %lvl; sub level_sets { my $L = shift; $lvl{$L} //= [ map { [ sort { $a <=> $b } map { $_ + 1 } @$_ ] } grep { @$_ <= $rmax } parts_of($L) ]; $lvl{$L} }
sub cnt { my %c; $c{$_}++ for @_; \%c }
sub submultisets { my @l = @_; my $c = cnt(@l); my @vals = sort { $a <=> $b } keys %$c; my @s = ([]);
  for my $v (@vals) { my @n; for my $x (@s) { for my $i (0..$c->{$v}) { push @n, [@$x, ($v) x $i]; } } @s = @n; } @s }
sub minus { my ($l, $s) = @_; my $c = cnt(@$s); my @r; for my $v (@$l) { if (($c->{$v} // 0) > 0) { $c->{$v}--; } else { push @r, $v; } } \@r }
sub distribute { my ($R, $q) = @_; return ([[@$R]]) if $q == 1; my @o; for my $s (submultisets(@$R)) { my $rest = minus($R, $s); push @o, [$s, @$_] for distribute($rest, $q - 1); } @o }
my ($nrel, $nsurv, $ncirc, %why) = (0, 0, 0);
my %l0zero;
sub classify_relation { my ($G, $R, $D) = @_; my $n = 2*$G; my ($target, $surv, $circ) = (0, 0, 0);
  for my $mono (@$P) { my ($c, $alpha) = @$mono; my $q = @$alpha;
    for my $sp (R2rel::setparts($q)) { my $nb = @$sp; my @Bi = map { [ map { $alpha->[$_] } @$_ ] } @$sp;
      for my $dist (distribute($R, $nb)) {
        my @blk; for my $ib (0..$nb-1) { my $B = $Bi[$ib]; my @U = sort { $a <=> $b } (@{$dist->[$ib]}, grep { $_ >= 2 } @$B); my $v = grep { $_ == 1 } @$B;
          my $J = 0; $J += $_ for @U; push @blk, { B => $B, U => \@U, v => $v, DU => $J - scalar(@U), odd => @$B % 2 }; }
        my $rec; $rec = sub { my ($ib, $budget, $tags) = @_;
          if ($ib == $nb) { return unless $budget == 0;
            my %t = map { $_ => 1 } @$tags;
            if ($t{TARGET}) { $target++; return; }
            for my $k (qw(FF L0 H1 H2)) { if ($t{$k}) { $why{$k}++; return; } }
            $surv++; return; }
          my $x = $blk[$ib];
          for my $g (0..int(($budget + 2) / 2)) { for my $k (($x->{odd} ? 1 : 0)..($budget + 2 - 2*$g - scalar(@{$x->{B}}))) {
              my $wt = 2*$g - 2 + scalar(@{$x->{B}}) + $k; next if $wt < 0 || $wt > $budget; my $tag = 'NONZERO';
              if (@{$x->{U}}) { my $delta = $x->{DU} - 2*$g + 2;
                if ($delta >= 0 && $delta < $k) { $tag = 'FF'; }
                elsif ($g == $G) { if ($x->{DU} > $D) { $circ++; $tag = 'CIRC'; }
                  elsif ($x->{DU} < $D) { $tag = 'H2'; }
                  else { # level D at genus G: must be the target
                    if ($nb == 1 && @$alpha == 1 && $alpha->[0] == $kP && $k == 1 && join(',', @{$x->{U}}) eq join(',', sort { $a <=> $b } (@$R, $kP))) { $tag = 'TARGET'; } else { $tag = 'NONTARGET-LEVEL-D'; } } }
                elsif ($g >= 2 && $delta <= -1) { $tag = 'H1'; } }
              else { my $key = "$g|$k|$x->{v}"; $l0zero{$key} //= ($g >= 2 || R2rel::L0($Q, $sigtest, $g, $k, $x->{v}) == 0) ? 1 : 0; $tag = 'L0' if $l0zero{$key}; }
              $rec->($ib + 1, $budget - $wt, [@$tags, $tag]); } } };
        $rec->(0, $n, []); } } }
  ($target, $surv, $circ) }
# ---- (L1)
for my $G ($Gmin..$Gmax) { my ($rel, $sv, $ci, $bt) = (0, 0, 0, 0);
  for my $D (1..2*$G - 3) { my $LR = $D - $kP + 1; next if $LR < 0;
    for my $R (@{ level_sets($LR) }) { my ($t, $s, $c) = classify_relation($G, $R, $D); $rel++; $sv += $s; $ci += $c; $bt++ unless $t == 1;
      printf "  SURVIVOR/CIRC at G=%d D=%d R=(%s): target=%d survivors=%d circular=%d\n", $G, $D, join(',', @$R), $t, $s, $c if $s || $c || $t != 1; } }
  $nrel += $rel; $nsurv += $sv; $ncirc += $ci;
  printf "(L1) P_%d G=%d: %d relations (R at levels D-%d, D = 1..%d, r <= %s), survivors %d, circular %d, relations without exactly one target %d\n", $kP, $G, $rel, $kP - 1, 2*$G - 3, $rmax, $sv, $ci, $bt; }
printf "(L1) SUMMARY P_%d G=%d..%d: %d relations, %d surviving non-target terms, %d circular terms; zero-reasons: %s\n", $kP, $Gmin, $Gmax, $nrel, $nsurv, $ncirc, join(', ', map { "$_ $why{$_}" } sort keys %why);
# ---- (L2) Theorem B induction (proof (i)), only meaningful with k = 2
exit 0 unless $kP == 2;
my ($tb, $tbad) = (0, 0);
for my $G ($Gmin..$Gmax) { my $gaps = 0;
  for my $D (1..2*$G - 3) { my @TL = @{ level_sets($D) }; my %known;
    $known{join(',', @$_)} = 1 for grep { $_->[0] == 2 } @TL;   # I_2 (TL sorted ascending: part 2 iff first part 2)
    for my $T (sort { @$a <=> @$b || $a->[0] <=> $b->[0] } grep { $_->[0] >= 3 } @TL) { my $m = $T->[0]; my @R = sort { $a <=> $b } ($m - 1, @{$T}[1..$#$T]); $tb++;
      my $ok = 1; my @why;
      # genus-G terms of Theorem B at R: R+2 (has a 2), R_i^+ (raise one part), R_ik (merge two parts); level D each
      my @terms = ( [ sort { $a <=> $b } (@R, 2) ] );
      for my $i (0..$#R) { my @x = @R; $x[$i]++; push @terms, [ sort { $a <=> $b } @x ]; }
      for my $i (0..$#R) { for my $j ($i+1..$#R) { my @x = map { $R[$_] } grep { $_ != $i && $_ != $j } 0..$#R; push @terms, [ sort { $a <=> $b } (@x, $R[$i] + $R[$j]) ]; } }
      my $tk = join(',', @$T); my $seenT = 0;
      for my $u (@terms) { my $uk = join(',', @$u); my $J = 0; $J += $_ for @$u; die "level" unless $J - @$u == $D;
        if ($uk eq $tk) { $seenT++; next; } next if $known{$uk}; $ok = 0; push @why, "($uk) unknown"; }
      # T must occur exactly via raising the unique part m-1 (count occurrences)
      $ok = 0, push @why, "T occurs $seenT times" unless $seenT == 1;
      # cut terms Xi_{G-1}(al,be,R-j) and join terms Xi_g1(al,R_I) Xi_g2(be,R_J): excess bookkeeping (unit parts by string)
      for my $i (0..$#R) { my $j = $R[$i]; my @Rm = map { $R[$_] } grep { $_ != $i } 0..$#R; my $Dm = 0; $Dm += $_ - 1 for @Rm;
        for my $al (1..$j-1) { my $be = $j - $al; my $Dc = $Dm + ($al - 1) + ($be - 1);   # cut profile level
          my $exc = 2*($G-1) - 2 - $Dc; unless ($G - 1 >= 2 && $exc >= 1) { $ok = 0; push @why, "cut term not in range (G-1=" . ($G-1) . ", excess $exc)"; }
          for my $I (submultisets(@Rm)) { my $Jr = minus(\@Rm, $I); my $D1 = ($al - 1); $D1 += $_ - 1 for @$I; my $D2 = ($be - 1); $D2 += $_ - 1 for @$Jr;
            for my $g1 (0..$G) { my $g2 = $G - $g1; my $e1 = 2*$g1 - 2 - $D1; my $e2 = 2*$g2 - 2 - $D2;
              my $z1 = ($g1 >= 2 && $g1 < $G && $e1 >= 1) || ($g1 == $G && $D1 < $D && $e1 >= 1) || ($g1 >= 1 && $D1 == 0 && $al == 1 && !@$I && $g1 >= 2);
              my $z2 = ($g2 >= 2 && $g2 < $G && $e2 >= 1) || ($g2 == $G && $D2 < $D && $e2 >= 1) || ($g2 >= 1 && $D2 == 0 && $be == 1 && !@$Jr && $g2 >= 2);
              unless ($z1 || $z2) { $ok = 0; push @why, "join term g=($g1,$g2) levels ($D1,$D2) excess ($e1,$e2) not killed"; } } } } }
      unless ($ok) { $tbad++; $gaps++; printf "  (L2) GAP G=%d D=%d T=(%s) via R=(%s): %s\n", $G, $D, $tk, join(',', @R), join('; ', @why[0..($#why < 3 ? $#why : 3)]); }
      $known{$tk} = 1; }
    my @miss = grep { !$known{join(',', @$_)} } @TL; if (@miss) { $tbad++; $gaps++; printf "  (L2) level G=%d D=%d: %d profiles never reached\n", $G, $D, scalar @miss; } }
  printf "(L2) Theorem B step, G=%d: %s\n", $G, ($gaps ? "$gaps GAPS" : "every profile at every level D <= 2G-3 reached, all side terms killed by (H1)/(H2)"); }
printf "(L2) SUMMARY G=%d..%d: %d Theorem-B steps, %d gaps\n", $Gmin, $Gmax, $tb, $tbad;
