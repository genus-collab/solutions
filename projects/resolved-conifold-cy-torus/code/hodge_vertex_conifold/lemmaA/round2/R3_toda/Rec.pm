package Rec;
# The leading-coefficient recursion (REC) derived from the Toda equation for single Hurwitz numbers (report, section 4):
#   for g >= 1, R a multiset of integers >= 2 (r = |R|, J = sum R), e = 2g-2+r, delta = J - e  (= D(R) - 2g + 2):
#   (e+1) delta Xi_g(R) = sum'_{pi, (g_B,k_B)} prod_{B in pi} (2/(2k_B)!) (delta_B)_{(2k_B)} Xi_{g_B}(R_B),
#   pi a set partition of the r labelled parts, g_B >= 0, k_B >= 1, sum_B (g_B + k_B - 1) = g, the single term
#   (pi = {R}, (g_B,k_B) = (g,1)) excluded; delta_B = J_B - (2g_B - 2 + r_B); (x)_(n) falling factorial.
#   Genus 0 input: Xi_0(R) = (J-1)(J-2)...(J-r+3)  (r >= 3),  Xi_0(j1,j2) = 1/J,  Xi_0(j) = 1/(j(j+1)).
# Xi_g(R) here is Xi^(0)_g(R) = Int_{Mbar_{g,r}} Lambda^vee(1) K prod 1/(1 - j_i psi_i) (Hurwitz component).
use strict; use warnings; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
sub ff { my ($x, $n) = @_; my $p = BR(1); $p *= BR($x - $_) for 0..$n-1; return $p }   # falling factorial x(x-1)...(x-n+1)
my %MEMO;
sub xi0 { my ($R) = @_; my $r = @$R; my $J = 0; $J += $_ for @$R;
  return BR(1) / BR($J * ($J + 1)) if $r == 1;
  return BR(1) / BR($J) if $r == 2;
  return ff($J - 1, $r - 3) }
# set partitions of index list
sub set_partitions { my @idx = @_; return ([]) if !@idx; my $first = shift @idx; my @out;
  for my $p (set_partitions(@idx)) { push @out, [[$first], @$p];
    for my $i (0..$#$p) { my @q = map { [@$_] } @$p; unshift @{$q[$i]}, $first; push @out, \@q } }
  return @out }
# compositions of the genus budget g among blocks: each block gets (g_B, k_B) with m_B = g_B + k_B - 1 >= 0, sum m_B = g
sub xi { my ($g, $R) = @_; my @R = sort { $a <=> $b } @$R; my $key = "$g|" . join(",", @R);
  return $MEMO{$key} if exists $MEMO{$key};
  my $res;
  if ($g == 0) { $res = xi0(\@R) }
  else { my $r = @R; my $J = 0; $J += $_ for @R; my $e = 2*$g - 2 + $r; my $delta = $J - $e;
    if ($delta == 0) { $res = undef }
    else { my $tot = BR(0);
      for my $pi (set_partitions(0..$r-1)) { my @blocks = map { [map { $R[$_] } @$_] } @$pi;
        # distribute m_B >= 0 with sum g; for each block, (g_B,k_B) with g_B + k_B - 1 = m_B, k_B >= 1
        my $rec; $rec = sub { my ($bi, $left, $acc, $single) = @_;
          if ($bi == @blocks) { return if $left != 0; $tot += $acc; return }
          my $B = $blocks[$bi]; my $rB = @$B; my $JB = 0; $JB += $_ for @$B;
          for my $m (0..$left) { for my $gB (0..$m) { my $kB = $m - $gB + 1;
              next if @blocks == 1 && $gB == $g && $kB == 1;           # the excluded linear term
              my $dB = $JB - (2*$gB - 2 + $rB);
              next if $dB >= 0 && $dB <= 2*$kB - 1;                    # falling factorial vanishes
              my $x = xi($gB, $B); die "needed unknown Xi_$gB(@$B)" unless defined $x; next if $x == 0;
              my $f = BR(2) / BR(fact(2*$kB)) * ff($dB, 2*$kB) * $x;
              $rec->($bi + 1, $left - $m, $acc * $f) } } };
        $rec->(0, $g, BR(1)) }
      $res = $tot / BR(($e + 1) * $delta) } }
  $MEMO{$key} = $res; return $res }
my @F = (Math::BigInt->new(1)); sub fact { my $n = shift; while (@F <= $n) { push @F, $F[-1] * scalar(@F) } return $F[$n]->copy }
1;
