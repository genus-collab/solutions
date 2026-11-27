#!/usr/bin/perl
# a3_kp_nullity2_sub.pl -- variant of a3_kp_nullity2.pl with env MAXSMALL (polynomiality imposed only for T whose smallest part <= MAXSMALL) -- same question as a3_kp_nullity.pl (does the KP hierarchy at the dilaton point force Lemma A at genus G, given lower
# genera?) in a reduced formulation that scales to higher genus.
# Unknowns: only the one-leaf polynomials X_t(sigma) = sigma^G Xi_G(t) = sum_{k=1}^{G} c_{t,k} sigma^k, t = 2..2G-2 (lambda_g term = 0).
# At each sample tau (mod p) every X_T (T in the range, r >= 2) is computed as a linear form in the c's from ONE KP relation
# (merge the two largest parts), level by level; every other KP relation of the same level gives a linear equation; finally every X_T must be
# a polynomial in sigma of degree <= G with zero constant term (left kernel of the sample Vandermonde matrix) -> more linear equations.
# KP relation (R, al, be), al,be >= 2, D(R)+al+be-2 <= 2G-3, all factors genus G or 0 (lower genera vanish in the range by induction):
#   sum_{a<=al, b<=be} sum_{B subset R} Pi_a Pi_b Xi_G(R_B + {a,b}) E0(al-a, be-b; R - R_B) = Pi_{al+be-1} Xi_G(R + {1, al+be-1}),
# with E0 = coefficients of exp(sum z^-a w^-b sigma Pi_a Pi_b Xi_0(. + {a,b})) and Xi_0(T) = Gamma(J)/Gamma(J-r+3).
# Usage: perl a3_kp_nullity2.pl G [samples]
use strict; use warnings;
my ($G, $S, $rmax) = @ARGV; $S //= $G + 8; $S = $G + 8 if $S == 0; $rmax //= 1000;   # rmax: only multisets with <= rmax parts (restricted subsystem)
my $p = 2147483647;
sub md { my $x = shift; $x %= $p; $x += $p if $x < 0; return $x; }
sub mm { return md($_[0] * $_[1]); }
sub mi { my $a = md(shift); my ($t, $nt, $r, $nr) = (0, 1, $p, $a); while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q*$nt); ($r, $nr) = ($nr, $r - $q*$nr); } die "not invertible" if $r != 1; return md($t); }
sub parts { my ($n, $max) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
sub key { return join(',', sort { $b <=> $a } @_); }
my $Dmax = 2*$G - 3;
my @T; my %Tidx; my @byLevel;
for my $D (1..$Dmax) { my @q = parts($D, $D); @q = sort { scalar(@$a) <=> scalar(@$b) } @q;   # fewer parts first
  for my $q (@q) { my @t = map { $_ + 1 } @$q; $Tidx{key(@t)} = scalar @T; push @{$byLevel[$D]}, scalar @T; push @T, \@t; } }
my $nT = scalar @T; my $nC = ($Dmax) * $G;   # single parts t = 2..2G-2  <->  level D = t-1 = 1..Dmax
sub cidx { my ($t, $k) = @_; return ($t - 2) * $G + $k - 1; }
my @Rs = ([]); for my $D (1..$Dmax) { for my $q (parts($D, $D)) { push @Rs, [ map { $_ + 1 } @$q ]; } }
srand(271828);
my @eqs; my %samplesV;   # samplesV{Tindex} = [ [sig, vec], ... ]
my $zero = sub { return [ (0) x $nC ]; };
for my $s (1..$S) {
  my $tau = 3 + int(rand($p - 10)); my $sig = mm($tau, md($tau + 1));
  my %Pc; my $Pi = sub { my $j = shift; return 1 if $j <= 1; return $Pc{$j} //= do { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } my $f = 1; $f = mm($f, $_) for 1..$j-1; mm($v, mi($f)); }; };
  my $X0 = sub { my @L = @_; my $J = 0; $J += $_ for @L; my $r = scalar @L; if ($r >= 3) { my $v = 1; $v = mm($v, md($J - $_)) for 1..($r-3); return $v; } elsif ($r == 2) { return mi($J); } else { return mi(mm($J, $J+1)); } };
  my @V;   # V[Tindex] = vector (linear form in c) of X_T at this sample
  for my $D (1..$Dmax) { my $t = $D + 1; my $i = $Tidx{key($t)}; my $v = $zero->(); my $sp = 1; for my $k (1..$G) { $sp = mm($sp, $sig); $v->[cidx($t,$k)] = $sp; } $V[$i] = $v; }
  # X_G of a list with possible unit parts, as a vector (string equation); undef-safe
  my $XGv = sub { my @L = @_; my @Tp = grep { $_ >= 2 } @L; my $b = scalar(@L) - scalar(@Tp); return undef unless @Tp; my $J = 0; $J += $_ for @Tp;
    my $f = 1; $f = mm($f, md($J + $_)) for 0..$b-1; my $i = $Tidx{key(@Tp)}; die "out of range @L" unless defined $i; my $v = $V[$i]; return undef unless defined $v;
    return [ map { mm($f, $_) } @$v ]; };
  # build KP relations grouped by target level; process levels in order
  my %relsByLevel;
  for my $R (@Rs) { my @R = @$R; my $DR = 0; $DR += $_ - 1 for @R; next if scalar(@R) + 2 > $rmax; my $wmax = 2*$G - 1 - $DR; next if $wmax < 4;
    for my $al (2..$wmax) { for my $be ($al..$wmax) { next if $al + $be > $wmax; push @{$relsByLevel{$DR + $al + $be - 2}}, [\@R, $al, $be]; } } }
  for my $D (1..$Dmax) { next unless $relsByLevel{$D};
    # each relation: main term Pi_al Pi_be X(R+{al,be}) + lower(known) - Pi_c (J+c) X(R+c) = 0 ; solve with the canonical relation for each multi-part T
    my %relvec;   # evaluate relation as: coef_main * X_main + rest_vector (rest includes RHS), where rest uses known V's
    for my $rel (@{$relsByLevel{$D}}) { my ($Rr, $al, $be) = @$rel; my @R = @$Rr; my $nR = scalar @R;
      my %m; $m{$_}++ for @R; my @vv = sort { $a <=> $b } keys %m; my @nn = map { $m{$_} } @vv; my $nv = scalar @vv;
      my @subs = ([]); for my $ii (0..$nv-1) { my @ns; for my $ss (@subs) { push @ns, [@$ss, $_] for 0..$nn[$ii]; } @subs = @ns; }
      my %sidx; $sidx{join(",", @{$subs[$_]})} = $_ for 0..$#subs; my @list = map { my $ss = $_; [ map { ($vv[$_]) x $ss->[$_] } 0..$nv-1 ] } @subs;
      my $full = $sidx{join(",", @nn)}; my $empty = $sidx{join(",", (0) x $nv)};
      my $binom = sub { my ($n, $k) = @_; my $v = 1; for my $q (1..$k) { $v = mm($v, mm($n - $q + 1, mi($q))); } return $v; };
      my @pairs; for my $si (0..$#subs) { my $Sv = $subs[$si]; my @pl = ([]); for my $ii (0..$nv-1) { my @np; for my $pp (@pl) { push @np, [@$pp, $_] for 0..$Sv->[$ii]; } @pl = @np; }
        for my $kv (@pl) { my $w = 1; $w = mm($w, $binom->($Sv->[$_], $kv->[$_])) for 0..$nv-1; push @{$pairs[$si]}, [ $sidx{join(",", @$kv)}, $sidx{join(",", map { $Sv->[$_] - $kv->[$_] } 0..$nv-1)}, $w ]; } }
      my %om0; my $om0 = sub { my ($a, $b, $si) = @_; return $om0{"$a,$b,$si"} //= mm($sig, mm(mm($Pi->($a), $Pi->($b)), $X0->(@{$list[$si]}, $a, $b))); };
      my %Em; my $E; $E = sub { my ($a, $b, $si) = @_; my $kk = "$a,$b,$si"; return $Em{$kk} if exists $Em{$kk}; my $res = 0;
        if ($a == 0 && $b == 0) { $res = $si == $empty ? 1 : 0; } elsif ($a <= 0 || $b <= 0) { $res = 0; }
        else { for my $a1 (1..$a) { for my $b1 (1..$b) { for my $pr (@{$pairs[$si]}) { my ($ki, $ci, $w) = @$pr; my $e = $E->($a-$a1, $b-$b1, $ci); next unless $e; $res = md($res + mm(mm($w, mm($a1, $om0->($a1, $b1, $ki))), $e)); } } } $res = mm($res, mi($a)); }
        $Em{$kk} = $res; return $res; };
      my $mainT = $Tidx{key(@R, $al, $be)}; my $mainc = mm($Pi->($al), $Pi->($be));
      my $build = sub { my $rest = $zero->();
      for my $a (1..$al) { for my $b (1..$be) { for my $pr (@{$pairs[$full]}) { my ($ki, $ci, $w) = @$pr; next if $a == $al && $b == $be && $ki == $full;
            my $e = $E->($al-$a, $be-$b, $ci); next unless $e; my $v = $XGv->(@{$list[$ki]}, $a, $b); next unless defined $v;
            my $c = mm(mm($w, $e), mm($Pi->($a), $Pi->($b))); for my $q (0..$nC-1) { $rest->[$q] = md($rest->[$q] + mm($c, $v->[$q])) if $v->[$q]; } } } }
      my $cc = $al + $be - 1; my $J = 0; $J += $_ for @R; my $vr = $XGv->(@R, $cc); die "rhs unknown" unless defined $vr;   # X(R+{1,cc}) = (J+cc) X(R+cc)
      my $f = mm(md($J + $cc), $Pi->($cc)); for my $q (0..$nC-1) { $rest->[$q] = md($rest->[$q] - mm($f, $vr->[$q])) if $vr->[$q]; }
      return $rest; };
      push @{$relvec{$mainT}}, [ $mainc, $build ]; }
    # within the level: process T in order of increasing number of parts (RHS X(R+cc) has fewer parts and is defined first)
    for my $i (@{$byLevel[$D]}) { next if defined $V[$i] || scalar(@{$T[$i]}) > $rmax; my $L = $relvec{$i} or die "no relation for T @{$T[$i]}";
      my ($mc, $bld) = @{$L->[0]}; my $rest = $bld->(); my $inv = mi($mc); $V[$i] = [ map { md(-mm($inv, $_)) } @$rest ];
      for my $k (1..$#$L) { my ($mc2, $bld2) = @{$L->[$k]}; my $rest2 = $bld2->(); push @eqs, [ map { md(mm($mc2, $V[$i][$_]) + $rest2->[$_]) } 0..$nC-1 ]; } }
    # single-part targets: relations whose main term is a single part do not exist (al,be >= 2 gives r >= 2); but rest vectors were built
    # using V of lower levels and the RHS at this level, which may be a single part (defined) or a multi-part T with fewer parts (defined above).
  }
  for my $i (0..$nT-1) { push @{$samplesV{$i}}, [ $sig, $V[$i] ]; }
}
# polynomiality of X_T for multi-part T: vector (X_T(s))_s must lie in span{ sig_s^k : k = 1..G }
for my $i (0..$nT-1) { next if scalar(@{$T[$i]}) == 1 || scalar(@{$T[$i]}) > $rmax; if (defined $ENV{MAXSMALL}) { my @tt = sort { $a <=> $b } @{$T[$i]}; next if $tt[0] > $ENV{MAXSMALL}; } my $L = $samplesV{$i}; my $n = scalar @$L;
  # left kernel of the n x G matrix M[s][k] = sig_s^k: solve by Gaussian elimination on the transpose (G x n) to get kernel basis of M^T
  my @M = map { my $sg = $L->[$_][0]; my @row; my $pw = 1; for my $k (1..$G) { $pw = mm($pw, $sg); push @row, $pw; } \@row } 0..$n-1;
  # kernel of M^T: vectors w (length n) with sum_s w_s M[s][k] = 0 for all k. Row-reduce M^T (G x n).
  my @A = map { my $k = $_; [ map { $M[$_][$k] } 0..$n-1 ] } 0..$G-1; my @pivc; my $row = 0;
  for my $col (0..$n-1) { my $pr; for my $rr ($row..$G-1) { if ($A[$rr][$col]) { $pr = $rr; last; } } next unless defined $pr; @A[$row, $pr] = @A[$pr, $row];
    my $inv = mi($A[$row][$col]); $A[$row] = [ map { mm($inv, $_) } @{$A[$row]} ];
    for my $rr (0..$G-1) { next if $rr == $row || !$A[$rr][$col]; my $f = $A[$rr][$col]; $A[$rr] = [ map { md($A[$rr][$_] - mm($f, $A[$row][$_])) } 0..$n-1 ]; }
    push @pivc, $col; $row++; last if $row == $G; }
  my %isp = map { $_ => 1 } @pivc; my @freec = grep { !$isp{$_} } 0..$n-1;
  for my $fc (@freec) { my @w = (0) x $n; $w[$fc] = 1; for my $rr (0..$#pivc) { $w[$pivc[$rr]] = md(-$A[$rr][$fc]); }
    my @e = (0) x $nC; for my $s2 (0..$n-1) { next unless $w[$s2]; my $v = $L->[$s2][1]; for my $q (0..$nC-1) { $e[$q] = md($e[$q] + mm($w[$s2], $v->[$q])) if $v->[$q]; } }
    push @eqs, \@e; } }
# rank of eqs (dense, nC columns)
my @piv; my %pivrow; my $rank = 0;
for my $e (@eqs) { last if $rank == $nC; my @v = @$e;
  for my $col (0..$nC-1) { next unless $v[$col]; if (exists $pivrow{$col}) { my $f = $v[$col]; my $pr = $pivrow{$col}; for my $q ($col..$nC-1) { $v[$q] = md($v[$q] - mm($f, $pr->[$q])) if $pr->[$q]; } next; }
    my $inv = mi($v[$col]); my @nr = map { mm($inv, $_) } @v; $pivrow{$col} = \@nr; $rank++; last; } }
my @freeT; for my $t (2..2*$G-2) { my $all = 1; for my $k (1..$G) { $all = 0 unless exists $pivrow{cidx($t,$k)}; } push @freeT, $t unless $all; }
printf "genus %d (reduced KP system, %d samples): %d one-leaf unknowns, %d equations, rank %d, nullity %d%s\n", $G, $S, $nC, scalar(@eqs), $rank, $nC - $rank,
  ($nC - $rank ? "; one-leaf coefficients not all pivots for t = @freeT (indicative)" : "");
