# ref_n2_n4.pl : referee re-computation of the negative results N2 (single ch_l insertions, l = 3, 5) and N4 (higher Taylor coefficients),
# from REC data only (referee's bivariate Y_h of ref_k1.pl and X_h of Theorem B), i.e. WITHOUT any Hodge engine or moment table:
#   Q_l(h,j) := Int c(T(-log)) ch_l lam_j - (irreducible-boundary part)
#            = (B_{l+1}/(l+1)!) [ <kappa_l> - <psi^l> + Sum_{g1=1}^{h-1} Sum_{a=0}^{l-1} (-1)^a [y^a]Y_{g1}(j) X_{h-g1}[l-1-a] ],
#   <kappa_l> = Sum_{c>=l+1} [y^c] Y_h = (2h-1) X_h - Sum_{c=1}^{l} [y^c] Y_h,   <psi^l> = Sum_a j^a X_h[a+l],   X_g[b] = [j^b] X_g(j).
# Compared with the author's irr_test_out.txt.   N4: A^(c)_g := [y^c]Y_g(1); test B^(c)_g(j) = Sum_{g1} A^(c)_{g-g1} T_{g1}(j)/j at j = 2.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL;
sub Q { Math::BigRat->new(shift) }
my $Hm = shift // 6; my $dir = "$FindBin::Bin/out";
sub loadpoly { my ($f) = @_; open my $fh, '<', $f or die "no $f"; my $l = <$fh>; close $fh; chomp $l; my ($h, $r, $k, $rest) = split /\s*\|\s*/, $l, 4; my %c;
  for my $t (split /\s*;\s*/, $rest) { my ($m, $v) = split /:/, $t; $c{$m} = Q($v) } return \%c }
my (%X, %Y);
for my $g (1..$Hm) { my $x = loadpoly("$dir/X_ref_h$g.txt"); my @p; for my $a (keys %$x) { $p[$a] = $x->{$a} } $_ //= Q(0) for @p; $X{$g} = \@p;
  my $y = loadpoly("$dir/Y_ref_h$g.txt"); my @b; for my $m (keys %$y) { my ($a, $c) = split /,/, $m; $b[$c][$a] = $y->{$m} } for my $row (@b) { $row //= []; $_ //= Q(0) for @$row; } $Y{$g} = \@b; }
sub Xc { my ($g, $b) = @_; return $X{$g}[$b] // Q(0) }
sub Yc { my ($g, $c, $j) = @_; my $row = $Y{$g}[$c] or return Q(0); return RL::pev($row, $j) }   # [y^c] Y_g(j, y) at j
sub Xv { my ($g, $j) = @_; return RL::pev($X{$g}, $j) }
my %irr; if (open my $fh, '<', "$FindBin::Bin/../irr_test_out.txt") { while (my $l = <$fh>) { if ($l =~ /^h=(\d+) l=(\d+): (.*)$/) { my ($h, $ll, $rest) = ($1, $2, $3); for my $t (split /\s+/, $rest) { my ($jj, $v) = $t =~ /j=(\d+):(\S+)/ or next; $irr{"$h,$ll,$jj"} = Q($v) } } } close $fh }
for my $l (3, 5) { my $Bc = RL::bern($l + 1) / RL::fact($l + 1);
  for my $h (2..$Hm) { next if $l > $h && $l == 5 && $h < 5; my @out; my ($cmp, $bad) = (0, 0);
    for my $j (2..2*$h) { my $jq = Q($j);
      my $kap = (2*$h - 1) * Xv($h, $jq); $kap -= Yc($h, $_, $jq) for 1..$l;
      my $ps = Q(0); for my $a (0..3*$h) { $ps += $jq ** $a * Xc($h, $a + $l) }
      my $sep = Q(0); for my $g1 (1..$h-1) { my $g2 = $h - $g1; for my $a (0..$l-1) { $sep += Q((-1) ** $a) * Yc($g1, $a, $jq) * Xc($g2, $l - 1 - $a) } }
      my $q = $Bc * ($kap - $ps + $sep); push @out, "j=$j:" . ($q == 0 ? '0' : $q);
      if (exists $irr{"$h,$l,$j"}) { $cmp++; $bad++ if $irr{"$h,$l,$j"} != $q } }
    printf "h=%d l=%d: %s\n   vs author's irr_test_out: %d compared, %s\n", $h, $l, join('  ', @out), $cmp, ($bad ? "$bad MISMATCHES" : "all equal"); } }
# N4
$RL::NT = 2*$Hm + 4; my $one = RL::szero(); $one->[0] = Q(1); my $L = RL::logH(); my %T; for my $n (0..$Hm) { $T{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one, $L, 2*$n)) }
for my $c (1..4) { my %A; my @rep;
  for my $g (1..$Hm) { $A{$g} = Yc($g, $c, Q(1)); my @good;
    for my $j (1..2*$g) { my $jq = Q($j); my $s = Q(0); for my $g1 (0..$g-1) { $s += $A{$g - $g1} * RL::pev($T{$g1}, $jq) / $jq } push @good, $j if $s == Yc($g, $c, $jq) }
    push @rep, "g=$g: " . join(',', @good); }
  print "N4 c=$c: boundary-sum ansatz with A^(c)_m = [y^c]Y_m(1) holds at j = ", join(' | ', @rep), "\n"; }
