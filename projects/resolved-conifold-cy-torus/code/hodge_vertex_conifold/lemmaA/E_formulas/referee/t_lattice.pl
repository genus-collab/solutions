# t_lattice.pl h r [file] : lattice-point analysis of a component file (default: the referee's ref_engine_h<h>_r<r>.txt).
#  (A) Lemma A: every sorted R with parts >= 2 and D <= 2h-3: all components must vanish;
#  (M) Conjecture M (vanishing direction) on the whole lattice: parts >= 0, J >= 1, D <= 2h-3;
#  (C) controls D = 2h-2 (parts >= 2): list of non-zero components;
#  (H) Hurwitz component Xi^(0) on the layers D = 2h-2, 2h-1, 2h (parts >= 2): the author claims Xi^(0) = 0 exactly on D = 2h-1;
#  (S) string/dilaton cannot be tested within one file; the top-degree/KW test is in t_kw.pl.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($h, $r, $f) = @ARGV; $f //= "$FindBin::Bin/ref_engine_h${h}_r${r}.txt";
my @C;   # C[k] = list of [exps, coef]
open my $fh, '<', $f or die "no $f"; while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4; $rest //= '';
  for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $m], Q($cf)]; } } close $fh;
sub ev { my ($k, $p) = @_; my $s = Q(0); for my $t (@{$C[$k] // []}) { my ($e, $cf) = @$t; my $x = $cf->copy; for my $i (0..$r-1) { if ($e->[$i]) { $x *= $p->[$i] ** $e->[$i] } } $s += $x; } return $s; }
sub sorted_tuples { my ($n, $lo, $Jmax) = @_; my @o; my $rec; $rec = sub { my ($pre, $min, $J) = @_; if (@$pre == $n) { push @o, [@$pre]; return } for (my $j = $min; $J + $j <= $Jmax; $j++) { $rec->([@$pre, $j], $j, $J + $j) } }; $rec->([], $lo, 0); return @o; }
my $N = 2*$h - 3 + $r;
# (A) and (C)
my ($nA, $bA, $nC, $nzC) = (0, 0, 0, 0); my @ctrl;
for my $p (sorted_tuples($r, 2, 2*$h - 2 + $r)) { my $D = 0; $D += $_ - 1 for @$p; my @nz = grep { ev($_, [map { Q($_) } @$p]) != 0 } 0..$h;
  if ($D <= 2*$h - 3) { $nA++; if (@nz) { $bA++; print "VIOLATION R = {@$p}: non-zero components @nz\n" } } else { $nC++; $nzC++ if @nz; push @ctrl, "{" . join(',', @$p) . "}:" . join('', @nz == $h + 1 ? ('all') : @nz); } }
printf "(A) genus %d r = %d: %d Lemma-A points (parts >= 2, D <= 2h-3), %d violations\n", $h, $r, $nA, $bA;
printf "(C) genus %d r = %d: %d control points (D = 2h-2), %d non-zero; non-zero components: %s\n", $h, $r, $nC, $nzC, join(' ', @ctrl);
# (M)
my ($nM, $bM) = (0, 0);
for my $p (sorted_tuples($r, 0, $N)) { my $J = 0; $J += $_ for @$p; next if $J < 1; my @nz = grep { ev($_, [map { Q($_) } @$p]) != 0 } 0..$h; $nM++; if (@nz) { $bM++; print "M-VIOLATION R = {@$p}: components @nz\n" } }
printf "(M) genus %d r = %d: %d lattice points with parts >= 0, 1 <= J, D <= 2h-3: %d violations\n", $h, $r, $nM, $bM;
# (H)
my %H;
for my $p (sorted_tuples($r, 2, 2*$h + $r)) { my $D = 0; $D += $_ - 1 for @$p; next if $D < 2*$h - 2; my $z = ev(0, [map { Q($_) } @$p]) == 0 ? 1 : 0; $H{"$D|$z"}++; }
for my $D (2*$h - 2 .. 2*$h) { printf "(H) genus %d r = %d, D = %d: Hurwitz component zero at %d points, non-zero at %d\n", $h, $r, $D, $H{"$D|1"} // 0, $H{"$D|0"} // 0; }
