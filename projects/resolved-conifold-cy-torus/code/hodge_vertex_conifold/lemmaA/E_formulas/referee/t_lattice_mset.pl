# t_lattice_mset.pl h r : as t_lattice.pl, for the compact multiset files ref_mset_h<h>_r<r>.txt written by COMPACT=1 perl t_moments.pl h r
# (value at a point j = Sum_A coef_A m_A(j), m_A = monomial symmetric polynomial of the exponent multiset A).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my ($h, $r) = @ARGV; my $f = "$FindBin::Bin/ref_mset_h${h}_r${r}.txt";
my @C; open my $fh, '<', $f or die "no $f"; while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $k, $t) = split /\s*\|\s*/, $line, 4; my ($A, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $A], Q($cf)]; } close $fh;
my %mAc;
sub ev { my ($k, $p) = @_; my $s = Q(0); my $pk = join(',', @$p); for my $t (@{$C[$k] // []}) { my ($A, $cf) = @$t; my $key = $pk . '|' . join(',', @$A);
    my $w = $mAc{$key} //= RefHodge::mA([map { Q($_) } @$p], $A); $s += $cf * $w if $w != 0; } return $s; }
sub sorted_tuples { my ($n, $lo, $Jmax) = @_; my @o; my $rec; $rec = sub { my ($pre, $min, $J) = @_; if (@$pre == $n) { push @o, [@$pre]; return } for (my $j = $min; $J + $j <= $Jmax; $j++) { $rec->([@$pre, $j], $j, $J + $j) } }; $rec->([], $lo, 0); return @o; }
my $N = 2*$h - 3 + $r;
my ($nA, $bA, $nC, $nzC) = (0, 0, 0, 0); my @ctrl;
for my $p (sorted_tuples($r, 2, 2*$h - 2 + $r)) { my $D = 0; $D += $_ - 1 for @$p; my @nz = grep { ev($_, $p) != 0 } 0..$h;
  if ($D <= 2*$h - 3) { $nA++; if (@nz) { $bA++; print "VIOLATION R = {@$p}: non-zero components @nz\n" } } else { $nC++; $nzC++ if @nz; push @ctrl, "{" . join(',', @$p) . "}:" . (@nz == $h + 1 ? 'all' : join('', @nz)); } }
printf "(A) genus %d r = %d: %d Lemma-A points (parts >= 2, D <= 2h-3), %d violations\n", $h, $r, $nA, $bA;
printf "(C) genus %d r = %d: %d control points (D = 2h-2), %d non-zero; non-zero components: %s\n", $h, $r, $nC, $nzC, join(' ', @ctrl);
my ($nM, $bM) = (0, 0);
if (!$ENV{NOM}) { for my $p (sorted_tuples($r, 0, $N)) { my $J = 0; $J += $_ for @$p; next if $J < 1; my @nz = grep { ev($_, $p) != 0 } 0..$h; $nM++; if (@nz) { $bM++; print "M-VIOLATION R = {@$p}: components @nz\n" } }
  printf "(M) genus %d r = %d: %d lattice points with parts >= 0, 1 <= J, D <= 2h-3: %d violations\n", $h, $r, $nM, $bM; }
my %H;
for my $p (sorted_tuples($r, 2, 2*$h + $r)) { my $D = 0; $D += $_ - 1 for @$p; next if $D < 2*$h - 2; my $z = ev(0, $p) == 0 ? 1 : 0; $H{"$D|$z"}++; }
for my $D (2*$h - 2 .. 2*$h) { printf "(H) genus %d r = %d, D = %d: Hurwitz component zero at %d points, non-zero at %d\n", $h, $r, $D, $H{"$D|1"} // 0, $H{"$D|0"} // 0; }
