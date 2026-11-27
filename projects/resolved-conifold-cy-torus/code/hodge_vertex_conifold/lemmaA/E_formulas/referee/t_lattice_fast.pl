# t_lattice_fast.pl h r : as t_lattice_mset.pl (compact files ref_mset_h<h>_r<r>.txt), with the monomial symmetric weights computed by
# grouping equal leaves: m_A(R) = Sum over splittings A = B_1 + ... + B_G (B_g of size c_g for the g-th distinct leaf value u_g)
# of prod_g (c_g! / prod_v mult_{B_g}(v)!) u_g^{|B_g|}.  Same checks: (A) Lemma A, (C) controls, (M) full lattice (skip with NOM=1),
# (H) Hurwitz layers.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($h, $r) = @ARGV; my $f = "$FindBin::Bin/ref_mset_h${h}_r${r}.txt";
my @C; open my $fh, '<', $f or die "no $f"; while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $k, $t) = split /\s*\|\s*/, $line, 4; my ($A, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $A], Q($cf)]; } close $fh;
my @F = (Q(1)); for my $n (1..40) { push @F, $F[-1] * $n; }
sub mAg { my ($groups, $A) = @_;   # groups: list of [u (BigRat), c]; A sorted list
  my %M; $M{$_}++ for @$A; my @vals = sort { $a <=> $b } keys %M;
  my $rec; $rec = sub { my ($g, $rem) = @_;   # rem: hash value -> remaining multiplicity
    my ($u, $c) = @{$groups->[$g]};
    if ($g == $#$groups) { my $n = 0; my $den = Q(1); my $deg = 0; for my $v (keys %$rem) { my $m = $rem->{$v}; next unless $m; $n += $m; $den *= $F[$m]; $deg += $v * $m; }
      return Q(0) unless $n == $c; return Q(0) if $u == 0 && $deg > 0; return $F[$c] / $den * ($deg ? $u ** $deg : Q(1)); }
    # choose a sub-multiset of size c from rem
    my @vs = grep { $rem->{$_} } sort { $a <=> $b } keys %$rem; my $tot = Q(0);
    my $choose; $choose = sub { my ($i, $left, $pick) = @_;
      if ($left == 0) { my $den = Q(1); my $deg = 0; my %r2 = %$rem; for my $v (keys %$pick) { $den *= $F[$pick->{$v}]; $deg += $v * $pick->{$v}; $r2{$v} -= $pick->{$v}; }
        return if $u == 0 && $deg > 0; my $w = $F[$c] / $den * ($deg ? $u ** $deg : Q(1)); my $rest = $rec->($g + 1, \%r2); $tot += $w * $rest if $rest != 0; return; }
      return if $i > $#vs; my $v = $vs[$i];
      for my $k (0..($rem->{$v} < $left ? $rem->{$v} : $left)) { my %p2 = %$pick; $p2{$v} = $k if $k; $choose->($i + 1, $left - $k, \%p2); } };
    $choose->(0, $c, {}); return $tot; };
  return $rec->(0, \%M); }
sub groups_of { my $p = shift; my %c; $c{$_}++ for @$p; return [ map { [Q($_), $c{$_}] } sort { $a <=> $b } keys %c ]; }
sub ev { my ($k, $p) = @_; my $g = groups_of($p); my $s = Q(0); for my $t (@{$C[$k] // []}) { my ($A, $cf) = @$t; my $w = mAg($g, $A); $s += $cf * $w if $w != 0; } return $s; }
sub sorted_tuples { my ($n, $lo, $Jmax) = @_; my @o; my $rec; $rec = sub { my ($pre, $min, $J) = @_; if (@$pre == $n) { push @o, [@$pre]; return } for (my $j = $min; $J + $j <= $Jmax; $j++) { $rec->([@$pre, $j], $j, $J + $j) } }; $rec->([], $lo, 0); return @o; }
my $N = 2*$h - 3 + $r;
my ($nA, $bA, $nC, $nzC) = (0, 0, 0, 0); my @ctrl; my %H;
for my $p (sorted_tuples($r, 2, ($ENV{ONLYA} ? 2*$h - 2 : 2*$h) + $r)) { my $D = 0; $D += $_ - 1 for @$p;   # ONLYA=1: only Lemma-A and control points
  if ($D <= 2*$h - 2) { my @nz = grep { ev($_, $p) != 0 } 0..$h;
    if ($D <= 2*$h - 3) { $nA++; if (@nz) { $bA++; print "VIOLATION R = {@$p}: non-zero components @nz\n" } } else { $nC++; $nzC++ if @nz; push @ctrl, "{" . join(',', @$p) . "}:" . (@nz == $h + 1 ? 'all' : join('', @nz)); } }
  if ($D >= 2*$h - 2 && !$ENV{ONLYA}) { my $z = ev(0, $p) == 0 ? 1 : 0; $H{"$D|$z"}++; } }
printf "(A) genus %d r = %d: %d Lemma-A points (parts >= 2, D <= 2h-3), %d violations\n", $h, $r, $nA, $bA;
printf "(C) genus %d r = %d: %d control points (D = 2h-2), %d non-zero; non-zero components: %s\n", $h, $r, $nC, $nzC, join(' ', @ctrl);
for my $D (2*$h - 2 .. 2*$h) { printf "(H) genus %d r = %d, D = %d: Hurwitz component zero at %d points, non-zero at %d\n", $h, $r, $D, $H{"$D|1"} // 0, $H{"$D|0"} // 0; }
unless ($ENV{NOM}) { my ($nM, $bM) = (0, 0);
  for my $p (sorted_tuples($r, 0, $N)) { my $J = 0; $J += $_ for @$p; next if $J < 1; my @nz = grep { ev($_, $p) != 0 } 0..$h; $nM++; if (@nz) { $bM++; print "M-VIOLATION R = {@$p}: components @nz\n" } }
  printf "(M) genus %d r = %d: %d lattice points with parts >= 0, 1 <= J, D <= 2h-3: %d violations\n", $h, $r, $nM, $bM; }
