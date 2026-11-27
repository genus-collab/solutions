# irr_test.pl hmax : "irreducible-boundary localization" of single Chern-character insertions for the Hurwitz class.
# For odd l, Mumford: ch_l(E) = B_{l+1}/(l+1)! [ kappa_l - psi^l + (1/2) iota_*( Sum_a (-1)^a psi'^a psi''^{l-1-a} ) ].
# Q_l(h,j) := Int_{Mbar_{h,1}} c(T(-log)) ch_l lam_j  -  (B_{l+1}/(l+1)!) (1/2) Sum_a (-1)^a [y2^a y3^{l-1-a}] Xi^(0)_{h-1}(j,y2,y3)
# (the second term is the irreducible-boundary part).  For l = 1, Q_1 = 0 on 2 <= j <= 2h-2 is the lambda_1 theorem.
# Prints Q_l(h,j) for 2 <= j <= 2h.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Pair; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my $hmax = shift // 6;
use File::Basename (); use File::Spec ();
my $E = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../E_formulas";
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
sub load3 { my $g = shift; my $f = -e "$E/engine_h${g}_r3.txt" ? "$E/engine_h${g}_r3.txt" : "$E/referee/ref_engine_h${g}_r3.txt"; open my $fh, '<', $f or die "no $f"; my %c;
  while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 0; for my $t (split /\s*;\s*/, $rest) { my ($m, $v) = split /:/, $t; $c{$m} = Q($v); } } close $fh; return \%c; }
for my $h (2..$hmax) {
  my $L1 = Pair::lamvee($h, 1); my $T3 = load3($h - 1);
  for my $l (grep { $_ % 2 } 1..(2*$h - 1)) { next unless grep { $_ == $l } @RefHodge::CHL; next if $l > $h && $l > 1;   # ch_l with l <= h available in the moments
    my $P = Pair::jpoly($h, Pair::mul($h, $L1, Pair::chm($l)));
    my $Bc = PJbern($l + 1) / fact($l + 1);
    my @out;
    for my $j (2..2*$h) { my $v = Pair::peval($P, Q($j)); my $irr = Q(0);
      for my $a (0..$l-1) { my $b = $l - 1 - $a;
        for my $m (keys %$T3) { my ($e1, $e2, $e3) = split /,/, $m; next unless $e2 == $a && $e3 == $b; $irr += (($a % 2) ? -1 : 1) * $T3->{$m} * Q($j) ** $e1; } }
      my $q = $v - $Bc * $irr / 2; push @out, "j=$j:" . ($q == 0 ? '0' : $q); }
    printf "h=%d l=%d: %s\n", $h, $l, join('  ', @out); } }
sub PJbern { my $n = shift; my @B = (Q(1)); for my $m (1..$n) { my $s = Q(0); for my $k (0..$m-1) { $s += binom($m + 1, $k) * $B[$k] } push @B, -$s / ($m + 1) } return $B[$n] }
sub binom { my ($n, $k) = @_; my $r = Q(1); for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r }
