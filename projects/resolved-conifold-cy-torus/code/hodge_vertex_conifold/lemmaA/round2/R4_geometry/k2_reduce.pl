# k2_reduce.pl hmax : Mumford reduction of the k = 2 component (report s.6), checked against the component tables.
#   e_2(alpha+alpha^2) = lambda_2 - 6 ch_3 = lambda_1^2/2 - 6 ch_3,
#   Xi^(2)_h(j) = (1/24) < lambda_1 (kappa_1 - psi + delta) >  +  (1/120) [ <kappa_3> - <psi^3> + (1/2) < iota_*(psi'^2 - psi' psi'' + psi''^2) > ],
#   <.> = Int_{Mbar_{h,1}} c(T(-log)) (.) lam_j.  Ingredients:
#   <lambda_1 kappa_1> = (2h-1) Xi^(1)_h(j) - B1_h(j),   B1_h(j) := [y^1] Xi^(1)_h(j,y)            (two-leaf table, k = 1)
#   <kappa_3> = (2h-1) X_h - Bc_1 - Bc_2 - Bc_3,   Bc_c := [y^c] Xi^(0)_h(j,y)                   (two-leaf table, k = 0)
#   irreducible boundary of ch_3: (1/2)( 2 [y^2 z^0] - [y z] ) Xi^(0)_{h-1}(j,y,z)              (three-leaf table, k = 0)
#   separating: Sum_{g1} ( [y^2]Xi^(0)_{g1}(j,y) X_{g2}[0] - [y^1] Xi^(0)_{g1}(j,y) X_{g2}[1] + [y^0] Xi^(0)_{g1}(j,y) X_{g2}[2] ),
#   X_g[b] = Int_{Mbar_{g,1}} c(T(-log)) psi^b.
# Also checks the two-leaf lambda_1 reduction:  12 B1_h = [y^1] < kappa_1 - psi_1 - psi_2 + delta >_{(j,y)} on Mbar_{h,2}.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $hmax = shift // 5;
use File::Basename (); use File::Spec ();
my $E = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../E_formulas";
sub loadr { my ($g, $r) = @_; my $f = "$E/engine_h${g}_r${r}.txt"; $f = "$E/referee/ref_engine_h${g}_r${r}.txt" unless -e $f; $f = "$FindBin::Bin/refdata/ref_engine_h${g}_r1.txt" if $r == 1;
  open my $fh, '<', $f or die "no $f"; my @C; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; $rest //= ''; my %c;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; $c{$m} = Q($v); } $C[$k] = \%c; } close $fh; return \@C; }
my (%R1, %R2, %R3); for my $g (1..$hmax) { $R1{$g} = loadr($g, 1); $R2{$g} = loadr($g, 2); $R3{$g} = loadr($g, 3) if $g <= $hmax - 1; }
# evaluators
sub e1 { my ($g, $k, $j) = @_; my $c = $R1{$g}[$k] or return Q(0); my $s = Q(0); for my $m (keys %$c) { $s += $c->{$m} * Q($j) ** $m } return $s }      # Xi^(k)_g(j)
sub e1c { my ($g, $k, $b) = @_; my $c = $R1{$g}[$k] or return Q(0); return $c->{$b} // Q(0) }                                                       # [j^b] Xi^(k)_g
sub e2y { my ($g, $k, $j, $c) = @_; my $C = $R2{$g}[$k] or return Q(0); my $s = Q(0); for my $m (keys %$C) { my ($a, $b) = split /,/, $m; next unless $b == $c; $s += $C->{$m} * Q($j) ** $a } return $s }   # [y^c] Xi^(k)_g(j,y)
sub e3yz { my ($g, $j, $cy, $cz) = @_; my $C = $R3{$g}[0] or return Q(0); my $s = Q(0); for my $m (keys %$C) { my ($a, $b, $d) = split /,/, $m; next unless $b == $cy && $d == $cz; $s += $C->{$m} * Q($j) ** $a } return $s }
my $bad = 0; my $n = 0; my $bad2 = 0;
for my $h (2..$hmax) {
  for my $j (1..3*$h) { my $jq = Q($j);
    my $X = sub { my ($g, $x) = @_; return e1($g, 0, $x) };
    # (a) < lambda_1 (kappa_1 - psi + delta) >
    my $A = (2*$h - 1) * e1($h, 1, $j) - e2y($h, 1, $j, 1) - (e1($h, 1, $j) - e1($h, 1, 0)) / $jq;
    $A += Q(1)/2 * ($jq - 2*$h + 2) * ($jq - 2*$h + 3) * e1($h - 1, 1, $j) if $h >= 2 && $h - 1 >= 1;
    for my $g1 (1..$h-1) { my $g2 = $h - $g1; $A += ($jq - 2*$g1 + 1) * (e1($g1, 1, $j) * e1($g2, 0, 0) + e1($g1, 0, $j) * e1($g2, 1, 0)); }
    # (b) [ <kappa_3> - <psi^3> + boundary ]
    my $k3 = (2*$h - 1) * $X->($h, $j) - e2y($h, 0, $j, 1) - e2y($h, 0, $j, 2) - e2y($h, 0, $j, 3);
    my $p3 = Q(0); for my $a (3..3*$h) { $p3 += e1c($h, 0, $a) * $jq ** ($a - 3); }
    my $irr = 0; if ($h - 1 >= 1) { $irr = Q(1)/2 * (2 * e3yz($h - 1, $j, 2, 0) - e3yz($h - 1, $j, 1, 1)); }
    my $sep = Q(0); for my $g1 (1..$h-1) { my $g2 = $h - $g1; $sep += e2y($g1, 0, $j, 2) * e1c($g2, 0, 0) - e2y($g1, 0, $j, 1) * e1c($g2, 0, 1) + e2y($g1, 0, $j, 0) * e1c($g2, 0, 2); }
    my $Bk = $k3 - $p3 + $irr + $sep;
    my $pred = $A / 24 + $Bk / 120;
    $n++; $bad++ if $pred != e1($h, 2, $j);
    # two-leaf lambda_1 reduction: 12 B1 = [y^1] of < kappa_1 - psi_1 - psi_2 + delta >_(j,y)
    my $k1 = 2*$h * e2y($h, 0, $j, 1) - ($h <= $hmax ? e3yz($h, $j, 1, 1) : 0);   # [y^1]( 2h Xi(j,y) - [z^1] Xi(j,y,z) )
    my $ps1 = (e2y($h, 0, $j, 1) - (Q(1) * e1($h, 0, 0) + (1 - 2*$h) * e1c($h, 0, 1))) / $jq;   # [y^1] (Xi(j,y) - Xi(0,y))/j, Xi(0,y) = (y-2h+1) X_h(y)
    my $ps2 = e2y($h, 0, $j, 2);
    my $dl = Q(0);   # [y^1] delta-terms
    if ($h - 1 >= 1) { # (1/2)(j+y-2h+1)(j+y-2h+2) Xi_{h-1}(j,y)
      my ($a0, $a1) = (e2y($h - 1, 0, $j, 0), e2y($h - 1, 0, $j, 1)); my ($u, $v) = ($jq - 2*$h + 1, $jq - 2*$h + 2);
      $dl += Q(1)/2 * ( ($u + $v) * $a0 + $u * $v * $a1 ); }
    $dl += e1($h, 0, 0) * 0;   # rational tail with both points (g1 = 0): constant in y -> no y^1 term
    for my $g1 (1..$h-1) { my $g2 = $h - $g1;   # both points on the g1 side: (j+y-2g1) Xi_{g1}(j,y) X_{g2}(0)
      $dl += e1($g2, 0, 0) * ( e2y($g1, 0, $j, 0) + ($jq - 2*$g1) * e2y($g1, 0, $j, 1) ); }
    for my $g1 (1..$h-1) { my $g2 = $h - $g1;   # x1 on g1, x2 on g2: (j-2g1+1) X_{g1}(j) * (y-2g2+1) X_{g2}(y)  -> [y^1]
      $dl += ($jq - 2*$g1 + 1) * e1($g1, 0, $j) * ( e1($g2, 0, 0) + (1 - 2*$g2) * e1c($g2, 0, 1) ); }
    my $pB1 = ($k1 - $ps1 - $ps2 + $dl) / 12;
    $bad2++ if $h <= $hmax - 1 && $pB1 != e2y($h, 1, $j, 1); }
  printf "genus %d: k=2 Mumford reduction vs table: %s;  two-leaf lambda_1 reduction (B1): %s\n", $h, ($bad ? "$bad mismatches so far" : "ok"), ($h <= $hmax - 1 ? ($bad2 ? "$bad2 mismatches so far" : "ok") : "skipped (no r=3 table)"); }
