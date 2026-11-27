#!/usr/bin/perl
# Xi_g(j_1..j_r) and its Hodge components from the psi table psi_g<g>.txt (no q table needed):
#   Xi_g(R) = Sum_m 1/m! < Omega prod_i lam_{j_i}(psi_i) prod_k T(psi_{r+k}) >_{g, r+m}   (c = 1 normalisation)
# lam_j = 1/(1 - j psi), T = -psi^2/(1 - psi). The table gives <Omega_c prod psi^a>_{g,n} in the c-normalisation;
# <Omega_{(1,tau,-1-tau)} prod psi^a> = c^(3g-3+n-|a|) <Omega_c prod psi^a>. Components as in components.pl.
# Output comp_g<g>_r<r>_psi.txt (same format as comp_g<g>_r<r>.txt).  Usage: perl dilsum.pl g r
use strict; use warnings;
use FindBin; use lib $FindBin::Bin;
use polylib; use Math::BigRat;
my ($g, $r) = @ARGV; my $dir = "$FindBin::Bin/../data";
open(my $fh, '<', "$dir/psi_g$g.txt") or die; my @lines = <$fh>; close $fh;
my %tab; my $nmax = 0;
for (@lines) { chomp; next unless /\|/; my @f = split /\s*\|\s*/; my $n = $f[1]; my $astr = $f[2]; $astr =~ s/[{}\s]//g;
  my @a = $n ? (sort { $a <=> $b } split /,/, $astr) : (); $tab{join(',', @a)} = $f[3]; $nmax = $n if $n > $nmax; }
my $dim = 3 * $g - 3 + $r; my $mmax = $dim; # T has degree >= 2: need 2m <= dim + m
my @xs = map { Math::BigRat->new($_) } ('2', '3', '5', '7', '11', '13', '17', '19');
my $np = $g + 2; my (@t, @P);
my @jv = map { "j$_" } 1 .. $r;
for my $ix (0 .. $np - 1) {
  my $x = $xs[$ix]; my $c = $x + 1; my $sigma = -$x / ($c * $c);
  my %val; # cache of <Omega_{c=1} psi^a> at this x
  my $get = sub { my @a = sort { $a <=> $b } @_; my $k = join(',', @a); my $n = scalar @a; my $A = 0; $A += $_ for @a;
    return Math::BigRat->new(0) if $A != 3 * $g - 3 + $n - 0 && 0; # (Omega inhomogeneous: any A <= dim allowed)
    return $val{$k} if exists $val{$k}; die "missing table entry {$k} (n = $n > nmax = $nmax?)" unless exists $tab{$k};
    my $p = parse_expr($tab{$k}, { x => $x }); my $v = pcoef($p, ''); $v *= $c ** (3 * $g - 3 + $n - $A); return $val{$k} = $v; };
  my $res = {};
  # enumerate m, then the psi-exponents of the dilaton leaves d_k >= 2, then the probe exponents a_i >= 0
  for my $m (0 .. $mmax) { last if 2 * $m > $dim + $m; my $n = $r + $m; my $D = 3 * $g - 3 + $n;
    my $fact = Math::BigRat->new(1); $fact *= $_ for 1 .. $m;
    my @ds = ([]); for (1 .. $m) { my @nd; for my $d (@ds) { my $s = 0; $s += $_ for @$d; for (my $e = ($d->[-1] // 2); $s + $e <= $D; $e++) { push @nd, [@$d, $e]; } } @ds = @nd; }
    for my $d (@ds) { my $sd = 0; $sd += $_ for @$d; # multiplicity of this sorted d-tuple among ordered tuples
      my %cnt; $cnt{$_}++ for @$d; my $mult = Math::BigRat->new($fact); $mult /= Math::BigRat->new(1) * do { my $q = 1; $q *= f($_) for values %cnt; $q } ;
      my $wT = Math::BigRat->new((-1) ** $m);  # coefficient of psi^d in T = -psi^2/(1-psi) is -1 for d >= 2
      # probes: all a with sum a <= D - sd
      my @as = ([]); for (1 .. $r) { my @na; for my $a (@as) { my $s = 0; $s += $_ for @$a; for my $e (0 .. $D - $sd - $s) { push @na, [@$a, $e]; } } @as = @na; }
      for my $a (@as) { my @all = (@$a, @$d); my $v = $get->(@all); next if $v->is_zero; my %e; @e{@jv} = @$a;
        my $mono = mono_key(\%e); my $cf = $v * $wT * $mult / $fact;
        $res = padd($res, { $mono => $cf }); } } }
  push @t, Math::BigRat->new(-1) / $sigma; push @P, $res;
}
sub f { my $n = shift; my $q = 1; $q *= $_ for 1 .. $n; return $q; }
my $n = $g + 1; my @A; my @B;
for my $i (0 .. $n - 1) { $A[$i] = [ map { $t[$i] ** $_ } 0 .. $g ]; $B[$i] = $P[$i]; }
for my $col (0 .. $n - 1) { my $piv = $col; $piv++ while $A[$piv][$col]->is_zero; @A[$col, $piv] = @A[$piv, $col]; @B[$col, $piv] = @B[$piv, $col];
  my $inv = Math::BigRat->new(1) / $A[$col][$col]; $A[$col] = [ map { $_ * $inv } @{$A[$col]} ]; $B[$col] = pscale($B[$col], $inv);
  for my $row (0 .. $n - 1) { next if $row == $col; my $fct = $A[$row][$col]; next if $fct->is_zero;
    $A[$row] = [ map { $A[$row][$_] - $fct * $A[$col][$_] } 0 .. $n - 1 ]; $B[$row] = psub($B[$row], pscale($B[$col], $fct)); } }
my $chk = {}; for my $k (0 .. $g) { $chk = padd($chk, pscale($B[$k], $t[$n] ** $k)); }
print STDERR "genus $g r $r (psi table): consistency check at extra x: ", (pis_zero(psub($chk, $P[$n])) ? "OK" : "FAIL"), "\n";
open(my $out, '>', "$FindBin::Bin/comp_g${g}_r${r}_psi.txt") or die;
for my $k (0 .. $g) { print $out "$g | $r | $k | ", pstr($B[$k]), "\n"; } close $out;
