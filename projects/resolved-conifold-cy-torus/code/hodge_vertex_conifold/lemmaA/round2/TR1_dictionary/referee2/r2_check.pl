#!/usr/bin/perl
# r2_check.pl MODE TAU S GMAX [E0] [XMAX]     (referee 2; uses only R2TR.pm, RefHodge.pm/RefWK.pm and the E_formulas tables)
#  MODE p2   : y = -1/w. Compares (a) every xi-coefficient with sigma^(1-g) X_{g,n}[k] (Route E engine + Route E referee engine
#              tables, genus-0 formula) and (b) END-TO-END in x-form, [x^mu]F_{g,n} (Lagrange inversion of the global forms E_l)
#              against sigma^(1-g) prod Pi_{mu_i} Xi_g(mu), for all mu with 1 <= mu_i and |mu| <= XMAX (default 2g-1+n).
#  MODE bkmp : y = log(1+w). Compares every xi-coefficient with sigma^(g-1+n) M_{g,n}[k], M = Int_{Mbar_{g,n}} Omega prod psi^k,
#              computed with the Route E referee's TR-free engine RefHodge at this numeric tau (genus 0: (n-3)!/prod k_i!).
#  MODE var  : dual numbers, y = log(1+w) + (E0 + eps)(-1/w). Tests the variational formula of Lemma 2 + the pairing of Lemma 3:
#              d/de C_{g,n}[k] = Sum_l C_{g,n+1}[k,l] v_l,  v_0 = 0, v_l = -1 (l >= 1), for every (g,n) with (g,n+1) computed.
#  MODE vbase: dual numbers at the base curve itself, y = (-1/w) + eps(-1/w) (homogeneity test d/de = -(2g-2+n) omega and the
#              dilaton-type identity Sum_l C_{g,n+1}[k,l] v_l = -(2g-2+n) C_{g,n}[k] at the dilaton point).
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use R2TR;
sub R { Math::BigRat->new(@_) }
$| = 1;
my ($mode, $tau, $S, $GMAX, $E0, $XMAX) = @ARGV; die "usage: r2_check.pl p2|bkmp|var|vbase tau S GMAX [E0] [XMAX]\n" unless defined $GMAX;
$E0 //= 0; my $t0 = time;
my $lmax = 1; for my $g (0..$GMAX) { my $n = $S - $g; next if $n < 1; my $d = 3*$g - 3 + $n; $lmax = $d if $d > $lmax; }
my $dual = ($mode eq 'var' || $mode eq 'vbase') ? 1 : 0;
R2TR::setup(tau => $tau, lmax => $lmax, dual => $dual);
my $sig = $R2TR::sigma; my $tauR = R($tau);
printf "mode %s, tau = %s, sigma = %s, S = %d, GMAX = %d, lmax = %d%s; setup %.1f s\n", $mode, $tauR->bstr, $sig->bstr, $S, $GMAX, $lmax,
  ($mode eq 'var' ? ", e0 = $E0" : ""), time - $t0;
# basis sanity: [x^m] xi_k = m^k Pi_m(tau), via xi_k = Sum_l T_kl Prim(E_l)
{ my $bad = 0; my $cnt = 0; for my $k (0..($lmax < 6 ? $lmax : 6)) { for my $m (1..7) { my $s = R(0); $s += $R2TR::T[$k][$_] * R2TR::PLm($_, $m) for 0..$k;
      $cnt++; $bad++ unless $s == R($m) ** $k * R2TR::Pi_mu($m); } }
  print "basis check [x^m] xi_k = m^k Pi_m: $cnt values, $bad failures\n"; }
my (@y, @eta);
if ($mode eq 'p2') { @y = R2TR::yseries('base'); }
elsif ($mode eq 'bkmp') { @y = R2TR::yseries('bkmp'); }
elsif ($mode eq 'var') { @y = R2TR::yseries('mix', $E0); @eta = R2TR::yseries('base'); }
elsif ($mode eq 'vbase') { @y = R2TR::yseries('base'); @eta = R2TR::yseries('base'); }
else { die "unknown mode" }
R2TR::set_y(\@y, \@eta);
my @done = R2TR::run($S, $GMAX);
printf "recursion: %d correlators, %.1f s\n", scalar(@done), time - $t0;

my $EDIR = "$FindBin::Bin/../../../E_formulas";
sub load_table { my ($f) = @_; return undef unless -e $f; open my $fh, '<', $f or die "open $f"; my %comp;
  while (my $line = <$fh>) { chomp $line; $line =~ s/\r//g; next unless $line =~ /\S/; my ($hh, $rr, $c, $body) = split /\s*\|\s*/, $line, 4;
    $comp{$c} //= {}; next unless defined $body && $body =~ /\S/;
    for my $t (split /\s*;\s*/, $body) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; $comp{$c}{$m} = R($v); } }
  close $fh; return \%comp; }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; return $r; }
sub tuples { my ($n, $maxsum) = @_; my @out; my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x (0..$left) { $rec->([@$pre, $x], $left - $x); } }; $rec->([], $maxsum); return @out; }
sub comps { my ($n, $tot) = @_; return () if $n == 0; return ([$tot]) if $n == 1; my @o; for my $a (1..$tot-$n+1) { push @o, [$a, @$_] for comps($n-1, $tot-$a); } return @o; }
sub g0X { my ($n, $k) = @_; my $K = 0; $K += $_ for @$k; return R(0) if $n < 3 || $K > $n - 3;
  my @c = (R(1)); for my $i (1..$n-3) { my @d = map { R(0) } 0..@c; for my $e (0..$#c) { $d[$e+1] += $c[$e]; $d[$e] -= $i * $c[$e]; } @c = @d; }
  my $m = fact($K); $m /= fact($_) for @$k; return ($c[$K] // R(0)) * $m; }
my ($TOT, $BAD) = (0, 0);
my %XI;   # xi-coefficients (real part)
for my $gn (keys %R2TR::W) { $XI{$gn} = R2TR::to_xi($R2TR::W{$gn}, 0); }

if ($mode eq 'p2') {
  for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %XI) {
    my ($g, $n) = split /,/, $gn; my $Xc = $XI{$gn}; my $dim = 3*$g - 3 + $n; my @lines;
    my ($symbad, $outside) = (0, 0); my %s;
    for my $key (keys %$Xc) { my @kk = split /,/, $key; my $sum = 0; $sum += $_ for @kk; $outside++ if $sum > $dim;
      push @{$s{join(",", sort { $a <=> $b } @kk)}}, $Xc->{$key}; }
    for my $sk (keys %s) { my @v = @{$s{$sk}}; for (@v) { $symbad++ if $_ != $v[0]; } }
    my @src; if ($g == 0) { push @src, ['genus0', undef]; }
    else { for my $f (["engine", "$EDIR/engine_h${g}_r${n}.txt"], ["refengine", "$EDIR/referee/ref_engine_h${g}_r${n}.txt"]) { my $T = load_table($f->[1]); push @src, [$f->[0], $T] if $T; } }
    my %Xhodge;
    for my $src (@src) { my ($name, $T) = @$src; my ($cmp, $bad) = (0, 0);
      for my $k (tuples($n, $dim)) { my $key = join(",", @$k); my $X;
        if ($name eq 'genus0') { $X = g0X($n, $k); } else { $X = R(0); my $ms = R(-1) / $sig; for my $c (keys %$T) { my $v = $T->{$c}{$key}; $X += $v * $ms ** $c if defined $v; } }
        $Xhodge{$key} //= $X; my $got = $Xc->{$key} // R(0); $cmp++; $bad++ if $got != $sig ** (1 - $g) * $X; }
      push @lines, "xi-form vs $name: $cmp compared, $bad mismatches"; $TOT += $cmp; $BAD += $bad; }
    if (%Xhodge) {   # end-to-end x-form
      my $xm = defined $XMAX ? $XMAX : 2*$g + 1 + $n; my ($cmp, $bad) = (0, 0); my $Wc = $R2TR::W{$gn};
      for my $tot ($n..$xm) { for my $mu (comps($n, $tot)) {
          my $lhs = R(0); for my $key (keys %$Wc) { my @Lv = split /,/, $key; my $t = $Wc->{$key}; for my $i (0..$n-1) { $t *= R2TR::PLm($Lv[$i], $mu->[$i]); } $lhs += $t; }
          my $xi = R(0); for my $key (keys %Xhodge) { my @kk = split /,/, $key; my $t = $Xhodge{$key}; next if $t->is_zero; for my $i (0..$n-1) { $t *= R($mu->[$i]) ** $kk[$i]; } $xi += $t; }
          my $rhs = $sig ** (1 - $g) * $xi; $rhs *= R2TR::Pi_mu($_) for @$mu;
          $cmp++; $bad++ if $lhs != $rhs; } }
      push @lines, "x-form [x^mu]F vs sigma^(1-g) prod Pi Xi_g(mu), |mu| <= $xm: $cmp compared, $bad mismatches"; $TOT += $cmp; $BAD += $bad; }
    printf "(%d,%d): %d xi-coefficients, symmetry violations %d, outside %d; %s\n", $g, $n, scalar(keys %$Xc), $symbad, $outside, join("; ", @lines);
    $BAD += $symbad + $outside; } }
elsif ($mode eq 'bkmp') {
  require RefHodge;
  my %OMC;
  for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %XI) {
    my ($g, $n) = split /,/, $gn; my $Xc = $XI{$gn}; my $dim = 3*$g - 3 + $n; my ($cmp, $bad) = (0, 0); my $tg = time;
    my $om; if ($g >= 1) { $om = $OMC{$g} //= RefHodge::omega($g, [R(1), $tauR, -1 - $tauR], ($g == 1 ? 1 : 3*$g - 3)); }
    for my $k (tuples($n, $dim)) { my $key = join(",", @$k); my $M;
      if ($g == 0) { my $K = 0; $K += $_ for @$k; $M = R(0); if ($K == $n - 3) { $M = fact($n - 3); $M /= fact($_) for @$k; } }
      else { $M = R(0); for my $mu (keys %$om) { my @m = split /,/, $mu; my $v = RefHodge::I($g, \@m, $k); $M += $om->{$mu} * $v unless $v->is_zero; } }
      my $got = $Xc->{$key} // R(0); $cmp++; $bad++ if $got != $sig ** ($g - 1 + $n) * $M; }
    printf "(%d,%d): BKMP vs RefHodge M_{g,n}[k]: %d compared, %d mismatches (%.1f s)\n", $g, $n, $cmp, $bad, time - $tg; $TOT += $cmp; $BAD += $bad; } }
else {   # var / vbase
  my %XE; for my $gn (keys %R2TR::W) { $XE{$gn} = R2TR::to_xi($R2TR::W{$gn}, 1); }
  for my $gn (sort { my ($ga,$na) = split /,/, $a; my ($gb,$nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %XI) {
    my ($g, $n) = split /,/, $gn; my $up = $g . "," . ($n + 1); next unless exists $XI{$up};
    my %pair; for my $key (keys %{$XI{$up}}) { my @kk = split /,/, $key; my $l = pop @kk; next if $l == 0; my $k = join(",", @kk);
      if (exists $pair{$k}) { $pair{$k} -= $XI{$up}{$key}; } else { $pair{$k} = -$XI{$up}{$key}; } }
    my %keys = map { $_ => 1 } (keys %pair, keys %{$XE{$gn}}); my ($cmp, $bad) = (0, 0); my ($cmp2, $bad2) = (0, 0);
    for my $k (keys %keys) { my $d = $XE{$gn}{$k} // R(0); my $p = $pair{$k} // R(0); $cmp++; $bad++ if $d != $p;
      if ($mode eq 'vbase') { my $c = $XI{$gn}{$k} // R(0); $cmp2++; $bad2++ if $d != -(2*$g - 2 + $n) * $c; } }
    my $extra = $mode eq 'vbase' ? sprintf("; homogeneity d/de = -(2g-2+n) omega: %d/%d mismatches", $bad2, $cmp2) : "";
    printf "(%d,%d): d/de C[k] vs Int Lambda omega_{g,n+1}: %d coefficients, %d mismatches%s\n", $g, $n, $cmp, $bad, $extra;
    $TOT += $cmp + $cmp2; $BAD += $bad + $bad2; } }
printf "TOTAL: %d comparisons, %d failures; time %.1f s\n", $TOT, $BAD, time - $t0;
