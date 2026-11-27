package HodgeSide;
# Hodge-side descendant correlators of the CohFT Omega = Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau)/(-sigma)^g (c = 1),
#   (a) at the dilaton point q_0 (dilaton leaves T = -psi^2/(1-psi)):  <Omega prod psi^{k_i}>^{q_0}_{g,n}
#       = Sum_c (-1/sigma)^c [j^k] Xi^(c)_g   from Route E's engine tables (E_formulas/engine_h<h>_r<r>.txt; Mumford GRR + DVV, no TR);
#       genus 0: [j^k] prod_{i=1}^{n-3} (J - i)  (Route A, Lemma 3.1: Xi_0(T) = Gamma(J)/Gamma(J-r+3));
#   (b) without dilaton leaves (Marino-Vafa integrals) and along the family K_s (Route A's Q_h(s;R)) from the Mathematica
#       tables lemmaA/data/mv_g<g>_l<l>.txt, q_g<g>_r<r>.txt (weights (c,-x,-1), c = 1+x = 1/(1+tau), x = -tau/(1+tau);
#       conversion to c = 1: multiply by c^(3g-3+2r)).
use strict; use warnings;
use FindBin;
use lib $FindBin::Bin;
use Math::BigRat;
use MParse; use Poly;
sub R { Math::BigRat->new(@_) }
our $EDIR = "$FindBin::Bin/../../E_formulas";
our $DDIR = "$FindBin::Bin/../../data";
my %eng;
sub engine { my ($h, $r) = @_; my $key = "$h,$r"; return $eng{$key} if exists $eng{$key};
  my $f = "$EDIR/engine_h${h}_r${r}.txt"; return undef unless -e $f;
  open my $fh, '<', $f or die "open $f"; my %comp;
  while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $c, $body) = split /\s*\|\s*/, $line, 4;
    die "bad line in $f" unless $hh == $h && $rr == $r; $comp{$c} = {};
    next unless defined $body && $body =~ /\S/;
    for my $t (split /\s*;\s*/, $body) { my ($m, $v) = split /:/, $t; $comp{$c}{$m} = R($v); } }
  close $fh; $eng{$key} = \%comp; return \%comp; }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; return $r; }
# <Omega prod psi^{k_i}>^{q_0}_{g,n} at the given sigma
sub hodge_q0 { my ($g, $n, $k, $sigma) = @_; my $K = 0; $K += $_ for @$k;
  if ($g == 0) { return R(0) if $n < 3; return R(0) if $K > $n - 3;
    # prod_{i=1}^{n-3} (J - i) = Sum_d c_d J^d
    my @c = (R(1)); for my $i (1..$n-3) { my @d = map { R(0) } 0..@c; for my $e (0..$#c) { $d[$e+1] += $c[$e]; $d[$e] -= $i * $c[$e]; } @c = @d; }
    my $mult = fact($K); $mult /= fact($_) for @$k; return ($c[$K] // R(0)) * $mult; }
  my $E = engine($g, $n) or die "no engine table for g=$g n=$n";
  my $m = join(',', @$k); my $s = R(0); my $ms = R(-1) / $sigma;
  for my $c (sort { $a <=> $b } keys %$E) { my $v = $E->{$c}{$m}; $s += $v * ($ms ** $c) if defined $v; }
  return $s; }
sub have_q0 { my ($g, $n) = @_; return 1 if $g == 0; return defined engine($g, $n); }
# Mathematica data tables: expression for (g, l) as a Poly after substituting x
my %dat;
sub datapoly { my ($kind, $g, $l, $tau) = @_;   # kind = 'mv' or 'q'
  my $key = "$kind,$g,$l,$tau"; return $dat{$key} if exists $dat{$key};
  my $f = $kind eq 'mv' ? "$DDIR/mv_g${g}_l${l}.txt" : "$DDIR/q_g${g}_r${l}.txt"; return undef unless -e $f;
  open my $fh, '<', $f or die; local $/; my $txt = <$fh>; close $fh; $txt =~ s/\r//g;
  my @lines = grep { /\S/ } split /\n/, $txt; my $expr;
  for my $ln (@lines) { my ($gg, $ll, $e) = split /\s*\|\s*/, $ln, 3; if ($gg == $g && $ll == $l) { $expr = $e; last; } }
  die "no expression in $f" unless defined $expr;
  my $tauR = R($tau); my $x = -$tauR / (1 + $tauR); my $c = 1 + $x;
  my $P = MParse::parse($expr, { x => $x });
  $P = Poly::scale($P, $c ** (3*$g - 3 + 2*$l));
  $dat{$key} = $P; return $P; }
# coefficient of prod v_i^{k_i} * s^m in a Poly (variables named <pre>1.. and s)
sub polycoef { my ($P, $pre, $k, $m) = @_; my %e; for my $i (0..$#$k) { $e{$pre . ($i+1)} = $k->[$i] if $k->[$i]; } $e{s} = $m if $m;
  my $mk = Poly::mono_key(\%e); return exists $P->{$mk} ? $P->{$mk} : R(0); }
1;
