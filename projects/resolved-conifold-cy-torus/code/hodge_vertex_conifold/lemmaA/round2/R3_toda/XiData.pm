package XiData;
# Exact values X_g(T) = Xi_g(T) (c = 1 normalisation of the brief) as v-polynomials, v = -1/sigma:
#   Xi_g(T) = Sum_k v^k Xi^(k)_g(T).
# Sources (read-only): Route E engine tables E_formulas/engine_h<h>_r<r>.txt (components as polynomials in the leaves),
# loaded lazily (only r <= $RTAB), fallback: the Route E point engine (copies Hodge.pm, HodgePt.pm, WK.pm in this folder).
# Genus 0: Xi_0(T) = (J-1)(J-2)...(J-r+3) (r >= 3), 1/J (r = 2), 1/(j(j+1)) (r = 1)  [Route A Lemma 3.1, Route B (U3),(U4)].
use strict; use warnings;
use Math::BigRat;
use FindBin;
use lib $FindBin::Bin;
use Hodge; use HodgePt;
our $EDIR = "$FindBin::Bin/../../E_formulas";
our $RTAB = $ENV{XI_RTAB} // 5;      # use tables only for r <= RTAB
our $NOTAB = $ENV{XI_NOTAB} // 0;    # 1: never use tables (engine only)
sub R { Math::BigRat->new(@_) }
our %tab;      # $tab{"$h|$r"} = [ [k, [exps], coef], ... ]  or 0 if absent
sub load_table { my ($h, $r) = @_; my $key = "$h|$r";
  return $tab{$key} if exists $tab{$key};
  my $f = "$EDIR/engine_h${h}_r${r}.txt";
  if ($NOTAB || $r > $RTAB || !-e $f) { $tab{$key} = 0; return 0; }
  open(my $fh, '<', $f) or die "cannot open $f";
  my @rows;
  while (my $line = <$fh>) {
    chomp $line; $line =~ s/\r//g; next unless $line =~ /\S/;
    my @p = split /\s*\|\s*/, $line, 4;
    next unless @p == 4;
    my ($hh, $rr, $k, $terms) = @p;
    for my $t (split /\s*;\s*/, $terms) {
      next unless $t =~ /:/;
      my ($e, $c) = split /:/, $t;
      push @rows, [$k, [split /,/, $e], R($c)];
    }
  }
  close $fh;
  $tab{$key} = \@rows; return \@rows;
}
my %cache;
my $curgen = -1;
sub engine_v {   # v-polynomial from the point engine
  my ($g, $T) = @_;
  if ($curgen != $g) { Hodge::setup($g); %HodgePt::memo = (); $curgen = $g; }
  my $v = HodgePt::WP([], [@$T]);
  return [ map { R("$_") } @$v ];
}
sub xi_v {   # v-polynomial of Xi_g(T), T = array ref of parts
  my ($g, $T) = @_;
  my @T = sort { $a <=> $b } @$T; my $r = scalar @T;
  die "xi_v: genus 0 has no v-poly here" if $g == 0;
  die "xi_v: empty profile not supported" if $r == 0;
  my $key = "$g|" . join(',', @T);
  return $cache{$key} if exists $cache{$key};
  my $res; my $tb = load_table($g, $r);
  if ($tb) {
    my @c;
    for my $row (@$tb) {
      my ($k, $e, $c) = @$row;
      my $m = $c->copy;
      for my $i (0..$r-1) { $m *= R($T[$i]) ** $e->[$i] if $e->[$i]; }
      $c[$k] = (defined $c[$k] ? $c[$k] : R(0)) + $m;
    }
    $res = [ map { defined $_ ? $_ : R(0) } @c ];
  } else {
    $res = engine_v($g, \@T);
  }
  $cache{$key} = $res; return $res;
}
sub xi0 {   # genus 0
  my ($T) = @_; my $r = scalar @$T; my $J = 0; $J += $_ for @$T;
  die "xi0: r = 0" if $r == 0;
  return R(1) / R($J * ($J + 1)) if $r == 1;
  return R(1) / R($J) if $r == 2;
  my $p = R(1); for my $i (1..$r-3) { $p *= ($J - $i); } return $p;
}
sub xi_at {   # value at tau (BigRat)
  my ($g, $T, $tau) = @_;
  return xi0($T) if $g == 0;
  my $v = xi_v($g, $T); my $sig = $tau * ($tau + 1); my $vv = R(-1) / $sig;
  my $s = R(0); my $p = R(1);
  for my $k (0..$#$v) { $s += $v->[$k] * $p; $p *= $vv; }
  return $s;
}
1;
