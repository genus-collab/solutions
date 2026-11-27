package ref_data;
# Referee data access (own parser): Xi_g(T) components from Route E's exact tables E_formulas/engine_h<g>_r<r>.txt,
# fallback Route E's point engine (read-only use of E_formulas/Hodge.pm, HodgePt.pm), extra values from files (env EXTRA).
# xi_num(g, T, tau): numeric value at rational tau (c = 1), unit parts handled by the string equation, genus 0 closed form.
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/../../../E_formulas"; use Math::BigRat; use Hodge; use HodgePt;
sub R { Math::BigRat->new(@_) }
my $EDIR = "$FindBin::Bin/../../../E_formulas";
our %RTAB = (1 => 3, 2 => 4, 3 => 4, 4 => 5, 5 => 5, 6 => 3);
my %EXTRA;
if (defined $ENV{EXTRA}) { for my $f (split /;/, $ENV{EXTRA}) { open(my $fh, '<', $f) or die "$f"; while (my $l = <$fh>) { $l =~ s/\r//g; chomp $l; next if $l =~ /^#/;
      my @p = split /\s*\|\s*/, $l; next unless @p == 3; $EXTRA{"$p[0]|" . join(',', sort { $a <=> $b } split /,/, $p[1])} = [ map { R($_) } split /\s*;\s*/, $p[2] ]; } close $fh; } }
my %tab; my %xc; my $eng_g = -1;
sub load { my ($g, $r) = @_; my $key = "$g|$r"; return $tab{$key} if exists $tab{$key};
  my $f = "$EDIR/engine_h${g}_r${r}.txt"; if (!-e $f || $r > ($RTAB{$g} // 0)) { return $tab{$key} = 0; }
  open(my $fh, '<', $f) or die $f; my @rows;
  while (my $l = <$fh>) { $l =~ s/\r//g; chomp $l; my @p = split /\s*\|\s*/, $l, 4; next unless @p == 4;
    for my $t (split /\s*;\s*/, $p[3]) { next unless $t =~ /:/; my ($e, $c) = split /:/, $t; push @rows, [$p[2], [split /,/, $e], R($c)]; } }
  close $fh; return $tab{$key} = \@rows; }
sub comps { my ($g, $T) = @_; my @T = sort { $a <=> $b } @$T; my $key = "$g|" . join(',', @T); return $xc{$key} if exists $xc{$key};
  return $xc{$key} = $EXTRA{$key} if exists $EXTRA{$key};
  my $r = scalar @T; my $tb = load($g, $r); my @c;
  if ($tb) { for my $row (@$tb) { my ($k, $e, $c) = @$row; my $m = $c->copy; for my $i (0..$r-1) { $m *= R($T[$i]) ** $e->[$i] if $e->[$i]; } $c[$k] = ($c[$k] // R(0)) + $m; } @c = map { $_ // R(0) } @c; }
  else { die "engine needed at genus $g for (@T): refused" if $g >= 7; if ($eng_g != $g) { Hodge::setup($g); %HodgePt::memo = (); $eng_g = $g; } @c = map { R("$_") } @{ HodgePt::WP([], [@T]) }; }
  return $xc{$key} = \@c; }
sub xi0 { my @T = @_; my $r = @T; my $J = 0; $J += $_ for @T; return R(1, $J*($J+1)) if $r == 1; return R(1, $J) if $r == 2; my $p = R(1); $p *= ($J - $_) for 1..$r-3; $p }
sub xi_num { my ($g, $T, $tau) = @_; my @T = sort { $a <=> $b } @$T; my $sig = $tau * ($tau + 1);
  return R(0) if $g < 0;
  if ($g == 0) { die "xi0 of empty" unless @T; return xi0(@T); }
  if (@T && $T[0] == 1) { shift @T; my $Js = 0; $Js += $_ for @T; if (!@T) { return $g == 1 ? -(1 + 1/$sig)/24 : R(0); } return $Js * xi_num($g, \@T, $tau); }
  die "xi_num: empty profile at genus $g" unless @T;
  my $c = comps($g, \@T); my $v = R(-1) / $sig; my $s = R(0); my $p = R(1); for my $k (0..$#$c) { $s += $c->[$k] * $p; $p *= $v; } return $s; }
1;
