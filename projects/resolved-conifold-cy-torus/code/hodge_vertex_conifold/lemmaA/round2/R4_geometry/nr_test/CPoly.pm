package CPoly;
# Load engine component files as Poly objects in variables j1..jr; helpers for structural tests.
use strict; use warnings; use Math::BigRat; use Poly;
sub R { Math::BigRat->new(shift) }
# returns arrayref [Xi^(0), ..., Xi^(h)] of Poly in variables j1..jr
sub load { my ($dir, $h, $r) = @_; my $f = "$dir/engine_h${h}_r${r}.txt"; open my $fh, '<', $f or die "no $f"; my @C;
  while (<$fh>) { chomp; next unless /\S/; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $_, 4; $rest //= ''; my %p;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; my @e = split /,/, $m;
      my $key = join('*', map { "j$_^$e[$_-1]" } grep { $e[$_-1] } 1..$r); $p{$key} = R($cf); }
    $C[$k] = \%p; }
  return \@C; }
# substitute j_i -> linear forms: $map = { j2 => Poly } etc.
sub subst_poly { my ($p, $map) = @_; my $res = {};
  for my $k (keys %$p) { my $e = Poly::mono_parse($k); my $t = Poly::const($p->{$k});
    for my $v (keys %$e) { my $f = exists $map->{$v} ? $map->{$v} : Poly::var($v); $t = Poly::mul($t, Poly::pow($f, $e->{$v})); }
    $res = Poly::add($res, $t); }
  return $res; }
# total polynomial in sigma-form: sigma^h Xi = Sum_k (-1)^k sigma^(h-k) Xi^(k)  -> Poly with extra variable 's' (= sigma)
sub sigma_form { my ($C, $h) = @_; my $res = {};
  for my $k (0..$h) { my $t = Poly::mul($C->[$k], Poly::pow(Poly::var('s'), $h - $k)); $t = Poly::scale($t, (($k % 2) ? -1 : 1)); $res = Poly::add($res, $t); }
  return $res; }
1;
