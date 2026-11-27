package mkp_eqs;
# The u_k-linear members of the modified-KP hierarchy of the 2D Toda lattice for (tau_0, tau_1) = (Z(p), Z(e^beta o p)):
#   P_k(D) tau_0 . tau_1 = 0,   P_k(D) = D_{t_k} + p_k(-D~),   D~ = (D_{t1}, D_{t2}/2, D_{t3}/3, ...),   t_k = p_k/k
# (u_k-coefficient of the bilinear identity Res_z tau_0(t+u-[1/z]) tau_1(t-u+[1/z]) e^{xi(2u,z)} dz/z = tau_1(t+u) tau_0(t-u)),
# p_k = complete homogeneous Schur polynomial, exp(Sum t_i z^i) = Sum p_k(t) z^k.
# Translated to F = log tau_0, A = log tau_1 - log tau_0 = Sum_{m>=1} beta^m E^m F / m!, S = 2F + A by
#   P(D) tau_0.tau_1/(tau_0 tau_1) = P(d_s) exp(Phi)|_{s=0},   Phi(s) = F(t+s) + F1(t-s) - F(t) - F1(t),
#   d_s^B Phi|_0 = d^B S (|B| even),  -d^B A (|B| odd)   (set-partition (Faa di Bruno) expansion of exp).
# Output: list of terms [coefficient (BigRat), [factor...]], factor = [type 'A'|'S', [inserted parts >= 2], u = number of d_{p1}],
# with d_{t_i} = i d_{p_i}.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub partitions { my ($n, $max) = @_; $max //= $n; return ([]) if $n == 0; my @o;
  for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for partitions($n - $k, $k); } @o }
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @out;
  for my $sp (setparts($n-1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n-1; push @out, \@c; }
    my @c = map { [@$_] } @$sp; push @c, [$n-1]; push @out, \@c; } @out }
sub canon { my ($fac) = @_; my @f = map { [$_->[0], [sort { $a <=> $b } @{$_->[1]}], $_->[2]] } @$fac;
  return join(';', sort map { $_->[0] . ':' . join(',', @{$_->[1]}) . ':' . $_->[2] } @f); }
sub eq_k { my ($k) = @_;
  # P_k as a polynomial in D_{t_i}: hash "sorted index list" -> coefficient
  my %P; $P{"$k"} = R(1);
  for my $lam (partitions($k)) { my %m; $m{$_}++ for @$lam; my $c = R(1);
    for my $i (keys %m) { $c *= (R(-1) / $i) ** $m{$i} / fact($m{$i}); }
    my $key = join(',', sort { $a <=> $b } @$lam); $P{$key} = ($P{$key} // R(0)) + $c; }
  my %T; my %TF;
  for my $key (keys %P) { my $c = $P{$key}; next if $c == 0; my @idx = split /,/, $key; my $q = scalar @idx;
    for my $sp (setparts($q)) { my $coef = $c->copy; my @fac;
      for my $B (@$sp) { my @ind = map { $idx[$_] } @$B; my $type = (@ind % 2) ? 'A' : 'S'; $coef = -$coef if @ind % 2;
        my $pw = R(1); $pw *= $_ for @ind;   # d_{t_i} = i d_{p_i}
        $coef *= $pw; push @fac, [$type, [grep { $_ >= 2 } @ind], scalar(grep { $_ == 1 } @ind)]; }
      my $cf = canon(\@fac); $T{$cf} = ($T{$cf} // R(0)) + $coef; $TF{$cf} = \@fac; } }
  my @terms; for my $cf (sort keys %T) { next if $T{$cf} == 0; push @terms, [$T{$cf}, $TF{$cf}]; }
  return \@terms; }
sub show { my ($terms) = @_; join(' + ', map { my ($c, $f) = @$_; "($c)" . join('', map { "[$_->[0];" . join(',', @{$_->[1]}) . ";u=$_->[2]]" } @$f) } @$terms); }
1;
