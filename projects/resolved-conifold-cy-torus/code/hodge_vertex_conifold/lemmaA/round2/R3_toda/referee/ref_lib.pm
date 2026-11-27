package ref_lib;
# Referee library (written from scratch, independent of mvlib.pm / mkp_eqs.pm):
#   partitions, Murnaghan-Nakayama characters (beta-sets), z_mu, content sums,
#   sparse polynomials in the power sums p_1, p_2, ... (keys = parts in decreasing order, "" = constant),
#   derivatives d/dt_k = k d/dp_k, products truncated by degree, Hirota bilinear operators.
use strict; use warnings; use Math::BigRat;
our @EXPORT_OK;
sub R { Math::BigRat->new(@_) }
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }

# partitions of n (decreasing parts), parts <= max
sub partitions { my ($n, $max) = @_; $max //= $n; return ([]) if $n == 0; my @o;
  for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for partitions($n - $k, $k); } return @o; }

# z_mu = prod_i i^{m_i} m_i!
sub zee { my @mu = @_; my %m; $m{$_}++ for @mu; my $z = R(1); for my $i (keys %m) { $z *= R($i) ** $m{$i} * fact($m{$i}); } $z }

# content sum f_2(nu) = sum over boxes (col - row)
sub f2 { my @nu = @_; my $s = 0; for my $i (0..$#nu) { for my $j (0..$nu[$i]-1) { $s += $j - $i; } } $s }

# Murnaghan-Nakayama: chi^nu(mu) via beta-sets
my %chimemo;
sub chi { my ($nu, $mu) = @_; my $L = scalar @$nu; my @beta = map { $nu->[$_] + $L - 1 - $_ } 0..$L-1;   # distinct, decreasing
  return _chi([sort { $a <=> $b } @beta], [@$mu], $L); }
sub _chi { my ($beta, $mu, $L) = @_;
  if (!@$mu) { for my $i (0..$L-1) { return 0 if $beta->[$i] != $i; } return 1; }
  my $key = join(',', @$beta) . '|' . join(',', @$mu); return $chimemo{$key} if exists $chimemo{$key};
  my @m = @$mu; my $k = shift @m; my %in = map { $_ => 1 } @$beta; my $tot = 0;
  for my $bb (@$beta) { my $nb = $bb - $k; next if $nb < 0 || $in{$nb};
    my $cnt = 0; for my $c (@$beta) { $cnt++ if $c > $nb && $c < $bb; }
    my @new = sort { $a <=> $b } ((grep { $_ != $bb } @$beta), $nb);
    $tot += (($cnt % 2) ? -1 : 1) * _chi(\@new, \@m, $L); }
  $chimemo{$key} = $tot; return $tot; }

# ---- sparse polynomials in p ----
sub mkey { my @p = sort { $b <=> $a } @_; return join(',', @p); }
sub parts { my $k = shift; return () if $k eq ''; return split /,/, $k; }
sub deg { my $k = shift; my $d = 0; $d += $_ for parts($k); $d }
sub padd { my ($a, $b, $c) = @_; my %r = %$a; $c //= R(1);
  for my $k (keys %$b) { my $v = ($r{$k} // R(0)) + $c * $b->{$k}; if ($v == 0) { delete $r{$k} } else { $r{$k} = $v } } return \%r; }
sub pscale { my ($a, $c) = @_; my %r; return \%r if $c == 0; $r{$_} = $a->{$_} * $c for keys %$a; return \%r; }
sub pmul { my ($a, $b, $dmax) = @_; my %r;
  for my $k1 (keys %$a) { my $d1 = deg($k1); my @p1 = parts($k1);
    for my $k2 (keys %$b) { my $d2 = deg($k2); next if defined $dmax && $d1 + $d2 > $dmax;
      my $k = mkey(@p1, parts($k2)); my $v = ($r{$k} // R(0)) + $a->{$k1} * $b->{$k2}; if ($v == 0) { delete $r{$k} } else { $r{$k} = $v } } }
  return \%r; }
# d/dt_k = k d/dp_k
sub dt { my ($a, $k) = @_; my %r;
  for my $key (keys %$a) { my @p = parts($key); my $m = grep { $_ == $k } @p; next unless $m;
    my $removed = 0; my @q = grep { ($_ == $k && !$removed++) ? 0 : 1 } @p;
    my $nk = mkey(@q); $r{$nk} = ($r{$nk} // R(0)) + $a->{$key} * $m * $k; }
  for (keys %r) { delete $r{$_} if $r{$_} == 0 } return \%r; }
sub dts { my ($a, @idx) = @_; my $r = $a; $r = dt($r, $_) for @idx; return $r; }
# Hirota D_{t_{i1}} ... D_{t_{im}} f.g, output truncated at degree <= dmax
sub hirota { my ($f, $g, $idx, $dmax) = @_; my @idx = @$idx; my $m = scalar @idx; my $res = {};
  for my $mask (0..(1<<$m)-1) { my (@a, @b); for my $i (0..$m-1) { if ($mask & (1<<$i)) { push @a, $idx[$i] } else { push @b, $idx[$i] } }
    my $sign = (scalar(@b) % 2) ? R(-1) : R(1);
    $res = padd($res, pmul(dts($f, @a), dts($g, @b), $dmax), $sign); }
  return $res; }
# elementary Schur polynomial p_k(-D~) as a list of [coef, [indices]],  D~ = (D_1, D_2/2, D_3/3, ...)
sub pk_minus_Dtilde { my ($k) = @_; my @out;
  for my $lam (partitions($k)) { my %m; $m{$_}++ for @$lam; my $c = R(1);
    for my $i (keys %m) { $c *= (R(-1) / $i) ** $m{$i} / fact($m{$i}); } push @out, [$c, [@$lam]]; }
  return @out; }
1;
