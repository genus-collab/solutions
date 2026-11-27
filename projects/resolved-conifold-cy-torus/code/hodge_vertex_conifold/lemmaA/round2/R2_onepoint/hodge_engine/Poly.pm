package Poly;
# Sparse multivariate polynomials with Math::BigRat coefficients.
# A polynomial is a hashref { monomial_key => BigRat }, monomial key = "v1^e1*v2^e2" sorted by variable name ("" = constant).
use strict; use warnings;
use Math::BigRat;
sub R { Math::BigRat->new(shift) }

sub mono_parse { my $k = shift; my %e; return \%e if $k eq ''; for my $f (split /\*/, $k) { my ($v, $x) = split /\^/, $f; $e{$v} = $x; } return \%e; }
sub mono_key { my $e = shift; return join('*', map { "$_^$e->{$_}" } grep { $e->{$_} } sort keys %$e); }
sub mono_mul { my ($a, $b) = @_; return $b if $a eq ''; return $a if $b eq '';
  my $ea = mono_parse($a); my $eb = mono_parse($b); $ea->{$_} += $eb->{$_} for keys %$eb; return mono_key($ea); }

sub const { my $c = shift; $c = R($c) unless ref $c; return {} if $c == 0; return { '' => $c }; }
sub var { my $v = shift; return { "$v^1" => R(1) }; }
sub clean { my $p = shift; for (keys %$p) { delete $p->{$_} if $p->{$_} == 0 } return $p; }
sub add { my ($p, $q) = @_; my %r = map { $_ => $p->{$_}->copy } keys %$p;
  for (keys %$q) { if (exists $r{$_}) { $r{$_} += $q->{$_} } else { $r{$_} = $q->{$_}->copy } } return clean(\%r); }
sub neg { my $p = shift; my %r = map { $_ => -$p->{$_} } keys %$p; return \%r; }
sub sub_ { my ($p, $q) = @_; return add($p, neg($q)); }
sub scale { my ($p, $c) = @_; $c = R($c) unless ref $c; return {} if $c == 0; my %r = map { $_ => $p->{$_} * $c } keys %$p; return \%r; }
sub mul { my ($p, $q) = @_; my %r;
  for my $a (keys %$p) { for my $b (keys %$q) { my $k = mono_mul($a, $b); my $c = $p->{$a} * $q->{$b};
      if (exists $r{$k}) { $r{$k} += $c } else { $r{$k} = $c } } }
  return clean(\%r); }
sub pow { my ($p, $n) = @_; die "negative power of polynomial" if $n < 0; my $r = const(1); my $b = $p;
  while ($n > 0) { $r = mul($r, $b) if $n & 1; $n >>= 1; $b = mul($b, $b) if $n; } return $r; }
sub is_const { my $p = shift; for (keys %$p) { return 0 if $_ ne '' } return 1; }
sub const_val { my $p = shift; return exists $p->{''} ? $p->{''} : R(0); }
sub is_zero { my $p = shift; return scalar(keys %$p) == 0; }
# evaluate (some) variables: $vals = { var => BigRat }
sub subst { my ($p, $vals) = @_; my %r;
  for my $k (keys %$p) { my $e = mono_parse($k); my $c = $p->{$k}->copy; my %rest;
    for my $v (keys %$e) { if (exists $vals->{$v}) { $c *= $vals->{$v} ** $e->{$v} } else { $rest{$v} = $e->{$v} } }
    my $nk = mono_key(\%rest); if (exists $r{$nk}) { $r{$nk} += $c } else { $r{$nk} = $c } }
  return clean(\%r); }
# coefficient of a monomial given as hash { var => exp } (other variables kept)
sub coeff_of { my ($p, $v, $e) = @_; my %r;   # coefficient of v^e
  for my $k (keys %$p) { my $m = mono_parse($k); my $ex = $m->{$v} // 0; next unless $ex == $e; delete $m->{$v};
    my $nk = mono_key($m); if (exists $r{$nk}) { $r{$nk} += $p->{$k} } else { $r{$nk} = $p->{$k}->copy } }
  return clean(\%r); }
sub degree_in { my ($p, $v) = @_; my $d = -1; for my $k (keys %$p) { my $m = mono_parse($k); my $e = $m->{$v} // 0; $d = $e if $e > $d } return $d; }
sub total_degree { my ($p, @vars) = @_; my $d = -1; my %s = map { $_ => 1 } @vars;
  for my $k (keys %$p) { my $m = mono_parse($k); my $t = 0; for (keys %$m) { $t += $m->{$_} if !@vars || $s{$_} } $d = $t if $t > $d } return $d; }
sub to_string { my $p = shift; my @k = sort keys %$p; return '0' unless @k;
  return join(' + ', map { ($_ eq '' ? '' : "") . "(" . $p->{$_} . ")" . ($_ eq '' ? '' : "*$_") } @k); }
sub equal { my ($p, $q) = @_; return is_zero(sub_($p, $q)); }
1;
