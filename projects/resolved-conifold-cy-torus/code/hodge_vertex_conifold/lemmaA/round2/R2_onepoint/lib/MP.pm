package MP;
# Minimal multivariate polynomials over Q (Math::BigRat), keys = comma-joined exponent vectors.
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(Q mp_const mp_var mp_add mp_sub mp_mul mp_scal mp_pow mp_iszero mp_str mp_coef mp_eval mp_deg mp_subst1);
sub Q { Math::BigRat->new(@_) }
sub mp_const { my ($c, $nv) = @_; return {} if Q($c)->is_zero; return { join(',', (0) x $nv) => Q($c) }; }
sub mp_var { my ($i, $nv) = @_; my @e = (0) x $nv; $e[$i] = 1; return { join(',', @e) => Q(1) }; }
sub mp_add { my ($a, $b) = @_; my %c = map { $_ => $a->{$_}->copy } keys %$a;
  for my $k (keys %$b) { if (exists $c{$k}) { $c{$k} += $b->{$k}; delete $c{$k} if $c{$k}->is_zero; } else { $c{$k} = $b->{$k}->copy; } } return \%c; }
sub mp_scal { my ($s, $a) = @_; $s = Q($s) unless ref $s; return {} if $s->is_zero; my %c = map { $_ => $a->{$_} * $s } keys %$a; return \%c; }
sub mp_sub { my ($a, $b) = @_; return mp_add($a, mp_scal(Q(-1), $b)); }
sub mp_mul { my ($a, $b) = @_; my %c;
  for my $k1 (keys %$a) { my @e1 = split /,/, $k1; for my $k2 (keys %$b) { my @e2 = split /,/, $k2;
      my $k = join(',', map { $e1[$_] + $e2[$_] } 0..$#e1); my $v = $a->{$k1} * $b->{$k2};
      if (exists $c{$k}) { $c{$k} += $v; } else { $c{$k} = $v; } } }
  for my $k (keys %c) { delete $c{$k} if $c{$k}->is_zero; } return \%c; }
sub mp_pow { my ($a, $n, $nv) = @_; my $r = mp_const(1, $nv); $r = mp_mul($r, $a) for 1..$n; return $r; }
sub mp_iszero { my $a = shift; return scalar(keys %$a) == 0; }
sub mp_coef { my ($a, $key) = @_; return exists $a->{$key} ? $a->{$key}->copy : Q(0); }
sub mp_deg { my ($a, $i) = @_; my $d = -1; for my $k (keys %$a) { my @e = split /,/, $k; $d = $e[$i] if $e[$i] > $d; } return $d; }
sub mp_eval { my ($a, $vals) = @_; my $s = Q(0);
  for my $k (keys %$a) { my @e = split /,/, $k; my $t = $a->{$k}->copy; for my $i (0..$#e) { $t *= Q($vals->[$i]) ** $e[$i] if $e[$i]; } $s += $t; } return $s; }
# substitute variable i by a polynomial p (same variable set)
sub mp_subst1 { my ($a, $i, $p, $nv) = @_; my $res = {}; my %pp;
  for my $k (keys %$a) { my @e = split /,/, $k; my $n = $e[$i]; $e[$i] = 0; my $mono = { join(',', @e) => $a->{$k}->copy };
    $pp{$n} //= mp_pow($p, $n, $nv); $res = mp_add($res, mp_mul($mono, $pp{$n})); } return $res; }
sub mp_str { my ($pp, $names) = @_; my @terms;
  for my $k (sort { my @x = split /,/, $a; my @y = split /,/, $b; for my $i (0..$#x) { return $x[$i] <=> $y[$i] if $x[$i] != $y[$i]; } 0 } keys %$pp) {
    my @e = split /,/, $k; my $m = join('*', map { $e[$_] ? ($e[$_] == 1 ? $names->[$_] : "$names->[$_]^$e[$_]") : () } 0..$#e);
    my $c = $pp->{$k}; push @terms, ($m eq '' ? "$c" : "($c)*$m"); }
  return @terms ? join(' + ', @terms) : '0'; }
1;
