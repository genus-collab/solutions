# ref_lib.pl -- referee's own helpers (written independently of hodge_lib.pl / xi_lib.pl / mparse.pl).
#   rq(...)              Math::BigRat constructor
#   ev_expr($s, \%v)     evaluate a Mathematica InputForm rational expression at rational values (own parser)
#   parse_mlist($s, \%v) parse a nested Mathematica list {..} into nested Perl arrays of BigRat (True -> 1, False -> 0)
#   series helpers on array refs (orders 0..N): s_add, s_scale, s_mul, s_inv, s_exp, s_pow_rat (1+u)^alpha, s_int
# ASCII only. No TR of any kind.
use strict; use warnings;
use Math::BigRat;
sub rq { Math::BigRat->new(@_) }
sub rfact { my $n = shift; my $r = rq(1); $r *= $_ for 2..$n; return $r; }
sub rbinom { my ($n, $k) = @_; return rq(0) if $k < 0 || $k > $n; my $r = rq(1); for my $i (1..$k) { $r *= ($n - $k + $i); $r /= $i; } return $r; }
sub gbinom { my ($al, $k) = @_; my $r = rq(1); for my $i (0..$k-1) { $r *= ($al - $i); } return $r / rfact($k); }   # alpha choose k, alpha rational
sub poch { my ($a, $n) = @_; my $r = rq(1); for my $i (0..$n-1) { $r *= ($a + $i); } return $r; }                 # rising factorial

# ---------- expression evaluation (tokenise, then precedence climbing) ----------
sub _tok { my $s = shift; $s =~ s/\s+//g; my @t;
  while (length $s) {
    if ($s =~ s/^(\d+)//) { push @t, ['n', $1]; }
    elsif ($s =~ s/^([A-Za-z][A-Za-z0-9]*)//) { push @t, ['v', $1]; }
    elsif ($s =~ s/^([-+*\/^(){},])//) { push @t, ['o', $1]; }
    else { die "ev_expr: cannot tokenise at '" . substr($s, 0, 20) . "'"; } }
  return \@t; }
sub _ipow { my ($b, $e) = @_; die "non-integer exponent $e" unless $e->is_int(); my $n = $e->numify();
  return $n >= 0 ? $b->copy->bpow($n) : rq(1) / $b->copy->bpow(-$n); }
sub _parse_expr { my ($T, $ip, $V) = @_;
  my $v = _parse_term($T, $ip, $V);
  while ($$ip < @$T && $T->[$$ip][0] eq 'o' && ($T->[$$ip][1] eq '+' || $T->[$$ip][1] eq '-')) {
    my $op = $T->[$$ip++][1]; my $w = _parse_term($T, $ip, $V); $v = $op eq '+' ? $v + $w : $v - $w; }
  return $v; }
sub _parse_term { my ($T, $ip, $V) = @_;
  my $v = _parse_unary($T, $ip, $V);
  while ($$ip < @$T && $T->[$$ip][0] eq 'o' && ($T->[$$ip][1] eq '*' || $T->[$$ip][1] eq '/')) {
    my $op = $T->[$$ip++][1]; my $w = _parse_unary($T, $ip, $V); $v = $op eq '*' ? $v * $w : $v / $w; }
  return $v; }
sub _parse_unary { my ($T, $ip, $V) = @_;
  if ($T->[$$ip][0] eq 'o' && $T->[$$ip][1] eq '-') { $$ip++; return -_parse_unary($T, $ip, $V); }
  if ($T->[$$ip][0] eq 'o' && $T->[$$ip][1] eq '+') { $$ip++; return _parse_unary($T, $ip, $V); }
  return _parse_pow($T, $ip, $V); }
sub _parse_pow { my ($T, $ip, $V) = @_;
  my $b = _parse_atom($T, $ip, $V);
  if ($$ip < @$T && $T->[$$ip][0] eq 'o' && $T->[$$ip][1] eq '^') { $$ip++; my $e = _parse_unary($T, $ip, $V); return _ipow($b, $e); }
  return $b; }
sub _parse_atom { my ($T, $ip, $V) = @_;
  my $t = $T->[$$ip++]; die "unexpected end" unless defined $t;
  if ($t->[0] eq 'n') { return rq($t->[1]); }
  if ($t->[0] eq 'v') { return rq(1) if $t->[1] eq 'True'; return rq(0) if $t->[1] eq 'False';
    die "unknown symbol $t->[1]" unless exists $V->{$t->[1]}; return rq($V->{$t->[1]}); }
  if ($t->[1] eq '(') { my $v = _parse_expr($T, $ip, $V); die "missing )" unless $T->[$$ip++][1] eq ')'; return $v; }
  die "unexpected token $t->[1]"; }
sub ev_expr { my ($s, $V) = @_; my $T = _tok($s); my $i = 0; my $v = _parse_expr($T, \$i, $V); die "trailing tokens" if $i < @$T; return $v; }
# nested lists
sub _parse_item { my ($T, $ip, $V) = @_;
  if ($T->[$$ip][0] eq 'o' && $T->[$$ip][1] eq '{') { $$ip++; my @L;
    if ($T->[$$ip][0] eq 'o' && $T->[$$ip][1] eq '}') { $$ip++; return \@L; }
    while (1) { push @L, _parse_item($T, $ip, $V); my $t = $T->[$$ip++];
      last if $t->[1] eq '}'; die "expected , or }" unless $t->[1] eq ','; }
    return \@L; }
  return _parse_expr($T, $ip, $V); }
sub parse_mlist { my ($s, $V) = @_; my $T = _tok($s); my $i = 0; my $v = _parse_item($T, \$i, $V); return $v; }

# ---------- truncated power series (array refs, orders 0..N) ----------
sub s_zero { my $N = shift; return [ map { rq(0) } 0..$N ]; }
sub s_const { my ($a, $N) = @_; my $s = s_zero($N); $s->[0] = rq($a); return $s; }
sub s_add { my ($A, $B) = @_; return [ map { $A->[$_] + $B->[$_] } 0..$#$A ]; }
sub s_sub { my ($A, $B) = @_; return [ map { $A->[$_] - $B->[$_] } 0..$#$A ]; }
sub s_scale { my ($A, $c) = @_; return [ map { $_ * $c } @$A ]; }
sub s_mul { my ($A, $B) = @_; my $N = $#$A; my $C = s_zero($N);
  for my $i (0..$N) { next if $A->[$i] == 0; for my $j (0..$N-$i) { $C->[$i+$j] += $A->[$i] * $B->[$j]; } } return $C; }
sub s_inv { my ($A) = @_; my $N = $#$A; die "s_inv: zero constant" if $A->[0] == 0; my $B = s_zero($N); $B->[0] = rq(1) / $A->[0];
  for my $n (1..$N) { my $s = rq(0); $s += $A->[$_] * $B->[$n-$_] for 1..$n; $B->[$n] = -$s / $A->[0]; } return $B; }
sub s_powint { my ($A, $e) = @_; my $N = $#$A; my $R = s_const(1, $N); $R = s_mul($R, $A) for 1..$e; return $R; }
sub s_exp { my ($A) = @_; my $N = $#$A; die "s_exp: nonzero constant" if $A->[0] != 0;
  my $E = s_const(1, $N); my $term = s_const(1, $N);
  for my $k (1..$N) { $term = s_scale(s_mul($term, $A), rq(1, $k)); $E = s_add($E, $term); } return $E; }
sub s_pow_rat { my ($U, $al) = @_; my $N = $#$U; die "s_pow_rat: u must have zero constant" if $U->[0] != 0;   # (1 + U)^alpha
  my $R = s_const(1, $N); my $P = s_const(1, $N);
  for my $k (1..$N) { $P = s_mul($P, $U); $R = s_add($R, s_scale($P, gbinom($al, $k))); } return $R; }
sub s_int { my ($A) = @_; my $N = $#$A; my $B = s_zero($N); $B->[$_] = $A->[$_-1] / $_ for 1..$N; return $B; }    # integral from 0
sub s_shift_down { my ($A) = @_; my $N = $#$A; die "s_shift_down: nonzero constant" if $A->[0] != 0;
  my $B = s_zero($N); $B->[$_] = $A->[$_+1] for 0..$N-1; return $B; }   # A/t (last coefficient lost -> set 0; caller keeps N+1 margin)
1;
