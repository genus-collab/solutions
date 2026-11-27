package rpoly;
# Referee's own exact polynomial toolkit (independent of the author's polylib.pm).
# Polynomials in the fixed variable list (x, s, j1, j2, j3) with Math::BigInt coefficients (integer numerators);
# key = "e_x,e_s,e_j1,e_j2,e_j3". Rational polynomials (for final components) use Math::BigRat values.
# Parser: Mathematica InputForm numerators (+ - * ^, parentheses, integers, variables x, s, j1..j3, jj[1]..jj[3], m1..m3).
use strict; use warnings;
use Math::BigInt; use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(NV vidx zkey kadd ip_add ip_scale ip_mul ip_pow ip_const ip_var parse_num split_frac
                 ip_divx1 ip_lowx ip_shiftx ip_subst_s rp_from_ip rp_add rp_scale rp_mul rp_eval rp_str parse_rp rp_iszero);
our $NV = 5; sub NV { $NV }
my %VI = (x => 0, s => 1, j1 => 2, j2 => 3, j3 => 4, 'jj[1]' => 2, 'jj[2]' => 3, 'jj[3]' => 4, m1 => 2, m2 => 3, m3 => 4);
sub vidx { my $v = shift; die "unknown variable $v" unless exists $VI{$v}; return $VI{$v}; }
sub zkey { return join(',', (0) x $NV); }
sub kadd { my ($a, $b) = @_; my @a = split /,/, $a; my @b = split /,/, $b; return join(',', map { $a[$_] + $b[$_] } 0 .. $NV - 1); }
# ---- integer polynomials ----
sub ip_const { my $c = Math::BigInt->new(shift); return $c->is_zero ? {} : { zkey() => $c }; }
sub ip_var { my $i = vidx(shift); my @e = (0) x $NV; $e[$i] = 1; return { join(',', @e) => Math::BigInt->new(1) }; }
sub ip_add { my ($a, $b) = @_; my %r = map { $_ => $a->{$_}->copy } keys %$a;
  for my $k (keys %$b) { if (exists $r{$k}) { $r{$k}->badd($b->{$k}); delete $r{$k} if $r{$k}->is_zero; } else { $r{$k} = $b->{$k}->copy; } } return \%r; }
sub ip_scale { my ($a, $c) = @_; $c = Math::BigInt->new($c) unless ref $c; return {} if $c->is_zero; my %r = map { $_ => $a->{$_}->copy->bmul($c) } keys %$a; return \%r; }
sub ip_mul { my ($a, $b) = @_; my %r;
  for my $k1 (keys %$a) { for my $k2 (keys %$b) { my $k = kadd($k1, $k2); my $p = $a->{$k1}->copy->bmul($b->{$k2});
      if (exists $r{$k}) { $r{$k}->badd($p); } else { $r{$k} = $p; } } }
  for (keys %r) { delete $r{$_} if $r{$_}->is_zero; } return \%r; }
sub ip_pow { my ($a, $n) = @_; my $r = ip_const(1); for (1 .. $n) { $r = ip_mul($r, $a); } return $r; }
# ---- parser for numerators ----
my ($S, $P);
sub peek { return $P < length($S) ? substr($S, $P, 1) : ''; }
sub parse_num { my $str = shift; $S = $str; $S =~ s/\s+//g; $P = 0; my $r = pe_sum(); die "parse error at $P: " . substr($S, $P, 40) if $P < length($S); return $r; }
sub pe_sum { my $r; my $c = peek(); if ($c eq '-') { $P++; $r = ip_scale(pe_prod(), -1); } elsif ($c eq '+') { $P++; $r = pe_prod(); } else { $r = pe_prod(); }
  while (1) { $c = peek(); if ($c eq '+') { $P++; $r = ip_add($r, pe_prod()); } elsif ($c eq '-') { $P++; $r = ip_add($r, ip_scale(pe_prod(), -1)); } else { last; } } return $r; }
sub pe_prod { my $r = pe_pow(); while (peek() eq '*') { $P++; $r = ip_mul($r, pe_pow()); } die "division inside numerator at $P" if peek() eq '/'; return $r; }
sub pe_pow { my $b = pe_atom(); if (peek() eq '^') { $P++; my $n = ''; while (peek() =~ /\d/) { $n .= peek(); $P++; } die "bad exponent" if $n eq ''; return ip_pow($b, $n); } return $b; }
sub pe_atom { my $c = peek();
  if ($c eq '(') { $P++; my $r = pe_sum(); die "expected ) at $P" unless peek() eq ')'; $P++; return $r; }
  if ($c eq '-') { $P++; return ip_scale(pe_pow(), -1); }
  if ($c =~ /\d/) { my $n = ''; while (peek() =~ /\d/) { $n .= peek(); $P++; } return ip_const($n); }
  if ($c =~ /[A-Za-z]/) { my $v = ''; while (peek() =~ /[A-Za-z0-9]/) { $v .= peek(); $P++; } if (peek() eq '[') { while (peek() ne ']') { $v .= peek(); $P++; } $v .= ']'; $P++; } return ip_var($v); }
  die "unexpected '$c' at $P"; }
# split "num/(den)" at the last top-level '/'
sub split_frac { my $str = shift; $str =~ s/\s+//g; my $depth = 0; my $pos = -1;
  for my $i (0 .. length($str) - 1) { my $ch = substr($str, $i, 1); $depth++ if $ch eq '('; $depth-- if $ch eq ')'; $pos = $i if ($ch eq '/' && $depth == 0); }
  return ($str, '1') if $pos < 0; return (substr($str, 0, $pos), substr($str, $pos + 1)); }
# lowest power of x present
sub ip_lowx { my $a = shift; my $m; for (keys %$a) { my ($e) = split /,/; $m = $e if (!defined $m || $e < $m); } return $m // 0; }
# multiply by x^d (d may be negative: Laurent shift; exponents of x may become negative)
sub ip_shiftx { my ($a, $d) = @_; my %r; for (keys %$a) { my @e = split /,/; $e[0] += $d; $r{join(',', @e)} = $a->{$_}->copy; } return \%r; }
# exact division by (1 + x) in the variable x (other variables as coefficients); dies if not exact
sub ip_divx1 { my $a = shift; my %byrest; # group by the other exponents
  for my $k (keys %$a) { my @e = split /,/, $k; my $ex = shift @e; push @{$byrest{join(',', @e)}}, [$ex, $a->{$k}]; }
  my %r;
  for my $rest (keys %byrest) { my %c; my ($lo, $hi); for my $t (@{$byrest{$rest}}) { $c{$t->[0]} = $t->[1]; $lo = $t->[0] if (!defined $lo || $t->[0] < $lo); $hi = $t->[0] if (!defined $hi || $t->[0] > $hi); }
    # divide sum_{e=lo..hi} c_e x^e by (1+x): q_{e} for e = lo..hi-1, from the top: q_{hi-1} = c_hi, q_{e-1} = c_e - q_e
    my %q; my $carry = Math::BigInt->new(0);
    for (my $e = $hi; $e > $lo; $e--) { my $ce = exists $c{$e} ? $c{$e}->copy : Math::BigInt->new(0); my $qe1 = $ce->bsub($carry); $q{$e - 1} = $qe1; $carry = $qe1->copy; }
    my $clo = exists $c{$lo} ? $c{$lo}->copy : Math::BigInt->new(0); die "division by (1+x) not exact" unless $clo->bcmp($carry) == 0;
    for my $e (keys %q) { next if $q{$e}->is_zero; $r{join(',', $e, $rest)} = $q{$e}; } }
  return \%r; }
# substitute s = value (integer)
sub ip_subst_s { my ($a, $v) = @_; my %r; my $V = Math::BigInt->new($v);
  for my $k (keys %$a) { my @e = split /,/, $k; my $es = $e[1]; $e[1] = 0; my $nk = join(',', @e); my $c = $a->{$k}->copy->bmul($V->copy->bpow($es));
    if (exists $r{$nk}) { $r{$nk}->badd($c); } else { $r{$nk} = $c; } }
  for (keys %r) { delete $r{$_} if $r{$_}->is_zero; } return \%r; }
# ---- rational polynomials in (j1, j2, j3) only: key "a,b,c" ----
sub rp_from_ip { my ($a, $den) = @_; my %r; for my $k (keys %$a) { my @e = split /,/, $k; die "x or s left in rp_from_ip" if ($e[0] || $e[1]);
    my $nk = join(',', @e[2 .. 4]); $r{$nk} = Math::BigRat->new($a->{$k}->bstr . "/" . $den); } return \%r; }
sub rp_add { my ($a, $b) = @_; my %r = map { $_ => $a->{$_}->copy } keys %$a; for my $k (keys %$b) { if (exists $r{$k}) { $r{$k} += $b->{$k}; delete $r{$k} if $r{$k}->is_zero; } else { $r{$k} = $b->{$k}->copy; } } return \%r; }
sub rp_scale { my ($a, $c) = @_; $c = Math::BigRat->new($c) unless ref $c; return {} if $c->is_zero; my %r = map { $_ => $a->{$_} * $c } keys %$a; return \%r; }
sub rp_mul { my ($a, $b) = @_; my %r; for my $k1 (keys %$a) { for my $k2 (keys %$b) { my @x = split /,/, $k1; my @y = split /,/, $k2; my $k = join(',', map { $x[$_] + $y[$_] } 0 .. 2);
      my $p = $a->{$k1} * $b->{$k2}; if (exists $r{$k}) { $r{$k} += $p; } else { $r{$k} = $p; } } } for (keys %r) { delete $r{$_} if $r{$_}->is_zero; } return \%r; }
sub rp_iszero { my $a = shift; return !%$a; }
sub rp_eval { my ($a, @v) = @_; my $s = Math::BigRat->new(0); for my $k (keys %$a) { my @e = split /,/, $k; my $t = $a->{$k}->copy;
    for my $i (0 .. 2) { next unless $e[$i]; $t *= Math::BigRat->new($v[$i] // 0) ** $e[$i]; } $s += $t; } return $s; }
sub rp_str { my $a = shift; return '0' unless %$a; my @t; for my $k (sort keys %$a) { my @e = split /,/, $k; my @m; for my $i (0 .. 2) { next unless $e[$i]; push @m, "j" . ($i + 1) . ($e[$i] > 1 ? "^$e[$i]" : ''); }
    push @t, "(" . $a->{$k}->bstr . ")" . (@m ? "*" . join('*', @m) : ''); } return join(' + ', @t); }
# parse the author's component format "(c)*j1^a*j2^b + (c) + ..." into an rp
sub parse_rp { my $str = shift; $str =~ s/\s+//g; return {} if $str eq '0'; my %r;
  while (length $str) { $str =~ s/^\+//; $str =~ s/^\(([-0-9\/]+)\)// or die "bad coefficient in '$str'"; my $c = Math::BigRat->new($1); my @e = (0, 0, 0);
    while ($str =~ s/^\*j(\d)(?:\^(\d+))?//) { $e[$1 - 1] += (defined $2 ? $2 : 1); } my $k = join(',', @e); $r{$k} = exists $r{$k} ? $r{$k} + $c : $c; }
  for (keys %r) { delete $r{$_} if $r{$_}->is_zero; } return \%r; }
1;
