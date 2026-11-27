package MParse;
# Parser for Mathematica InputForm polynomial / rational expressions (as stored in the HVC data files).
# parse($string, \%subst) returns a Poly (hashref) in the remaining symbols; divisions must be by expressions
# that are constant after substitution. jj[i] -> variable "jj<i>"; wk[g, {a1,...}] -> Witten-Kontsevich number.
use strict; use warnings;
use Math::BigRat;
use Poly; use WK;
my @tok; my $pos; my $sub;
sub R { Math::BigRat->new(shift) }
sub tokenize { my $s = shift; my @t;
  while (length $s) {
    if ($s =~ s/^\s+//) { next }
    if ($s =~ s/^(\d+)//) { push @t, ['num', $1]; next }
    if ($s =~ s/^([A-Za-z\$][A-Za-z0-9\$]*)//) { push @t, ['id', $1]; next }
    if ($s =~ s/^([-+*\/^()\[\]{},])//) { push @t, ['op', $1]; next }
    die "tokenize: cannot parse near '" . substr($s, 0, 30) . "'";
  }
  return @t; }
sub peek { return $pos < @tok ? $tok[$pos] : ['eof', '']; }
sub take { my $t = peek(); $pos++; return $t; }
sub expect { my $v = shift; my $t = take(); die "expected $v, got $t->[1]" unless $t->[0] eq 'op' && $t->[1] eq $v; }
sub parse { my ($s, $subst) = @_; $sub = $subst || {}; @tok = tokenize($s); $pos = 0;
  my $r = p_expr(); die "trailing tokens at $pos" if $pos < @tok; return $r; }
sub p_expr {
  my $sign = 1; my $t = peek();
  if ($t->[0] eq 'op' && ($t->[1] eq '-' || $t->[1] eq '+')) { take(); $sign = -1 if $t->[1] eq '-'; }
  my $r = p_term(); $r = Poly::neg($r) if $sign < 0;
  while (1) { $t = peek(); last unless $t->[0] eq 'op' && ($t->[1] eq '+' || $t->[1] eq '-'); take();
    my $q = p_term(); $r = $t->[1] eq '+' ? Poly::add($r, $q) : Poly::sub_($r, $q); }
  return $r; }
sub p_term {
  my $r = p_factor();
  while (1) { my $t = peek(); last unless $t->[0] eq 'op' && ($t->[1] eq '*' || $t->[1] eq '/'); take();
    my $q = p_factor();
    if ($t->[1] eq '*') { $r = Poly::mul($r, $q) }
    else { die "division by non-constant" unless Poly::is_const($q); my $c = Poly::const_val($q); die "division by zero" if $c == 0; $r = Poly::scale($r, R(1)/$c); } }
  return $r; }
sub p_factor { my $t = peek(); if ($t->[0] eq 'op' && $t->[1] eq '-') { take(); return Poly::neg(p_factor()); } return p_power(); }
sub p_power { my $b = p_atom(); my $t = peek();
  if ($t->[0] eq 'op' && $t->[1] eq '^') { take(); my $neg = 0; my $u = peek(); if ($u->[0] eq 'op' && $u->[1] eq '-') { take(); $neg = 1; }
    my $e = p_atom(); die "non-constant exponent" unless Poly::is_const($e); my $n = Poly::const_val($e); die "non-integer exponent" unless $n->is_int();
    $n = $n->numify(); $n = -$n if $neg;
    if ($n >= 0) { return Poly::pow($b, $n) }
    die "negative power of non-constant" unless Poly::is_const($b); my $c = Poly::const_val($b); return Poly::const((R(1)/$c) ** (-$n)); }
  return $b; }
sub p_atom { my $t = take();
  if ($t->[0] eq 'num') { return Poly::const(R($t->[1])) }
  if ($t->[0] eq 'op' && $t->[1] eq '(') { my $r = p_expr(); expect(')'); return $r; }
  if ($t->[0] eq 'id') { my $name = $t->[1];
    if ($name eq 'wk') { expect('['); my $g = take()->[1]; expect(','); expect('{'); my @a;
      while (1) { my $u = peek(); if ($u->[0] eq 'op' && $u->[1] eq '}') { take(); last } if ($u->[0] eq 'op' && $u->[1] eq ',') { take(); next } push @a, take()->[1]; }
      expect(']'); return Poly::const(WK::wk($g, @a)); }
    my $u = peek();
    if ($u->[0] eq 'op' && $u->[1] eq '[') { take(); my $i = take()->[1]; expect(']'); $name = "$name$i"; }
    if (exists $sub->{$name}) { my $v = $sub->{$name}; return ref($v) eq 'HASH' ? $v : Poly::const($v); }
    return Poly::var($name); }
  die "unexpected token $t->[1]"; }
1;
