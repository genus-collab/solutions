package RefData;
# Referee's own reader for the Mathematica InputForm tables in lemmaA\data (no code shared with HB.pm / Data.pm).
# Grammar: expr = ['+'|'-'] term {('+'|'-') term}; term = factor {('*'|'/') factor}; factor = unary ['^' factor];
# unary = '-' unary | atom; atom = integer | name | name '[' integer ']' | '(' expr ')'.   (jj[1] -> variable "jj1")
# k = 0 component: a table value E is in the c-normalisation (c = 1 + x); F := c^(3g-3+2r) E is a polynomial of degree <= g in
# t = c^2/x (t = -1/sigma); F(t = 0) is the single-class (Lambda^vee(1)) component.  The referee uses the NEGATIVE sample points
# x = -2, -3, ..., which the author did not use (he used x = 1, 2, ...), and checks the degree with two extra points.
use strict; use warnings; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }

sub tokens { my $s = shift; my @t;
  while (length $s) { $s =~ s/^\s+//; last unless length $s;
    if ($s =~ s/^(\d+)//) { push @t, ['n', $1] }
    elsif ($s =~ s/^([A-Za-z][A-Za-z0-9]*)\[(\d+)\]//) { push @t, ['v', "$1$2"] }
    elsif ($s =~ s/^([A-Za-z][A-Za-z0-9]*)//) { push @t, ['v', $1] }
    elsif ($s =~ s/^([-+*\/^()])//) { push @t, ['o', $1] }
    else { die "RefData: cannot tokenize at: " . substr($s, 0, 30) } }
  return \@t }
# AST: ['c', BigRat] | ['v', name] | ['add', [s, node]...] | ['mul', [e, node]...] | ['pow', node, int]
sub parse { my $s = shift; my $t = tokens($s); my $i = 0; my $n = p_expr($t, \$i); die "RefData: trailing tokens" if $i < @$t; return $n }
sub peek { my ($t, $i) = @_; return $$i < @$t ? $t->[$$i] : ['eof', ''] }
sub p_expr { my ($t, $i) = @_; my @terms; my $sg = 1;
  my $k = peek($t, $i); if ($k->[0] eq 'o' && ($k->[1] eq '-' || $k->[1] eq '+')) { $sg = $k->[1] eq '-' ? -1 : 1; $$i++ }
  push @terms, [$sg, p_term($t, $i)];
  while (1) { my $k2 = peek($t, $i); last unless $k2->[0] eq 'o' && ($k2->[1] eq '+' || $k2->[1] eq '-'); $$i++;
    push @terms, [($k2->[1] eq '-' ? -1 : 1), p_term($t, $i)] }
  return ['add', @terms] }
sub p_term { my ($t, $i) = @_; my @f = ([1, p_factor($t, $i)]);
  while (1) { my $k = peek($t, $i); last unless $k->[0] eq 'o' && ($k->[1] eq '*' || $k->[1] eq '/'); $$i++;
    push @f, [($k->[1] eq '/' ? -1 : 1), p_factor($t, $i)] }
  return ['mul', @f] }
sub p_factor { my ($t, $i) = @_; my $b = p_unary($t, $i); my $k = peek($t, $i);
  if ($k->[0] eq 'o' && $k->[1] eq '^') { $$i++; my $neg = 0; my $k2 = peek($t, $i);
    if ($k2->[0] eq 'o' && $k2->[1] eq '-') { $neg = 1; $$i++ }
    my $k3 = peek($t, $i); die "RefData: non-integer exponent" unless $k3->[0] eq 'n'; $$i++;
    return ['pow', $b, ($neg ? -1 : 1) * $k3->[1]] }
  return $b }
sub p_unary { my ($t, $i) = @_; my $k = peek($t, $i);
  if ($k->[0] eq 'o' && $k->[1] eq '-') { $$i++; return ['mul', [1, ['c', BR(-1)]], [1, p_unary($t, $i)]] }
  return p_atom($t, $i) }
sub p_atom { my ($t, $i) = @_; my $k = peek($t, $i);
  if ($k->[0] eq 'n') { $$i++; return ['c', BR($k->[1])] }
  if ($k->[0] eq 'v') { $$i++; return ['v', $k->[1]] }
  if ($k->[0] eq 'o' && $k->[1] eq '(') { $$i++; my $e = p_expr($t, $i); my $c = peek($t, $i);
    die "RefData: missing )" unless $c->[0] eq 'o' && $c->[1] eq ')'; $$i++; return $e }
  die "RefData: unexpected token $k->[1]" }
sub ev { my ($n, $env) = @_; my $k = $n->[0];
  return $n->[1]->copy if $k eq 'c';
  if ($k eq 'v') { die "RefData: unbound $n->[1]" unless defined $env->{$n->[1]}; return BR($env->{$n->[1]}) }
  if ($k eq 'add') { my $s = BR(0); for my $j (1..$#$n) { my $v = ev($n->[$j][1], $env); if ($n->[$j][0] > 0) { $s += $v } else { $s -= $v } } return $s }
  if ($k eq 'mul') { my ($nu, $de) = (BR(1), BR(1)); for my $j (1..$#$n) { my $v = ev($n->[$j][1], $env); if ($n->[$j][0] > 0) { $nu *= $v } else { $de *= $v } }
    die "RefData: division by zero" if $de->is_zero; return $nu / $de }
  if ($k eq 'pow') { my $b = ev($n->[1], $env); my $e = $n->[2]; my $r = BR(1); my $bb = $e < 0 ? BR(1) / $b : $b; $r *= $bb for 1..abs($e); return $r }
  die "RefData: bad node $k" }
sub read_table { my $file = shift; open my $fh, '<', $file or die "RefData: cannot open $file"; my @rows;
  while (my $l = <$fh>) { chomp $l; next unless $l =~ /\S/; my @f = split /\s*\|\s*/, $l; push @rows, { fields => \@f, ast => parse($f[-1]) } }
  close $fh; return \@rows }
sub lagr0 { my ($ts, $vs, $t0) = @_; my $s = BR(0);
  for my $i (0..$#$ts) { my $w = $vs->[$i]->copy; for my $j (0..$#$ts) { next if $j == $i; $w *= ($t0 - $ts->[$j]) / ($ts->[$i] - $ts->[$j]) } $s += $w }
  return $s }
# k = 0 component of a table value at fixed other variables; g = genus, r = number of probe leaves
sub k0 { my ($ast, $g, $r, $env) = @_; my (@t, @v);
  for my $j (0..$g + 2) { my $x = BR(-2 - $j); my $c = BR(1) + $x; my %e = (%$env, x => $x);
    push @t, $c * $c / $x; push @v, ev($ast, \%e) * ($c ** (3 * $g - 3 + 2 * $r)) }
  my @t0 = @t[0..$g]; my @v0 = @v[0..$g];
  for my $j ($g + 1, $g + 2) { die "RefData: degree-in-t check failed" unless lagr0(\@t0, \@v0, $t[$j]) == $v[$j] }
  return lagr0(\@t0, \@v0, BR(0)) }
1;
