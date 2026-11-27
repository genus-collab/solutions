package RefParse;
# Referee's own evaluator for Mathematica InputForm rational expressions (+ - * / ^, parentheses, integers, symbols, sym[n]).
# evalx($string, { sym => BigRat, ... }) returns a BigRat.  (Independent of ../MParse.pm.)
use strict; use warnings; use Math::BigRat;
my ($S, $P, $V);
sub evalx { my ($str, $vals) = @_; ($S = $str) =~ s/\s+//g; $P = 0; $V = $vals; my $r = expr(); die "parse error at $P in '$S'" if $P < length $S; return $r; }
sub peek { return substr($S, $P, 1); }
sub expr { my $v = term(); while ($P < length $S) { my $c = peek(); if ($c eq '+') { $P++; $v = $v + term(); } elsif ($c eq '-') { $P++; $v = $v - term(); } else { last } } return $v; }
sub term { my $v = unary(); while ($P < length $S) { my $c = peek(); if ($c eq '*') { $P++; $v = $v * unary(); } elsif ($c eq '/') { $P++; $v = $v / unary(); } else { last } } return $v; }
sub unary { my $c = peek(); if ($c eq '-') { $P++; return -unary(); } if ($c eq '+') { $P++; return unary(); } return power(); }
sub power { my $b = atom(); if ($P < length $S && peek() eq '^') { $P++; my $e = unary(); die "non-integer exponent" unless $e->is_int(); my $n = $e->numify(); return $n >= 0 ? $b ** $n : Math::BigRat->new(1) / ($b ** (-$n)); } return $b; }
sub atom { my $c = peek();
  if ($c eq '(') { $P++; my $v = expr(); die "missing ) at $P" unless peek() eq ')'; $P++; return $v; }
  if ($c =~ /\d/) { my ($num) = substr($S, $P) =~ /^(\d+)/; $P += length $num; return Math::BigRat->new($num); }
  if ($c =~ /[A-Za-z]/) { my ($name) = substr($S, $P) =~ /^([A-Za-z]\w*)/; $P += length $name;
    if ($P < length $S && peek() eq '[') { my ($idx) = substr($S, $P) =~ /^\[(\d+)\]/; $P += length($idx) + 2; $name .= "[$idx]"; }
    die "unknown symbol $name" unless exists $V->{$name}; return Math::BigRat->new($V->{$name}); }
  die "unexpected '$c' at $P"; }
1;
