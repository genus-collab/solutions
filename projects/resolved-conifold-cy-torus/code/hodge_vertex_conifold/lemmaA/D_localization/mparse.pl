# mparse.pl -- evaluate a Mathematica InputForm rational expression (numbers, symbols, + - * / ^, parentheses)
# at given rational values of the symbols, exactly (Math::BigRat). Usage: meval($string, {x => Q(2), m1 => Q(3)}).
use strict; use warnings;
use Math::BigRat;
sub meval {
  my ($s, $vals) = @_; $s =~ s/\s+//g;
  my @tok = ($s =~ /(\d+|[A-Za-z][A-Za-z0-9]*|[-+*\/^()])/g);
  my $i = 0;
  my ($expr, $term, $factor, $power, $atom);
  $expr = sub { my $v = $term->(); while ($i < @tok && ($tok[$i] eq '+' || $tok[$i] eq '-')) { my $op = $tok[$i++]; my $w = $term->(); $v = $op eq '+' ? $v + $w : $v - $w; } return $v; };
  $term = sub { my $v = $factor->(); while ($i < @tok && ($tok[$i] eq '*' || $tok[$i] eq '/')) { my $op = $tok[$i++]; my $w = $factor->(); $v = $op eq '*' ? $v * $w : $v / $w; } return $v; };
  $factor = sub { if ($tok[$i] eq '-') { $i++; return -$factor->(); } if ($tok[$i] eq '+') { $i++; return $factor->(); } return $power->(); };
  $power = sub { my $b = $atom->(); if ($i < @tok && $tok[$i] eq '^') { $i++; my $e = $factor->(); die "non-integer exponent" unless $e->is_int(); my $n = $e->numify(); return $n >= 0 ? $b**$n : Math::BigRat->new(1) / $b**(-$n); } return $b; };
  $atom = sub { my $t = $tok[$i++];
    if ($t eq '(') { my $v = $expr->(); die "missing )" unless $tok[$i++] eq ')'; return $v; }
    if ($t =~ /^\d+$/) { return Math::BigRat->new($t); }
    die "unknown symbol $t" unless exists $vals->{$t}; return Math::BigRat->new($vals->{$t}); };
  my $v = $expr->(); die "trailing tokens" if $i < @tok; return $v;
}
1;
