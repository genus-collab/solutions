package RParse;
# Referee's evaluator for Mathematica InputForm rational expressions in symbols (j1, j2, s, x, ...), exact (Math::BigRat).
# Tokenizer + recursive descent; integer literals are read as BigInt strings (no floating point).
use strict; use warnings; use Math::BigRat;
my (@tok, $pos, $env);
sub tokenize { my $s = shift; $s =~ s/\s+//g; my @t;
  while (length $s) {
    if ($s =~ s/^(\d+)//) { push @t, ['n', $1] }
    elsif ($s =~ s/^([A-Za-z]\w*)//) { push @t, ['v', $1] }
    elsif ($s =~ s/^([-+*\/^()])//) { push @t, ['o', $1] }
    else { die "RParse: bad character at '" . substr($s, 0, 10) . "'" } }
  return @t }
my %TC;
sub evalx { my ($str, $vals) = @_; my $tc = $TC{$str} //= [ tokenize($str) ]; @tok = @$tc; $pos = 0; $env = $vals; my $r = sum_(); die "RParse: trailing tokens" if $pos < @tok; return $r }
sub peek { return $pos < @tok ? $tok[$pos] : ['e', ''] }
sub sum_ { my $v = prod_(); while (1) { my $t = peek(); last unless $t->[0] eq 'o' && ($t->[1] eq '+' || $t->[1] eq '-'); $pos++; my $w = prod_(); $v = $t->[1] eq '+' ? $v + $w : $v - $w } return $v }
sub prod_ { my $v = un_(); while (1) { my $t = peek(); last unless $t->[0] eq 'o' && ($t->[1] eq '*' || $t->[1] eq '/'); $pos++; my $w = un_(); $v = $t->[1] eq '*' ? $v * $w : $v / $w } return $v }
sub un_ { my $t = peek(); if ($t->[0] eq 'o' && $t->[1] eq '-') { $pos++; return -un_() } if ($t->[0] eq 'o' && $t->[1] eq '+') { $pos++; return un_() } return pow_() }
sub pow_ { my $b = at_(); my $t = peek(); if ($t->[0] eq 'o' && $t->[1] eq '^') { $pos++; my $e = un_(); die "RParse: non-integer exponent" unless $e->is_int();
    my $n = $e->numify(); my $r = Math::BigRat->new(1); my $base = $n >= 0 ? $b : Math::BigRat->new(1) / $b; $r *= $base for 1..abs($n); return $r } return $b }
sub at_ { my $t = peek(); $pos++;
  if ($t->[0] eq 'n') { return Math::BigRat->new($t->[1]) }
  if ($t->[0] eq 'v') { die "RParse: unknown symbol $t->[1]" unless exists $env->{$t->[1]}; return Math::BigRat->new($env->{$t->[1]}) }
  if ($t->[0] eq 'o' && $t->[1] eq '(') { my $v = sum_(); my $c = peek(); die "RParse: missing )" unless $c->[0] eq 'o' && $c->[1] eq ')'; $pos++; return $v }
  die "RParse: unexpected token $t->[1]" }
# read a q table file "g | r | expr"; returns the expression string
sub qexpr { my $f = shift; open my $fh, '<', $f or die "no $f"; local $/; my $s = <$fh>; close $fh; my ($g, $r, $e) = split /\s*\|\s*/, $s, 3; $e =~ s/\s+$//; return ($g, $r, $e) }
1;
