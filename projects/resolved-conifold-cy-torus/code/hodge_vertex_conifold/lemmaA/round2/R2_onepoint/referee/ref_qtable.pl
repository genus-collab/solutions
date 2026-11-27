#!/usr/bin/perl
# ref_qtable.pl TAU G R JMAX [PIDX] : evaluate lemmaA/data/q_g<G>_r<R>.txt (Mathematica, Mumford's formula; c-normalisation, symbolic x)
# directly at x = -tau/(1+tau), s = -1, and convert to c = 1:  Xi^(c=1) = c^{3G-3+2R} Q(x; s=-1; j),  c = 1 + x.
# Prints "j1,..,jR value" lines (mod p) for all descending tuples with parts <= JMAX; used by ref_cmp_q.sh against ref_xy_mod.pl.
use strict; use warnings;
my ($TAUS, $G, $R, $JMAX, $PIDX) = @ARGV; $PIDX //= 0;
sub isprime { my $q = shift; for (my $d = 2; $d * $d <= $q; $d++) { return 0 if $q % $d == 0; } return 1; }
my @PR; { my $q = 2147483647; while (@PR < 4) { push @PR, $q if isprime($q); $q -= 2; } }
my $P = $PR[$PIDX];
sub mm { return ($_[0] * $_[1]) % $P; }
sub ad { my $s = $_[0] + $_[1]; return $s >= $P ? $s - $P : $s; }
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub iv { my $a = shift; die "division by zero mod p" if $a % $P == 0; return pw($a, $P - 2); }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; $s =~ s/\s+//g; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), iv(decmod($2))); } return decmod($s); }
sub im { my $k = shift; return $k >= 0 ? $k % $P : ng((-$k) % $P); }
my $file = "../../../data/q_g${G}_r${R}.txt";
open my $fh, '<', $file or die "no $file"; local $/; my $txt = <$fh>; close $fh;
$txt =~ s/\s+//g; my @f = split /\|/, $txt; die "bad header" unless $f[0] == $G && $f[1] == $R; my $expr = $f[2];
# tokenizer + recursive descent: expr := term (('+'|'-') term)* ; term := factor (('*'|'/') factor)* ; factor := unary ('^' integer)?
my @tok = ($expr =~ /(\d+|[A-Za-z]\w*|[-+*\/^()])/g);
my %VAL; my $pos;
sub peek { return $pos < @tok ? $tok[$pos] : ''; }
sub expr_ { my $v = term_(); while (peek() eq '+' || peek() eq '-') { my $op = $tok[$pos++]; my $w = term_(); $v = $op eq '+' ? ad($v, $w) : ad($v, ng($w)); } return $v; }
sub term_ { my $v = factor_(); while (peek() eq '*' || peek() eq '/') { my $op = $tok[$pos++]; my $w = factor_(); $v = $op eq '*' ? mm($v, $w) : mm($v, iv($w)); } return $v; }
sub factor_ { my $t = $tok[$pos++]; my $v;
  if ($t eq '-') { return ng(factor_()); }
  if ($t eq '(') { $v = expr_(); die "missing )" unless $tok[$pos++] eq ')'; }
  elsif ($t =~ /^\d+$/) { $v = decmod($t); }
  elsif (exists $VAL{$t}) { $v = $VAL{$t}; } else { die "unknown symbol $t"; }
  if (peek() eq '^') { $pos++; my $e = $tok[$pos++]; die "non-integer power" unless $e =~ /^\d+$/; $v = pw($v, $e); }
  return $v; }
my $tau = ratmod($TAUS); my $x = mm(ng($tau), iv(ad($tau, 1))); my $c = ad(1, $x);
my @tuples; { my $rec; $rec = sub { my ($pre, $max) = @_; if (@$pre == $R) { push @tuples, [@$pre]; return; } for my $y (reverse 1..$max) { $rec->([@$pre, $y], $y); } }; $rec->([], $JMAX); }
for my $t (@tuples) { %VAL = (x => $x, s => ng(1)); $VAL{"j$_"} = im($t->[$_-1]) for 1..$R; $pos = 0; my $v = expr_(); die "trailing tokens" if $pos != @tok;
  $v = mm($v, pw($c, 3*$G - 3 + 2*$R)); print join(',', @$t), " $v\n"; }
