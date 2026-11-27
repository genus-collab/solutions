package rparse;
# rparse.pm -- referee's own evaluator for the Mathematica InputForm tables in lemmaA/data (written independently of mvlib.pm).
# Evaluates an expression exactly (Math::BigRat) at given rational values of its variables; or modulo a prime.
# Grammar: sum := ['+'|'-'] prod (('+'|'-') prod)* ; prod := pow (('*'|'/') pow)* ; pow := atom ['^' ['-'] int | '^' '(' ['-'] int ')'] ;
#          atom := int | name | '(' sum ')' | '-' pow
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(read_table eval_expr eval_expr_mod);

sub read_table {   # returns list of [field1, field2, ..., expr] split on " | "
  my ($file) = @_; open my $fh, '<', $file or die "cannot open $file: $!"; local $/; my $t = <$fh>; close $fh;
  my @rows; for my $line (split /\r?\n/, $t) { next unless $line =~ /\S/; push @rows, [ split /\s*\|\s*/, $line ]; } return @rows; }

sub toks { my $s = shift; my @t; pos($s) = 0;
  while (pos($s) < length($s)) { if ($s =~ /\G\s+/gc) { next; } if ($s =~ /\G(\d+)/gc) { push @t, ['n', $1]; next; }
    if ($s =~ /\G([A-Za-z]\w*)/gc) { push @t, ['v', $1]; next; } if ($s =~ /\G([-+*\/^()])/gc) { push @t, ['o', $1]; next; }
    die "rparse: bad character at " . pos($s) . ": '" . substr($s, pos($s), 20) . "'"; }
  return \@t; }

# generic evaluator with arithmetic callbacks
sub _eval { my ($expr, $vals, $A) = @_; my $t = toks($expr); my $i = 0;
  my ($sum, $prod, $pow, $atom);
  my $peek = sub { return $i < @$t ? $t->[$i] : ['e', '']; };
  $sum = sub { my $v; my $tk = $peek->();
    if ($tk->[0] eq 'o' && ($tk->[1] eq '+' || $tk->[1] eq '-')) { $i++; my $p = $prod->(); $v = $tk->[1] eq '-' ? $A->{neg}->($p) : $p; } else { $v = $prod->(); }
    while (1) { $tk = $peek->(); last unless $tk->[0] eq 'o' && ($tk->[1] eq '+' || $tk->[1] eq '-'); $i++; my $p = $prod->();
      $v = $tk->[1] eq '+' ? $A->{add}->($v, $p) : $A->{add}->($v, $A->{neg}->($p)); }
    return $v; };
  $prod = sub { my $v = $pow->();
    while (1) { my $tk = $peek->(); last unless $tk->[0] eq 'o' && ($tk->[1] eq '*' || $tk->[1] eq '/'); $i++; my $p = $pow->();
      $v = $tk->[1] eq '*' ? $A->{mul}->($v, $p) : $A->{mul}->($v, $A->{inv}->($p)); }
    return $v; };
  $pow = sub { my $b = $atom->(); my $tk = $peek->();
    if ($tk->[0] eq 'o' && $tk->[1] eq '^') { $i++; my $sg = 1; my $par = 0; $tk = $peek->();
      if ($tk->[0] eq 'o' && $tk->[1] eq '(') { $par = 1; $i++; $tk = $peek->(); }
      if ($tk->[0] eq 'o' && $tk->[1] eq '-') { $sg = -1; $i++; $tk = $peek->(); }
      die "rparse: exponent must be an integer" unless $tk->[0] eq 'n'; my $e = $tk->[1]; $i++;
      if ($par) { $tk = $peek->(); die "rparse: missing )" unless $tk->[0] eq 'o' && $tk->[1] eq ')'; $i++; }
      my $r = $A->{one}->(); for (1..$e) { $r = $A->{mul}->($r, $b); } $r = $A->{inv}->($r) if $sg < 0; return $r; }
    return $b; };
  $atom = sub { my $tk = $peek->(); $i++;
    if ($tk->[0] eq 'n') { return $A->{num}->($tk->[1]); }
    if ($tk->[0] eq 'v') { die "rparse: no value for variable $tk->[1]" unless exists $vals->{$tk->[1]}; return $A->{num}->($vals->{$tk->[1]}); }
    if ($tk->[0] eq 'o' && $tk->[1] eq '(') { my $v = $sum->(); my $c = $peek->(); die "rparse: expected )" unless $c->[0] eq 'o' && $c->[1] eq ')'; $i++; return $v; }
    if ($tk->[0] eq 'o' && $tk->[1] eq '-') { return $A->{neg}->($pow->()); }
    die "rparse: unexpected token $tk->[1]"; };
  my $v = $sum->(); die "rparse: trailing input" if $i < @$t; return $v; }

my %BR = ( num => sub { Math::BigRat->new($_[0]) }, one => sub { Math::BigRat->new(1) },
  add => sub { $_[0] + $_[1] }, mul => sub { $_[0] * $_[1] }, neg => sub { -$_[0] }, inv => sub { die "rparse: division by zero" if $_[0]->is_zero; Math::BigRat->new(1) / $_[0] } );
sub eval_expr { my ($e, $vals) = @_; return _eval($e, $vals, \%BR); }

sub eval_expr_mod { my ($e, $vals, $p) = @_;   # values: integers mod p (or strings of integers)
  my $md = sub { my $x = $_[0] % $p; $x += $p if $x < 0; $x };
  my $inv = sub { my $a = $md->($_[0]); die "rparse: non-invertible mod p" if $a == 0; my ($t, $nt, $r, $nr) = (0, 1, $p, $a);
    while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q * $nt); ($r, $nr) = ($nr, $r - $q * $nr); } return $md->($t); };
  my %A = ( num => sub { my $s = "$_[0]"; if (length($s) > 15) { return Math::BigInt->new($s)->bmod($p)->numify; } return $md->($s); }, one => sub { 1 },
    add => sub { $md->($_[0] + $_[1]) }, mul => sub { ($_[0] * $_[1]) % $p }, neg => sub { $md->(-$_[0]) }, inv => $inv );
  return _eval($e, $vals, \%A); }
1;
