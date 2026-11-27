package HodgeRef;
# Referee's own readers for the exact Hodge-side tables (no code shared with the author's HodgeSide.pm / MParse.pm).
#  - Mathematica tables lemmaA/data/q_g<g>_r<r>.txt  (Q(s; j) = Sum_k s^k/k! <Omega prod 1/(c - j_i psi_i) Tt^k>, Tt = psi^2/(c-psi),
#    weights (c,-x,-1), c = 1+x) and mv_g<g>_l<l>.txt (P(m) = <Omega prod 1/(c - m_i psi_i)>, no dilaton leaves).
#    With x = -tau/(1+tau) (c = 1/(1+tau)) the c = 1 values (weights (1,tau,-1-tau)) are c^(3g-3+2r) times the table values.
#  - Perl engine tables (Route E author: E_formulas/engine_h*_r*.txt; Route E referee: E_formulas/referee/ref_engine_h*_r*.txt),
#    lines "h | r | c | k1,..,kr:coef ; ...", Xi = Sum_c (-1/sigma)^c Xi^(c).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
use File::Basename (); use File::Spec ();
our $ROOT = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../..";

# ---------- tiny polynomial algebra: { monokey => BigRat }, monokey = "v1^e1*v2^e2" sorted, "" = constant ----------
sub mk { my %e = @_; return join('*', map { "$_^$e{$_}" } sort grep { $e{$_} } keys %e); }
sub mparse { my $k = shift; my %e; return %e if $k eq ''; for my $t (split /\*/, $k) { my ($v, $x) = split /\^/, $t; $e{$v} = $x; } return %e; }
sub padd { my ($A, $B, $sgn) = @_; my %C = %$A; for my $k (keys %$B) { my $v = $sgn ? -$B->{$k} : $B->{$k}; $C{$k} = exists $C{$k} ? $C{$k} + $v : $v; delete $C{$k} if $C{$k}->is_zero; } return \%C; }
sub pmul { my ($A, $B) = @_; my %C;
  for my $ka (keys %$A) { my %ea = mparse($ka);
    for my $kb (keys %$B) { my %e = %ea; my %eb = mparse($kb); $e{$_} = ($e{$_} // 0) + $eb{$_} for keys %eb;
      my $k = mk(%e); my $v = $A->{$ka} * $B->{$kb}; $C{$k} = exists $C{$k} ? $C{$k} + $v : $v; delete $C{$k} if $C{$k}->is_zero; } }
  return \%C; }
sub pconst { my $c = shift; return $c->is_zero ? {} : { '' => $c }; }
sub isconst { my $A = shift; for (keys %$A) { return 0 if $_ ne ''; } return 1; }
sub cval { my $A = shift; return $A->{''} // R(0); }

# ---------- recursive-descent parser for Mathematica InputForm (numbers, symbols, + - * / ^, parentheses) ----------
my (@tok, $pos, %subst);
sub tokenize { my $s = shift; $s =~ s/\s+//g; my @t;
  while (length $s) { if ($s =~ s/^(\d+)//) { push @t, ['n', $1]; } elsif ($s =~ s/^([A-Za-z][A-Za-z0-9]*)//) { push @t, ['v', $1]; }
    elsif ($s =~ s/^([-+*\/^()])//) { push @t, ['o', $1]; } else { die "tokenize: cannot parse near '" . substr($s, 0, 20) . "'"; } }
  return @t; }
sub peek { return $pos <= $#tok ? $tok[$pos] : ['e', '']; }
sub expr { my $v; my $p = peek();
  if ($p->[0] eq 'o' && ($p->[1] eq '-' || $p->[1] eq '+')) { $pos++; $v = term(); $v = padd({}, $v, 1) if $p->[1] eq '-'; } else { $v = term(); }
  while (1) { my $q = peek(); last unless $q->[0] eq 'o' && ($q->[1] eq '+' || $q->[1] eq '-'); $pos++; my $t = term(); $v = padd($v, $t, $q->[1] eq '-'); }
  return $v; }
sub term { my $v = factor();
  while (1) { my $q = peek(); last unless $q->[0] eq 'o' && ($q->[1] eq '*' || $q->[1] eq '/'); $pos++; my $f = factor();
    if ($q->[1] eq '*') { $v = pmul($v, $f); } else { die "division by non-constant" unless isconst($f); my $c = cval($f); die "division by 0" if $c->is_zero;
      my %C; $C{$_} = $v->{$_} / $c for keys %$v; $v = \%C; } }
  return $v; }
sub factor { my $q = peek(); if ($q->[0] eq 'o' && $q->[1] eq '-') { $pos++; return padd({}, factor(), 1); } return power(); }
sub power { my $b = atom(); my $q = peek();
  if ($q->[0] eq 'o' && $q->[1] eq '^') { $pos++; my $e = factor(); die "non-constant exponent" unless isconst($e); my $n = cval($e); die "non-integer exponent" unless $n->is_int;
    $n = $n->numify;
    if ($n < 0) { die "negative power of non-constant" unless isconst($b); return pconst(cval($b) ** $n); }
    my $r = pconst(R(1)); $r = pmul($r, $b) for 1..$n; return $r; }
  return $b; }
sub atom { my $q = $tok[$pos++] // die "unexpected end";
  if ($q->[0] eq 'n') { return pconst(R($q->[1])); }
  if ($q->[0] eq 'v') { return exists $subst{$q->[1]} ? pconst($subst{$q->[1]}) : { mk($q->[1] => 1) => R(1) }; }
  if ($q->[0] eq 'o' && $q->[1] eq '(') { my $v = expr(); my $c = $tok[$pos++]; die "expected )" unless $c && $c->[1] eq ')'; return $v; }
  die "unexpected token $q->[1]"; }
sub parse_poly { my ($str, $sub) = @_; @tok = tokenize($str); $pos = 0; %subst = %$sub; my $v = expr(); die "trailing tokens" if $pos <= $#tok; return $v; }
sub peval { my ($P, $vals) = @_; my $tot = R(0);
  for my $k (keys %$P) { my %e = mparse($k); my $t = $P->{$k}->copy; for my $v (keys %e) { die "no value for $v" unless exists $vals->{$v}; $t *= $vals->{$v} ** $e{$v}; } $tot += $t; }
  return $tot; }

# ---------- table access ----------
my %cache;
sub table_expr { my ($file, $g, $l) = @_; open my $fh, '<', $file or return undef; local $/; my $txt = <$fh>; close $fh; $txt =~ s/\r//g;
  for my $ln (split /\n/, $txt) { next unless $ln =~ /\S/; my ($gg, $ll, $e) = split /\s*\|\s*/, $ln, 3; return $e if $gg == $g && $ll == $l; }
  return undef; }
# c = 1 polynomial in j1..jr (kind 'q', at the given s) or m1..ml (kind 'mv')
sub data_poly { my ($kind, $g, $l, $tau, $s) = @_; my $key = join('|', $kind, $g, $l, "$tau", defined $s ? "$s" : '');
  return $cache{$key} if exists $cache{$key};
  my $file = $kind eq 'q' ? "$ROOT/data/q_g${g}_r${l}.txt" : "$ROOT/data/mv_g${g}_l${l}.txt";
  my $e = table_expr($file, $g, $l); return ($cache{$key} = undef) unless defined $e;
  my $T = R($tau); my $x = -$T / (1 + $T); my $c = 1 + $x;
  my %sub = (x => $x); $sub{s} = R($s) if $kind eq 'q';
  my $P = parse_poly($e, \%sub); my $f = $c ** (3*$g - 3 + 2*$l); my %Q; $Q{$_} = $P->{$_} * $f for keys %$P;
  $cache{$key} = \%Q; return \%Q; }
sub eval_data { my ($kind, $g, $mu, $tau, $s) = @_; my $l = scalar @$mu; my $P = data_poly($kind, $g, $l, $tau, $s); return undef unless $P;
  my $pre = $kind eq 'q' ? 'j' : 'm'; my %v; $v{"$pre" . ($_ + 1)} = R($mu->[$_]) for 0..$#$mu; $v{s} = R($s) if $kind eq 'q' && 0;
  return peval($P, \%v); }
# engine tables: returns { c => { "k1,..,kr" => coef } }
sub engine { my ($h, $r, $which) = @_; my $key = "eng|$h|$r|$which"; return $cache{$key} if exists $cache{$key};
  my $file = $which eq 'ref' ? "$ROOT/E_formulas/referee/ref_engine_h${h}_r${r}.txt" : "$ROOT/E_formulas/engine_h${h}_r${r}.txt";
  open my $fh, '<', $file or return ($cache{$key} = undef); my %comp;
  while (my $ln = <$fh>) { $ln =~ s/\r//g; chomp $ln; next unless $ln =~ /\S/; my ($hh, $rr, $c, $body) = split /\s*\|\s*/, $ln, 4;
    $comp{$c} //= {}; next unless defined $body && $body =~ /\S/;
    for my $t (split /\s*;\s*/, $body) { my ($k, $v) = split /:/, $t; $comp{$c}{$k} = R($v); } }
  close $fh; $cache{$key} = \%comp; return \%comp; }
sub eval_engine { my ($h, $mu, $sigma, $which) = @_; my $E = engine($h, scalar(@$mu), $which); return undef unless $E;
  my $tot = R(0); my $ms = R(-1) / $sigma;
  for my $c (keys %$E) { my $sub = R(0);
    for my $k (keys %{$E->{$c}}) { my @k = split /,/, $k; my $t = $E->{$c}{$k}->copy; $t *= R($mu->[$_]) ** $k[$_] for 0..$#k; $sub += $t; }
    $tot += $sub * ($ms ** $c); }
  return $tot; }
1;
