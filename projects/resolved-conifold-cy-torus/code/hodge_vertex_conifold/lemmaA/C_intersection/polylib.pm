package polylib;
# Minimal exact multivariate polynomial library (Math::BigRat coefficients) and a parser for
# Mathematica InputForm expressions (+ - * / ^, integers, variables, parentheses).
# Variables listed in %$env are substituted by the given rational values while parsing;
# all other variables stay symbolic. Division is allowed only by constants (after substitution).
# Polynomial = hash ref { monomial_key => Math::BigRat }, monomial_key = "v1^e1*v2^e2" (sorted), "" = 1.
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(pconst padd psub pmul pscale ppow pis_zero pcoef pstr parse_expr mono_parse mono_key pvars peval_at);

sub R { return Math::BigRat->new($_[0]); }
sub mono_parse { my ($k) = @_; my %e; return \%e if $k eq ''; for my $f (split /\*/, $k) { my ($v,$p) = split /\^/, $f; $e{$v} += (defined $p ? $p : 1); } return \%e; }
sub mono_key { my ($e) = @_; my @f; for my $v (sort keys %$e) { next unless $e->{$v}; push @f, ($e->{$v} == 1 ? $v : "$v^$e->{$v}"); } return join('*', @f); }
sub pconst { my ($c) = @_; my $r = R($c); return $r->is_zero ? {} : { '' => $r }; }
sub pvar { my ($v) = @_; return { $v => R(1) }; }
sub padd { my ($a, $b) = @_; my %r = %$a; for my $k (keys %$b) { my $s = (exists $r{$k} ? $r{$k} + $b->{$k} : $b->{$k}->copy); if ($s->is_zero) { delete $r{$k}; } else { $r{$k} = $s; } } return \%r; }
sub pscale { my ($a, $c) = @_; $c = R($c) unless ref $c; return {} if $c->is_zero; my %r; for my $k (keys %$a) { $r{$k} = $a->{$k} * $c; } return \%r; }
sub psub { my ($a, $b) = @_; return padd($a, pscale($b, -1)); }
sub mkey_mul { my ($k1, $k2) = @_; return $k2 if $k1 eq ''; return $k1 if $k2 eq ''; my $e1 = mono_parse($k1); my $e2 = mono_parse($k2); for my $v (keys %$e2) { $e1->{$v} += $e2->{$v}; } return mono_key($e1); }
my %mcache;
sub pmul { my ($a, $b) = @_; my %r; for my $k1 (keys %$a) { for my $k2 (keys %$b) { my $ck = "$k1|$k2"; my $k = $mcache{$ck} //= mkey_mul($k1, $k2); my $p = $a->{$k1} * $b->{$k2}; if (exists $r{$k}) { $r{$k} += $p; delete $r{$k} if $r{$k}->is_zero; } else { $r{$k} = $p; } } } return \%r; }
sub ppow { my ($a, $n) = @_; die "negative power of polynomial" if $n < 0; my $r = pconst(1); for (1..$n) { $r = pmul($r, $a); } return $r; }
sub pis_zero { my ($a) = @_; return !%$a; }
sub pcoef { my ($a, $k) = @_; return exists $a->{$k} ? $a->{$k} : R(0); }
sub pvars { my ($a) = @_; my %v; for my $k (keys %$a) { my $e = mono_parse($k); $v{$_} = 1 for keys %$e; } return sort keys %v; }
sub pstr { my ($a) = @_; return '0' unless %$a; return join(' + ', map { "(" . $a->{$_}->bstr . ")" . ($_ eq '' ? '' : "*$_") } sort keys %$a); }
# evaluate polynomial at values for some variables (hash), returns polynomial in the others
sub peval_at { my ($a, $vals) = @_; my %r; for my $k (keys %$a) { my $e = mono_parse($k); my $c = $a->{$k}->copy; my %rest; for my $v (keys %$e) { if (exists $vals->{$v}) { $c *= R($vals->{$v}) ** $e->{$v}; } else { $rest{$v} = $e->{$v}; } } next if $c->is_zero; my $nk = mono_key(\%rest); if (exists $r{$nk}) { $r{$nk} += $c; delete $r{$nk} if $r{$nk}->is_zero; } else { $r{$nk} = $c; } } return \%r; }

# ---------------- parser ----------------
my ($S, $P, $ENV);
sub parse_expr { my ($str, $env) = @_; $S = $str; $S =~ s/\s+//g; $P = 0; $ENV = $env || {}; my $r = p_expr(); die "parse error at $P: " . substr($S, $P, 30) if $P < length($S); return $r; }
sub peek { return $P < length($S) ? substr($S, $P, 1) : ''; }
sub p_expr { my $r; if (peek() eq '-') { $P++; $r = pscale(p_term(), -1); } elsif (peek() eq '+') { $P++; $r = p_term(); } else { $r = p_term(); } while (1) { my $c = peek(); if ($c eq '+') { $P++; $r = padd($r, p_term()); } elsif ($c eq '-') { $P++; $r = psub($r, p_term()); } else { last; } } return $r; }
sub p_term { my $r = p_unary(); while (1) { my $c = peek(); if ($c eq '*') { $P++; $r = pmul($r, p_unary()); } elsif ($c eq '/') { $P++; my $d = p_unary(); my @k = keys %$d; die "division by non-constant" unless (@k == 1 && $k[0] eq ''); $r = pscale($r, R(1) / $d->{''}); } else { last; } } return $r; }
sub p_unary { if (peek() eq '-') { $P++; return pscale(p_unary(), -1); } return p_power(); }
sub p_power { my $b = p_base(); if (peek() eq '^') { $P++; my $neg = 0; my $par = 0; if (peek() eq '(') { $P++; $par = 1; } if (peek() eq '-') { $neg = 1; $P++; } my $n = ''; while (peek() =~ /\d/) { $n .= peek(); $P++; } if ($par) { die "expected )" unless peek() eq ')'; $P++; } if ($neg) { my @k = keys %$b; die "negative power of non-constant" unless (@k == 1 && $k[0] eq ''); return pconst(R(1) / ($b->{''} ** $n)); } return ppow($b, $n); } return $b; }
sub p_base { my $c = peek(); if ($c eq '(') { $P++; my $r = p_expr(); die "expected ) at $P" unless peek() eq ')'; $P++; return $r; } if ($c =~ /\d/) { my $n = ''; while (peek() =~ /\d/) { $n .= peek(); $P++; } return pconst($n); } if ($c =~ /[A-Za-z]/) { my $v = ''; while (peek() =~ /[A-Za-z0-9]/) { $v .= peek(); $P++; } if (peek() eq '[') { $P++; my $i = ''; while (peek() ne ']') { $i .= peek(); $P++; } $P++; $v .= $i; } if (exists $ENV->{$v}) { return pconst($ENV->{$v}); } return pvar($v); } die "unexpected char '$c' at $P"; }
1;
