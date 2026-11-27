package HB;
# Library for route B (Hurwitz component of Lemma A).  Exact arithmetic only (Math::BigInt / Math::BigRat).
#  - KW intersection numbers <tau_a1 ... tau_an>_g (DVV recursion), used to evaluate the wk[g,{..}] symbols of the data files
#  - a small evaluator for Mathematica InputForm expressions (numbers, symbols, + - * / ^, wk[g,{..}])
#  - partitions, dimensions (hook formula), content sums, Murnaghan-Nakayama characters chi^lambda(nu u 1^m)
#  - disconnected and connected Hurwitz counts with labelled non-unit cycles (exponential formula)
#  - ELSV-normalised integrals A_h(R,1^c) = Int_{Mbar_{h,r+c}} Lambda^vee(1) prod 1/(1-j_i psi_i) prod_{k<=c} 1/(1-psi)
#  - the finite (Laguerre) formula for Xi^(0)_h(R)
use strict; use warnings;
use Math::BigInt; use Math::BigRat;
our @EXPORT_OK;

sub BR { Math::BigRat->new(shift) }
sub BI { Math::BigInt->new(shift) }

my @FACT = (BI(1));
sub fact { my $n = shift; die "fact neg" if $n < 0; while (@FACT <= $n) { push @FACT, $FACT[-1] * scalar(@FACT) } return $FACT[$n]->copy }
my %BINOM;
sub binom { my ($n, $k) = @_; return BI(0) if $k < 0 || $k > $n || $n < 0;
  my $key = "$n,$k"; return $BINOM{$key}->copy if exists $BINOM{$key};
  my $r = fact($n) / (fact($k) * fact($n - $k)); $BINOM{$key} = $r; return $r->copy }

# ---------------------------------------------------------------- KW numbers
my %KW;
sub dfact { my $n = shift; my $r = BI(1); for (my $i = $n; $i > 1; $i -= 2) { $r *= $i } return $r }
sub kw { my ($g, @a) = @_;
  for (@a) { return BR(0) if $_ < 0 }
  my $n = scalar @a;
  return BR(0) if $g < 0 || 2*$g - 2 + $n <= 0;
  my $s = 0; $s += $_ for @a;
  return BR(0) if $s != 3*$g - 3 + $n;
  @a = sort { $a <=> $b } @a;
  my $key = "$g:" . join(",", @a);
  return $KW{$key}->copy if exists $KW{$key};
  my $res;
  if ($g == 0 && $n == 3) { $res = BR(1) }
  elsif ($g == 1 && $n == 1) { $res = BR("1/24") }
  elsif ($a[0] == 0) {                       # string equation
    shift @a; $res = BR(0);
    for my $i (0..$#a) { next if $a[$i] == 0; my @b = @a; $b[$i]--; $res += kw($g, @b) }
  }
  elsif ($a[0] == 1) {                       # dilaton equation
    shift @a; $res = BR(2*$g - 2 + $n - 1) * kw($g, @a);
  }
  else {                                     # DVV on the largest index tau_{k+1}
    my $kp1 = pop @a; my $k = $kp1 - 1; my $m = scalar @a;
    my $tot = BR(0);
    for my $i (0..$m-1) { my @b = @a; my $aj = $b[$i]; $b[$i] = $aj + $k;
      $tot += BR(dfact(2*$k + 2*$aj + 1)) / BR(dfact(2*$aj - 1)) * kw($g, @b) }
    for my $r (0..$k-1) { my $s2 = $k - 1 - $r;
      my $c = BR(dfact(2*$r + 1) * dfact(2*$s2 + 1)) / 2;
      $tot += $c * kw($g - 1, $r, $s2, @a) if $g >= 1;
      for my $mask (0..(1 << $m) - 1) { my (@I, @J);
        for my $t (0..$m-1) { if ($mask & (1 << $t)) { push @I, $a[$t] } else { push @J, $a[$t] } }
        for my $g1 (0..$g) { my $v1 = kw($g1, $r, @I); next if $v1 == 0; $tot += $c * $v1 * kw($g - $g1, $s2, @J) } }
    }
    $res = $tot / BR(dfact(2*$k + 3));
  }
  $KW{$key} = $res; return $res->copy;
}

# ---------------------------------------------------------------- InputForm parser / evaluator
# AST nodes: ['n',BigRat] ['v',name] ['+',[sign,node]...] ['*',[pow(+1/-1),node]...] ['^',base,exp] ['wk',g,[list]] ['{',[nodes]]
sub tokenize { my $s = shift; my @t;
  while ($s =~ /\G\s*(\d+|[A-Za-z][A-Za-z0-9]*|[-+*\/^()\[\]{},])/gc) { push @t, $1 }
  $s =~ /\G\s*/gc; die "tokenize: trailing garbage at " . pos($s) . ": " . substr($s, pos($s) // 0, 40) if (pos($s) // 0) < length($s);
  return \@t }
sub parse { my $s = shift; my $t = tokenize($s); my $p = 0; my $ast = p_expr($t, \$p); die "parse: leftover tokens at $p ($t->[$p])" if $p < @$t; return $ast }
sub p_expr { my ($t, $p) = @_; my @terms; my $sign = 1;
  if ($$p < @$t && ($t->[$$p] eq '-' || $t->[$$p] eq '+')) { $sign = $t->[$$p] eq '-' ? -1 : 1; $$p++ }
  push @terms, [$sign, p_term($t, $p)];
  while ($$p < @$t && ($t->[$$p] eq '+' || $t->[$$p] eq '-')) { my $sg = $t->[$$p] eq '-' ? -1 : 1; $$p++; push @terms, [$sg, p_term($t, $p)] }
  return @terms == 1 && $terms[0][0] == 1 ? $terms[0][1] : ['+', @terms] }
sub p_term { my ($t, $p) = @_; my @f = ([1, p_unary($t, $p)]);
  while ($$p < @$t && ($t->[$$p] eq '*' || $t->[$$p] eq '/')) { my $pw = $t->[$$p] eq '/' ? -1 : 1; $$p++; push @f, [$pw, p_unary($t, $p)] }
  return @f == 1 && $f[0][0] == 1 ? $f[0][1] : ['*', @f] }
sub p_unary { my ($t, $p) = @_;
  if ($t->[$$p] eq '-') { $$p++; return ['*', [1, ['n', BR(-1)]], [1, p_unary($t, $p)]] }
  return p_pow($t, $p) }
sub p_pow { my ($t, $p) = @_; my $b = p_prim($t, $p);
  if ($$p < @$t && $t->[$$p] eq '^') { $$p++; my $e = p_unary($t, $p); return ['^', $b, $e] }
  return $b }
sub p_prim { my ($t, $p) = @_; my $tok = $t->[$$p];
  die "parse: unexpected end" unless defined $tok;
  if ($tok =~ /^\d+$/) { $$p++; return ['n', BR($tok)] }
  if ($tok eq '(') { $$p++; my $e = p_expr($t, $p); die "parse: expected )" unless $t->[$$p] eq ')'; $$p++; return $e }
  if ($tok eq '{') { $$p++; my @l; if ($t->[$$p] ne '}') { push @l, p_expr($t, $p); while ($t->[$$p] eq ',') { $$p++; push @l, p_expr($t, $p) } }
    die "parse: expected }" unless $t->[$$p] eq '}'; $$p++; return ['{', @l] }
  if ($tok =~ /^[A-Za-z]/) { $$p++;
    if ($$p < @$t && $t->[$$p] eq '[') { $$p++; my @args = (p_expr($t, $p)); while ($t->[$$p] eq ',') { $$p++; push @args, p_expr($t, $p) }
      die "parse: expected ]" unless $t->[$$p] eq ']'; $$p++;
      if ($tok eq 'wk') { my $g = ev_const($args[0]); my @l = map { ev_const($_) } @{$args[1]}[1..$#{$args[1]}];
        return ['n', kw($g, @l)] }
      die "parse: unknown function $tok" }
    return ['v', $tok] }
  die "parse: unexpected token $tok" }
sub ev_const { my $n = shift; my $v = evaluate($n, {}); die "not integer" unless $v->is_int; return $v->numify }
sub evaluate { my ($n, $env) = @_; my $k = $n->[0];
  if ($k eq 'n') { return $n->[1]->copy }
  if ($k eq 'v') { die "unbound variable $n->[1]" unless exists $env->{$n->[1]}; return BR($env->{$n->[1]}) }
  if ($k eq '+') { my $s = BR(0); for my $i (1..$#$n) { my $v = evaluate($n->[$i][1], $env); if ($n->[$i][0] > 0) { $s += $v } else { $s -= $v } } return $s }
  if ($k eq '*') { my $num = BR(1); my $den = BR(1);
    for my $i (1..$#$n) { my $v = evaluate($n->[$i][1], $env); if ($n->[$i][0] > 0) { $num *= $v } else { $den *= $v } }
    die "division by zero" if $den == 0; return $num / $den }
  if ($k eq '^') { my $b = evaluate($n->[1], $env); my $e = evaluate($n->[2], $env); die "non-integer power" unless $e->is_int; my $ei = $e->numify;
    my $r = BR(1); my $bb = $ei < 0 ? BR(1) / $b : $b; for (1..abs($ei)) { $r *= $bb } return $r }
  die "evaluate: bad node $k" }
# simplify: fold constant subtrees (after wk evaluation) -- speeds up repeated evaluation
sub fold { my $n = shift; my $k = $n->[0];
  return $n if $k eq 'n' || $k eq 'v';
  if ($k eq '+' || $k eq '*') { my @c = map { [$_->[0], fold($_->[1])] } @{$n}[1..$#$n];
    if (!grep { $_->[1][0] ne 'n' } @c) { return ['n', evaluate([$k, @c], {})] }
    return [$k, @c] }
  if ($k eq '^') { my $b = fold($n->[1]); my $e = fold($n->[2]); return $b->[0] eq 'n' && $e->[0] eq 'n' ? ['n', evaluate(['^', $b, $e], {})] : ['^', $b, $e] }
  return $n }
sub read_record { my $file = shift; open my $fh, '<', $file or die "cannot open $file"; local $/; my $s = <$fh>; close $fh;
  my @recs; for my $line (split /\n/, $s) { next unless $line =~ /\S/; my @f = split /\s*\|\s*/, $line; push @recs, \@f } return \@recs }

# ---------------------------------------------------------------- partitions and characters
my %PARTS;
sub partitions { my $n = shift; return $PARTS{$n} if $PARTS{$n};
  my @out; my @a;
  my $rec; $rec = sub { my ($rem, $max) = @_; if ($rem == 0) { push @out, [@a]; return }
    for (my $k = ($max < $rem ? $max : $rem); $k >= 1; $k--) { push @a, $k; $rec->($rem - $k, $k); pop @a } };
  $rec->($n, $n); $PARTS{$n} = \@out; return \@out }
my %DIM;
sub dimlam { my $l = shift; my $key = join(",", @$l); return $DIM{$key}->copy if exists $DIM{$key};
  my $n = 0; $n += $_ for @$l; my @conj; for my $i (0..$#$l) { for my $j (0..$l->[$i]-1) { $conj[$j]++ } }
  my $hp = BI(1); for my $i (0..$#$l) { for my $j (0..$l->[$i]-1) { $hp *= ($l->[$i] - $j - 1) + ($conj[$j] - $i - 1) + 1 } }
  my $d = fact($n) / $hp; $DIM{$key} = $d; return $d->copy }
sub content_sum { my $l = shift; my $s = 0; for my $i (0..$#$l) { $s += $l->[$i] * ($l->[$i] - 1) / 2 - $i * $l->[$i] } return $s }
# chi^lambda(nu u 1^m) via Murnaghan-Nakayama on beta-sets; nu = list of hook sizes to remove (any order), then dim of the rest
my %CHI;
sub chi { my ($l, $nu) = @_; return dimlam($l) if !@$nu;
  my $key = join(",", @$l) . "|" . join(",", @$nu); return $CHI{$key}->copy if exists $CHI{$key};
  my @nu = @$nu; my $k = shift @nu;
  my $L = scalar @$l; my @beta = map { $l->[$_] + ($L - 1 - $_) } 0..$L-1;   # distinct, decreasing
  my %in = map { $_ => 1 } @beta; my $tot = BI(0);
  for my $i (0..$L-1) { my $nb = $beta[$i] - $k; next if $nb < 0 || $in{$nb};
    my $between = 0; for my $b (@beta) { $between++ if $b > $nb && $b < $beta[$i] }
    my @nbeta = @beta; $nbeta[$i] = $nb; @nbeta = sort { $b <=> $a } @nbeta;
    my @nl; for my $t (0..$L-1) { my $part = $nbeta[$t] - ($L - 1 - $t); push @nl, $part if $part > 0 }
    my $v = chi(\@nl, \@nu); if ($between % 2) { $tot -= $v } else { $tot += $v } }
  $CHI{$key} = $tot; return $tot->copy }

# T(nu; n, b) = #{(t_1..t_b; c_1..c_k): t_i transpositions in S_n, c_i disjoint cycles of lengths nu_i (ordered), t_1...t_b c_1...c_k = 1}
#            = sum_{lambda |- n} dim(lambda) chi^lambda(nu u 1^{n-|nu|}) f2(lambda)^b / (prod nu_i * (n-|nu|)!)
# returns hashref {b => BigInt} for the requested list of b
sub Tcount { my ($nu, $n, $blist) = @_;
  my $sz = 0; $sz += $_ for @$nu; return { map { $_ => BI(0) } @$blist } if $n < $sz;
  my %acc = map { $_ => BI(0) } @$blist;
  for my $l (@{partitions($n)}) { my $c = chi($l, $nu); next if $c->is_zero; my $w = dimlam($l) * $c; my $f2 = content_sum($l);
    for my $b (@$blist) { $acc{$b} += $w * (BI($f2) ** $b) } }
  my $den = fact($n - $sz); $den *= $_ for @$nu;
  for my $b (@$blist) { my ($q, $r) = $acc{$b}->copy->bdiv($den); die "Tcount not integral" unless $r->is_zero; $acc{$b} = $q }
  return \%acc }

1;
