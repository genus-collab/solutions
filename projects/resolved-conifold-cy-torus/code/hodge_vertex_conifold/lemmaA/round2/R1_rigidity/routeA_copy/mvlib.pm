package mvlib;
# mvlib.pm -- exact (Math::BigRat) tools for route A (Marino-Vafa formula), lemmaA/A_mv.
# No topological recursion anywhere: Witten-Kontsevich numbers by string + dilaton from the classical base table,
# characters by Murnaghan-Nakayama, Schur side of the Marino-Vafa formula as Laurent series in x = sqrt(-1)*lambda.
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(R fact binom bern wk partitions chi hooks kappa zee
                 ser_mul ser_add ser_scal ser_exp ser_zero
                 parse_file parse_expr eval_poly poly_eval_at
                 mv_conn mv_unit lag);

sub R { Math::BigRat->new(@_) }
my @factc = (R(1));
sub fact { my $n = shift; for my $k (scalar(@factc)..$n) { $factc[$k] = $factc[$k-1] * $k; } return $factc[$n]->copy; }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; return fact($n) / (fact($k) * fact($n - $k)); }

# Bernoulli numbers B_n (B_1 = -1/2)
my @bernc = (R(1));
sub bern { my $n = shift;
  for my $m (scalar(@bernc)..$n) { my $s = R(0); for my $k (0..$m-1) { $s += binom($m+1, $k) * $bernc[$k]; } $bernc[$m] = -$s / ($m + 1); }
  return $bernc[$n]->copy; }

# ---------------- Witten-Kontsevich numbers <tau_a1 ... tau_an>_g ----------------
# base values (all a_i >= 2, plus the genus 0 and 1 seeds): Witten 1991 (genus <= 2), Faber's tables (genus 3).
my %wkbase = (
 '0:0,0,0' => '1', '1:1' => '1/24',
 '2:4' => '1/1152', '2:2,3' => '29/5760', '2:2,2,2' => '7/240',
 '3:7' => '1/82944', '3:2,6' => '77/414720', '3:3,5' => '503/1451520', '3:4,4' => '607/1451520',
 '3:2,2,5' => '17/5760', '3:2,3,4' => '1121/241920', '3:3,3,3' => '583/96768',
 '3:2,2,2,4' => '7/240', '3:2,2,3,3' => '205/3456', '3:2,2,2,2,3' => '193/288', '3:2,2,2,2,2,2' => '1225/144',
);
my %wkmemo;
sub wk { my ($g, @a) = @_; @a = sort { $a <=> $b } @a; my $n = scalar @a;
  my $key = "$g:" . join(',', @a);
  return $wkmemo{$key} if exists $wkmemo{$key};
  my $sum = 0; $sum += $_ for @a; my $res;
  if ($g < 0 || 2*$g - 2 + $n <= 0 || $sum != 3*$g - 3 + $n || (grep { $_ < 0 } @a)) { $res = R(0); }
  elsif (exists $wkbase{$key}) { $res = R($wkbase{$key}); }
  elsif ($a[0] == 0) { my @rest = @a[1..$#a]; $res = R(0);
     for my $i (0..$#rest) { next if $rest[$i] == 0; my @b = @rest; $b[$i]--; $res += wk($g, @b); } }
  elsif ($a[0] == 1) { my @rest = @a[1..$#a]; $res = (2*$g - 2 + $n - 1) * wk($g, @rest); }
  else { die "wk: base value missing for $key\n"; }
  $wkmemo{$key} = $res; return $res; }

# ---------------- partitions, characters, hooks ----------------
my %partc;
sub partitions { my ($n, $max) = @_; $max = $n unless defined $max; $max = $n if $max > $n;
  my $key = "$n:$max"; return @{$partc{$key}} if exists $partc{$key};
  my @out; if ($n == 0) { @out = ([]); }
  else { for my $k (reverse 1..$max) { for my $p (partitions($n - $k, $k)) { push @out, [$k, @$p]; } } }
  $partc{$key} = \@out; return @out; }
# Murnaghan-Nakayama: chi_lambda(mu), lambda and mu arrays (mu any order)
my %chic;
sub chi { my ($lam, $mu) = @_;
  my @mu = sort { $b <=> $a } @$mu; my $key = join(',', @$lam) . '|' . join(',', @mu);
  return $chic{$key} if exists $chic{$key};
  my $n = 0; $n += $_ for @$lam; my $m = 0; $m += $_ for @mu; die "chi: size mismatch" unless $n == $m;
  my $res;
  if ($n == 0) { $res = 1; }
  else { my $k = shift @mu; my $l = scalar @$lam;
    my @beta = map { $lam->[$_] + $l - 1 - $_ } 0..$l-1;   # strictly decreasing beta numbers
    my %occ = map { $_ => 1 } @beta; $res = 0;
    for my $i (0..$l-1) { my $bt = $beta[$i]; my $t = $bt - $k; next if $t < 0 || $occ{$t};
      my $ht = 0; for my $bb (@beta) { $ht++ if $bb > $t && $bb < $bt; }
      my @nb = sort { $b <=> $a } (map { $_ == $bt ? $t : $_ } @beta);
      my $L = scalar @nb; my @nl = map { $nb[$_] - ($L - 1 - $_) } 0..$L-1; @nl = grep { $_ > 0 } @nl;
      $res += (($ht % 2) ? -1 : 1) * chi(\@nl, \@mu); } }
  $chic{$key} = $res; return $res; }
sub hooks { my ($lam) = @_; my @l = @$lam; my @conj; for my $i (0..$#l) { for my $j (0..$l[$i]-1) { $conj[$j]++; } }
  my @h; for my $i (0..$#l) { for my $j (0..$l[$i]-1) { push @h, ($l[$i] - $j - 1) + ($conj[$j] - $i - 1) + 1; } } return @h; }
sub kappa { my ($lam) = @_; my $k = 0; for my $i (0..$#$lam) { $k += $lam->[$i] * ($lam->[$i] - 2*($i+1) + 1); } return $k; }
sub zee { my @mu = @_; my %m; $m{$_}++ for @mu; my $z = R(1); for my $k (keys %m) { $z *= R($k) ** $m{$k}; $z *= fact($m{$k}); } return $z; }

# ---------------- truncated power series (arrays, fixed length P+1) ----------------
sub ser_zero { my $P = shift; return [ map { R(0) } 0..$P ]; }
sub ser_mul { my ($a, $b, $P) = @_; my @c = map { R(0) } 0..$P;
  for my $i (0..$P) { next if $a->[$i]->is_zero; for my $j (0..$P-$i) { next if $b->[$j]->is_zero; $c[$i+$j] += $a->[$i] * $b->[$j]; } } return \@c; }
sub ser_add { my ($a, $b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$#$a ]; }
sub ser_scal { my ($s, $a) = @_; return [ map { $a->[$_] * $s } 0..$#$a ]; }
sub ser_exp { my ($e, $P) = @_; # e[0] must be 0
  die "ser_exp: constant term" unless $e->[0]->is_zero; my @f = (R(1));
  for my $n (1..$P) { my $s = R(0); for my $m (1..$n) { next if $e->[$m]->is_zero; $s += $m * $e->[$m] * $f[$n-$m]; } $f[$n] = $s / $n; } return \@f; }

# ---------------- Marino-Vafa Schur side ----------------
# W_nu(x) = exp((tau+1/2) kappa_nu x/2) prod_boxes 1/(2 sinh(h x/2)) = x^{-|nu|} wnorm(nu), wnorm a power series (length P+1)
# Integer form: W_nu(x) = x^{-|nu|} (1/H_nu) sum_n F_n x^n / (d^n n!), F_n integers, d = dcommon(tau, P):
#   e_1 = (tau+1/2) kappa/2, e_{2k} = c_k P_{2k}(nu), c_k = -B_{2k}/(2k (2k)!), E_m = d e_m integers,
#   F_0 = 1, F_n = sum_{m=1}^n m (n-1)!/(n-m)! E_m d^{m-1} F_{n-m}.
my @logSc;   # c_k
my %dcom;
sub dcommon { my ($tau, $P) = @_; my $key = "$tau|$P"; return $dcom{$key} if exists $dcom{$key};
  my $d = Math::BigInt->new(4) * Math::BigInt->new($tau->denominator);
  for (my $k = 1; 2*$k <= $P; $k++) { $logSc[$k] //= -bern(2*$k) / (2*$k * fact(2*$k)); my $dk = Math::BigInt->new($logSc[$k]->denominator);
    $d = $d * $dk / Math::BigInt::bgcd($d, $dk); }
  $dcom{$key} = $d; return $d; }
my %Wc;
sub conjp { my ($nu) = @_; my @c; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { $c[$j]++; } } return \@c; }
sub WnuInt { my ($nu, $tau, $P) = @_; my $key = join(",", @$nu) . "|$tau|$P"; return $Wc{$key} if exists $Wc{$key};
  # conjugate: kappa -> -kappa, same hooks, so F_n(nu') = (-1)^n F_n(nu)
  my $ckey = join(",", @{conjp($nu)}) . "|$tau|$P";
  if (exists $Wc{$ckey}) { my $Fc = $Wc{$ckey}; my @F = map { ($_ % 2) ? -$Fc->[$_] : $Fc->[$_]->copy } 0..$#$Fc; $Wc{$key} = \@F; return \@F; }
  my $d = dcommon($tau, $P); my @h = hooks($nu); my $kap = kappa($nu);
  my @E = map { Math::BigInt->new(0) } 0..$P;
  if ($P >= 1) { my $e1 = ($tau + R("1/2")) * $kap / 2 * Math::BigRat->new($d); die "E1" unless $e1->is_int; $E[1] = Math::BigInt->new("$e1"); }
  my %hm; $hm{$_}++ for @h;
  for (my $k = 1; 2*$k <= $P; $k++) { my $p2k = Math::BigInt->new(0); for my $hh (keys %hm) { $p2k += Math::BigInt->new($hh) ** (2*$k) * $hm{$hh}; }
    my $ek = $logSc[$k] * Math::BigRat->new($d) * Math::BigRat->new($p2k); die "E2k" unless $ek->is_int; $E[2*$k] = Math::BigInt->new("$ek"); }
  my @dp = (Math::BigInt->new(1)); for my $m (1..$P) { $dp[$m] = $dp[$m-1] * $d; }
  my @F = (Math::BigInt->new(1));
  for my $n (1..$P) { my $s = Math::BigInt->new(0); my $ff = Math::BigInt->new(1);   # ff = (n-1)!/(n-m)!
    for my $m (1..$n) { $ff *= ($n - $m + 1) if $m > 1; next if $E[$m]->is_zero; $s += $ff * $m * $E[$m] * $dp[$m-1] * $F[$n-$m]; }
    $F[$n] = $s; }
  $Wc{$key} = \@F; return \@F; }
sub flam { my ($nu) = @_; my @h = hooks($nu); my $n = 0; $n += $_ for @$nu; my $H = Math::BigInt->new(1); $H *= $_ for @h; return Math::BigInt->new($n)->bfac() / $H; }
my %flc; sub flamc { my $k = join(",", @{$_[0]}); $flc{$k} //= flam($_[0]); return $flc{$k}; }
# [p_mu] R (mu = B + 1^na) = x^{-|mu|} * series(P):  coefficient_n = (1/(z_mu |mu|! d^n n!)) sum_nu chi_nu(mu) f^nu F_n(nu)
my %Rc;
sub Rcoef { my ($B, $na, $tau, $P) = @_; my @mu = sort { $b <=> $a } (@$B, (1) x $na); my $key = join(",", @mu) . "|$tau|$P";
  return $Rc{$key} if exists $Rc{$key};
  my $n = 0; $n += $_ for @mu; my $d = dcommon($tau, $P);
  my @S = map { Math::BigInt->new(0) } 0..$P;
  for my $nu (partitions($n)) { my $ch = chi($nu, \@mu); next unless $ch; my $F = WnuInt($nu, $tau, $P); my $w = flamc($nu) * $ch;
    for my $i (0..$P) { $S[$i] += $w * $F->[$i] unless $F->[$i]->is_zero; } }
  my $den0 = zee(@mu) * fact($n); my @s; my $dpow = Math::BigInt->new(1);
  for my $i (0..$P) { $s[$i] = Math::BigRat->new($S[$i]) / ($den0 * Math::BigRat->new($dpow) * fact($i)); $dpow *= $d; }
  $Rc{$key} = \@s; return \@s; }

sub autf { my @mu = @_; my %m; $m{$_}++ for @mu; my $a = R(1); $a *= fact($m{$_}) for keys %m; return $a; }
# set partitions of (0..n-1)
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @out;
  for my $sp (setparts($n-1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n-1; push @out, \@c; }
    my @c = map { [@$_] } @$sp; push @c, [$n-1]; push @out, \@c; } return @out; }

# Connected coefficient [p_mu] log R (no unit-leg variable), returns series s with [p_mu] log R = x^{-|mu|} s(x).
sub mv_conn { my ($mu, $tau, $P) = @_; my @mu = @$mu; my $l = scalar @mu; my $tot = ser_zero($P);
  for my $sp (setparts($l)) { my $k = scalar @$sp; my $prod;
    for my $blk (@$sp) { my @B = map { $mu[$_] } @$blk; my $f = ser_scal(autf(@B), Rcoef(\@B, 0, $tau, $P));
      $prod = defined $prod ? ser_mul($prod, $f, $P) : $f; }
    my $c = R((-1) ** ($k - 1)) * fact($k - 1); $tot = ser_add($tot, ser_scal($c, $prod)); }
  return ser_scal(R(1) / autf(@mu), $tot); }

# Unit legs: for R (parts >= 2), returns M[g][a] = <Omega lam_R lam_1^a>_g (c = 1, Omega = Lambda(1)Lambda(tau)Lambda(-1-tau)/(-sigma)^g)
# for g = 0..gmax, a = 0..A.  Uses p_1 = zeta*x and the set-partition (cumulant) formula over the labelled parts of R.
sub mv_unit { my ($Rp, $tau, $gmax, $A) = @_; my @Rp = @$Rp; my $r = scalar @Rp; my $J = 0; $J += $_ for @Rp;
  my $P = 2*$gmax - 2 + $r + $J + 2*$A; $P = 0 if $P < 0;
  # F_S[a] : series with valuation -|B_S| (x^a * [p_B p_1^a] R * |Aut B|)
  my %F;
  my @subsets = (0..(2**$r - 1));
  for my $S (@subsets) { my @B = map { $Rp[$_] } grep { ($S >> $_) & 1 } 0..$r-1; my $au = autf(@B);
    $F{$S} = [ map { ser_scal($au, Rcoef(\@B, $_, $tau, $P)) } 0..$A ]; }
  # 2D helpers: arrays over a of series
  my $mul2 = sub { my ($X, $Y) = @_; my @Z; for my $c (0..$A) { my $s = ser_zero($P); for my $a (0..$c) { $s = ser_add($s, ser_mul($X->[$a], $Y->[$c-$a], $P)); } $Z[$c] = $s; } return \@Z; };
  # inverse of F_empty (valuation 0, constant 1)
  my $F0 = $F{0}; my @G; $G[0] = ser_zero($P); $G[0][0] = R(1);
  die "F0[0] != 1" unless $F0->[0][0] == 1 && (grep { !$_->is_zero } @{$F0->[0]}[1..$P]) == 0;
  for my $c (1..$A) { my $s = ser_zero($P); for my $b (1..$c) { $s = ser_add($s, ser_mul($F0->[$b], $G[$c-$b], $P)); } $G[$c] = ser_scal(R(-1), $s); }
  my $Ginv = \@G;
  my $full = 2**$r - 1; my @tot = map { ser_zero($P) } 0..$A;
  for my $sp (setparts($r)) { my $k = scalar @$sp; my $prod;
    for my $blk (@$sp) { my $S = 0; $S |= (1 << $_) for @$blk; my $f = $mul2->($F{$S}, $Ginv); $prod = defined $prod ? $mul2->($prod, $f) : $f; }
    my $c = R((-1) ** ($k - 1)) * fact($k - 1);
    for my $a (0..$A) { $tot[$a] = ser_add($tot[$a], ser_scal($c, $prod->[$a])); } }
  # tot[a] = x^{-J} * series = |Aut R| x^a [p_R p_1^a] log R;   [x^{2g-2+r+2a}] tot[a] = sigma^{g+r+a-1} Pi_R M_g(R+1^a)/a!
  my $sig = $tau * (1 + $tau); my $Pi = R(1);
  for my $j (@Rp) { for my $b (1..$j-1) { $Pi *= ($j * $tau + $b); } $Pi /= fact($j - 1); }
  my @M;
  for my $g (0..$gmax) { for my $a (0..$A) { my $idx = 2*$g - 2 + $r + 2*$a + $J;
      if ($idx < 0) { $M[$g][$a] = undef; next; }
      die "precision" if $idx > $P;
      $M[$g][$a] = $tot[$a][$idx] * fact($a) / ($sig ** ($g + $r + $a - 1) * $Pi); } }
  return \@M; }

# Laguerre polynomial L_n^{(al)}(x)
sub lag { my ($n, $al, $x) = @_; my $s = R(0); for my $i (0..$n) { $s += R((-1) ** $i) * binom($n + $al, $n - $i) * ($x ** $i) / fact($i); } return $s; }

# ---------------- parser for Mathematica InputForm (subset) ----------------
# AST: ['n', BigRat] | ['v', name] | ['wk', g, [a..]] | ['+', [[sign, node]..]] | ['*', [[exp(+1/-1), node]..]] | ['^', node, int]
sub tokenize { my $s = shift; my @t; while ($s =~ /\G\s*(\d+|[A-Za-z][A-Za-z0-9]*|[-+*\/^()\[\]{},])/gc) { push @t, $1; }
  $s =~ /\G\s*/gc; die "tokenize: junk at " . pos($s) . ": " . substr($s, pos($s), 30) if pos($s) < length($s); return \@t; }
sub parse_expr { my $s = shift; my $t = tokenize($s); my $i = 0; my $node = p_sum($t, \$i); die "parse: trailing tokens at $i" if $i < @$t; return $node; }
sub p_sum { my ($t, $i) = @_; my @terms; my $sign = 1;
  if ($$i < @$t && ($t->[$$i] eq '-' || $t->[$$i] eq '+')) { $sign = $t->[$$i] eq '-' ? -1 : 1; $$i++; }
  push @terms, [$sign, p_prod($t, $i)];
  while ($$i < @$t && ($t->[$$i] eq '+' || $t->[$$i] eq '-')) { my $sg = $t->[$$i] eq '-' ? -1 : 1; $$i++; push @terms, [$sg, p_prod($t, $i)]; }
  return @terms == 1 && $terms[0][0] == 1 ? $terms[0][1] : ['+', \@terms]; }
sub p_prod { my ($t, $i) = @_; my @f = ([1, p_pow($t, $i)]);
  while ($$i < @$t && ($t->[$$i] eq '*' || $t->[$$i] eq '/')) { my $e = $t->[$$i] eq '/' ? -1 : 1; $$i++; push @f, [$e, p_pow($t, $i)]; }
  return @f == 1 && $f[0][0] == 1 ? $f[0][1] : ['*', \@f]; }
sub p_pow { my ($t, $i) = @_; my $b = p_atom($t, $i);
  if ($$i < @$t && $t->[$$i] eq '^') { $$i++; my $sg = 1; if ($t->[$$i] eq '(') { $$i++; if ($t->[$$i] eq '-') { $sg = -1; $$i++; } my $e = $t->[$$i]; $$i++; die "pow )" unless $t->[$$i] eq ')'; $$i++; return ['^', $b, $sg * $e]; }
    if ($t->[$$i] eq '-') { $sg = -1; $$i++; } my $e = $t->[$$i]; die "pow exponent $e" unless $e =~ /^\d+$/; $$i++; return ['^', $b, $sg * $e]; }
  return $b; }
sub p_atom { my ($t, $i) = @_; my $tok = $t->[$$i];
  if ($tok eq '(') { $$i++; my $n = p_sum($t, $i); die "expected )" unless $t->[$$i] eq ')'; $$i++; return $n; }
  if ($tok eq '-') { $$i++; my $n = p_pow($t, $i); return ['*', [[1, ['n', R(-1)]], [1, $n]]]; }
  if ($tok =~ /^\d+$/) { $$i++; return ['n', R($tok)]; }
  if ($tok eq 'wk') { $$i++; die "wk[" unless $t->[$$i] eq '['; $$i++; my $g = $t->[$$i]; $$i++; die "wk," unless $t->[$$i] eq ','; $$i++;
    die "wk{" unless $t->[$$i] eq '{'; $$i++; my @a; while ($t->[$$i] ne '}') { push @a, $t->[$$i] if $t->[$$i] =~ /^\d+$/; $$i++; } $$i++;
    die "wk]" unless $t->[$$i] eq ']'; $$i++; return ['n', wk($g, @a)]; }
  if ($tok eq 'jj') { $$i++; die "jj[" unless $t->[$$i] eq '['; $$i++; my $k = $t->[$$i]; $$i++; die "jj]" unless $t->[$$i] eq ']'; $$i++; return ['v', "jj$k"]; }
  if ($tok =~ /^[A-Za-z]/) { $$i++; return ['v', $tok]; }
  die "parse: unexpected token $tok at $$i"; }
# polynomials in variables @V: hash key "e1,e2,..." -> BigRat
sub pconst { my ($c, $nv) = @_; return { join(',', (0) x $nv) => $c }; }
sub padd { my ($a, $b, $sg) = @_; my %c = map { $_ => $a->{$_}->copy } keys %$a;
  for my $k (keys %$b) { my $v = $sg < 0 ? -$b->{$k} : $b->{$k}; $c{$k} = exists $c{$k} ? $c{$k} + $v : $v->copy; delete $c{$k} if $c{$k}->is_zero; } return \%c; }
sub pmul { my ($a, $b) = @_; my %c;
  for my $k1 (keys %$a) { my @e1 = split /,/, $k1; for my $k2 (keys %$b) { my @e2 = split /,/, $k2;
      my $k = join(',', map { $e1[$_] + $e2[$_] } 0..$#e1); my $v = $a->{$k1} * $b->{$k2}; $c{$k} = exists $c{$k} ? $c{$k} + $v : $v; } }
  for my $k (keys %c) { delete $c{$k} if $c{$k}->is_zero; } return \%c; }
sub pisconst { my ($a) = @_; for my $k (keys %$a) { return 0 if $k =~ /[1-9]/; } return 1; }
sub pconstval { my ($a) = @_; for my $k (keys %$a) { return $a->{$k}; } return R(0); }
# evaluate AST to a polynomial in @$V with numeric values %$num for other variables
sub eval_poly { my ($node, $V, $num) = @_; my $nv = scalar @$V; my %vi = map { $V->[$_] => $_ } 0..$nv-1;
  my $rec; $rec = sub { my $n = shift; my $ty = $n->[0];
    if ($ty eq 'n') { return pconst($n->[1]->copy, $nv); }
    if ($ty eq 'v') { my $name = $n->[1];
      if (exists $vi{$name}) { my @e = (0) x $nv; $e[$vi{$name}] = 1; return { join(',', @e) => R(1) }; }
      die "eval: no value for $name" unless exists $num->{$name}; return pconst($num->{$name}->copy, $nv); }
    if ($ty eq '+') { my $acc = {}; for my $t (@{$n->[1]}) { $acc = padd($acc, $rec->($t->[1]), $t->[0]); } return $acc; }
    if ($ty eq '*') { my $acc = pconst(R(1), $nv); my $den = R(1);
      for my $f (@{$n->[1]}) { my $v = $rec->($f->[1]); if ($f->[0] > 0) { $acc = pmul($acc, $v); } else { die "division by non-constant" unless pisconst($v); $den *= pconstval($v); } }
      if ($den != 1) { for my $k (keys %$acc) { $acc->{$k} /= $den; } } return $acc; }
    if ($ty eq '^') { my $b = $rec->($n->[1]); my $e = $n->[2];
      if ($e < 0) { die "negative power of non-constant" unless pisconst($b); return pconst(pconstval($b) ** $e, $nv); }
      my $acc = pconst(R(1), $nv); $acc = pmul($acc, $b) for 1..$e; return $acc; }
    die "eval: bad node $ty"; };
  return $rec->($node); }
sub poly_eval_at { my ($p, $vals) = @_; my $s = R(0);
  for my $k (keys %$p) { my @e = split /,/, $k; my $t = $p->{$k}->copy; for my $i (0..$#e) { $t *= $vals->[$i] ** $e[$i] if $e[$i]; } $s += $t; } return $s; }
sub parse_file { my ($file) = @_; open my $fh, '<', $file or die "open $file: $!"; local $/; my $txt = <$fh>; close $fh;
  my @recs; for my $line (split /\n/, $txt) { next unless $line =~ /\S/; my @f = split /\s*\|\s*/, $line; push @recs, \@f; } return \@recs; }
1;
