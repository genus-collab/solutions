#!/usr/bin/perl
# Referee: fully independent computation of the Hurwitz component
#     Xi^(0)_h(R) = Sum_m 1/m! Int_{Mbar_{h,r+m}} Lambda^vee(1) prod_i 1/(1 - j_i psi_i) prod_k T(psi_{r+k}),  T = -psi^2/(1-psi)
# (= Int c(T(-log)) prod lam_{j_i}, Lemma 1 of proof_M_residue.md), for parts j_i >= 1, from HURWITZ NUMBERS ONLY:
#  (1) disconnected Hurwitz numbers by the Burnside/Frobenius character formula
#        H*(mu, b) = (1/d!)(1/z_mu) Sum_{lambda |- d} dim(lambda) chi_lambda(mu) f2(lambda)^b ,  f2 = content sum,
#      characters by Murnaghan-Nakayama (beta-numbers), dimensions by the hook formula;
#  (2) connected numbers by the exponential formula Z = exp(F) (graded recursion d F_d = d Z_d - Sum_k k F_k Z_{d-k}),
#      truncated to the monomials p_nu p_1^s with nu a sub-multiset of R_{>=2} (exact for the coefficients needed);
#  (3) ELSV (Ekedahl-Lando-Shapiro-Vainshtein 2001):  H(mu) = (b!/|Aut mu|) prod mu_i^mu_i/mu_i! Int Lambda^vee(1)/prod(1 - mu_i psi_i),
#      b = 2h - 2 + n + |mu| (normalisation checked below on genus 0 and genus 1 values);
#  (4) dilaton leaves: T = (1 + psi) - lam_1; the (1+psi) and 1 legs are removed by the string and dilaton equations
#      (Lambda^vee(1) is a pull-back), recursion G(t,u) below; the m-sum stops at m = 3h-3+r (T has psi-degree >= 2).
# No Mumford formula, no Witten-Kontsevich numbers, no tables, no topological recursion.
# Usage: perl ref_elsv.pl h j1 [j2 ...]      (prints Xi^(0)_h(R) exactly)
#        perl ref_elsv.pl selftest
use strict; use warnings; use Math::BigInt; use Math::BigRat;
my %fact; sub fct { my $n = shift; return $fact{$n} //= Math::BigInt->new($n)->bfac(); }
sub binom { my ($n, $k) = @_; return Math::BigInt->new(0) if ($k < 0 || $k > $n); return fct($n)->copy->bdiv(fct($k))->bdiv(fct($n - $k)); }
# ---------- partitions, dimensions, contents, characters ----------
my %parts_cache;
sub partitions { my $n = shift; return $parts_cache{$n} if exists $parts_cache{$n}; my @res; my @cur;
  my $rec; $rec = sub { my ($rem, $max) = @_; if ($rem == 0) { push @res, [@cur]; return; } for (my $p = ($rem < $max ? $rem : $max); $p >= 1; $p--) { push @cur, $p; $rec->($rem - $p, $p); pop @cur; } };
  $rec->($n, $n); return $parts_cache{$n} = \@res; }
my %dim_cache;
sub dimp { my @l = grep { $_ > 0 } @_; my $key = join(',', @l); return $dim_cache{$key} if exists $dim_cache{$key};
  my $n = 0; $n += $_ for @l; my $num = fct($n)->copy; my $den = Math::BigInt->new(1);
  my @conj; for my $i (0 .. $#l) { for my $c (0 .. $l[$i] - 1) { $conj[$c]++; } }
  for my $i (0 .. $#l) { for my $c (0 .. $l[$i] - 1) { my $hook = ($l[$i] - $c - 1) + ($conj[$c] - $i - 1) + 1; $den->bmul($hook); } }
  return $dim_cache{$key} = $num->bdiv($den); }
sub f2 { my @l = @_; my $s = 0; for my $i (0 .. $#l) { for my $c (0 .. $l[$i] - 1) { $s += $c - $i; } } return $s; }
# chi_lambda(nu u 1^rest) by Murnaghan-Nakayama: remove rim hooks of sizes nu (beta-numbers), then dim of the rest
sub chi_nu { my ($lam, $nu) = @_; my $L = scalar(@$lam) + 1; my @beta = map { ($lam->[$_] // 0) + ($L - 1 - $_) } 0 .. $L - 1; # distinct, decreasing
  return mn_rec(\@beta, $nu, 0); }
sub mn_rec { my ($beta, $nu, $idx) = @_;
  if ($idx >= @$nu) { my @b = sort { $b <=> $a } @$beta; my $L = scalar @b; my @lam = map { $b[$_] - ($L - 1 - $_) } 0 .. $L - 1; return dimp(@lam); }
  my $k = $nu->[$idx]; my %occ = map { $_ => 1 } @$beta; my $tot = Math::BigInt->new(0);
  for my $i (0 .. $#$beta) { my $bb = $beta->[$i]; my $nb = $bb - $k; next if $nb < 0 || $occ{$nb};
    my $between = 0; for my $x (@$beta) { $between++ if ($x > $nb && $x < $bb); }
    my @nbeta = @$beta; $nbeta[$i] = $nb; my $v = mn_rec(\@nbeta, $nu, $idx + 1); $v = $v->copy->bneg() if $between % 2; $tot->badd($v); }
  return $tot; }
# ---------- main computation ----------
sub xi0 { my ($h, @R) = @_;
  my @R2 = sort { $a <=> $b } grep { $_ >= 2 } @R; my $r1 = scalar(grep { $_ == 1 } @R); die "parts must be >= 1" if grep { $_ < 1 } @R;
  my $r = scalar @R; my $J = 0; $J += $_ for @R;
  my $M = 3 * $h - 3 + $r; my $Smax = $M; my $S1max = $r1 + $Smax;   # number of 1-parts in the ELSV multiset
  my $Bmax = 2 * $h - 2 + ($r + $Smax) + ($J + $Smax);
  # sub-multisets of R2, encoded as multiplicity vectors over distinct parts
  my @dist; my %mult; for (@R2) { push @dist, $_ unless $mult{$_}++; }
  my @nus = ([]); for my $p (@dist) { my @n2; for my $v (@nus) { for my $m (0 .. $mult{$p}) { push @n2, [@$v, $m]; } } @nus = @n2; }
  my $nkey = sub { join(',', @{$_[0]}) }; my %nuidx; $nuidx{$nkey->($nus[$_])} = $_ for 0 .. $#nus;
  my $nusize = sub { my $v = shift; my $s = 0; $s += $dist[$_] * $v->[$_] for 0 .. $#dist; $s };
  my $nulist = sub { my $v = shift; my @l; for my $i (0 .. $#dist) { push @l, ($dist[$i]) x $v->[$i]; } @l };
  # Z[nu][s][b] = H*(nu u 1^s, b)  (BigRat)
  my (@Z, @F);
  for my $iv (0 .. $#nus) { my $nu = $nus[$iv]; my @nl = $nulist->($nu); my $dn = $nusize->($nu);
    for my $s (0 .. $S1max) { my $d = $dn + $s; my @acc = map { Math::BigInt->new(0) } 0 .. $Bmax;
      if ($d == 0) { $acc[0] = Math::BigInt->new(1); }
      else { for my $lam (@{partitions($d)}) { my $w = dimp(@$lam)->copy->bmul(chi_nu($lam, \@nl)); next if $w->is_zero; my $f = Math::BigInt->new(f2(@$lam));
          my $pw = $w->copy; for my $b (0 .. $Bmax) { $acc[$b]->badd($pw); $pw = $pw->copy->bmul($f); } } }
      # z_mu = prod_k k^{m_k} m_k!
      my %m; $m{$_}++ for (@nl, (1) x $s); my $z = Math::BigInt->new(1); for my $k (keys %m) { $z->bmul(Math::BigInt->new($k)->bpow($m{$k}))->bmul(fct($m{$k})); }
      my $den = fct($d)->copy->bmul($z);
      $Z[$iv][$s] = [ map { $acc[$_]->is_zero ? Math::BigRat->new(0) : Math::BigRat->new($acc[$_]->bstr . "/" . $den->bstr) } 0 .. $Bmax ]; } }
  # connected: process (nu, s) by total size
  my @order; for my $iv (0 .. $#nus) { for my $s (0 .. $S1max) { push @order, [$iv, $s, $nusize->($nus[$iv]) + $s]; } } @order = sort { $a->[2] <=> $b->[2] } @order;
  # work with ordinary coefficients in beta (EGF coefficient / b!) so that products are plain convolutions
  my @bf = map { Math::BigRat->new(fct($_)->bstr) } 0 .. $Bmax;
  my (@z, @f); # z[iv][s] = list of [b, Z_b / b!] (non-zero only)
  for my $iv (0 .. $#nus) { for my $s (0 .. $S1max) { $z[$iv][$s] = [ map { [$_, $Z[$iv][$s][$_] / $bf[$_]] } grep { !$Z[$iv][$s][$_]->is_zero } 0 .. $Bmax ]; } }
  for my $o (@order) { my ($iv, $s, $d) = @$o; if ($d == 0) { $f[$iv][$s] = []; next; }
    my %acc = map { $_->[0] => $_->[1]->copy } @{$z[$iv][$s]}; my $nu = $nus[$iv];
    for my $o1 (@order) { my ($iv1, $s1, $d1) = @$o1; last if $d1 >= $d; next if $d1 == 0; next if $s1 > $s; my $nu1 = $nus[$iv1];
      my $ok = 1; my @nu2; for my $i (0 .. $#dist) { my $x = $nu->[$i] - $nu1->[$i]; if ($x < 0) { $ok = 0; last; } push @nu2, $x; } next unless $ok;
      my $iv2 = $nuidx{join(',', @nu2)}; my $s2 = $s - $s1; my $A = $f[$iv1][$s1]; my $B = $z[$iv2][$s2]; next unless (@$A && @$B);
      my %tmp; for my $pa (@$A) { my ($b1, $va) = @$pa; for my $pb (@$B) { my $bb = $b1 + $pb->[0]; next if $bb > $Bmax; my $t = $va * $pb->[1]; if (exists $tmp{$bb}) { $tmp{$bb} += $t; } else { $tmp{$bb} = $t; } } }
      my $fac = Math::BigRat->new("$d1/$d"); for my $bb (keys %tmp) { my $t = $tmp{$bb} * $fac; if (exists $acc{$bb}) { $acc{$bb} -= $t; } else { $acc{$bb} = -$t; } } }
    $f[$iv][$s] = [ map { [$_, $acc{$_}] } grep { !$acc{$_}->is_zero } sort { $a <=> $b } keys %acc ]; }
  for my $iv (0 .. $#nus) { for my $s (0 .. $S1max) { my @arr = map { Math::BigRat->new(0) } 0 .. $Bmax; for my $p (@{$f[$iv][$s] // []}) { $arr[$p->[0]] = $p->[1] * $bf[$p->[0]]; } $F[$iv][$s] = \@arr; } }
  # ELSV integrals E(s) = Int_{Mbar_{h, r+s}} Lambda^vee(1) prod lam_R prod lam_1^s
  my $ivR = $nuidx{join(',', map { $mult{$_} } @dist)}; my @E;
  for my $s (0 .. $Smax) { my $n = $r + $s; my $mu = $J + $s; my $b = 2 * $h - 2 + $n + $mu; die "Bmax too small" if $b > $Bmax;
    my $Hc = $F[$ivR][$r1 + $s][$b]; my %m; $m{$_}++ for (@R2, (1) x ($r1 + $s)); my $aut = Math::BigInt->new(1); $aut->bmul(fct($m{$_})) for keys %m;
    my $pp = Math::BigRat->new(1); for my $p (@R2) { $pp *= Math::BigRat->new(Math::BigInt->new($p)->bpow($p)->bstr) / Math::BigRat->new(fct($p)->bstr); }
    $E[$s] = $Hc * Math::BigRat->new($aut->bstr) / (Math::BigRat->new(fct($b)->bstr) * $pp); }
  # dilaton-leaf sum: G_s(t,u) = Int Lambda^vee(1) lam_R lam_1^s (1+psi)^t 1^u
  my $total = Math::BigRat->new(0);
  for my $s (0 .. $Smax) { my %G; my $mu = $J + $s; my $n0 = $r + $s;
    my $Gf; $Gf = sub { my ($t, $u) = @_; my $k = "$t,$u"; return $G{$k} if exists $G{$k}; my $v;
      if ($t == 0 && $u == 0) { $v = $E[$s]->copy; }
      elsif ($u >= 1) { $v = Math::BigRat->new($mu) * $Gf->($t, $u - 1); $v += Math::BigRat->new($t) * $Gf->($t - 1, $u) if $t >= 1; }
      else { $v = Math::BigRat->new($mu + 2 * $h - 2 + $n0 + $t - 1) * $Gf->($t - 1, 0); $v += Math::BigRat->new($t - 1) * $Gf->($t - 2, 1) if $t >= 2; }
      return $G{$k} = $v; };
    for my $m ($s .. $M) { my $c = Math::BigRat->new(binom($m, $s)->bstr) / Math::BigRat->new(fct($m)->bstr); $c = -$c if $s % 2; $total += $c * $Gf->($m - $s, 0); } }
  return ($total, \@E); }
unless (caller) {
if (@ARGV && $ARGV[0] eq q(selftest)) {
  # genus 0: Int_{Mbar_{0,n}} prod 1/(1 - mu_i psi_i) = |mu|^(n-3);  genus 1: Int_{Mbar_{1,1}} Lambda^vee(1)/(1 - j psi) = (j-1)/24
  for my $t ([0, 1, 1, 1], [0, 2, 1, 1], [0, 3, 2, 1, 1], [0, 2, 2, 2, 1], [1, 1], [1, 2], [1, 5], [1, 3, 2]) { my ($h, @R) = @$t; my (undef, $E) = xi0($h, @R);
    my $exp; if ($h == 0) { my $J = 0; $J += $_ for @R; $exp = Math::BigRat->new($J) ** (scalar(@R) - 3); } elsif (@R == 1) { $exp = Math::BigRat->new($R[0] - 1) / 24; }
    else { $exp = Math::BigRat->new(($R[0] + $R[1] - 1) * ($R[0] + $R[1]) - 0) / 24 - Math::BigRat->new($R[0]*$R[1]) / 24 ; }
    # genus 1, two legs: Int (1 - lambda_1) lam_a lam_b = <tau_2> (a^2+b^2)... computed directly below
    if ($h == 1 && @R == 2) { my ($a, $b) = @R; $exp = (Math::BigRat->new($a*$a + $b*$b + $a*$b) - Math::BigRat->new($a + $b)) / 24; }
    printf "selftest h=%d R=(%s): ELSV E(0) = %s, expected %s  %s\n", $h, join(',', @R), $E->[0]->bstr, $exp->bstr, ($E->[0] == $exp ? 'OK' : 'MISMATCH'); }
  exit; }
my ($h, @R) = @ARGV; my ($v) = xi0($h, @R); print "Xi^(0)_$h(", join(',', @R), ") = ", $v->bstr, "\n";
}
1;
