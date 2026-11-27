package HBconn;
# Connected Hurwitz counts with labelled non-unit cycles R and free unit parts, by the exponential formula,
# and the ELSV-normalised integrals A_h(R,1^c), and the finite formula (B2) for Xi^(0)_h(R).
#
# Conventions.  For a list nu of cycle lengths (all >= 2), n >= |nu|, b >= 0:
#   T(nu;n,b)      = #{(t_1..t_b; c_1..c_k) : t_i transpositions of S_n, c_i disjoint cycles of length nu_i (ordered tuple),
#                      t_1...t_b c_1...c_k = 1}                                    (HB::Tcount, Frobenius/Burnside formula)
#   Tc(nu;n,b)     = the same count restricted to tuples generating a transitive subgroup.
# EGF in (n,b):  sum_B y^B sum_{n,b} T(nu_B;n,b) x^n z^b/(n! b!) = exp( sum_B y^B sum Tc(nu_B;n,b) x^n z^b/(n! b!) ).
# ELSV connected Hurwitz number (Graber-Vakil convention, 1/|Aut| included):  H_{g,mu} = Tc(nu;n,b) / (n! |Aut nu|),
#   mu = nu u 1^(n-|nu|), b = 2g-2+l(mu)+|mu|.
# ELSV:  H_{g,mu} = (b!/|Aut mu|) prod mu_i^mu_i/mu_i!  Int_{Mbar_{g,l}} Lambda^vee(1)/prod(1-mu_i psi_i), hence
#   A_g(R,1^c) := Int Lambda^vee(1) prod_i 1/(1-j_i psi_i) prod_{k<=c} 1/(1-psi_{r+k}) = Tc(R;J+c,b) c!/((J+c)! b!) prod_i j_i!/j_i^j_i.
use strict; use warnings;
use Math::BigInt; use Math::BigRat;
use HB;

sub BR { Math::BigRat->new(shift) }
sub BI { Math::BigInt->new(shift) }

# arrays are hashes: $F->{n}{b} = BigInt
sub egf_mul { my ($F, $G, $N, $Bmax) = @_; my %P;
  for my $n1 (keys %$F) { for my $b1 (keys %{$F->{$n1}}) { my $f = $F->{$n1}{$b1}; next if $f->is_zero;
    for my $n2 (keys %$G) { my $n = $n1 + $n2; next if $n > $N; my $cn = HB::binom($n, $n1);
      for my $b2 (keys %{$G->{$n2}}) { my $b = $b1 + $b2; next if $b > $Bmax; my $g = $G->{$n2}{$b2}; next if $g->is_zero;
        my $t = $cn * HB::binom($b, $b1) * $f * $g;
        if (exists $P{$n}{$b}) { $P{$n}{$b} += $t } else { $P{$n}{$b} = $t } } } } }
  return \%P }
sub egf_sub { my ($F, $G) = @_; my %P; for my $n (keys %$F) { for my $b (keys %{$F->{$n}}) { $P{$n}{$b} = $F->{$n}{$b}->copy } }
  for my $n (keys %$G) { for my $b (keys %{$G->{$n}}) { if (exists $P{$n}{$b}) { $P{$n}{$b} -= $G->{$n}{$b} } else { $P{$n}{$b} = -$G->{$n}{$b} } } }
  return \%P }
sub egf_scale { my ($F, $c) = @_; my %P; for my $n (keys %$F) { for my $b (keys %{$F->{$n}}) { $P{$n}{$b} = $F->{$n}{$b} * $c } } return \%P }
# Q = F / Z with Z[0][0] = 1, Z supported on n >= 1 otherwise (Z = Z_empty); Q supported on n >= nmin
sub egf_div { my ($F, $Z, $N, $Bmax, $nmin) = @_; my %Q;
  for my $n ($nmin..$N) { for my $b (0..$Bmax) { my $v = exists $F->{$n}{$b} ? $F->{$n}{$b}->copy : BI(0);
      for my $n1 (1..$n - $nmin) { next unless exists $Z->{$n1}; my $cn = HB::binom($n, $n1);
        for my $b1 (keys %{$Z->{$n1}}) { next if $b1 > $b; my $q = $Q{$n - $n1}{$b - $b1}; next unless defined $q; next if $q->is_zero;
          $v -= $cn * HB::binom($b, $b1) * $Z->{$n1}{$b1} * $q } }
      $Q{$n}{$b} = $v unless $v->is_zero } }
  return \%Q }

# disconnected counts T(nu;n,b) on the grid n in [|nu|, N], b in [0,Bmax] (parity forced)
sub Zgrid { my ($nu, $N, $Bmax) = @_; my $s = 0; $s += $_ for @$nu; my $l = scalar @$nu; my %Z;
  for my $n ($s..$N) { my @bl = grep { ($_ - ($s - $l)) % 2 == 0 } 0..$Bmax; next unless @bl;
    my $T = HB::Tcount($nu, $n, \@bl); for my $b (@bl) { $Z{$n}{$b} = $T->{$b} unless $T->{$b}->is_zero } }
  return \%Z }

# main routine: connected labelled counts Tc(R; n, b) for all n <= N, b <= Bmax, via moment-cumulant over types
my %ZEMPTY;
sub conn_counts { my ($R, $N, $Bmax) = @_;
  my @R = sort { $a <=> $b } @$R; my %mult; $mult{$_}++ for @R; my @vals = sort { $a <=> $b } keys %mult; my $q = scalar @vals;
  my $zkey = "$N:$Bmax"; # Z_empty on n <= N
  my $Z0 = $ZEMPTY{$zkey} //= Zgrid([], $N, $Bmax);
  # enumerate types
  my @types = ([]); for my $i (0..$q-1) { my @nt; for my $t (@types) { for my $k (0..$mult{$vals[$i]}) { push @nt, [@$t, $k] } } @types = @nt }
  @types = sort { (eval join('+', 0, @$a)) <=> (eval join('+', 0, @$b)) } @types;
  my (%M, %C);
  for my $t (@types) { my $key = join(",", @$t); my @nu; my $s = 0;
    for my $i (0..$q-1) { push @nu, ($vals[$i]) x $t->[$i]; $s += $vals[$i] * $t->[$i] }
    if (!@nu) { $M{$key} = { 0 => { 0 => BI(1) } }; next }
    my $Zk = Zgrid(\@nu, $N, $Bmax);
    $M{$key} = egf_div($Zk, $Z0, $N, $Bmax, $s);
    # cumulant recursion with distinguished label of value index i0 (last nonzero)
    my $i0; for my $i (0..$q-1) { $i0 = $i if $t->[$i] > 0 }
    my $acc = $M{$key};
    for my $tp (@types) { my $kp = join(",", @$tp); next if $kp eq $key; next unless $tp->[$i0] >= 1;
      my $ok = 1; for my $i (0..$q-1) { $ok = 0 if $tp->[$i] > $t->[$i] } next unless $ok;
      my $coef = BI(1); for my $i (0..$q-1) { my $d = $i == $i0 ? 1 : 0; $coef *= HB::binom($t->[$i] - $d, $tp->[$i] - $d) }
      my @rest = map { $t->[$_] - $tp->[$_] } 0..$q-1; my $rk = join(",", @rest);
      die "missing C $kp" unless $C{$kp};
      my $prod = egf_mul($C{$kp}, $M{$rk}, $N, $Bmax);
      $acc = egf_sub($acc, egf_scale($prod, $coef)) }
    $C{$key} = $acc }
  my $full = join(",", map { $mult{$_} } @vals);
  return log_unit($Z0, $N, $Bmax) if !@R;
  return $C{$full} }
# connected unit-only counts: H = log Z (EGF in x and z); recursion on the component containing sheet 1:
#   Z[n][b] = sum_{m=1}^{n} sum_{b1} C(n-1,m-1) C(b,b1) H[m][b1] Z[n-m][b-b1]
sub log_unit { my ($Z, $N, $Bmax) = @_; my %H;
  for my $n (1..$N) { for my $b (0..$Bmax) { my $v = exists $Z->{$n}{$b} ? $Z->{$n}{$b}->copy : BI(0);
      for my $m (1..$n-1) { next unless exists $H{$m}; my $cn = HB::binom($n - 1, $m - 1);
        for my $b1 (keys %{$H{$m}}) { next if $b1 > $b; my $z = $Z->{$n - $m}{$b - $b1}; next unless defined $z;
          $v -= $cn * HB::binom($b, $b1) * $H{$m}{$b1} * $z } }
      $H{$n}{$b} = $v unless $v->is_zero } }
  return \%H }

# A_h(R,1^c) for c = 0..cmax  (BigRat), from connected counts
sub A_list { my ($R, $h, $cmax) = @_; my $r = scalar @$R; my $J = 0; $J += $_ for @$R;
  my $N = $J + $cmax; my $Bmax = 2*$h - 2 + $r + $J + 2*$cmax;
  my $C = conn_counts($R, $N, $Bmax);
  my $kap = BR(1); for my $j (@$R) { $kap *= BR(HB::fact($j)) / BR(BI($j) ** $j) }
  my @A;
  for my $c (0..$cmax) { my $n = $J + $c; my $b = 2*$h - 2 + $r + $J + 2*$c;
    my $tc = exists $C->{$n}{$b} ? $C->{$n}{$b} : BI(0);
    push @A, BR($tc) * BR(HB::fact($c)) / (BR(HB::fact($n)) * BR(HB::fact($b))) * $kap }
  return \@A }

# finite formula (B2):  Xi^(0)_h(R) = sum_{c=0}^{n0} (-1)^c (A_h(R,1^c)/c!) sum_{i=0}^{n0-c} C(n0+A0, n0-c-i) (J+c)^i/i!,
#   A0 = 2h-2+r, n0 = 3h-3+r.
sub xi0_from_A { my ($A, $h, $r, $J) = @_; my $A0 = 2*$h - 2 + $r; my $n0 = 3*$h - 3 + $r; my $tot = BR(0);
  for my $c (0..$n0) { my $inner = BR(0);
    for my $i (0..$n0 - $c) { $inner += BR(HB::binom($n0 + $A0, $n0 - $c - $i)) * BR(BI($J + $c) ** $i) / BR(HB::fact($i)) }
    my $term = $A->[$c] / BR(HB::fact($c)) * $inner; if ($c % 2) { $tot -= $term } else { $tot += $term } }
  return $tot }
sub xi0 { my ($R, $h) = @_; my $r = scalar @$R; my $J = 0; $J += $_ for @$R; my $n0 = 3*$h - 3 + $r;
  my $A = A_list($R, $h, $n0); return xi0_from_A($A, $h, $r, $J) }
# the polynomial Q(s) = sum_m q_m s^m (coefficients q_0..q_n0) of the k=0 component:  Q(s) = <Lambda^vee(1) K_s lam_R>
#   q_m = sum_{c<=m} (A_c/c!) [s^(m-c)] (1+s)^(-(A0+c)) exp(-s (J+c)/(1+s))
sub Q_from_A { my ($A, $h, $r, $J) = @_; my $A0 = 2*$h - 2 + $r; my $n0 = 3*$h - 3 + $r; my @q = map { BR(0) } 0..$n0;
  for my $c (0..$n0) { my $a = $A->[$c] / BR(HB::fact($c)); next if $a == 0;
    # series f(s) = (1+s)^(-(A0+c)) exp(-x s/(1+s)), x = J+c, up to s^(n0-c):  exp(-x s/(1+s)) = sum_i (-x)^i s^i (1+s)^(-i)/i!
    my $x = $J + $c; my @f = map { BR(0) } 0..$n0 - $c;
    for my $i (0..$n0 - $c) { my $ci = BR(BI(-$x) ** $i) / BR(HB::fact($i));
      # (1+s)^(-(A0+c+i)) = sum_k C(-(A0+c+i), k) s^k = sum_k (-1)^k C(A0+c+i+k-1, k) s^k
      for my $k (0..$n0 - $c - $i) { my $bk = BR(HB::binom($A0 + $c + $i + $k - 1, $k)); $bk = -$bk if $k % 2; $f[$i + $k] += $ci * $bk } }
    for my $m ($c..$n0) { $q[$m] += $a * $f[$m - $c] } }
  return \@q }

1;
