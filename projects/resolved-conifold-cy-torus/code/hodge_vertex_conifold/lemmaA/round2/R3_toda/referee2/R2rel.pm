package R2rel;
# Second referee: independent leading-order (dilaton point) evaluation of ANY Hirota equation P(D) tau_0 . tau_1 = 0,
# tau_0 = Z(p), tau_1 = Z(e^{tau x} o p), derived from scratch (does not use mkp_lead.pl / mkp_general_lead.pl / mkp_eqs.pm):
#   P(D) tau_0.tau_1/(tau_0 tau_1) = P(d_s) exp(Phi)|_0,  Phi = F(t+s) + F1(t-s) - F - F1,  d^B Phi = d^B S (|B| even), -d^B A (|B| odd),
#   A = F1 - F = sum_{k>=1} (tau x)^k E^k F/k!,  S = F + F1,  F = sum_g hbar^(g-1) F_g(q),  hbar = sigma x^2,  q_i = x sigma Pi_i p_i,
#   d_{t_i} = i d_{p_i} = i x sigma Pi_i d_{q_i}.  Block (indices B, profile T_B, genus g, Euler power k):
#   x-weight 2g-2+|B|+k;  factor (-1 if |B| odd)(2 if |B| even, k = 0) prod_{i in B}(i Pi_i) sigma^(|B|+g-1) tau^k/k! * L,
#   L = (J_U)^(v) (delta_g(U))_k X_g(U),  U = T_B + {i in B: i >= 2},  v = #{i in B: i = 1},  delta_g(U) = D(U) - 2g + 2,
#   L = L0(g,k,v) for U empty (computed here from G^(0)_0 = u - 3u^2/4 + u^3/6 and y d_y G^(1)_0 = (w-1)(w-1-1/sigma)/24;
#   L0 = 0 for g >= 2 since G^(g)_0 lies in w^(2g-2) Q[w]).
# The relation at (x^n, q_R) is the sum over monomials, set partitions of the monomial's index positions, labelled assignments of
# the parts of R to blocks, and (g_B, k_B) with sum (2g_B - 2 + |B| + k_B) = n.
use strict; use warnings; use Math::BigRat;

# ---------- rings ----------
sub ring_Q { +{ name => 'Q', z => Math::BigRat->new(0), o => Math::BigRat->new(1),
    add => sub { $_[0] + $_[1] }, mul => sub { $_[0] * $_[1] }, neg => sub { -$_[0] }, inv => sub { Math::BigRat->new(1) / $_[0] },
    int => sub { Math::BigRat->new($_[0]) }, rat => sub { Math::BigRat->new("$_[0]") }, iszero => sub { $_[0] == 0 }, str => sub { "$_[0]" } } }
sub ring_P { my $p = shift; my $md = sub { my $x = $_[0] % $p; $x += $p if $x < 0; $x };
  my $inv = sub { my $a = $md->($_[0]); die "ring_P inv 0" unless $a; my ($t, $nt, $r, $nr) = (0, 1, $p, $a);
    while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q * $nt); ($r, $nr) = ($nr, $r - $q * $nr); } $md->($t) };
  +{ name => "Fp($p)", p => $p, z => 0, o => 1, add => sub { ($_[0] + $_[1]) % $p }, mul => sub { $_[0] * $_[1] % $p }, neg => sub { $md->(-$_[0]) },
     inv => $inv, int => sub { $md->($_[0]) },
     rat => sub { my $q = Math::BigRat->new("$_[0]"); my $n = $q->numerator->copy->bmod($p)->numify; my $d = $q->denominator->copy->bmod($p)->numify; $md->($n) * $inv->($d) % $p },
     iszero => sub { $_[0] == 0 }, str => sub { "$_[0]" } } }

# ---------- Hirota members ----------
sub partitions { my ($n, $mx) = @_; $mx //= $n; return ([]) if $n == 0; my @o; for my $k (reverse 1..($mx < $n ? $mx : $n)) { push @o, [$k, @$_] for partitions($n - $k, $k); } @o }
sub bfact { my $n = shift; my $f = Math::BigRat->new(1); $f *= $_ for 2..$n; $f }
# P_k = D_{t_k} + p_k(-Dtilde),  p_k: exp(sum t_i z^i) = sum p_k z^k,  Dtilde = (D_1, D_2/2, D_3/3, ...)
sub hirota_Pk { my $k = shift; my %P; $P{"$k"} = Math::BigRat->new(1);
  for my $lam (partitions($k)) { my %m; $m{$_}++ for @$lam; my $c = Math::BigRat->new(1);
    for my $i (keys %m) { $c *= (Math::BigRat->new(-1) / $i) ** $m{$i} / bfact($m{$i}); }
    my $key = join(',', sort { $a <=> $b } @$lam); $P{$key} = ($P{$key} // Math::BigRat->new(0)) + $c; }
  [ map { [ $P{$_}, [ split /,/, $_ ] ] } grep { $P{$_} != 0 } sort keys %P ] }
sub show_P { my $P = shift; join(' + ', map { "($_->[0]) D[" . join(',', @{$_->[1]}) . "]" } @$P) }
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @o; for my $sp (setparts($n - 1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n - 1; push @o, \@c; } my @c = map { [@$_] } @$sp; push @c, [$n - 1]; push @o, \@c; } @o }

# ---------- empty-profile leading coefficients ----------
# Laurent polynomials in w: hash exponent -> ring element.  Dop_J: w^m -> (J-m) w^{m+1} + m w^{m+2};  dy_J: w^m -> J w^{m+1} + m w^{m+2} (sector J -> J+1).
sub lp_add { my ($R, $h, $e, $c) = @_; $h->{$e} = exists $h->{$e} ? $R->{add}->($h->{$e}, $c) : $c; }
sub lp_D { my ($R, $f, $J) = @_; my %g; for my $m (keys %$f) { my $c = $f->{$m}; lp_add($R, \%g, $m + 1, $R->{mul}->($c, $R->{int}->($J - $m))); lp_add($R, \%g, $m + 2, $R->{mul}->($c, $R->{int}->($m))); } \%g }
sub lp_dy { my ($R, $f, $J) = @_; my %g; for my $m (keys %$f) { my $c = $f->{$m}; lp_add($R, \%g, $m + 1, $R->{mul}->($c, $R->{int}->($J))); lp_add($R, \%g, $m + 2, $R->{mul}->($c, $R->{int}->($m))); } \%g }
my %L0cache;
sub L0 { my ($R, $sig, $g, $k, $v) = @_; die "L0: v >= 1 needed" unless $v >= 1;
  return $R->{z} if $g >= 2;
  my $key = "$R->{name}|$sig|$g|$k|$v"; return $L0cache{$key} if exists $L0cache{$key};
  my ($f, $J);
  if ($g == 0) { $f = { 0 => $R->{rat}->('5/12'), -2 => $R->{rat}->('-1/4'), -3 => $R->{rat}->('-1/6') }; $J = 0;   # u - 3u^2/4 + u^3/6, u = 1 - 1/w
    $f = lp_D($R, $f, 0) for 1..$k; for (1..$v) { $f = lp_dy($R, $f, $J); $J++; } }
  else { my $is = $R->{inv}->($sig); my $c24 = $R->{rat}->('1/24');
    if ($k >= 1) { $f = { 2 => $c24, 1 => $R->{mul}->($c24, $R->{neg}->($R->{add}->($R->{int}->(2), $is))), 0 => $R->{mul}->($c24, $R->{add}->($R->{o}, $is)) }; $J = 0;   # D G1
      $f = lp_D($R, $f, 0) for 2..$k; for (1..$v) { $f = lp_dy($R, $f, $J); $J++; } }
    else { $f = { 2 => $c24, 1 => $R->{mul}->($c24, $R->{neg}->($R->{add}->($R->{o}, $is))) }; $J = 1;   # e^{-u} d_y G1
      for (2..$v) { $f = lp_dy($R, $f, $J); $J++; } } }
  my $nom = 2*$g - 2 + $k + $v;
  for my $m (keys %$f) { die "L0: term w^$m below nominal order $nom (g=$g k=$k v=$v)" if $m < $nom && !$R->{iszero}->($f->{$m}); }
  $L0cache{$key} = $f->{$nom} // $R->{z} }

# ---------- the relation ----------
sub Pi { my ($R, $j, $tau) = @_; my $v = $R->{o}; return $v if $j <= 1; for my $a (1..$j-1) { $v = $R->{mul}->($v, $R->{add}->($R->{mul}->($R->{int}->($j), $tau), $R->{int}->($a))); } $R->{mul}->($v, $R->{inv}->($R->{rat}->(bfact($j - 1)))) }
sub rpow { my ($R, $b, $e) = @_; my $r = $R->{o}; $r = $R->{mul}->($r, $b) for 1..$e; $r }
# evaluate: P = monomial list, G, eps, Rp = parts of R, R = ring, tau = ring element, X = sub(g, [U sorted asc]) -> ring element (g >= 0)
sub relation { my (%a) = @_; my ($P, $G, $eps, $Rp, $R, $tau, $X) = @a{qw(P G eps Rparts ring tau X)};
  my $n = 2*$G + $eps; my $sig = $R->{mul}->($tau, $R->{add}->($tau, $R->{o})); my $r = scalar @$Rp; my $tot = $R->{z};
  my %Pic; my $PiF = sub { $Pic{$_[0]} //= Pi($R, $_[0], $tau) };
  my @tk = ($R->{o}); for my $i (1..$n + 4) { $tk[$i] = $R->{mul}->($tk[$i-1], $tau); }
  my @sk = ($R->{o}); for my $i (1..$n + 8) { $sk[$i] = $R->{mul}->($sk[$i-1], $sig); }
  my $nterms = 0;
  for my $mono (@$P) { my ($c, $alpha) = @$mono; my $q = scalar @$alpha; my $cR = $R->{rat}->($c);
    for my $sp (setparts($q)) { my $nb = scalar @$sp; my @Bi = map { [ map { $alpha->[$_] } @$_ ] } @$sp;
      my @pre = map { my $B = $_; my $f = $R->{o}; for my $i (@$B) { $f = $R->{mul}->($f, $R->{mul}->($R->{int}->($i), $PiF->($i))); } $f } @Bi;
      for my $asg (0..$nb ** $r - 1) { my @T = map { [] } 1..$nb; my $x = $asg; for my $j (0..$r-1) { push @{$T[$x % $nb]}, $Rp->[$j]; $x = int($x / $nb); }
        # enumerate (g_b, k_b)
        my $rec; $rec = sub { my ($bi, $budget, $acc) = @_;
          if ($bi == $nb) { return unless $budget == 0; $tot = $R->{add}->($tot, $R->{mul}->($cR, $acc)); $nterms++; return; }
          my $B = $Bi[$bi]; my $odd = @$B % 2; my @U = sort { $a <=> $b } (@{$T[$bi]}, grep { $_ >= 2 } @$B); my $v = grep { $_ == 1 } @$B;
          my $JU = 0; $JU += $_ for @U; my $DU = $JU - scalar(@U);
          for my $g (0..int(($budget + 2) / 2)) { for my $k (($odd ? 1 : 0)..($budget + 2 - 2*$g - scalar(@$B))) {
              my $wt = 2*$g - 2 + scalar(@$B) + $k; next if $wt < 0 || $wt > $budget;
              my $L;
              if (@U) { my $delta = $DU - 2*$g + 2; next if $delta >= 0 && $delta < $k;   # falling factorial (delta)_k = 0
                my $cf = $R->{o}; for my $i (0..$k-1) { $cf = $R->{mul}->($cf, $R->{int}->($delta - $i)); } for my $i (0..$v-1) { $cf = $R->{mul}->($cf, $R->{int}->($JU + $i)); }
                $L = $R->{mul}->($cf, $X->($g, \@U)); }
              else { $L = L0($R, $sig, $g, $k, $v); }
              next if $R->{iszero}->($L);
              my $f = $pre[$bi]; $f = $R->{neg}->($f) if $odd; $f = $R->{mul}->($f, $R->{int}->(2)) if !$odd && $k == 0;
              $f = $R->{mul}->($f, $sk[scalar(@$B) + $g - 1]); $f = $R->{mul}->($f, $tk[$k]); $f = $R->{mul}->($f, $R->{inv}->($R->{rat}->(bfact($k))));
              $rec->($bi + 1, $budget - $wt, $R->{mul}->($acc, $R->{mul}->($f, $L))); } } };
        $rec->(0, $n, $R->{o}); } } }
  wantarray ? ($tot, $nterms) : $tot }
1;
