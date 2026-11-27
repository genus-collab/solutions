package RefHur;
# Referee's independent Hurwitz-number engine for the Route B report (no code shared with HB.pm / HBconn.pm).
#  - disconnected counts by TWO methods: the cut-and-join recursion on class types (no characters), and the
#    Frobenius character sum (own Murnaghan-Nakayama code); cross-checked against each other;
#  - connected counts by the integer "E-logarithm" recursion  E Z = Z . E H  on sub-monomials of p_R p_1^C
#    (the author used a moment-cumulant recursion over types with a distinguished label);
#  - ELSV-normalised integrals A_h(R,1^c) = Int Lambda^vee(1) prod 1/(1-j psi) prod 1/(1-psi);
#  - p_0 = Xi^(0)_h(R) by fitting the Lambert form e^{-uJ} G(u e^{-u}) = sum_k p_k (1-u)^{-(e+k)} (overdetermined).
# Normalisations.  Zhat(mu)[b] := d! b! [p_mu beta^b] Z = |C_mu| N_b(mu), N_b(mu) = #{(t_1..t_b) transpositions : t_1...t_b = pi},
#   pi a fixed permutation of type mu.  Hhat(mu)[b] := d! b! [p_mu beta^b] log Z  (integer).
#   Connected Hurwitz number H_{g,mu} = Hhat(mu)[b]/d!  (weight 1/|Aut|), b = 2g-2+l(mu)+|mu|.
#   ELSV: H_{g,mu} = (b!/|Aut mu|) prod mu_i^mu_i/mu_i! Int Lambda^vee(1)/prod(1 - mu_i psi_i).
use strict; use warnings;
use Math::BigInt; use Math::BigRat;

sub BI { Math::BigInt->new(shift) }
sub BR { Math::BigRat->new(shift) }

my @FAC = (BI(1));
sub fac { my $n = shift; die "fac($n)" if $n < 0; while (@FAC <= $n) { push @FAC, $FAC[-1] * BI(scalar @FAC) } return $FAC[$n]->copy }
my %BIN;
sub bin { my ($n, $k) = @_; return BI(0) if $k < 0 || $n < 0 || $k > $n; my $key = "$n:$k";
  $BIN{$key} = fac($n) / (fac($k) * fac($n - $k)) unless exists $BIN{$key}; return $BIN{$key}->copy }

# ------------------------------------------------------------------ partitions (descending)
my %PART;
sub parts_of { my $d = shift; return $PART{$d} if $PART{$d};
  my @out; my $gen; $gen = sub { my ($pre, $rem, $max) = @_;
    if ($rem == 0) { push @out, [@$pre]; return }
    my $top = $max < $rem ? $max : $rem; for (my $k = $top; $k >= 1; $k--) { $gen->([@$pre, $k], $rem - $k, $k) } };
  $gen->([], $d, $d); $PART{$d} = \@out; return \@out }
sub pkey { join ",", @{$_[0]} }
sub mkey { my $m = shift; join ",", sort { $b <=> $a } map { ($_) x $m->{$_} } grep { $m->{$_} > 0 } keys %$m }
sub zmu { my $mu = shift; my %m; $m{$_}++ for @$mu; my $z = BI(1); for my $k (keys %m) { $z *= BI($k) ** $m{$k}; $z *= fac($m{$k}) } return $z }
sub classsize { my $mu = shift; my $d = 0; $d += $_ for @$mu; return fac($d) / zmu($mu) }

# ------------------------------------------------------------------ method 1: cut-and-join on class types
# N_b(pi) = sum_t N_{b-1}(pi t); pi t joins two cycles (lengths a,b: a*b transpositions) or cuts one (length m into k, m-k:
# m transpositions if k != m-k, m/2 if k = m-k).  Returns { key => [N_0..N_B] } for all mu |- d.
my %CJ;
sub cj_table { my ($d, $B) = @_;
  my $ck = "$d:$B"; return $CJ{$ck} if $CJ{$ck};
  my $P = parts_of($d); my %idx; $idx{pkey($P->[$_])} = $_ for 0..$#$P;
  my @moves;
  for my $i (0..$#$P) { my %m; $m{$_}++ for @{$P->[$i]}; my @vals = sort { $a <=> $b } keys %m; my %acc; my $tot = 0;
    for my $x (0..$#vals) { for my $y ($x..$#vals) { my ($a, $b) = ($vals[$x], $vals[$y]);
        my $cnt = $x == $y ? ($m{$a} * ($m{$a} - 1) / 2) * $a * $a : $m{$a} * $m{$b} * $a * $b; next unless $cnt;
        my %n = %m; $n{$a}--; $n{$b}--; $n{$a + $b}++; $acc{mkey(\%n)} += $cnt; $tot += $cnt } }
    for my $v (@vals) { next if $v < 2; for my $k (1..int($v / 2)) { my $cnt = $m{$v} * ($k == $v - $k ? $v / 2 : $v);
        my %n = %m; $n{$v}--; $n{$k}++; $n{$v - $k}++; $acc{mkey(\%n)} += $cnt; $tot += $cnt } }
    die "cj: move count $tot != C($d,2)" unless $tot == $d * ($d - 1) / 2;
    $moves[$i] = [ map { [ $idx{$_}, BI($acc{$_}) ] } keys %acc ] }
  my @cur = map { BI(0) } @$P; my $id = $d == 0 ? 0 : $idx{join(",", (1) x $d)}; $cur[$id] = BI(1);
  my %res; for my $i (0..$#$P) { $res{pkey($P->[$i])} = [ $cur[$i]->copy ] }
  for my $b (1..$B) { my @new;
    for my $i (0..$#$P) { my $s = BI(0); for my $mv (@{$moves[$i]}) { my $v = $cur[$mv->[0]]; next if $v->is_zero; $s += $mv->[1] * $v } $new[$i] = $s }
    @cur = @new; for my $i (0..$#$P) { push @{$res{pkey($P->[$i])}}, $cur[$i]->copy } }
  $CJ{$ck} = \%res; return \%res }
sub zhat_cj { my ($nu, $c, $B) = @_; my @mu = (sort({ $b <=> $a } @$nu), (1) x $c); my $d = 0; $d += $_ for @mu;
  my $T = cj_table($d, $B); my $row = $T->{join(",", @mu)}; my $cs = classsize(\@mu);
  return [ map { $row->[$_] * $cs } 0..$B ] }

# ------------------------------------------------------------------ method 2: characters (own Murnaghan-Nakayama)
my %DIM;
sub dimlam { my $l = shift; my $k = pkey($l); return $DIM{$k}->copy if exists $DIM{$k};
  my $n = 0; $n += $_ for @$l; my @conj; for my $i (0..$#$l) { $conj[$_]++ for 0..$l->[$i] - 1 }
  my $hp = BI(1); for my $i (0..$#$l) { for my $j (0..$l->[$i] - 1) { $hp *= BI($l->[$i] - $j + $conj[$j] - $i - 1) } }
  my $r = fac($n) / $hp; $DIM{$k} = $r; return $r->copy }
sub contents { my $l = shift; my $s = 0; for my $i (0..$#$l) { $s += $l->[$i] * ($l->[$i] - 1) / 2 - $i * $l->[$i] } return $s }
my %CHI;
sub chi { my ($l, $nu) = @_;            # chi^l at a permutation of type nu u 1^(|l|-|nu|); nu removed as rim hooks in the given order
  return dimlam($l) unless @$nu;
  my $key = pkey($l) . "|" . pkey($nu); return $CHI{$key}->copy if exists $CHI{$key};
  my ($k, @rest) = @$nu; my $L = scalar @$l;
  my @beta = map { $l->[$_] + $L - 1 - $_ } 0..$L - 1; my %has = map { $_ => 1 } @beta; my $tot = BI(0);
  for my $i (0..$L - 1) { my $nb = $beta[$i] - $k; next if $nb < 0 || $has{$nb};
    my $ht = grep { $_ > $nb && $_ < $beta[$i] } @beta;
    my @nbeta = @beta; $nbeta[$i] = $nb; @nbeta = sort { $b <=> $a } @nbeta;
    my @nl = grep { $_ > 0 } map { $nbeta[$_] - ($L - 1 - $_) } 0..$L - 1;
    my $v = chi(\@nl, \@rest); if ($ht % 2) { $tot -= $v } else { $tot += $v } }
  $CHI{$key} = $tot; return $tot->copy }
sub zhat_char { my ($nu, $c, $B) = @_; my @nus = sort { $b <=> $a } @$nu; my $d = $c; $d += $_ for @nus;
  my @mu = (@nus, (1) x $c); my $z = zmu(\@mu); my @acc = map { BI(0) } 0..$B;
  for my $l (@{parts_of($d)}) { my $ch = chi($l, \@nus); next if $ch->is_zero; my $w = dimlam($l) * $ch; my $f = BI(contents($l)); my $p = BI(1);
    for my $b (0..$B) { $acc[$b] += $w * $p; $p *= $f } }
  for my $b (0..$B) { my ($q, $rem) = $acc[$b]->copy->bdiv($z); die "zhat_char: not integral" unless $rem->is_zero; $acc[$b] = $q }
  return \@acc }

# ------------------------------------------------------------------ connected counts on the lattice of sub-monomials of p_R p_1^C
# Hhat(mu)[b] = Zhat(mu)[b] - sum_{m1 proper non-empty sub-monomial} C(d-1, d1-1) sum_{b1} C(b,b1) Hhat(m1)[b1] Zhat(mu/m1)[b-b1]
# Returns a hash  "a1,..,aq;c" => [Hhat[0..B]] with a = multiplicity vector over the distinct values of R.
sub connected { my ($R, $C, $B, $method) = @_; $method //= 'char';
  my %mult; $mult{$_}++ for @$R; my @vals = sort { $a <=> $b } keys %mult; my $q = @vals;
  my @subs = ([]); for my $i (0..$q - 1) { my @n; for my $s (@subs) { push @n, [@$s, $_] for 0..$mult{$vals[$i]} } @subs = @n }
  my @mons; for my $s (@subs) { for my $c (0..$C) { push @mons, [$s, $c] } }
  my $deg = sub { my $m = shift; my $d = $m->[1]; $d += $m->[0][$_] * $vals[$_] for 0..$q - 1; $d };
  my $mk = sub { my $m = shift; join(",", @{$m->[0]}) . ";" . $m->[1] };
  @mons = sort { $deg->($a) <=> $deg->($b) } @mons;
  my (%Z, %H);
  for my $m (@mons) { my @nu; for my $i (0..$q - 1) { push @nu, ($vals[$i]) x $m->[0][$i] }
    $Z{$mk->($m)} = ($method eq 'cj') ? zhat_cj(\@nu, $m->[1], $B) : zhat_char(\@nu, $m->[1], $B) }
  for my $m (@mons) { my $key = $mk->($m); my $d = $deg->($m); next if $d == 0;
    my @h = map { $_->copy } @{$Z{$key}};
    # proper non-empty sub-monomials m1 of m
    my @s1 = ([]); for my $i (0..$q - 1) { my @n; for my $s (@s1) { push @n, [@$s, $_] for 0..$m->[0][$i] } @s1 = @n }
    for my $a (@s1) { for my $c1 (0..$m->[1]) { my $m1 = [$a, $c1]; my $d1 = $deg->($m1); next if $d1 == 0 || $d1 == $d;
        my $m2 = [[map { $m->[0][$_] - $a->[$_] } 0..$q - 1], $m->[1] - $c1];
        my $H1 = $H{$mk->($m1)}; my $Z2 = $Z{$mk->($m2)}; my $cd = bin($d - 1, $d1 - 1);
        for my $b1 (0..$B) { next if $H1->[$b1]->is_zero; for my $b (($b1)..$B) { my $z2 = $Z2->[$b - $b1]; next if $z2->is_zero;
            $h[$b] -= $cd * bin($b, $b1) * $H1->[$b1] * $z2 } } } }
    $H{$key} = \@h }
  return (\%H, \@vals) }

# A_h(R, 1^c) for c = 0..C  (BigRat)
sub A_list { my ($R, $h, $C, $method) = @_; my $r = @$R; my $J = 0; $J += $_ for @$R;
  my $B = 2 * $h - 2 + $r + $J + 2 * $C; $B = 0 if $B < 0;
  my ($H, $vals) = connected($R, $C, $B, $method); my %mult; $mult{$_}++ for @$R;
  my $aut = BI(1); $aut *= fac($_) for values %mult; my $kap = BR(1); for my $j (@$R) { $kap *= BR(fac($j)) / BR(BI($j) ** $j) }
  my $full = join(",", map { $mult{$_} } @$vals); my @A;
  for my $c (0..$C) { my $d = $J + $c; my $b = 2 * $h - 2 + $r + $J + 2 * $c;
    if ($b < 0 || $d == 0) { push @A, BR(0); next }      # d = 0: empty cover; Int_{Mbar_h} Lambda^vee(1) = 0 for h >= 2
    my $hh = $H->{"$full;$c"}[$b];
    push @A, BR($hh) / (BR(fac($d)) * BR(fac($b))) * BR($aut * fac($c)) * $kap }
  return \@A }

# Lambert fit: F(u) = sum_c (A_c/c!) u^c e^{-(J+c)u} = sum_{k=0}^{n0} p_k (1-u)^{-(e+k)}.
# P(u) := (1-u)^{e+n0} F(u) must be a polynomial of degree <= n0; p_0 = (-1)^{n0} [u^{n0}] P.  Returns (p0, ok_structure, [p_k]).
sub lambert_fit { my ($A, $h, $r, $J, $N) = @_;       # A must have entries 0..N
  my $e = 2 * $h - 2 + $r; my $n0 = 3 * $h - 3 + $r;
  my @F = map { BR(0) } 0..$N;
  for my $n (0..$N) { for my $c (0..$n) { my $t = $A->[$c] / BR(fac($c)) * BR(BI(-($J + $c)) ** ($n - $c)) / BR(fac($n - $c)); $F[$n] += $t } }
  my $m = $e + $n0; my @P = map { BR(0) } 0..$N;     # (1-u)^m, m >= 0 here
  for my $n (0..$N) { for my $i (0..$n) { next if $i > $m; my $t = BR(bin($m, $i)) * $F[$n - $i]; if ($i % 2) { $P[$n] -= $t } else { $P[$n] += $t } } }
  my $ok = 1; for my $n ($n0 + 1..$N) { $ok = 0 unless $P[$n]->is_zero }
  # recover p_k from P(u) = sum_k p_k (1-u)^{n0-k}: substitute v = 1-u
  my @Pv = map { BR(0) } 0..$n0;   # P as polynomial in v
  for my $n (0..$n0) { # u^n = (1-v)^n
    for my $i (0..$n) { my $t = $P[$n] * BR(bin($n, $i)); if ($i % 2) { $Pv[$i] -= $t } else { $Pv[$i] += $t } } }
  my @p = map { $Pv[$n0 - $_] } 0..$n0;
  return ($p[0], $ok, \@p) }

# Proposition D (finite formula), as an alternative extraction of p_0
sub propD { my ($A, $h, $r, $J) = @_; my $e = 2 * $h - 2 + $r; my $n0 = 3 * $h - 3 + $r; my $tot = BR(0);
  for my $c (0..$n0) { my $in = BR(0);
    for my $i (0..$n0 - $c) { $in += BR(bin($n0 + $e, $n0 - $c - $i)) * BR(BI($J + $c) ** $i) / BR(fac($i)) }
    my $t = $A->[$c] / BR(fac($c)) * $in; if ($c % 2) { $tot -= $t } else { $tot += $t } }
  return $tot }

sub xi0 { my ($R, $h, $extra, $method) = @_; $extra //= 2; my $r = @$R; my $J = 0; $J += $_ for @$R; my $n0 = 3 * $h - 3 + $r;
  my $A = A_list($R, $h, $n0 + $extra, $method);
  my ($p0, $ok, $p) = lambert_fit($A, $h, $r, $J, $n0 + $extra);
  my $pd = propD($A, $h, $r, $J);
  return ($p0, $ok, $pd, $p, $A) }
1;
