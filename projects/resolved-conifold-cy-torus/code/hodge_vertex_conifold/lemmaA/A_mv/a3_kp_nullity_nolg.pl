#!/usr/bin/perl
# a3_kp_nullity_nolg.pl -- variant WITHOUT the lambda_g input (sigma^0 coefficients free) of a3_kp_nullity.pl -- what do the KP hierarchy (dispersionless Hirota at the dilaton point) and the cut-and-join 2-insertion relations
# determine at genus G, assuming Lemma A in all genera < G (induction)?   (report A3.5-A3.6)
# Unknowns: sigma^G Xi_G(T) = sum_{k=1}^{G} c_{T,k} sigma^k for all T (parts >= 2) with D(T) <= 2G-3 (lambda_g term = 0 in the range).
# Relations (exact, linear in the genus-G unknowns, genus-0 factors explicit: Xi_0(T) = Gamma(J)/Gamma(J-r+3)):
#  KP:  for every R and (al,be) with D(R)+al+be-2 <= 2G-3:  [z^-al w^-be p_R] ( exp(Omega) )_lin = om^(G)_{1,al+be-1}(R),
#       Omega = sum_{a,b} z^-a w^-b sum_g om^(g)_ab,  om^(g)_ab(R) = sigma^{g+1} Pi_a Pi_b Xi_g(R+{a,b})  (Pi_1 = 1, unit parts by string),
#       "_lin" = terms with exactly one genus-G factor and genus-0 factors otherwise.
#  CJ:  the 2-insertion relation of a3_relation.pl with its genus-(0,G) product terms (optional, flag C).
# Evaluated at random tau mod p; prints rank and nullity, and which T are not forced to vanish.
# Usage: perl a3_kp_nullity.pl G [C] [Z=t1.t2 ...]
use strict; use warnings;
my ($G, @opt) = @ARGV; my %opt = map { $_ => 1 } @opt;
my $p = 2147483647;
sub md { my $x = shift; $x %= $p; $x += $p if $x < 0; return $x; }
sub mm { return md($_[0] * $_[1]); }
sub mi { my $a = md(shift); my ($t, $nt, $r, $nr) = (0, 1, $p, $a); while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q*$nt); ($r, $nr) = ($nr, $r - $q*$nr); } die "not invertible" if $r != 1; return md($t); }
sub parts { my ($n, $max) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
sub key { return join(',', sort { $b <=> $a } @_); }
my @T; my %Tidx; for my $D (1..2*$G-3) { for my $q (parts($D, $D)) { my @t = map { $_ + 1 } @$q; $Tidx{key(@t)} = scalar @T; push @T, \@t; } }
my $nT = scalar @T; my $W = $G + 1; my $nU = $nT * $W;
my @Rs = ([]); for my $D (1..2*$G-3) { for my $q (parts($D, $D)) { push @Rs, [ map { $_ + 1 } @$q ]; } }
srand(4711);
my @rows;
for my $sample (1..(4*$G + 12)) {
  my $tau = 3 + int(rand($p - 10)); my $sig = mm($tau, md($tau + 1)); my $isig = mi($sig); my $sp = md(2*$tau + 1);
  my %Pc; my $Pi = sub { my $j = shift; return 1 if $j <= 1; return $Pc{$j} //= do { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } my $f = 1; $f = mm($f, $_) for 1..$j-1; mm($v, mi($f)); }; };
  my @sp = (1); push @sp, mm($sp[-1], $sig) for 1..2*$G+4;
  # genus-0 value Xi_0(list) mod p
  my $X0 = sub { my @L = @_; my $J = 0; $J += $_ for @L; my $r = scalar @L;
    if ($r >= 3) { my $v = 1; $v = mm($v, md($J - $_)) for 1..($r-3); return $v; } elsif ($r == 2) { return mi($J); } else { return mi(mm($J, $J+1)); } };
  # genus-G value as a linear form (hash col => coef): Xi_G(list) = string factor * sigma^{-G} sum_k c_{T,k} sigma^k
  my $XG = sub { my @L = @_; my @Tp = grep { $_ >= 2 } @L; my $b = scalar(@L) - scalar(@Tp); my $J = 0; $J += $_ for @Tp;
    return {} unless @Tp;   # Xi_G(1^b) = 0 for G >= 2
    my $f = 1; $f = mm($f, md($J + $_)) for 0..$b-1;
    my $i = $Tidx{key(@Tp)}; die "XG out of range: @L" unless defined $i;
    my $isg = 1; $isg = mm($isg, $isig) for 1..$G; my %v; for my $k (0..$G) { $v{$i*$W + $k} = mm($f, mm($isg, $sp[$k])); } return \%v; };
  my $addf = sub { my ($acc, $v, $c) = @_; for my $k (keys %$v) { $acc->{$k} = md(($acc->{$k} // 0) + mm($c, $v->{$k})); } };
  # ---- KP relations ----
  for my $R (@Rs) { my @R = @$R; my $r = scalar @R; my $DR = 0; $DR += $_ - 1 for @R;
    my $wmax = 2*$G - 1 - $DR; next if $wmax < 4;   # need al, be >= 2
    # multiset structure of R: distinct values @vv with multiplicities @nn; sub-multisets encoded as vectors -> index
    my %m; $m{$_}++ for @R; my @vv = sort { $a <=> $b } keys %m; my @nn = map { $m{$_} } @vv; my $nv = scalar @vv;
    my @subs = ([]); for my $i (0..$nv-1) { my @ns; for my $s (@subs) { for my $k (0..$nn[$i]) { push @ns, [@$s, $k]; } } @subs = @ns; }
    my %sidx; for my $i (0..$#subs) { $sidx{join(",", @{$subs[$i]})} = $i; }
    my @list = map { my $s = $_; [ map { ($vv[$_]) x $s->[$_] } 0..$nv-1 ] } @subs;   # parts of each sub-multiset
    my $fullk = join(",", @nn); my $full = $sidx{$fullk};
    # binomial weights and complements
    my %bin; my $binom = sub { my ($n, $k) = @_; my $v = 1; for my $i (1..$k) { $v = mm($v, mm($n - $i + 1, mi($i))); } return $v; };
    my @om0; for my $si (0..$#subs) { for my $a (1..$wmax) { for my $b (1..$wmax-$a) { $om0[$si][$a][$b] = mm($sig, mm(mm($Pi->($a), $Pi->($b)), $X0->(@{$list[$si]}, $a, $b))); } } }
    # sub-multiset pairs: for S (index), list of (k index, complement index, weight)
    my @pairs; for my $si (0..$#subs) { my $S = $subs[$si]; my @pl = ([]);
      for my $i (0..$nv-1) { my @np; for my $p (@pl) { for my $k (0..$S->[$i]) { push @np, [@$p, $k]; } } @pl = @np; }
      for my $kv (@pl) { my $w = 1; for my $i (0..$nv-1) { $w = mm($w, $binom->($S->[$i], $kv->[$i])); }
        my @cv = map { $S->[$_] - $kv->[$_] } 0..$nv-1; push @{$pairs[$si]}, [ $sidx{join(",", @$kv)}, $sidx{join(",", @cv)}, $w ]; } }
    my %Em; my $Erec; $Erec = sub { my ($a, $b, $si) = @_; my $kk = "$a,$b,$si"; return $Em{$kk} if exists $Em{$kk}; my $res = 0;
      if ($a == 0 && $b == 0) { $res = ($si == $sidx{join(",", (0) x $nv)}) ? 1 : 0; }
      elsif ($a <= 0 || $b <= 0) { $res = 0; }
      else { for my $a1 (1..$a) { for my $b1 (1..$b) { for my $pr (@{$pairs[$si]}) { my ($ki, $ci, $w) = @$pr; my $e = $Erec->($a-$a1, $b-$b1, $ci); next unless $e;
              $res = md($res + mm(mm($w, mm($a1, $om0[$ki][$a1][$b1])), $e)); } } } $res = mm($res, mi($a)); }
      $Em{$kk} = $res; return $res; };
    for my $al (2..$wmax) { for my $be ($al..$wmax) { next if $al + $be > $wmax;
      my %row;
      for my $a (1..$al) { for my $b (1..$be) { for my $pr (@{$pairs[$full]}) { my ($ki, $ci, $w) = @$pr; my $e = $Erec->($al-$a, $be-$b, $ci); next unless $e;
            my $v = $XG->(@{$list[$ki]}, $a, $b); my $c = mm(mm($w, $e), mm($sp[$G+1], mm($Pi->($a), $Pi->($b)))); $addf->(\%row, $v, $c); } } }
      my $cc = $al + $be - 1; my $v = $XG->(@R, 1, $cc); $addf->(\%row, $v, md(-mm($sp[$G+1], $Pi->($cc))));
      push @rows, \%row; } } }
  # ---- cut-and-join 2-insertion relations (optional) ----
  if ($opt{C}) { for my $R (@Rs) { my @R = @$R; my $r = scalar @R; my $DR = 0; $DR += $_ - 1 for @R; next if $DR + 1 > 2*$G - 3;
      my %row; $addf->(\%row, $XG->(@R, 2), 1);
      for my $i (0..$r-1) { my $j = $R[$i]; my @Rp = @R; $Rp[$i] = $j+1; my $gam = mm(mm($j+1, $Pi->($j+1)), mi(mm($Pi->($j), $sp))); $addf->(\%row, $XG->(@Rp), md(-$gam)); }
      for my $i (0..$r-1) { for my $k ($i+1..$r-1) { my $jj = $R[$i] + $R[$k]; my @Rm = ((map { $R[$_] } grep { $_ != $i && $_ != $k } 0..$r-1), $jj);
          my $gam = mm(mm($jj, $Pi->($jj)), mi(mm(mm($Pi->($R[$i]), $Pi->($R[$k])), $sp))); $addf->(\%row, $XG->(@Rm), $gam); } }
      for my $i (0..$r-1) { my $j = $R[$i]; my @rest = map { $R[$_] } grep { $_ != $i } 0..$r-1; my $nr = scalar @rest;
        for my $a (1..$j-1) { my $b = $j - $a; my $cf = mm(mm($sig, mi(mm(2, $sp))), mm(mm($a, $b), mm(mm($Pi->($a), $Pi->($b)), mi($Pi->($j)))));
          for my $m (0..(1 << $nr) - 1) { my (@I, @Jx); for my $q (0..$nr-1) { if (($m >> $q) & 1) { push @I, $rest[$q]; } else { push @Jx, $rest[$q]; } }
            # (g1,g2) = (0,G): Xi_0(a,I) Xi_G(b,J);  (G,0): Xi_G(a,I) Xi_0(b,J)
            $addf->(\%row, $XG->($b, @Jx), mm($cf, $X0->($a, @I)));
            $addf->(\%row, $XG->($a, @I), mm($cf, $X0->($b, @Jx))); } } }
      push @rows, \%row; } }
}
# Z=t1.t2: impose full vanishing of T; H=t..: impose its Hurwitz component (coefficient of sigma^G); L=t..: its lambda_1 component (sigma^(G-1))
for my $o (@opt) { if ($o =~ /^([ZHL])=([\d.]+)$/) { my ($ty, $tl) = ($1, $2); my $i = $Tidx{key(split /\./, $tl)}; die "no T $tl" unless defined $i;
    if ($ty eq 'Z') { push @rows, { $i*$G + $_ - 1 => 1 } for 1..$G; } elsif ($ty eq 'H') { push @rows, { $i*$G + $G - 1 => 1 }; } else { push @rows, { $i*$G + $G - 2 => 1 }; } } }
my $rank = 0; my %piv;
for my $row (@rows) { last if $rank == $nU; my %v = %$row; for my $c (keys %v) { delete $v{$c} unless $v{$c}; }
  for my $col (sort { $a <=> $b } keys %piv) { next unless $v{$col}; my $f = $v{$col}; my $pr = $piv{$col}; for my $c2 (keys %$pr) { $v{$c2} = md(($v{$c2} // 0) - mm($f, $pr->{$c2})); } delete $v{$col}; }
  for my $c (keys %v) { delete $v{$c} unless $v{$c}; } next unless %v;
  my $c0 = (sort { $a <=> $b } keys %v)[0]; my $inv = mi($v{$c0}); my %nr = map { $_ => mm($v{$_}, $inv) } keys %v;
  for my $pc (keys %piv) { my $pr = $piv{$pc}; if ($pr->{$c0}) { my $f = $pr->{$c0}; for my $c2 (keys %nr) { $pr->{$c2} = md(($pr->{$c2} // 0) - mm($f, $nr{$c2})); } delete $pr->{$c0}; for my $c (keys %$pr) { delete $pr->{$c} unless $pr->{$c}; } } }
  $piv{$c0} = \%nr; $rank++; }
my @free; for my $i (0..$nT-1) { my $forced = 1; for my $k (0..$G) { my $col = $i*$W + $k; if (!exists $piv{$col}) { $forced = 0; last; } my @o = grep { $_ != $col } keys %{$piv{$col}}; if (@o) { $forced = 0; last; } } push @free, "(" . join(',', @{$T[$i]}) . ")" unless $forced; }
printf "genus %d: %d multisets in the range, %d unknowns, %d equations (KP%s%s); rank %d, nullity %d\n   not forced to vanish: %s\n",
  $G, $nT, $nU, scalar(@rows), ($opt{C} ? " + cut-and-join" : ""), (join("", map { " $_" } grep { /^[ZHL]=/ } @opt)), $rank, $nU - $rank, (@free ? join(' ', @free) : 'none');
