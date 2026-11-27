#!/usr/bin/perl
# t_generator.pl -- the generator of the leaf-curve evolution as an explicit vector field (dispersionless formula):
#   v_t(lambda) = -h(lambda) ( [d_t P(lambda;t) / h(lambda)]_+ + alpha'(t) ) / ( alpha(t) - [d_lambda P(lambda;t)]_+ ),
#   h = u/u' = -lambda(lambda-1)/(lambda-1-tau),  P(mu;t) = sum_j Phat_j(t) u(mu)^{-j},  Phat_j = P_j alpha^j,  [.]_+ = polynomial part.
# (1) checks the branch-point condition  [d_t P / h]_+ (1+tau) = -alpha'(t)  (needed for v_t to be a polynomial);
# (2) compares v_t with the generator d_t phi / d_lambda phi of the flow fitted from the ODE state (file of xi_m).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $file) = @ARGV; my $tau = R($taus);
my $NEG = 3*$M + 6;
# leaf data
my $N = $M + 3; my $Tt = T_of_t_series($N); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N), $Tt, $N);
my %Pj; for my $j (2..$M+3) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N), $Tt, $N); my $p = [R(1), map { R(0) } 1..$N]; $p = series_mul($p, $h1, $N) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N), $N);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$N ]; }
my @al = (R(1), map { R(0) } 1..$N);
for my $it (0..$N+1) { my @nw = (R(1), map { R(0) } 1..$N); my $wp = [@al]; for my $j (2..$M+3) { $wp = series_mul($wp, \@al, $N); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $N); $nw[$_] = $nw[$_] + $q->[$_] for 0..$N; } @al = @nw; }
my @dal = map { ($_+1) * ($al[$_+1] // R(0)) } 0..$N-1;
my %Ph; for my $j (2..$M+3) { my $ap = [R(1), map { R(0) } 1..$N]; $ap = series_mul($ap, \@al, $N) for 1..$j; $Ph{$j} = series_mul($Pj{$j}, $ap, $N); }
# Laurent series in mu (hash power -> coeff), truncated below -NEG
sub lmul { my ($a, $b) = @_; my %c; for my $i (keys %$a) { for my $j (keys %$b) { my $k = $i + $j; next if $k < -$NEG; $c{$k} = ($c{$k} // R(0)) + $a->{$i}*$b->{$j}; } } return \%c; }
sub ladd { my ($a, $b, $s) = @_; my %c = %$a; for my $k (keys %$b) { $c{$k} = ($c{$k} // R(0)) + ($s // R(1))*$b->{$k}; } return \%c; }
sub lpos { my $a = shift; my %c; for my $k (keys %$a) { $c{$k} = $a->{$k} if $k >= 0; } return \%c; }
sub leval { my ($a, $x) = @_; my $v = R(0); for my $k (keys %$a) { $v = $v + $a->{$k} * $x**$k; } return $v; }
# u^{-1}(mu) = mu (1 - 1/mu)^{-tau} ; u'/u = -1/mu * (1 - (1+tau)/mu)/(1 - 1/mu)
my %uinv; for my $i (0..$NEG+1) { my $b = R(1); for my $q (0..$i-1) { $b = $b * ($tau + $q) / ($q + 1); } $uinv{1-$i} = $b; }
my %upu = (-1 => R(-1)); { my %s; for my $k (1..$NEG) { $upu{-1-$k} = ($upu{-1-$k} // R(0)) + $tau; } }   # -1/mu + tau/mu^2 + tau/mu^3 + ...
my @upow = ({0 => R(1)}); for my $j (1..$M+3) { $upow[$j] = lmul($upow[$j-1], \%uinv); }
# d_t P (mu) at each t-order m: sum_j Phat_j'(t)[m] u^{-j}
my $bad = 0;
for my $m (0..$M-1) { my %dP; for my $j (2..$M+3) { my $c = ($m+1) * ($Ph{$j}[$m+1] // R(0)); next if $c->is_zero; %dP = %{ ladd(\%dP, $upow[$j], $c) }; }
  my $q = lpos(lmul(\%dP, \%upu)); my $val = leval($q, 1 + $tau);
  my $ok = ($val == -$dal[$m]); $bad++ unless $ok;
  print "t^$m: [d_t P * u'/u]_+ (1+tau) = $val ; -alpha'(t) = ", -$dal[$m], "  ", ($ok ? "OK" : "FAIL"), "\n"; }
print "branch-point condition: ", ($bad ? "$bad failures" : "holds to order $M"), "\n";
# ---- (2) explicit generator v_t and comparison with the fitted flow
exit 0 unless $file;
my %xi; open my $fh, '<', $file or die; while (<$fh>) { next unless /^xi_(\d+) = (.*?)\s+\(rank (\d+)\/(\d+)/; my ($mm, $rest, $r, $n) = ($1, $2, $3, $4); next unless $r == $n; my @c; while ($rest =~ /\(([^)]*)\)l\^(\d+)/g) { $c[$2] = R($1); } $xi{$mm} = [ map { $_ // R(0) } @c ]; } close $fh;
my $MF = 0; $MF++ while exists $xi{$MF+1}; print "fitted flow known to order $MF\n";
# polynomials in lambda: arrays; t-series of polynomials: arrays of arrays
sub padd { my ($a, $b, $s) = @_; $s //= R(1); my $n = @$a > @$b ? @$a : @$b; return [ map { ($a->[$_] // R(0)) + $s*($b->[$_] // R(0)) } 0..$n-1 ]; }
sub pmulp { my ($a, $b) = @_; return [R(0)] unless @$a && @$b; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i]->is_zero; for my $j (0..$#$b) { $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
sub tmul { my ($A, $B, $K) = @_; my @C = map { [R(0)] } 0..$K; for my $i (0..$K) { for my $j (0..$K-$i) { $C[$i+$j] = padd($C[$i+$j], pmulp($A->[$i], $B->[$j])); } } return \@C; }
sub h2p { my $h = shift; my $mx = 0; for (keys %$h) { $mx = $_ if $_ > $mx; } return [ map { $h->{$_} // R(0) } 0..$mx ]; }
my $K = $MF - 1;
# N_t = [d_t P / h]_+ + alpha'  ;  [d_lambda P]_+  ;  D_t = alpha - [d P]_+
my (@Nt, @Dt);
for my $m (0..$K) { my %dP; for my $j (2..$M+3) { my $c = ($m+1) * ($Ph{$j}[$m+1] // R(0)); next if $c->is_zero; %dP = %{ ladd(\%dP, $upow[$j], $c) }; }
  my $q = h2p(lpos(lmul(\%dP, \%upu))); $q->[0] = $q->[0] + $dal[$m]; $Nt[$m] = $q;
  my %P; for my $j (2..$M+3) { my $c = $Ph{$j}[$m] // R(0); next if $c->is_zero; %P = %{ ladd(\%P, $upow[$j], $c) }; }
  my %dPl; for my $k (keys %P) { next if $k == 0; $dPl{$k-1} = ($dPl{$k-1} // R(0)) + R($k)*$P{$k}; }
  my $dp = h2p(lpos(\%dPl)); my $d = [ map { -$_ } @$dp ]; $d->[0] = ($d->[0] // R(0)) + $al[$m]; $Dt[$m] = $d; }
# divide N_t by (lambda - 1 - tau) exactly
my @Nq; for my $m (0..$K) { my @p = @{$Nt[$m]}; my @q; my $r = R(0); for (my $i = $#p; $i >= 1; $i--) { $q[$i-1] = $p[$i] + $r; $r = $q[$i-1]*(1+$tau); } my $rem = $p[0] + $r; die "N_t not divisible at t^$m (rem $rem)" unless $rem->is_zero; $Nq[$m] = [ map { $_ // R(0) } @q ]; }
# inverse of D_t as t-series (D_t[0] = 1)
my @Di = ([R(1)]); for my $m (1..$K) { my $s = [R(0)]; for my $k (1..$m) { $s = padd($s, pmulp($Dt[$k], $Di[$m-$k])); } $Di[$m] = [ map { -$_ } @$s ]; }
my $v = tmul(\@Nq, \@Di, $K); for my $m (0..$K) { $v->[$m] = pmulp($v->[$m], [R(0), R(-1), R(1)]); }   # times lambda(lambda-1)
# fitted generator: d_t phi / d_lambda phi
my @dtphi = map { pmulp([ R($_+1) ], $xi{$_+1}) } 0..$K; my @dlphi = ([R(1)]); for my $m (1..$K) { my $p = $xi{$m}; $dlphi[$m] = [ map { R($_+1) * ($p->[$_+1] // R(0)) } 0..$#$p-1 ]; }
my @dli = ([R(1)]); for my $m (1..$K) { my $s = [R(0)]; for my $k (1..$m) { $s = padd($s, pmulp($dlphi[$k], $dli[$m-$k])); } $dli[$m] = [ map { -$_ } @$s ]; }
my $vf = tmul(\@dtphi, \@dli, $K);
my $badv = 0; for my $m (0..$K) { my $d = padd($v->[$m], $vf->[$m], R(-1)); my $nz = grep { !$_->is_zero } @$d; $badv++ if $nz;
  print "t^$m: v_t = ", join(' + ', map { "(" . ($v->[$m][$_] // 0) . ")l^$_" } 0..$#{$v->[$m]}), "   fitted: ", ($nz ? "DIFFERENT" : "same"), "\n"; }
print "generator check: ", ($badv ? "$badv orders differ" : "explicit generator = fitted flow generator to order $K"), "\n";
