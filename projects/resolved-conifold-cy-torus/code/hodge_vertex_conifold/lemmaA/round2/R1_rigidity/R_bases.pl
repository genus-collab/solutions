#!/usr/bin/perl
# R_bases.pl -- the normalised closed-form relations R_k = relation_k / kappa_k (k even) as polynomials in tau, written in three bases:
#   y-basis  (y_L = Pi_{L+1} Xi(L+1)),  F-basis (F_n = [y^n] F(y), F(y) = V~(u(y)), u = y(1-y)^tau),  f-basis (f_m = [X^m] V, X = y/(1+tau y)).
# Input: out/beta_poly_K14.txt (exact beta_{k,L}(tau)).  Output: coefficients factored as theta^e * poly(sigma) when possible.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $file = shift // "$FindBin::Bin/out/beta_poly_K14.txt"; my $Kmax = shift // 14;
my %B; open my $fh, '<', $file or die; while (<$fh>) { next unless /^beta (\d+) (\d+) : (.*)$/; $B{"$1,$2"} = [ map { R($_) } split ' ', $3 ]; } close $fh;
my @Bn = (R(1)); for my $m (1..40) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
# polynomial helpers (arrays of BigRat, index = power of tau)
sub padd { my ($a, $b) = @_; my $n = @$a > @$b ? @$a : @$b; return [ map { ($a->[$_] // R(0)) + ($b->[$_] // R(0)) } 0..$n-1 ]; }
sub pmul { my ($a, $b) = @_; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i]->is_zero; for my $j (0..$#$b) { $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
sub pscal { my ($a, $s) = @_; return [ map { $_ * $s } @$a ]; }
sub pdiv { my ($p, $d) = @_; my @p = @$p; my @q; my $dd = $#$d; while (@p - 1 >= $dd && grep { !$_->is_zero } @p) { my $c = $p[-1] / $d->[-1]; my $sh = $#p - $dd; $q[$sh] = $c; for my $i (0..$dd) { $p[$sh+$i] = $p[$sh+$i] - $c*$d->[$i]; } pop @p; } @q = map { $_ // R(0) } @q; return (\@q, \@p); }
sub sigform { my $p = shift; my @p = @$p; pop @p while @p && $p[-1]->is_zero; return "0" unless @p;
  my $pref = ''; my ($q, $r) = pdiv(\@p, [R('1/2'), R(1)]); if (!grep { !$_->is_zero } @$r) { # divisible by theta: check parity
    my $test = $q; my @f; my $ok = 1; my $cur = [@$test]; while (grep { !$_->is_zero } @$cur) { my ($qq, $rr) = pdiv($cur, [R(0), R(1), R(1)]); my @rr = @$rr; push @rr, R(0) while @rr < 2; if (!$rr[1]->is_zero) { $ok = 0; last; } push @f, $rr[0]; $cur = $qq; last unless @$qq; }
    return "theta*[" . join(', ', @f) . "]" if $ok; }
  my @f; my $cur = [@p]; while (grep { !$_->is_zero } @$cur) { my ($qq, $rr) = pdiv($cur, [R(0), R(1), R(1)]); my @rr = @$rr; push @rr, R(0) while @rr < 2; return "tau-poly[" . join(', ', @p) . "]" unless $rr[1]->is_zero; push @f, $rr[0]; $cur = $qq; last unless @$qq; }
  return "[" . join(', ', @f) . "]"; }
# binom(c tau + c - n - 1, m) as a polynomial in tau
sub gbin { my ($c, $n, $m) = @_; my $p = [R(1)]; for my $i (0..$m-1) { $p = pmul($p, [R($c - $n - 1 - $i), R($c)]); } return pscal($p, R(1)/fact($m)); }
for (my $k = 2; $k <= $Kmax; $k += 2) { my $kap = $Bn[$k] / ($k * fact($k));
  my %Q; for my $L (1..$k) { my $b = $B{"$k,$L"} // [R(0)]; $Q{$L} = pscal($b, R(1)/$kap); }
  print "R_$k (y-basis): ", join('  ', map { "y$_:" . sigform($Q{$_}) } grep { grep { !$_->is_zero } @{$Q{$_}} } 1..$k), "\n";
  # F-basis: coefficient of F_n = sum_{c>=n} Q_{k,c-1} (n/c^2) binom(c tau + c - n - 1, c - n)
  my @outF; for my $n (2..$k+1) { my $acc = [R(0)]; for my $c ($n..$k+1) { my $q = $Q{$c-1} or next; $acc = padd($acc, pscal(pmul($q, gbin($c, $n, $c - $n)), R($n)/R($c*$c))); }
    push @outF, "F$n:" . sigform($acc) if grep { !$_->is_zero } @$acc; }
  print "R_$k (F-basis): ", join('  ', @outF), "\n";
}
