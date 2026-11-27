# t_engine_cmp.pl h r tau1,tau2,... pt1 pt2 ... : the referee engine (dilaton leaves as explicit points, no kappa classes)
# against the author's engine polynomials ../engine_h<h>_r<r>.txt (components Xi^(k), Xi = Sum_k (-1/sigma)^k Xi^(k), c = 1).
# A point is a colon-separated list, e.g. 2:3 ; rational entries allowed (1/2).  Prints both values and the verdict.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($h, $r, $taus, @pts) = @ARGV; my $t0 = time;
my $f = "$FindBin::Bin/../engine_h${h}_r${r}.txt";
my %C;   # C{k} = [ [exps], coef ]
open my $fh, '<', $f or die "no $f"; while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4; $rest //= '';
  for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; push @{$C{$k}}, [[split /,/, $m], Q($cf)]; } } close $fh;
printf "loaded %s [%.1f s]\n", $f, time - $t0;
sub author { my ($p, $tau) = @_; my $sig = $tau * (1 + $tau); my $v = Q(-1) / $sig; my $s = Q(0);
  for my $k (0..$h) { my $c = Q(0); for my $t (@{$C{$k} // []}) { my ($e, $cf) = @$t; my $x = $cf->copy; for my $i (0..$r-1) { $x *= $p->[$i] ** $e->[$i] if $e->[$i]; } $c += $x; } $s += $c * $v ** $k; }
  return $s; }
my ($n, $bad) = (0, 0);
for my $tau (map { Q($_) } split /,/, $taus) { my $maxd = $h == 1 ? 1 : 3*$h - 3; %RefHodge::IM = () if 0;
  my $om = RefHodge::omega($h, [Q(1), $tau, -1 - $tau], $maxd);
  for my $ps (@pts) { my @p = map { Q($_) } split /:/, $ps; die "wrong r" unless @p == $r;
    my $q = RefHodge::xi($h, $om, [map { "$_" } @p]); my $mine = RefHodge::xi_value($q, Q(-1)); my $theirs = author(\@p, $tau);
    $n++; my $ok = $mine == $theirs; $bad++ unless $ok;
    printf "genus %d R = {%s} tau = %s: referee %s | author %s : %s [%.1f s]\n", $h, join(',', @p), $tau, ($mine == 0 ? "0" : "non-zero"), ($theirs == 0 ? "0" : "non-zero"), ($ok ? "AGREE" : "DISAGREE ($mine vs $theirs)"), time - $t0; } }
printf "SUMMARY genus %d r = %d: %d comparisons, %d disagreements, memo %d, %.1f s\n", $h, $r, $n, $bad, scalar(keys %RefHodge::IM), time - $t0;
