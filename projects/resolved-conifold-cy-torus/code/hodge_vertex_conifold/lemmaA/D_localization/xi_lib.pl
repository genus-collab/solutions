# xi_lib.pl -- the Lemma A functional Xi_h(R) and variants (c = 1 normalisation), on top of hodge_lib.pl
#   w = (1, tau, -1-tau), Omega = prod Lambda^vee(w_a)/w_a^h, leaves lam_j = 1/(1 - j psi), dilaton leaves T = -psi^2/(1 - psi)
#   Xi_h(R) = sum_m 1/m! Int_{Mbar_{h,r+m}} Omega prod lam_{j_i} prod_{k=1}^m T(psi_{r+k})
# Everything is first computed as tau-free "components" I_{k1,k3}(R) = sum_m 1/m! Int L^k1 C^k3 prod lam prod T
# (L = lambda_1, C = ch_3), then contracted with the coefficients of any Hodge class.
use strict; use warnings;
use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/hodge_lib.pl"; }
our %XC;
sub monomials { my $g = shift; my @m; for my $k3 (0..2) { for my $k1 (0..9) { push @m, [$k1, $k3] if $k1 + 3*$k3 <= dmax($g); } } return @m; }
sub group_leaves { my ($N, @R) = @_; my %c; $c{$_}++ for @R; return map { [ser_geom(Q($_), $N), $c{$_}] } sort { $a <=> $b } keys %c; }
# components with a dilaton-type series D (array ref, starting at psi^2) summed with 1/m!
sub components {
  my ($h, $Dser, $Dname, @R) = @_;
  my $key = join("|", $h, $Dname, sort { $a <=> $b } @R);
  return $XC{$key} if exists $XC{$key};
  my $r = scalar @R; my $N = 3*$h - 3 + $r + 2;
  my %I;
  for my $mo (monomials($h)) { my ($k1, $k3) = @$mo; my $tot = Q(0);
    for (my $m = 0; ; $m++) { my $dim = 3*$h - 3 + $r + $m; last if 2*$m > $dim - $k1 - 3*$k3; last if !defined $Dser && $m > 0;
      my @grp = group_leaves($N + $m, @R); push @grp, [$Dser, $m] if $m > 0;
      $tot += integral_groups($h, $k1, $k3, @grp) / fact($m); }
    $I{"$k1,$k3"} = $tot; }
  $XC{$key} = \%I; return \%I;
}
sub contract { my ($cls, $I) = @_; my $t = Q(0); for my $k (keys %$cls) { $t += $cls->{$k} * (exists $I->{$k} ? $I->{$k} : Q(0)); } return $t; }
sub Tseries { my $N = shift; return ser_T(Q(1), $N); }
sub xi_comp { my ($h, @R) = @_; my $N = 3*$h - 3 + scalar(@R) + 12; return components($h, Tseries($N), "T", @R); }
sub mv_comp { my ($h, @R) = @_; return components($h, undef, "none", @R); }
sub xi { my ($h, $tau, @R) = @_; return contract(omega_class($h, Q(1), $tau, -1 - $tau), xi_comp($h, @R)); }
sub mv { my ($h, $tau, @R) = @_; return contract(omega_class($h, Q(1), $tau, -1 - $tau), mv_comp($h, @R)); }
sub xi_cls { my ($h, $cls, @R) = @_; return contract($cls, xi_comp($h, @R)); }
# s-interpolated dilaton leaves s*psi^2/(1-psi) (= -s T): Q(s) of the brief, components for a given rational s
sub q_comp { my ($h, $s, @R) = @_; my $N = 3*$h - 3 + scalar(@R) + 12; my $T = ser_T(Q(1), $N); my @D = map { -$s * $_ } @$T;
  return components($h, \@D, "s=$s", @R); }
1;
