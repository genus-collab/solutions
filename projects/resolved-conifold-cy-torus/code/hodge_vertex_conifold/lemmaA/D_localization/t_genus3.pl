# t_genus3.pl -- genus-3 irreducible points of Lemma A and controls, exact, at two values of tau
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
my $t0 = time;
for my $R ([2],[3],[4],[5],[2,2],[2,3],[2,2,2],[3,3],[2,4],[2,2,3]) {
  my @v = map { xi(3, Q($_), @$R) } ("2/7", "-3/5");
  my $D = 0; $D += $_ - 1 for @$R;
  printf "h=3 R=(%s) D=%d  Xi(tau=2/7)=%s  Xi(tau=-3/5)=%s  [%d s]\n", join(",", @$R), $D, $v[0], $v[1], time - $t0;
}
