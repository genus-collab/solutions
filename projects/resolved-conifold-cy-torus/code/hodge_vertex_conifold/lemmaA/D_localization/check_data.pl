# check_data.pl -- cross-check of the Perl engine (hodge_lib/xi_lib) against the main session's Mathematica tables
#   mv_g<g>_l<l>.txt : <Omega prod lam_{m_i}>_g in the c-normalisation (w = (c,-x,-1), c = 1+x, leaves 1/(c - m psi))
#   q_g<g>_r<r>.txt  : Q(s; j) = sum_k s^k/k! <Omega prod lam_{j_i} Tt^k>, Tt = psi^2/(c - psi)
# c-normalisation vs c = 1: MV^(c) = c^{-(n + 3h - 3 + n)} MV^(c=1)(tau = -x/c) for n leaves; same for Q (dilaton leaves have degree 0).
use strict; use warnings; use File::Basename qw(dirname);
{ my $d = dirname(__FILE__); require "$d/xi_lib.pl"; require "$d/mparse.pl"; }
my $dd = dirname(__FILE__) . "/../data";
sub readval { my $f = shift; open my $fh, '<', $f or die "$f: $!"; local $/; my $s = <$fh>; close $fh; my @p = split /\s*\|\s*/, $s; return $p[-1]; }
my $x = Q(2); my $c = 1 + $x; my $tau = -$x / $c;
for my $g (1, 2, 3) { for my $l (1, 2) {
  my $e = readval("$dd/mv_g${g}_l$l.txt");
  for my $ms ([2,3], [1,4], [3,5]) { my @m = @$ms[0..$l-1];
    my %v = (x => $x); $v{"m$_"} = Q($m[$_-1]) for 1..$l;
    my $data = meval($e, \%v);
    my $eng = mv($g, $tau, @m) / $c**($l + 3*$g - 3 + $l);
    printf "mv g=%d l=%d m=(%s): data=%s engine=%s %s\n", $g, $l, join(",", @m), $data, $eng, ($data == $eng ? "ok" : "MISMATCH"); } } }
for my $g (2, 3) { for my $r (1, 2) {
  my $e = readval("$dd/q_g${g}_r$r.txt");
  for my $js ([2,3], [3,2]) { my @j = @$js[0..$r-1];
    my %v = (x => $x, s => Q(-1)); $v{"j$_"} = Q($j[$_-1]) for 1..$r;
    my $data = meval($e, \%v);
    my $eng = xi($g, $tau, @j) / $c**($r + 3*$g - 3 + $r);
    printf "Q(s=-1) = Xi, g=%d r=%d j=(%s): data=%s engine=%s %s\n", $g, $r, join(",", @j), $data, $eng, ($data == $eng ? "ok" : "MISMATCH"); } } }
