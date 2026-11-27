# ref_engine_g3.pl -- the author's engine was compared with the Mathematica q-tables only at genus-3 points where Xi = 0 (check_data.pl).
# Here the NON-ZERO genus-3 values printed by the author's runs (t_genus3.out, check_genus3_alltau.out; c = 1 normalisation at tau) are
# compared with the q-tables (c-normalisation, x symbolic): Xi(tau) = c^(3h-3+2r) Q(s = -1; j; x) with x = -tau/(1+tau), c = 1 + x.
# Also: all genus-3 and genus-4 one- and two-leaf irreducible points with D <= 2h-3 vanish in the tables, and the first points with
# D = 2h-2 do not (controls). No engine code is used.
use strict; use warnings; use FindBin; require "$FindBin::Bin/ref_lib.pl";
my $data = "$FindBin::Bin/../../data"; my $dir = "$FindBin::Bin/..";
my %QT; for my $g (3, 4) { for my $r (1, 2) { open my $fh, '<', "$data/q_g${g}_r$r.txt" or die $!; local $/; my $s = <$fh>; close $fh; my @p = split /\s*\|\s*/, $s; $QT{"$g|$r"} = $p[-1]; } }
sub xi_tab { my ($g, $tau, @R) = @_; my $x = -$tau / (1 + $tau); my $c = 1 + $x; my $r = scalar @R;
  my %v = (x => $x, s => -1); $v{"j$_"} = $R[$_-1] for 1..$r; return ev_expr($QT{"$g|$r"}, \%v) * $c**(3*$g - 3 + 2*$r); }
my $bad = 0; my $n = 0; sub ok { my ($nm, $a, $b) = @_; $n++; my $g = ($a == $b); $bad++ unless $g; printf "   %-66s %s\n", $nm, $g ? "ok" : "MISMATCH ($a vs $b)"; }
# t_genus3.out: "h=3 R=(5) D=4  Xi(tau=2/7)=...  Xi(tau=-3/5)=..."
open my $fh, '<', "$dir/t_genus3.out" or die $!;
while (my $l = <$fh>) { next unless $l =~ /R=\(([\d,]+)\).*Xi\(tau=2\/7\)=(\S+)\s+Xi\(tau=-3\/5\)=(\S+)/; my @R = split /,/, $1; next if @R > 2;
  ok("t_genus3.out  R=(@R) tau=2/7", rq($2), xi_tab(3, rq(2,7), @R)); ok("t_genus3.out  R=(@R) tau=-3/5", rq($3), xi_tab(3, rq(-3,5), @R)); }
close $fh;
# check_genus3_alltau.out: values at tau = 1, 2, -1/3, 3/4, 2/7
open $fh, '<', "$dir/check_genus3_alltau.out" or die $!;
my @taus = (rq(1), rq(2), rq(-1,3), rq(3,4), rq(2,7));
while (my $l = <$fh>) { next unless $l =~ /R=\(([\d,]+)\).*?3\/4,2\/7: (.*?) \|/; my @R = split /,/, $1; next if @R > 2; my @v = split ' ', $2;
  for my $i (0..4) { ok("check_genus3_alltau.out R=(@R) tau=$taus[$i]", rq($v[$i]), xi_tab(3, $taus[$i], @R)); } }
close $fh;
# M_irr in the tables (r <= 2) and controls, genus 3 and 4, at three values of tau
for my $g (3, 4) { for my $tau (rq(2,7), rq(-3,5), rq(5,3)) {
  my @in; my @ctl;
  for my $j (2..2*$g) { my $D = $j - 1; push @in, [$j] if $D <= 2*$g - 3; push @ctl, [$j] if $D == 2*$g - 2; }
  for my $j1 (2..2*$g) { for my $j2 ($j1..2*$g) { my $D = $j1 + $j2 - 2; push @in, [$j1, $j2] if $D <= 2*$g - 3; push @ctl, [$j1, $j2] if $D == 2*$g - 2; } }
  my $z = 0; for my $R (@in) { $z++ if xi_tab($g, $tau, @$R) == 0; }
  my $nz = 0; for my $R (@ctl) { $nz++ if xi_tab($g, $tau, @$R) != 0; }
  ok("genus $g tau=$tau: in-range points (r<=2) vanish (" . scalar(@in) . ")", rq($z), rq(scalar @in));
  ok("genus $g tau=$tau: D=2h-2 controls (r<=2) non-zero (" . scalar(@ctl) . ")", rq($nz), rq(scalar @ctl)); } }
print "$n checks: ", ($bad ? "SOME MISMATCH ($bad)" : "ALL OK"), "\n";
