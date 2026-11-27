# xi_engine.pl h r [outfile] -- exact components Xi^(k)_h(j_1..j_r), k = 0..h, of
#   Xi_h(R; tau) = Int_{Mbar_{h,r}} Omega K prod_i 1/(1 - j_i psi_i) = Sum_k (-1/sigma)^k Xi^(k)_h(R)    (c = 1)
# computed from scratch (Hodge.pm: Mumford's GRR formula, AC kappa classes, Witten-Kontsevich numbers; no TR).
# Output lines:  h | r | k | a1,..,ar:coef ; ...   (coefficient of j1^a1 ... jr^ar in Xi^(k))
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Hodge; use Time::HiRes qw(time);
my ($h, $r, $out) = @ARGV; $out //= "$FindBin::Bin/engine_h${h}_r${r}.txt";
my $t0 = time;
Hodge::setup($h);
my $N = 3*$h - 3 + $r;
# sorted (descending) tuples a with |a| <= N
my @tuples;
sub gen { my ($pre, $left, $max, $k) = @_; if ($k == 0) { push @tuples, [@$pre]; return; }
  for my $x (reverse 0..($max < $left ? $max : $left)) { gen([@$pre, $x], $left - $x, $x, $k - 1); } }
gen([], $N, $N, $r);
my %poly;   # key "k|e1,..,er" -> BigRat
for my $a (@tuples) {
  my $v = Hodge::W([], $a);
  # distinct permutations of a
  my %seen; my @perms = uniqperms(@$a);
  for my $p (@perms) { my $mk = join(',', @$p);
    for my $k (0..$#$v) { next if $v->[$k] == 0; my $key = "$k|$mk"; $poly{$key} = exists $poly{$key} ? $poly{$key} + $v->[$k] : $v->[$k]->copy; } }
}
open my $o, '>', $out or die;
for my $k (0..$h) { my @terms;
  for my $key (sort keys %poly) { my ($kk, $mk) = split /\|/, $key; next unless $kk == $k; next if $poly{$key} == 0; push @terms, "$mk:$poly{$key}"; }
  print $o "$h | $r | $k | ", join(' ; ', @terms), "\n"; }
close $o;
printf "genus %d r=%d: %d tuples, time %.1f s, memo HK %d HK0 %d WK %d W %d T %d\n", $h, $r, scalar(@tuples), time - $t0,
  scalar(keys %Hodge::memo), scalar(keys %Hodge::memo0), scalar(keys %WK::memo), scalar(keys %Hodge::Wmemo), scalar(keys %Hodge::Tmemo);
sub uniqperms { my @x = @_; my %s; my @out;
  my $rec; $rec = sub { my ($pre, $rest) = @_; if (!@$rest) { my $k = join(',', @$pre); push @out, [@$pre] unless $s{$k}++; return; }
    my %used; for my $i (0..$#$rest) { next if $used{$rest->[$i]}++; my @r2 = @$rest; splice(@r2, $i, 1); $rec->([@$pre, $rest->[$i]], \@r2); } };
  $rec->([], \@x); return @out; }
