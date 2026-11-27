# xi_engine_multi.pl h r1 r2 ... -- as xi_engine.pl, several numbers of leaves in one run (shared memo tables).
# Writes engine_h<h>_r<r>.txt for each r.
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Hodge; use Time::HiRes qw(time);
my ($h, @rs) = @ARGV;
my $t0 = time;
Hodge::setup($h);
for my $r (@rs) {
  my $N = 3*$h - 3 + $r; my @tuples;
  my $gen; $gen = sub { my ($pre, $left, $max, $k) = @_; if ($k == 0) { push @tuples, [@$pre]; return; }
    for my $x (reverse 0..($max < $left ? $max : $left)) { $gen->([@$pre, $x], $left - $x, $x, $k - 1); } };
  $gen->([], $N, $N, $r);
  my %poly;
  for my $a (@tuples) { my $v = Hodge::W([], $a);
    for my $p (uniqperms(@$a)) { my $mk = join(',', @$p);
      for my $k (0..$#$v) { next if $v->[$k] == 0; my $key = "$k|$mk"; $poly{$key} = exists $poly{$key} ? $poly{$key} + $v->[$k] : $v->[$k]->copy; } } }
  open my $o, '>', "$FindBin::Bin/engine_h${h}_r${r}.txt" or die;
  for my $k (0..$h) { my @terms;
    for my $key (sort keys %poly) { my ($kk, $mk) = split /\|/, $key; next unless $kk == $k; next if $poly{$key} == 0; push @terms, "$mk:$poly{$key}"; }
    print $o "$h | $r | $k | ", join(' ; ', @terms), "\n"; }
  close $o;
  printf "genus %d r=%d: %d tuples, cumulative time %.1f s\n", $h, $r, scalar(@tuples), time - $t0;
}
sub uniqperms { my @x = @_; my %s; my @out;
  my $rec; $rec = sub { my ($pre, $rest) = @_; if (!@$rest) { my $k = join(',', @$pre); push @out, [@$pre] unless $s{$k}++; return; }
    my %used; for my $i (0..$#$rest) { next if $used{$rest->[$i]}++; my @r2 = @$rest; splice(@r2, $i, 1); $rec->([@$pre, $rest->[$i]], \@r2); } };
  $rec->([], \@x); return @out; }
