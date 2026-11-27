# p2_engine_ref.pl (referee 2): the dictionary P2 tested against Route E's exact Hodge engine tables, multi-leaf.
#   Hodge side: Xi_h(mu; tau) = sum_k (-1/sigma)^k comp_k(mu), comp_k from E_formulas/engine_h<h>_r<r>.txt (polynomials in the j_i);
#   TR side (RefTR.pm, independent recursion): sigma^(h-1) [x^mu] F_{h,r} / prod_i C(mu_i(tau+1)-1, mu_i-1).
#   All sorted mu with mu_i >= 1 and |mu| <= 2h-2+r+X (the range |mu| <= 2h-3+r, the first layer and X-1 further layers).
# usage: perl p2_engine_ref.pl tn td h r [X] [prime]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefTR; use Math::BigInt;
my ($tn,$td,$h,$r,$X,$p) = @ARGV; $X //= 2; $p //= 2147483647;
my $file = "$FindBin::Bin/../../../E_formulas/engine_h${h}_r${r}.txt"; open my $fh, '<', $file or die "no $file";
my %C;   # C{k}{"e1,..,er"} = coef mod p
while (<$fh>) { chomp; next unless /\S/; my ($hh,$rr,$k,$rest) = split /\s*\|\s*/, $_, 4; next unless $hh == $h && $rr == $r;
  for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m,$cf) = split /:/, $t; my ($nu,$de) = split m{/}, $cf; $de //= 1;
    my $nm = Math::BigInt->new($nu) % $p; my $dm = Math::BigInt->new($de) % $p; $C{$k}{$m} = RefTR::rat($nm->numify, $dm->numify, $p); } }
close $fh; die "no components read" unless %C;
my $t0 = time; my $S = RefTR->new($tn,$td,3*$h-3+$r,$p); $S->run($h,$r);
my $sig = $S->{sig}; my $ms = RefTR::md(-RefTR::inv($sig,$p),$p);
my @tup; my $rec; $rec = sub { my ($pre,$left,$min)=@_; if (@$pre == $r) { push @tup, [@$pre]; return; }
  for my $m ($min..$left) { last if $m*($r-@$pre) > $left; $rec->([@$pre,$m],$left-$m,$m); } }; $rec->([], 2*$h-2+$r+$X, 1);
my ($n,$bad,$skip,$nzr) = (0,0,0,0);
for my $mu (@tup) { my $s = 0; $s += $_ for @$mu;
  my $hod = 0; for my $k (keys %C) { my $v = 0; for my $m (keys %{$C{$k}}) { my @e = split /,/, $m; my $t = $C{$k}{$m};
        for my $i (0..$r-1) { $t = $t * RefTR::pw($mu->[$i],$e[$i],$p) % $p if $e[$i]; } $v = ($v + $t) % $p; }
    $hod = ($hod + $v * RefTR::pw($ms,$k,$p)) % $p; }
  my $tr = $S->xi($h,@$mu); unless (defined $tr) { $skip++; next; }
  $n++; my $inr = $s <= 2*$h-3+$r; $nzr++ if $inr && $hod;
  if ($tr != $hod) { $bad++; print "MISMATCH tau=$tn/$td h=$h mu=(@$mu) ".($inr?"[range]":"")."\n"; } }
print "tau=$tn/$td h=$h r=$r: $n values of mu (|mu| <= ".(2*$h-2+$r+$X)."), $bad mismatches, $skip skipped (vanishing P2 normalisation); ".
  "Hodge side non-zero in the range: $nzr; ".(time-$t0)."s\n";
