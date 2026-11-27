#!/usr/bin/perl
# (B1) Normalisation check of ELSV against the exact Hodge data.
#  mv_g<g>_l<l>.txt gives P(m) = <Omega prod lam_{m_i}>_g (c-normalisation, CY triple class). Its k = 0 component (sigma -> oo) is
#  A_g(m) = Int Lambda^vee(1) prod 1/(1 - m_i psi_i).  We compare it with the ELSV value computed from the symmetric group:
#  A_g(m) = H_{g,m} |Aut m| / b! * prod m_i!/m_i^m_i,   H_{g,m} from the Frobenius formula + exponential formula (HBconn).
# Usage: perl b1_check.pl [maxg]
use strict; use warnings; use lib '.'; use HB; use HBconn; use Data; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my $maxg = shift // 3;
my ($n, $bad) = (0, 0);
for my $g (1..$maxg) { for my $l (1..3) {
  my $file = "mv_g${g}_l${l}.txt"; next unless -e "../data/$file";
  my $recs = Data::load($file); my $ast = $recs->[0]{ast};
  my $D = 3*$g - 3 + $l;        # total degree of the polynomial in m
  my $maxsum = ($g == 3 && $l >= 2) ? $l + 4 : $l + $D;     # unisolvent simplex sum(m_i - 1) <= D for g <= 2; partial for g = 3
  my @tuples; my $rec; $rec = sub { my ($k, $pre, $sum) = @_; if ($k == $l) { push @tuples, [@$pre] if $sum <= $maxsum; return }
      my $lo = @$pre ? $pre->[-1] : 1; for (my $m = $lo; $sum + $m + ($l - $k - 1) * $m <= $maxsum; $m++) { $rec->($k+1, [@$pre, $m], $sum + $m) } };
  $rec->(0, [], 0);
  for my $m (@tuples) {
    my %env; $env{"m$_"} = BR($m->[$_ - 1]) for 1..$l;
    my $comp = Data::components($ast, $g, $l, \%env); my $Adata = $comp->[0];
    my @R = grep { $_ >= 2 } @$m; my $c = grep { $_ == 1 } @$m;
    my $A = HBconn::A_list(\@R, $g, $c)->[$c];
    $n++; my $ok = ($A == $Adata) ? "ok" : "MISMATCH"; $bad++ if $ok ne "ok";
    printf "g=%d m=(%s): data %s, ELSV/Hurwitz %s %s\n", $g, join(",", @$m), $Adata, $A, $ok } } }
print "B1: $n comparisons, $bad mismatches\n";
