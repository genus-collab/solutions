# Referee check: independent EO recursion on C_tau (RefTR.pm) against Route E's exact Hodge tables through P2,
#   Xi_h(mu) = sigma^(h-1) [x^mu] F_{h,r} / prod C(mu_i(tau+1)-1, mu_i-1),
# for every multiset mu (parts >= 1) with r parts, |mu| <= 2h-2+r+X, all (h,r) with h <= G, r <= N and a table present.
# Also tests Lemma D directly in TR: [x^mu] F_{h,r} = 0 for |mu| <= 2h-3+r.  Mod NP primes.
# usage: perl cmp_hodge.pl tn td G N [X] [NP] [GMIN]
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use RefMod; use RefTR; use RefHodge;
my ($tn,$td,$G,$N,$X,$NP,$GMIN) = @ARGV; $X //= 2; $NP //= 2; $GMIN //= 1;
my @P = prime_list($NP); my %res; my $t0 = time;
sub parts { my ($r,$tot,$min) = @_; return ([]) if $r == 0 && $tot == 0; return () if $r == 0; my @o;
  for my $a ($min..$tot) { last if $a*$r > $tot; for my $rest (parts($r-1,$tot-$a,$a)) { push @o, [$a,@$rest]; } } @o }
my ($ncmp,$nbad,$nzero,$nviol,$nskip) = (0,0,0,0,0);
for my $p (@P) { setp($p);
  my $C = RefTR::new_curve($tn,$td,6*$G+2*$N-4); RefTR::run($C,$G,$N);
  for my $g ($GMIN..$G) { for my $r (1..$N) { next unless 2*$g-2+$r > 0; next if $g == 0;
      for my $tot ($r..2*$g-2+$r+$X) { for my $mu (parts($r,$tot,1)) {
          my $xc = RefTR::xcoef($C,$g,@$mu); my $xi = RefTR::xi_p2($C,$g,@$mu); my $hx = RefHodge::xi($g,$tn,$td,@$mu);
          push @{$res{"$g;@$mu"}}, [$xc, $xi, $hx]; } } } } }
for my $key (sort { my ($ga)=split /;/,$a; my ($gb)=split /;/,$b; $ga <=> $gb || $a cmp $b } keys %res) {
  my ($g,$mus) = split /;/, $key; my @mu = split / /, $mus; my $r = @mu; my $tot = 0; $tot += $_ for @mu;
  my @L = @{$res{$key}}; my $allzero = !grep { $_->[0] } @L;
  if ($tot <= 2*$g-3+$r) { if ($allzero) { $nzero++ } else { $nviol++; print "LEMMA D VIOLATION g=$g mu=(@mu)\n"; } }
  my $have = !grep { !defined $_->[2] } @L; my $norm = !grep { !defined $_->[1] } @L;
  if (!$have || !$norm) { $nskip++; next; }
  $ncmp++; my $ok = !grep { $_->[1] != $_->[2] } @L; $nbad++ unless $ok;
  print "g=$g mu=(@mu): ".($ok ? "agree" : "MISMATCH tr=".join("/",map {$_->[1]} @L)." hodge=".join("/",map {$_->[2]} @L)).(($L[0][2]//1) == 0 ? " (zero)" : "")."\n" if !$ok || $ENV{VERBOSE}; }
print "tau=$tn/$td G=$G N=$N X=$X primes=".scalar(@P).": $ncmp TR-vs-Hodge(P2) comparisons, $nbad mismatches, $nskip skipped (no table / zero normalisation); Lemma D: $nzero zeros in range, $nviol violations; ".(time-$t0)."s\n";
