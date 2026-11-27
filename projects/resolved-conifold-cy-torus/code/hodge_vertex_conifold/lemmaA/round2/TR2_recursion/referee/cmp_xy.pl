# Referee check: the integrated ABDKS formula (RefXY.pm) against the referee's EO recursion (RefTR.pm) and, through P2,
# against Route E's exact Hodge tables, for [x^mu] F_{g,n}, n = 1..N, all mu (parts >= 1) with |mu| <= 2g-2+n+X.  Mod NP primes.
# usage: perl cmp_xy.pl tn td G N [X] [NP] [GMIN]
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use RefMod; use RefTR; use RefXY; use RefHodge;
my ($tn,$td,$G,$N,$X,$NP,$GMIN) = @ARGV; $X //= 2; $NP //= 1; $GMIN //= 1;
my @P = prime_list($NP); my $t0 = time;
sub parts { my ($r,$tot,$min) = @_; return ([]) if $r == 0 && $tot == 0; return () if $r == 0; my @o;
  for my $a ($min..$tot) { last if $a*$r > $tot; for my $rest (parts($r-1,$tot-$a,$a)) { push @o, [$a,@$rest]; } } @o }
my ($nc,$nbad,$nh,$nhbad,$nlz,$nlv) = (0,0,0,0,0,0); my ($ttr,$txy) = (0,0);
for my $p (@P) { setp($p);
  my $t1 = time; my $C = RefTR::new_curve($tn,$td,6*$G+2*$N-4); RefTR::run($C,$G,$N); $ttr += time - $t1;
  my $XY = RefXY::new($tn,$td,$G, 2*(2*$G-2+$N+$X) + 8);
  for my $g ($GMIN..$G) { for my $n (1..$N) { next unless 2*$g-2+$n > 0;
      for my $tot ($n..2*$g-2+$n+$X) { for my $mu (parts($n,$tot,1)) {
          my $t2 = time; my $f = RefXY::xF($XY,$g,@$mu); $txy += time - $t2; my $t = RefTR::xcoef($C,$g,@$mu); $nc++;
          if ($f != $t) { $nbad++; print "MISMATCH formula vs TR: p=$p g=$g mu=(@$mu) formula=$f TR=$t\n"; }
          if ($tot <= 2*$g-3+$n) { if ($f) { $nlv++; print "Lemma D violation (formula) g=$g mu=(@$mu)\n"; } else { $nlz++; } }
          my $hx = RefHodge::xi($g,$tn,$td,@$mu); next unless defined $hx;
          my $v = mulm($f, powm($C->{sig},$g-1)); my $okn = 1;
          for my $m (@$mu) { my $b = binm(subm(mulm($m,$C->{tau1}),1),$m-1); if (!$b) { $okn = 0; last; } $v = mulm($v, invm($b)); }
          next unless $okn; $nh++; if ($v != $hx) { $nhbad++; print "MISMATCH formula vs Hodge(P2): p=$p g=$g mu=(@$mu)\n"; } } } } } }
print "tau=$tn/$td G=$G N=$N X=$X primes=".scalar(@P).": formula vs TR $nc comparisons, $nbad mismatches; formula vs Hodge(P2) $nh comparisons, $nhbad mismatches; Lemma D in formula: $nlz zeros, $nlv violations; time TR ${ttr}s, formula ${txy}s, total ".(time-$t0)."s\n";
