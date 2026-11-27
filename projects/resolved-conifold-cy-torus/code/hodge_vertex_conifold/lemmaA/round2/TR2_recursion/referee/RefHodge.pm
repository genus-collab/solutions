package RefHodge;
# Referee: reads Route E's exact Hodge component tables E_formulas/engine_h<h>_r<r>.txt (no TR involved in them):
#   line "h | r | k | e_1,..,e_r:coef ; ..."  means  Xi^(k)_h(j_1..j_r) = Sum coef prod j_i^(e_i),
#   Xi_h(R; tau) = Sum_k (-1/sigma)^k Xi^(k)_h(R)   (brief s.1, proof_status 1.1(e)).
# Evaluates Xi_h(R; tau) mod P.
use strict; use warnings; use Math::BigInt; use Math::BigRat;
use RefMod;
use File::Basename (); use File::Spec ();
my $DIR = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../../E_formulas";
my %cache;
sub load { my ($h,$r) = @_; my $key = "$h,$r"; return $cache{$key} if $cache{$key};
  my $f = "$DIR/engine_h${h}_r${r}.txt"; return undef unless -e $f;
  open(my $fh, '<', $f) or die "cannot open $f\n"; my @comp;
  while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my @F = split /\s*\|\s*/, $line, 4; next unless @F == 4;
    my ($hh,$rr,$k,$body) = @F; next unless $hh == $h && $rr == $r;
    for my $t (split /\s*;\s*/, $body) { next unless $t =~ /\S/; my ($ex,$c) = split /:/, $t; my @e = split /,/, $ex; die "bad arity in $f\n" unless @e == $r;
      push @{$comp[$k]}, [\@e, Math::BigRat->new($c)]; } }
  close $fh; $cache{$key} = \@comp; \@comp }
sub bigmod { my ($q) = @_; my $p = $RefMod::P; my $n = $q->numerator->copy->bmod($p)->numify; my $d = $q->denominator->copy->bmod($p)->numify; ratm($n,$d) }
# Xi_h(R) at tau = tn/td, mod P
sub xi { my ($h,$tn,$td,@R) = @_; my $comp = load($h,scalar @R) or return undef;
  my $tau = ratm($tn,$td); my $sig = mulm($tau, addm($tau,1)); my $ms = negm(invm($sig)); my $tot = 0;
  for my $k (0..$#$comp) { next unless $comp->[$k]; my $s = 0;
    for my $t (@{$comp->[$k]}) { my ($e,$c) = @$t; my $v = bigmod($c); for my $i (0..$#R) { $v = mulm($v, powm($R[$i],$e->[$i])); } $s = addm($s,$v); }
    $tot = addm($tot, mulm($s, powm($ms,$k))); }
  $tot }
1;
