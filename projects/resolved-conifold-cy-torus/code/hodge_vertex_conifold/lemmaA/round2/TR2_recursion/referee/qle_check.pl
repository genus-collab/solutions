# Referee check of Proposition 5 (local form): with Phi_g = Brk_g(v,v)/(2 dX^2) - v omega_{g,1}(v)/dX,
# Brk_g(v,v') = omega_{g-1,2}(v,v') + Sum_{h=1}^{g-1} omega_{h,1}(v) omega_{g-h,1}(v'),
# the function Phi_g(v) + Phi_g(sigma(v)) is holomorphic at v = tau (sigma = global involution near tau).
# Uses the referee's engine (RefTR.pm, coordinate e = v - tau, s(e) = sigma(v) - tau).  usage: perl qle_check.pl tn td G
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin;
use RefMod; use RefTR;
my ($tn,$td,$G) = @ARGV; my ($p) = prime_list(1); setp($p);
my $C = RefTR::new_curve($tn,$td,6*$G+2*2-4); RefTR::run($C,$G,2);
my $tau = $C->{tau}; my $tau1 = $C->{tau1}; my $NR = $C->{NR}; my $LO = 6*$G+2*2; # Laurent offset
# Laurent series: hash power -> coef
sub ladd { my ($A,$k,$v) = @_; $A->{$k} = addm($A->{$k} // 0, $v); }
my $te = psmul([$tau,1],[$tau1,1],$NR);            # (tau+e)(tau+1+e);  1/X'(tau+e) = -te/e
my $te2 = psmul($te,$te,$NR);
for my $g (2..$G) {
  my %phi;
  # -v omega_{g,1}/dX = -(tau+e) * Sum c_k e^(-k) * (-te/e) = (tau+e) te Sum c_k e^(-k-1)
  my $vte = psmul([$tau,1],$te,$NR);
  for my $key (keys %{$C->{W}{"$g,1"}}) { my $c = $C->{W}{"$g,1"}{$key}; for my $i (0..$NR) { next unless $vte->[$i]; ladd(\%phi, $i-$key-1, mulm($c,$vte->[$i])); } }
  # Brk(v,v)/(2 dX^2) = (1/2) te^2 e^(-2) Brk(e)
  my %brk;
  for my $key (keys %{$C->{W}{($g-1).",2"}}) { my ($a,$b) = split /,/, $key; ladd(\%brk, -$a-$b, $C->{W}{($g-1).",2"}{$key}); }
  for my $h (1..$g-1) { for my $k1 (keys %{$C->{W}{"$h,1"}}) { for my $k2 (keys %{$C->{W}{($g-$h).",1"}}) {
        ladd(\%brk, -$k1-$k2, mulm($C->{W}{"$h,1"}{$k1}, $C->{W}{($g-$h).",1"}{$k2})); } } }
  my $half = invm(2);
  for my $k (keys %brk) { for my $i (0..$NR) { next unless $te2->[$i]; ladd(\%phi, $k-2+$i, mulm(mulm($half,$brk{$k}),$te2->[$i])); } }
  # Phi(sigma v): substitute e -> s(e) = e r(e):  Sum_m phi_m e^m r^m
  my %tot = %phi; my $r = $C->{r}; my $minp = 0; for (keys %phi) { $minp = $_ if $_ < $minp; }
  for my $m (keys %phi) { next unless $phi{$m}; my $rm = pspow($r,$m,$NR); for my $i (0..$NR) { next unless $rm->[$i]; my $pw = $m+$i; last if $pw > 0; ladd(\%tot, $pw, mulm($phi{$m},$rm->[$i])); } }
  my @polar = sort { $a <=> $b } grep { $_ < 0 && $tot{$_} } keys %tot;
  my @phipol = grep { $_ < 0 && $phi{$_} } keys %phi;
  print "tau=$tn/$td g=$g: Phi_g has ".scalar(@phipol)." polar coefficients (lowest power $minp); Phi_g(v)+Phi_g(sigma v): ".(@polar ? "NON-ZERO polar part at powers @polar" : "polar part vanishes")."\n"; }
