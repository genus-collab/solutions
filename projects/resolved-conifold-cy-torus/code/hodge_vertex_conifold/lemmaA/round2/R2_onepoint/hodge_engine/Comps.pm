package Comps;
# load / evaluate component files written by xi_engine.pl ("h | r | k | e1,..,er:coef ; ...")
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
sub load { my $f = shift; my %c; open my $fh, '<', $f or die "cannot open $f";
  while (<$fh>) { chomp; next unless /\S/; my ($h, $r, $k, $rest) = split /\s*\|\s*/, $_, 4; $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; $c{$h}{$r}{$k}{$m} = R($cf); } }
  return \%c; }
# value of component k at the point j (array ref)
sub comp_at { my ($C, $h, $r, $k, $j) = @_; my $s = R(0); my $P = $C->{$h}{$r}{$k} or return $s;
  for my $m (keys %$P) { my @e = split /,/, $m; my $t = $P->{$m}->copy; for my $i (0..$#e) { $t *= R($j->[$i]) ** $e[$i] if $e[$i]; } $s += $t; }
  return $s; }
# Xi_h(j; tau) in the c = 1 normalisation
sub xi_at { my ($C, $h, $r, $j, $tau) = @_; $tau = R($tau) unless ref $tau; my $sig = $tau * (1 + $tau); my $v = R(-1) / $sig;
  my $s = R(0); for my $k (0..$h) { $s += comp_at($C, $h, $r, $k, $j) * $v ** $k; } return $s; }
1;
