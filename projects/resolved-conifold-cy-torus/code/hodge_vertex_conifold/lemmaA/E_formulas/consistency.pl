# consistency.pl h r : internal consistency checks of engine_h<h>_r<r>.txt (needs engine_h<h>_r<r-1>.txt):
#  (a) string and dilaton: Xi(1, R') = J' Xi(R'), Xi(0, R') = (D' - 2h + 2) Xi(R') for all components, at lattice points R';
#  (b) top j-degree part = Witten-Kontsevich numbers <tau_a1..tau_ar>_h  (component k = 0 only, degree 3h-3+r);
#  (c) lambda_g component = b_h (J-1)...(J-N) at lattice points with J <= N.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use WK;
sub R { Math::BigRat->new(shift) }
my ($h, $r) = @ARGV;
sub load { my ($h, $r) = @_; my %c; open my $fh, '<', "$FindBin::Bin/engine_h${h}_r${r}.txt" or die; while (<$fh>) { chomp; next unless /\S/;
    my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $_, 4; $rest //= ''; for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; $c{$k}{$m} = R($cf); } } return \%c; }
my $C = load($h, $r); my $C1 = $r > 1 ? load($h, $r - 1) : undef;
sub at { my ($c, $k, $p) = @_; my $s = R(0); my $P = $c->{$k} or return $s; for my $m (keys %$P) { my @e = split /,/, $m; my $x = $P->{$m}->copy; for my $i (0..$#e) { $x *= R($p->[$i]) ** $e[$i] if $e[$i]; } $s += $x; } return $s; }
# (a)
if ($C1) { my ($n, $bad) = (0, 0);
  my @pts = ([]); for (1..$r-1) { my @nn; for my $t (@pts) { for my $x (0..2*$h) { push @nn, [@$t, $x] } } @pts = @nn; }
  for my $p (@pts) { my $J = 0; $J += $_ for @$p; my $D = $J - ($r - 1);
    for my $k (0..$h) { my $x = at($C1, $k, $p); $n += 2;
      $bad++ if at($C, $k, [1, @$p]) != $J * $x; $bad++ if at($C, $k, [0, @$p]) != ($D - 2*$h + 2) * $x; } }
  print "genus $h r=$r: string/dilaton relations: $n checked, $bad failures\n"; }
# (b)
{ my ($n, $bad) = (0, 0); my $top = 3*$h - 3 + $r;
  for my $m (keys %{$C->{0}}) { my @e = split /,/, $m; my $s = 0; $s += $_ for @e; next unless $s == $top; $n++; $bad++ if $C->{0}{$m} != WK::wk($h, @e); }
  for my $k (1..$h) { for my $m (keys %{$C->{$k}}) { my @e = split /,/, $m; my $s = 0; $s += $_ for @e; $bad++ if $s > $top - $k; } }
  print "genus $h r=$r: top-degree = KW numbers: $n coefficients, $bad failures (incl. degree bound deg Xi^(k) <= 3h-3+r-k)\n"; }
