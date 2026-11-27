#!/usr/bin/perl
# a3_nullity.pl -- how much do the 2-insertion relations (Marino-Vafa cut-and-join at the dilaton point, report A3.3) determine?
# Genus g, range D <= 2g-3.  Unknowns: sigma^g Xi_g(T) = sum_{k=1}^{g} c_{T,k} sigma^k for every multiset T of parts >= 2 in the range
# (the sigma^0 coefficient is the lambda_g part, zero in the range).  Inside the range, by induction on the genus (and on the level for the
# products), the relation for every R with D(R) <= 2g-4 is homogeneous:
#     Xi(R+2) - sum_i gam_{j_i} Xi(R_i^+) + sum_{i<k} gam_{j_i j_k} Xi(R_ik) = 0   (identically in tau).
# The relations are imposed at many random tau (mod p), which is equivalent to the polynomial identity.  Prints the nullity
# (free parameters left) with optional extra conditions: H = Hurwitz components c_{T,g} = 0, L = lambda_1 components c_{T,g-1} = 0.
# Usage: perl a3_nullity.pl g [H] [L]
use strict; use warnings;
my ($g, @opt) = @ARGV; my %opt = map { $_ => 1 } @opt;
my $p = 2147483647;
sub md { my $x = shift; $x %= $p; $x += $p if $x < 0; return $x; }
sub mmul { return md($_[0] * $_[1]); }
sub minv { my $a = md(shift); my ($t, $nt, $r, $nr) = (0, 1, $p, $a); while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q*$nt); ($r, $nr) = ($nr, $r - $q*$nr); } die "not invertible" if $r != 1; return md($t); }
# partitions of D into parts >= 1, mapped to T with parts = part+1 (so D(T) = D)
sub parts { my ($n, $max) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
my @T; my %Tidx;
for my $D (1..2*$g-3) { for my $q (parts($D, $D)) { my @t = map { $_ + 1 } @$q; my $key = join(',', @t); $Tidx{$key} = scalar @T; push @T, \@t; } }
my $nT = scalar @T; my $nU = $nT * $g;
sub key { return join(',', sort { $b <=> $a } @_); }
my @Rs = ([]); for my $D (1..2*$g-4) { for my $q (parts($D, $D)) { push @Rs, [ map { $_ + 1 } @$q ]; } }
srand(12345);
my @rows;
for my $sample (1..(6*$g + 40)) {
  my $tau = 2 + int(rand($p - 10)); my $sig = mmul($tau, md($tau + 1)); my $sp = md(2*$tau + 1);
  my %Pi; my $Pi = sub { my $j = shift; return $Pi{$j} //= do { my $v = 1; for my $b (1..$j-1) { $v = mmul($v, md($j*$tau + $b)); } my $f = 1; $f = mmul($f, $_) for 1..$j-1; mmul($v, minv($f)); }; };
  my @sigp = (1); push @sigp, mmul($sigp[-1], $sig) for 1..$g;
  for my $R (@Rs) { my %row; my $add = sub { my ($t, $c) = @_; my $i = $Tidx{key(@$t)}; die "not in range: @$t" unless defined $i; for my $k (1..$g) { my $col = $i*$g + $k - 1; $row{$col} = md(($row{$col} // 0) + mmul($c, $sigp[$k])); } };
    my @R = @$R; my $r = scalar @R;
    $add->([@R, 2], 1);
    for my $i (0..$r-1) { my $j = $R[$i]; my @Rp = @R; $Rp[$i] = $j + 1; my $gam = mmul(mmul($j+1, $Pi->($j+1)), minv(mmul($Pi->($j), $sp))); $add->(\@Rp, md(-$gam)); }
    for my $i (0..$r-1) { for my $k ($i+1..$r-1) { my $jj = $R[$i] + $R[$k]; my @Rm = ((map { $R[$_] } grep { $_ != $i && $_ != $k } 0..$r-1), $jj);
        my $gam = mmul(mmul($jj, $Pi->($jj)), minv(mmul(mmul($Pi->($R[$i]), $Pi->($R[$k])), $sp))); $add->(\@Rm, $gam); } }
    push @rows, \%row; } }
# extra conditions
if ($opt{H}) { for my $i (0..$nT-1) { push @rows, { $i*$g + $g - 1 => 1 }; } }
# Z=t1.t2 : impose full vanishing of T = (t1,t2,..);  H=t : Hurwitz component of (t) only; L=t : lambda_1 component of (t) only
for my $o (@opt) { if ($o =~ /^([ZHL])=([\d.]+)$/) { my ($ty, $tl) = ($1, $2); my $i = $Tidx{key(split /\./, $tl)}; die "no T $tl" unless defined $i;
    if ($ty eq "Z") { push @rows, { $i*$g + $_ - 1 => 1 } for 1..$g; } elsif ($ty eq "H") { push @rows, { $i*$g + $g - 1 => 1 }; } else { push @rows, { $i*$g + $g - 2 => 1 }; } } }
if ($opt{L} && $g >= 2) { for my $i (0..$nT-1) { push @rows, { $i*$g + $g - 2 => 1 }; } }
# rank mod p (sparse-ish Gaussian elimination)
my @piv; my $rank = 0; my %pivrow;
for my $row (@rows) { my %v = %$row;
  for my $col (sort { $a <=> $b } keys %v) { next unless $v{$col}; if (exists $pivrow{$col}) { my $pr = $pivrow{$col}; my $f = $v{$col}; for my $c2 (keys %$pr) { $v{$c2} = md(($v{$c2} // 0) - mmul($f, $pr->{$c2})); } } }
  my @nz = sort { $a <=> $b } grep { $v{$_} } keys %v; next unless @nz; my $c0 = $nz[0]; my $inv = minv($v{$c0});
  my %nr = map { $_ => mmul($v{$_}, $inv) } @nz;
  # reduce existing pivot rows by the new one (keep reduced form)
  for my $pc (keys %pivrow) { my $pr = $pivrow{$pc}; if ($pr->{$c0}) { my $f = $pr->{$c0}; for my $c2 (keys %nr) { $pr->{$c2} = md(($pr->{$c2} // 0) - mmul($f, $nr{$c2})); } delete $pr->{$c0}; } }
  $pivrow{$c0} = \%nr; $rank++; }
my $null = $nU - $rank;
print "genus $g: range D <= ", 2*$g-3, ", ", $nT, " multisets, ", $nU, " unknown coefficients, ", scalar(@Rs), " relations; extra: ", (join(',', @opt) || 'none'), "; rank $rank, nullity $null\n";
# which unknowns are forced to vanish?
my @free; for my $i (0..$nT-1) { my $forced = 1; for my $k (1..$g) { my $col = $i*$g + $k - 1; if (!exists $pivrow{$col}) { $forced = 0; } else { my @o = grep { $_ != $col && $pivrow{$col}{$_} } keys %{$pivrow{$col}}; $forced = 0 if @o; } } push @free, "(" . join(',', @{$T[$i]}) . ")" unless $forced; }
print "   not forced to vanish: ", (@free ? join(' ', @free) : 'none'), "\n";
