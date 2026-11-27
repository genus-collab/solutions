# t_compact_check.pl h r : consistency of the compact multiset files with the full-format files (and of the two weight algorithms
# mA (leaf by leaf) and mAg (grouped, t_lattice_fast.pl)) at a few random rational points.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my ($h, $r) = @ARGV; srand(12345);
my (@F, @M);
open my $f1, '<', "$FindBin::Bin/ref_engine_h${h}_r${r}.txt" or die; while (my $line = <$f1>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
  for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; push @{$F[$k]}, [[split /,/, $m], Q($cf)]; } } close $f1;
open my $f2, '<', "$FindBin::Bin/ref_mset_h${h}_r${r}.txt" or die; while (my $line = <$f2>) { chomp $line; my ($hh, $rr, $k, $t) = split /\s*\|\s*/, $line, 4; my ($A, $cf) = split /:/, $t; push @{$M[$k]}, [[split /,/, $A], Q($cf)]; } close $f2;
my @Fa = (Q(1)); push @Fa, $Fa[-1] * $_ for 1..40;
sub mAg { my ($p, $A) = @_; my %c; $c{"$_"}++ for @$p; my @groups = map { [Q($_), $c{$_}] } keys %c; my %M; $M{$_}++ for @$A;
  my $rec; $rec = sub { my ($g, $rem) = @_; my ($u, $cc) = @{$groups[$g]};
    if ($g == $#groups) { my $n = 0; my $den = Q(1); my $deg = 0; for my $v (keys %$rem) { my $m = $rem->{$v}; next unless $m; $n += $m; $den *= $Fa[$m]; $deg += $v * $m; } return Q(0) unless $n == $cc; return $Fa[$cc] / $den * ($deg ? $u ** $deg : Q(1)); }
    my @vs = grep { $rem->{$_} } keys %$rem; my $tot = Q(0);
    my $ch; $ch = sub { my ($i, $left, $pick) = @_; if ($left == 0) { my $den = Q(1); my $deg = 0; my %r2 = %$rem; for my $v (keys %$pick) { $den *= $Fa[$pick->{$v}]; $deg += $v * $pick->{$v}; $r2{$v} -= $pick->{$v}; }
        my $w = $Fa[$cc] / $den * ($deg ? $u ** $deg : Q(1)); my $rest = $rec->($g + 1, \%r2); $tot += $w * $rest; return; }
      return if $i > $#vs; my $v = $vs[$i]; for my $k (0..($rem->{$v} < $left ? $rem->{$v} : $left)) { my %p2 = %$pick; $p2{$v} = $k if $k; $ch->($i + 1, $left - $k, \%p2); } };
    $ch->(0, $cc, {}); return $tot; };
  return $rec->(0, \%M); }
my ($n, $bad) = (0, 0);
for my $trial (1..3) { my @p = map { Q(int(rand(9)) - 2) / (1 + int(rand(3))) } 1..$r; $p[1] = $p[0] if $trial == 3;   # a repeated value too
  for my $k (0..$h) { my $a = Q(0); for my $t (@{$F[$k] // []}) { my ($e, $cf) = @$t; my $x = $cf->copy; for my $i (0..$r-1) { $x *= $p[$i] ** $e->[$i] if $e->[$i]; } $a += $x; }
    my $b = Q(0); my $c = Q(0); for my $t (@{$M[$k] // []}) { my ($A, $cf) = @$t; $b += $cf * RefHodge::mA(\@p, $A); $c += $cf * mAg(\@p, $A); }
    $n++; $bad++ if $a != $b || $a != $c; } }
printf "genus %d r = %d: full vs compact (mA and grouped weights) at 3 random points x %d components: %d comparisons, %d mismatches\n", $h, $r, $h + 1, $n, $bad;
