package R2tab;
# Second referee: reader for Route E's exact component tables E_formulas/engine_h<h>_r<r>.txt (read-only).
# Line format: h | r | k | e_1,..,e_r:coef ; ...   meaning Xi^(k)_h(j_1..j_r) = sum coef prod j_i^e_i  (tau-independent components);
# Xi_h(T) = sum_k (-1/sigma)^k Xi^(k)_h(T).  Genus 0: Xi_0(T) = Gamma(J)/Gamma(J-r+3) (Route A Lemma 3.1).
use strict; use warnings; use Math::BigRat; use FindBin;
our $EDIR = "$FindBin::Bin/../../../E_formulas";
my %tab;
sub load { my ($h, $r) = @_; my $key = "$h|$r"; return $tab{$key} if exists $tab{$key};
  my $f = "$EDIR/engine_h${h}_r${r}.txt"; unless (-e $f) { $tab{$key} = undef; return undef; }
  open(my $fh, '<', $f) or die "open $f"; my %rows;
  while (my $line = <$fh>) { $line =~ s/[\r\n]//g; next unless $line =~ /\S/; my @p = split /\s*\|\s*/, $line, 4; next unless @p == 4;
    my ($hh, $rr, $k, $terms) = @p; next unless $hh == $h && $rr == $r;
    for my $t (split /\s*;\s*/, $terms) { next unless $t =~ /:/; my ($e, $c) = split /:/, $t; push @{$rows{$k}}, [[split /,/, $e], Math::BigRat->new($c)]; } }
  close $fh; $tab{$key} = \%rows; \%rows }
sub available { my ($h, $r) = @_; defined load($h, $r) ? 1 : 0 }
# components (exact) at the parts T (order as given)
sub comps { my ($h, $T) = @_; my $r = scalar @$T; my $rows = load($h, $r); die "R2tab: no table for h=$h r=$r" unless $rows;
  my @c; for my $k (sort { $a <=> $b } keys %$rows) { my $s = Math::BigRat->new(0);
    for my $row (@{$rows->{$k}}) { my ($e, $c) = @$row; my $m = $c->copy; for my $i (0..$r-1) { $m *= Math::BigRat->new($T->[$i]) ** $e->[$i] if $e->[$i]; } $s += $m; }
    $c[$k] = $s; }
  [ map { $_ // Math::BigRat->new(0) } @c ] }
sub xi0 { my ($T) = @_; my $r = scalar @$T; my $J = 0; $J += $_ for @$T; die "xi0: empty" unless $r;
  return Math::BigRat->new(1) / ($J * ($J + 1)) if $r == 1; return Math::BigRat->new(1) / $J if $r == 2;
  my $v = Math::BigRat->new(1); $v *= ($J - $_) for 1..$r-3; $v }
# exact value at rational tau
sub xi_at { my ($h, $T, $tau) = @_; my @T = sort { $a <=> $b } @$T; return xi0(\@T) if $h == 0;
  my $c = comps($h, \@T); my $sig = $tau * ($tau + 1); my $v = Math::BigRat->new(-1) / $sig; my $s = Math::BigRat->new(0); my $q = Math::BigRat->new(1);
  for my $k (0..$#$c) { $s += $c->[$k] * $q; $q *= $v; } $s }
1;
