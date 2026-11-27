# check_lemmaA.pl : from the engine component files, check
#  (1) Lemma A = M_irr(h): Xi^(k)_h(R) = 0 for all k and every multiset R of integers >= 2 with D(R) <= 2h-3;
#  (2) control points D = 2h-2 (some component non-zero);
#  (3) the lambda_g component Xi^(h)(R) = b_h (J-1)...(J-N), N = 2h-3+r, as a polynomial identity (all lattice points j_i <= 3h);
#  (4) string/dilaton Xi(1,R) = J Xi(R), Xi(0,R) = (D - 2h + 2) Xi(R) at all lattice points with r+1 <= rmax;
#  (5) the genus-4 values of hvc_multiprobe_g4b_out.txt (c-normalisation) at x = 2 and x = 1/3.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Comps; use MParse; use Poly;
sub R { Math::BigRat->new(shift) }
my %BB = (2 => R("1/30"), 3 => R("1/42"), 4 => R("1/30"), 5 => R("5/66"), 6 => R("691/2730")); my %bh; for my $h (2..6) { my $f = R(1); $f *= $_ for 1..2*$h; $bh{$h} = (R(2)**(2*$h-1) - 1) * $BB{$h} / (R(2)**(2*$h-1) * $f); }
my %C;
opendir my $dh, $FindBin::Bin or die; my @files = grep { /^engine_h\d+_r\d+\.txt$/ } readdir $dh; closedir $dh;
for my $f (@files) { next unless $f =~ /engine_h(\d+)_r(\d+)\.txt$/; my ($h, $r) = ($1, $2);
  my $c = Comps::load("$FindBin::Bin/$f"); $C{$h}{$r} = $c; }
sub msets { my ($r, $lo, $maxD) = @_; my @out; my $rec; $rec = sub { my ($pre, $min, $D, $k) = @_; if ($k == 0) { push @out, [@$pre]; return }
    for (my $j = $min; ($j - 1) + $D <= $maxD; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1, $k - 1); } }; $rec->([], $lo, 0, $r); return @out; }
for my $h (sort { $a <=> $b } keys %C) {
  for my $r (sort { $a <=> $b } keys %{$C{$h}}) { my $c = $C{$h}{$r};
    # (1)+(2)
    my ($nz, $bad, $nctl, $ctlnz) = (0, 0, 0, 0);
    for my $R (msets($r, 2, 2*$h - 2)) { my $D = 0; $D += $_ - 1 for @$R;
      my $allzero = 1; for my $k (0..$h) { $allzero = 0 if Comps::comp_at($c, $h, $r, $k, $R) != 0; }
      if ($D <= 2*$h - 3) { $nz++; $bad++ unless $allzero; } else { $nctl++; $ctlnz++ unless $allzero; } }
    # (3)
    my ($n3, $b3) = (0, 0); my $N = 2*$h - 3 + $r;
    my @pts = msets_all($r, $N, $N);   # degree of Xi^(h) is <= N, so J <= N determines it
    for my $p (@pts) { $n3++; my $J = 0; $J += $_ for @$p; my $pred = $bh{$h}; $pred *= ($J - $_) for 1..$N;
      $b3++ if Comps::comp_at($c, $h, $r, $h, $p) != $pred; }
    # (4)
    my ($n4, $b4) = (0, 0);
    if (exists $C{$h}{$r+1} && $r <= 3) { my $c1 = $C{$h}{$r+1};
      for my $p (msets_all($r, 2*$h, 3*$h - 2 + $r)) { my $J = 0; $J += $_ for @$p; my $D = $J - $r;
        for my $k (0..$h) { my $x = Comps::comp_at($c, $h, $r, $k, $p);
          $n4++; $b4++ if Comps::comp_at($c1, $h, $r+1, $k, [1, @$p]) != $J * $x;
          $n4++; $b4++ if Comps::comp_at($c1, $h, $r+1, $k, [0, @$p]) != ($D - 2*$h + 2) * $x; } } }
    printf "genus %d r=%d: Lemma A points %d, non-vanishing %d | control D=2h-2 points %d, non-zero %d | lambda_g formula %d/%d ok | string+dilaton %d/%d ok\n",
      $h, $r, $nz, $bad, $nctl, $ctlnz, $n3 - $b3, $n3, $n4 - $b4, $n4;
  } }
# (5) genus-4 values from hvc_multiprobe_g4b_out.txt
open my $fh, '<', "$FindBin::Bin/../../hvc_multiprobe_g4b_out.txt" or die;
while (<$fh>) { next unless /^genus 4, R = \{([\d, ]+)\}, D = \d+ \([^)]*\): (.*?)\s+\[/; my @R = split /,\s*/, $1; my $val = $2; my $r = scalar @R;
  next unless exists $C{4}{$r};
  for my $x (R(2), R("1/3")) { my $tau = -$x / (1 + $x); my $mine = Comps::xi_at($C{4}, 4, $r, \@R, $tau) / (1 + $x) ** (9 + 2*$r);
    my $theirs = $val eq '0' ? R(0) : do { (my $e = $val) =~ s/^non-zero, //; Poly::const_val(MParse::parse($e, { x => $x })) };
    print "genus 4 R = {@R} at x = $x: engine ", ($mine == $theirs ? "AGREES" : "DISAGREES ($mine vs $theirs)"), " with hvc_multiprobe_g4b\n"; } }
sub msets_all { # non-decreasing r-tuples with entries in 0..jmax and total J <= $Jmax (default: no bound)
  my ($r, $jmax, $Jmax) = @_; $Jmax //= $r * $jmax; my @out; my $rec;
  $rec = sub { my ($pre, $min, $s, $k) = @_; if ($k == 0) { push @out, [@$pre]; return }
    for my $j ($min..$jmax) { last if $s + $j > $Jmax; $rec->([@$pre, $j], $j, $s + $j, $k - 1) } };
  $rec->([], 0, 0, $r); return @out; }
