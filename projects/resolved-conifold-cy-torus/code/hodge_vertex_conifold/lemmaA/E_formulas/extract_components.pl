# extract_components.pl -- read the symbolic Hodge data (hvc_multiprobe_out.txt: genus 2, 3 r = 2; hvc_multiprobe_g4_out.txt:
# genus 4 r = 1), convert to c = 1 and split into the tau-independent components Xi^(k) (Comp.pm).
# Output: components.txt, one line per (h, r, k):  h | r | k | a1,a2:coef ; ...   (coefficient of j1^a1 j2^a2)
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Poly; use MParse; use Comp;
my $HVC = "$FindBin::Bin/../..";
my @jobs;
open my $f, '<', "$HVC/hvc_multiprobe_out.txt" or die;
while (<$f>) { if (/^genus (\d) symbolic r=2: (.*)$/) { push @jobs, [$1, 2, $2]; } }
close $f;
open $f, '<', "$HVC/hvc_multiprobe_g4_out.txt" or die;
while (<$f>) { if (/^genus 4 r=1 symbolic: (.*)$/) { push @jobs, [4, 1, $1]; } }
close $f;
open my $o, '>', "$FindBin::Bin/components.txt" or die;
for my $jb (@jobs) {
  my ($h, $r, $e) = @$jb;
  my ($comp, $ok) = Comp::components($e, $h, $r);
  print "genus $h r=$r: consistency of the 1/sigma-polynomial structure (deg <= h) at extra points: ", ($ok ? "OK" : "FAILED"), "\n";
  for my $k (0..$h) { my $p = $comp->[$k];
    my @terms;
    for my $m (sort keys %$p) { my $ex = Poly::mono_parse($m); my @a = map { $ex->{"jj$_"} // 0 } 1..$r; push @terms, join(',', @a) . ":" . $p->{$m}; }
    print $o "$h | $r | $k | ", join(' ; ', @terms), "\n"; }
}
close $o;
print "written components.txt\n";
