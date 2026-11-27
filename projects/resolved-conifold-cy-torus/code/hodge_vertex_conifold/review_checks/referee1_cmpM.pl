# Compare TR Xi on C_1 with the Mathematica symbolic formulas at x = -1/2 (tau = 1).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
use File::Basename (); use File::Spec ();
my $dir = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/..";
sub getexpr { my ($file,$tag)=@_; open my $f, '<', "$dir/$file" or die; while (<$f>){ if (index($_,$tag)>=0){ my $e=substr($_, index($_,$tag)+length($tag)); chomp $e; $e =~ s/\r//g; return $e; } } die "no $tag" }
sub ev { my ($e,$x,@j)=@_; my $s=$e; $s =~ s/jj\[(\d)\]/\$J[$1]/g; $s =~ s/\bx\b/\$X/g; $s =~ s/\^/**/g;
  my @J=(undef, map { R($_) } @j); my $X=R($x); my $v = eval $s; die $@ if $@; $v }
my %F;
$F{"2,2"} = getexpr("hvc_multiprobe_out.txt","genus 2 symbolic r=2: ");
$F{"3,2"} = getexpr("hvc_multiprobe_out.txt","genus 3 symbolic r=2: ");
$F{"4,1"} = getexpr("hvc_multiprobe_g4_out.txt","genus 4 r=1 symbolic: ");
my %tr = ( "2,2" => [ [1,3,"-3/512"],[1,4,"-11/640"],[2,2,"-5/1536"],[2,3,"-23/1920"] ],
           "3,2" => [ [1,5,"-15/4096"],[1,6,"-309/14336"],[2,4,"-9/4096"],[2,5,"-619/43008"],[3,3,"-47/20480"],[3,4,"-943/71680"] ] );
for my $k (sort keys %tr){ for my $t (@{$tr{$k}}){ my ($a,$b,$v)=@$t; my $m = ev($F{$k}, "-1/2", $a, $b);
    print "($k) ($a,$b): TR = $v, Mathematica(x=-1/2) = ", $m->bstr, ", ratio TR/M = ", (R($v)/$m)->bstr, "\n"; } }
# zeros of the Mathematica formulas in range, at x = -1/2, for completeness
print "M genus 4 r=1 at x=-1/2: ", join("  ", map { "($_)=".ev($F{"4,1"},"-1/2",$_)->bstr } 1..10), "\n";
