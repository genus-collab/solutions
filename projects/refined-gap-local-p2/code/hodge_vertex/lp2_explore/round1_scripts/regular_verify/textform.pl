# textual check of the symbolic stored entries: denominator power of (e1+e2+3L), max L power elsewhere, Sqrt3 placement
# also writes the entries (with a,b,L,R substituted) one per line to ARGV[1]
use strict; local $/; open my $f,'<',$ARGV[0] or die; my $s=<$f>; $s=~s/\s+//g;
my $bs=chr(92); my $e1=$bs."[Epsilon][1]"; my $e2=$bs."[Epsilon][2]";
$s=~s/\Q$e1\E/a/g; $s=~s/\Q$e2\E/b/g; $s=~s/LLL\[0\]\[y\]/L/g; $s=~s/Sqrt\[3\]/R/g;
$s=substr($s,1,-1); my $d=0; my @c; my $cur="";
for my $ch (split //,$s){ $d++ if $ch=~/[\(\[\{]/; $d-- if $ch=~/[\)\]\}]/; if($ch eq ',' && $d==0){push @c,$cur;$cur=""} else {$cur.=$ch} } push @c,$cur;
print "entries: ",scalar(@c),"\n";
for my $k (0..$#c){ my $e=$c[$k]; my @m; while($e=~/\(a\+b\+3\*L\)(?:\^(\d+))?/g){push @m, defined $1?$1:1}
 (my $t=$e)=~s/\(a\+b\+3\*L\)(\^\d+)?//g; my $mx=0; while($t=~/L(?:\^(\d+))?/g){my $v=defined $1?$1:1; $mx=$v if $v>$mx}
 my $nR=()=$e=~/R/g; my %sy=map{$_=>1}($e=~/[A-Za-z]/g);
 print "k=$k denom (E+3L) powers: @m ; max L power elsewhere: $mx ; Sqrt3 occurrences: $nR ; symbols: ",join("",sort keys %sy),"\n"; }
open my $o,'>',"$ARGV[1]"; print $o join("\n",@c),"\n";
