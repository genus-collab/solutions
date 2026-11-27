use strict; use warnings; use feature 'current_sub'; use JSON::PP;
# usage: dumpres.pl journal lineno title > out.md   (dumps the "result" object of that journal line as markdown)
my ($file, $n, $title) = @ARGV;
open my $fh, '<', $file or die; my $l; while (<$fh>) { if ($. == $n) { $l = $_; last } } close $fh;
my $o = decode_json($l); my $r = $o->{result};
binmode STDOUT, ':encoding(UTF-8)';
print "# $title\n\n";
my $d; $d = sub { my ($x, $lev) = @_;
  if (ref $x eq 'HASH') { for my $k (sort keys %$x) { my $v = $x->{$k}; print '#' x $lev, " $k\n\n"; if (ref $v) { __SUB__->($v, $lev+1) } else { print "$v\n\n" } } }
  elsif (ref $x eq 'ARRAY') { my $i = 0; for my $v (@$x) { $i++; if (ref $v) { print '#' x $lev, " [$i]\n\n"; __SUB__->($v, $lev+1) } else { print "- $v\n" } } print "\n"; }
  else { print "$x\n\n" } };
$d->($r, 2);
