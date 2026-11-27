use strict; use warnings; use feature 'current_sub'; use JSON::PP;
# usage: save_r2.pl journal outdir_base label=folder/file.md ...   (writes the result of each labelled agent as markdown; report_text last)
my ($file, $base, @maps) = @ARGV;
my %map = map { my ($l, $f) = split /=/, $_, 2; ($l => $f) } @maps;
open my $fh, '<', $file or die; my (%lab, %res);
while (my $l = <$fh>) { my $o = eval { decode_json($l) } or next; my $t = $o->{type} // '';
  if ($t eq 'started') { $lab{$o->{key}} = $o->{label} } elsif ($t eq 'result') { $res{$lab{$o->{key}} // '?'} = $o->{result} } }
my $d; $d = sub { my ($x, $lev, $out) = @_;
  if (ref $x eq 'HASH') { for my $k (sort keys %$x) { my $v = $x->{$k}; $$out .= ('#' x $lev) . " $k\n\n"; if (ref $v) { __SUB__->($v, $lev + 1, $out) } else { $$out .= "$v\n\n" } } }
  elsif (ref $x eq 'ARRAY') { my $i = 0; for my $v (@$x) { $i++; if (ref $v) { $$out .= ('#' x $lev) . " [$i]\n\n"; __SUB__->($v, $lev + 1, $out) } else { $$out .= "- $v\n" } } $$out .= "\n" }
  else { $$out .= "$x\n\n" } };
for my $label (sort keys %map) {
  my $r = $res{$label}; unless ($r) { print "no result yet for '$label'\n"; next }
  my $path = "$base/$map{$label}"; die "refusing to overwrite $path\n" if -e $path;
  my %rest = %$r; my $rt = delete $rest{report_text};
  my $out = "# Lemma A round 2: $label (saved by the main session from the workflow journal, 2026-10-31)\n\n"
          . "Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.\n\n";
  $d->(\%rest, 2, \$out);
  $out .= "## report_text\n\n" . ($rt // '(none)') . "\n" if defined $rt;
  open my $o, '>:encoding(UTF-8)', $path or die "$path: $!"; print $o $out; close $o;
  printf "%-12s -> %s (%d bytes)\n", $label, $map{$label}, length $out; }
