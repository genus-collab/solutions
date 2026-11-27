#!/usr/bin/perl
# Referee check (TR1, round 2).  Compare the author's SAVED TR coefficients (runs/save_tau<t>_g7n1*.txt, written by
# run_xi.pl with SAVE=..., coefficients of prod d xi_k in omega_{g,n} for y = -1/w) with the Hodge side computed from the
# Route E REFEREE engine tables (E_formulas/referee/ref_engine_h<h>_r<r>.txt; Mumford GRR + DVV, no TR), which the
# author did not use (the author's report says "(7,1): no Hodge table").
#   P2:  C[k] = sigma^(1-g) X[k],   X[k] = Sum_c (-1/sigma)^c [j^k] Xi^(c)_g.
# Usage: perl cmp_saved_vs_refengine.pl TAU SAVEFILE
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my ($tau, $save) = @ARGV; die "usage: TAU SAVEFILE\n" unless defined $save;
my $T = R($tau); my $sig = $T * ($T + 1);
use File::Basename (); use File::Spec ();
my $REF = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../../E_formulas/referee";
my $ENG = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../../E_formulas";
sub load_tab { my ($f) = @_; open my $fh, '<', $f or return undef; my %comp;
  while (my $l = <$fh>) { $l =~ s/\r//g; chomp $l; next unless $l =~ /\S/; my ($h, $r, $c, $body) = split /\s*\|\s*/, $l, 4;
    $comp{$c} //= {}; next unless defined $body && $body =~ /\S/;
    for my $t (split /\s*;\s*/, $body) { my ($k, $v) = split /:/, $t; $comp{$c}{$k} = R($v); } }
  close $fh; return \%comp; }
sub X_of { my ($tab, $key) = @_; my $s = R(0); my $ms = R(-1) / $sig;
  for my $c (keys %$tab) { my $v = $tab->{$c}{$key}; $s += $v * ($ms ** $c) if defined $v; } return $s; }
open my $fh, '<', $save or die "no $save"; my ($tot, $bad) = (0, 0);
while (my $l = <$fh>) { next if $l =~ /^#/; $l =~ s/\r//g; chomp $l; next unless $l =~ /\S/;
  my ($g, $n, $body) = split /\s*\|\s*/, $l, 3; next if $g < 1;
  my %C; for my $t (split /\s*;\s*/, $body) { my ($k, $v) = split /:/, $t; $C{$k} = R($v); }
  my $tab = load_tab("$REF/ref_engine_h${g}_r${n}.txt"); my $src = 'ref_engine';
  unless ($tab) { $tab = load_tab("$ENG/engine_h${g}_r${n}.txt"); $src = 'author_engine'; }
  unless ($tab) { print "($g,$n): no Hodge table\n"; next; }
  # all keys with |k| <= 3g-3+n
  my $dim = 3*$g - 3 + $n; my @keys;
  my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @keys, join(',', @$pre); return; } for my $x (0..$left) { $rec->([@$pre, $x], $left - $x); } };
  $rec->([], $dim);
  my ($cnt, $b) = (0, 0); my @ex;
  for my $k (@keys) { my $tr = $C{$k} // R(0); my $hx = ($sig ** (1 - $g)) * X_of($tab, $k); $cnt++;
    if ($tr != $hx) { $b++; push @ex, "k=$k TR=$tr Hodge=$hx" if @ex < 3; } }
  for my $k (keys %C) { my $s = 0; $s += $_ for split /,/, $k; if ($s > $dim) { $b++; push @ex, "outside $k"; } }
  $tot += $cnt; $bad += $b;
  printf "tau=%s (%d,%d) [%s]: %d coefficients compared, %d mismatches%s\n", $T->bstr, $g, $n, $src, $cnt, $b, (@ex ? "  [" . join(' | ', @ex) . "]" : ''); }
print "TOTAL $tot compared, $bad mismatches\n";
