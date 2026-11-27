#!/usr/bin/perl
# run_w.pl TAU G N [ycase] [e] [src]
#   Referee's independent check: exact EO recursion on C_tau with TRw.pm (w-basis, involution series), then for every
#   (g,n) with g <= G, n <= N: symmetry, pole-order bound, and comparison of the x-expansion at w = 0 (Lagrange inversion)
#   at all lattice points mu in Z_{>=1}^n with |mu| <= (3g-3+n) + n + 2 (a unisolvent set for degree 3g-3+n, plus redundancy):
#     base : [x^mu]F = sigma^(1-g) prod Pi_mu Xi_g(mu)                                (P2)
#     bkmp : [x^mu]F = sigma^(g-1+n) prod Pi_mu M_g(mu),  M_g(mu) = <Omega prod 1/(1-mu_i psi_i)>      (BKMP, y = log(1+w))
#     zhou : (-1)^n times bkmp (y = -log(1+w), Zhou's eq. (6) literally)
#     mix  : y = log(1+w) - e/w:  [x^mu]F = sigma^(g-1+n) (1+e sigma)^(2-2g-n) prod Pi_mu Q_g(-e sigma/(1+e sigma); mu)
#   Hodge data: src = 'q' (default for base: Mathematica q-tables at s = -1 when present, else Route E engine tables),
#   'eng' (Route E author engine), 'ref' (Route E referee engine).  Genus 0 base: prod_{i=1}^{n-3}(J-i); bkmp: J^(n-3).
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use TRw; use HodgeRef;
sub R { Math::BigRat->new(@_) }
my ($tau, $G, $N, $ycase, $e, $src) = @ARGV; $ycase //= 'base'; $e //= 0; $src //= 'q';
die "usage: run_w.pl TAU G N [base|bkmp|zhou|mix] [e] [q|eng|ref]\n" unless defined $N;
my $t0 = time;
my $Pmax = 6*$G - 4 + 2*$N; $Pmax = 2*($G + $N) if $Pmax < 2*($G + $N);
my $L = 2*$Pmax + 8;
TRw::setup(tau => $tau, L => $L, y => $ycase, e => $e);
my $sig = $TRw::sigma; my $T = R($tau); my $eR = R($e);
printf "tau = %s (sigma = %s), y = %s%s, G = %d, N = %d, series order %d; setup %.1f s\n", $T->bstr, $sig->bstr, $ycase,
  ($ycase eq 'mix' ? " (e = " . $eR->bstr . ")" : ''), $G, $N, $L, time - $t0;
my @order = TRw::run_all($G + $N, $G);
printf "recursion: %d correlators, %.1f s\n", scalar(@order), time - $t0;
sub perms { my @x = @_; return ([]) unless @x; my @out; my %seen;
  for my $i (0..$#x) { next if $seen{$x[$i]}++; my @rest = @x; my ($h) = splice(@rest, $i, 1); push @out, [$h, @$_] for perms(@rest); } return @out; }
sub lattice { my ($n, $max) = @_; my @out; my $rec; $rec = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @out, [@$pre]; return; }
    for my $x (1..$left - ($n - @$pre - 1)) { $rec->([@$pre, $x], $left - $x); } }; $rec->([], $max); return @out; }
my ($TOT, $BAD) = (0, 0);
my @gn = sort { my ($ga, $na) = split /,/, $a; my ($gb, $nb) = split /,/, $b; $ga <=> $gb || $na <=> $nb } keys %TRw::W;
for my $key (@gn) { my ($g, $n) = split /,/, $key; my $NC = $ENV{NCMP} // $N; next if $g > $G || $n > $NC;
  my $Wh = $TRw::W{$key}; my ($sym, $pole) = (0, 0);
  for my $k (keys %$Wh) { my @m = split /,/, $k; my $sum = 0; $sum += $_ for @m; $pole++ if $sum > 6*$g - 6 + 4*$n || grep { $_ < 2 } @m;
    for my $p (perms(@m)) { my $pk = join(',', @$p); $sym++ unless exists $Wh->{$pk} && $Wh->{$pk} == $Wh->{$k}; } }
  my $d = 3*$g - 3 + $n; my @pts = lattice($n, $d + $n + 2);
  my ($cnt, $bad, $nz, $info) = (0, 0, 0, ''); my @ex;
  for my $mu (@pts) { my $pi = R(1); $pi *= TRw::Pi_mu($_) for @$mu; my $J = 0; $J += $_ for @$mu;
    my $tr = TRw::xcoef($Wh, @$mu); my $pred;
    if ($ycase eq 'base') { my $Xi;
      if ($g == 0) { $Xi = R(1); $Xi *= ($J - $_) for 1..$n-3; $info = 'genus-0 formula'; }
      else { if ($src eq 'q' && $g >= 2) { $Xi = HodgeRef::eval_data('q', $g, $mu, $tau, -1); $info = 'Mathematica q-table at s=-1' if defined $Xi; }
        if (!defined $Xi) { my $w = $src eq 'ref' ? 'ref' : 'eng'; $Xi = HodgeRef::eval_engine($g, $mu, $sig, $w);
          $Xi = HodgeRef::eval_engine($g, $mu, $sig, 'ref') if !defined $Xi && $w eq 'eng'; $info = "engine ($w)" if defined $Xi; } }
      if (!defined $Xi) { $info = 'no Hodge data'; last; }
      $pred = $sig ** (1 - $g) * $pi * $Xi; }
    elsif ($ycase eq 'bkmp' || $ycase eq 'zhou') { my $M;
      if ($g == 0) { $M = R($J) ** ($n - 3); $info = 'genus-0 formula'; } else { $M = HodgeRef::eval_data('mv', $g, $mu, $tau); $info = 'Mathematica mv-table'; }
      if (!defined $M) { $info = 'no Hodge data'; last; }
      $pred = $sig ** ($g - 1 + $n) * $pi * $M; $pred = -$pred if $ycase eq 'zhou' && $n % 2; }
    elsif ($ycase eq 'mix') { my $s = -$eR * $sig / (1 + $eR * $sig); my $Q = $g >= 2 ? HodgeRef::eval_data('q', $g, $mu, $tau, $s) : undef;
      if (!defined $Q) { $info = 'no Hodge data'; last; } $info = 'Mathematica q-table at s = -e sigma/(1+e sigma)';
      $pred = $sig ** ($g - 1 + $n) * (1 + $eR * $sig) ** (2 - 2*$g - $n) * $pi * $Q; }
    $cnt++; $nz++ unless $pi->is_zero;
    if ($tr != $pred) { $bad++; push @ex, "mu=(" . join(',', @$mu) . ") TR=$tr Hodge=$pred" if @ex < 2; } }
  my $line = sprintf "(%d,%d): %4d w-terms, symmetry viol. %d, pole-order viol. %d", $g, $n, scalar(keys %$Wh), $sym, $pole;
  if ($info ne 'no Hodge data' && $cnt) { $line .= sprintf "; lattice |mu| <= %d: %d points (%d with prod Pi != 0), mismatches %d [%s]", $d + $n + 2, $cnt, $nz, $bad, $info;
    $TOT += $cnt; $BAD += $bad; } else { $line .= "; (no Hodge data)"; }
  $BAD += $sym + $pole; $line .= "  {" . join(' | ', @ex) . "}" if @ex; print "$line\n"; }
printf "TOTAL: %d lattice comparisons, %d failures; time %.1f s\n", $TOT, $BAD, time - $t0;
if ($ENV{SHOW11}) { for my $k (sort keys %{$TRw::W{"1,1"}}) { print "  omega_{1,1} e_$k : ", $TRw::W{"1,1"}{$k}->bstr, "\n"; } }
