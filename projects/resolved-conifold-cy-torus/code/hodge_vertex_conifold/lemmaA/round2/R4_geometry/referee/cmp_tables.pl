# cmp_tables.pl fileA kA fileB kB [label] : compare component k_A of table A with component k_B of table B
# (format "h | r | k | e[,e2]:coef ; ..."), exactly, coefficient by coefficient.  Symmetrises r = 2 tables (j1 <-> j2) before comparing.
use strict; use warnings; use Math::BigRat;
my ($fa, $ka, $fb, $kb, $label) = @ARGV; $label //= '';
sub load { my ($f, $k) = @_; open my $fh, '<', $f or die "no $f"; my %c; my $found = 0;
  while (my $l = <$fh>) { chomp $l; next unless $l =~ /\|/; my ($h, $r, $kk, $rest) = split /\s*\|\s*/, $l, 4; next unless $kk == $k; $found = 1; $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; $m =~ s/\s//g; my $q = Math::BigRat->new($v); die "bad number '$v' in $f" if $q->is_nan(); $c{$m} = $q unless $q == 0; } }
  close $fh; die "component $k not found in $f" unless $found; return \%c; }
my $A = load($fa, $ka); my $B = load($fb, $kb);
my %keys = map { $_ => 1 } (keys %$A, keys %$B); my $bad = 0; my $n = 0;
for my $m (keys %keys) { $n++; my $x = $A->{$m} // Math::BigRat->new(0); my $y = $B->{$m} // Math::BigRat->new(0); $bad++ if $x != $y; }
printf "%-40s %s\n", $label, ($bad ? "$bad of $n coefficients DIFFER" : "equal ($n coefficients)");
