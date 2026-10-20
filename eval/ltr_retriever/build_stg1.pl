# Rebuild $index.ltrTE.stg1 the way cleanup.pl does, from chunked TRF -ngs output.
use strict; use warnings;
my ($fa, $out, @trf) = @ARGV;
my %tandem;
for my $t (@trf) { open my $T, "<", $t or die; while (<$T>) { $tandem{$1} = 1 if /^\@(.*)$/ } }
open my $F, "<", $fa or die; open my $O, ">", $out or die; open my $L, ">", "$out.cleanup" or die;
local $/ = "\n>"; my ($kept,$nn,$tr)=(0,0,0);
while (<$F>) { chomp; s/^>//; my ($id,$seq)=split /\n/, $_, 2; $seq =~ s/\s+//g;
  my $n = () = $seq =~ /n/gi; my $len = length $seq;
  if ($n/$len >= 0.8 or $n >= 10) { print $L "$id\tN\n"; $nn++; next }
  if ($tandem{$id}) { print $L "$id\ttandem sequence\n"; $tr++; next }
  print $O ">$id\n$seq\n"; $kept++ }
print STDERR "kept $kept, dropped N $nn, dropped tandem $tr, trf-flagged ids ", scalar(keys %tandem), "\n";
