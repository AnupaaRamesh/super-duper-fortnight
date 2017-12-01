#To create Phylogenetic Trees for the respective input sequences.
#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Bio::AlignIO;
my $in = Bio::AlignIO->new( -file   => "<seq.fasta.aln",
                                   -format => "clustalw" );
my $out = Bio::AlignIO->new(-file   => ">seq.fasta.phy" ,
                             -format => 'phylip');
  my $seq; 
  while($seq = $in->next_aln) { $out->write_aln($seq) }
  
#!/usr/bin/perl
use warnings;
use strict;
use Bio::Phylo::IO;
use Bio::Phylo::Treedrawer;

use Getopt::Long;
use Pod::Usage;
my $in  = '';
my $out = '';
my $usage      = "\n\n$0 [options] \n

Options:
       -$in        name the input file
       -$out       name the output file
       -help             Show this message
\n";
GetOptions(
'in=s'  => \$in,
'out=s' => \$out,
help           => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($in) {
die "\nProvide an input file, -in <string>", $usage;

}
unless ($out) {
die "\nProvide an output file, -out <string>", $usage;

}
open( SVG,    ">", $out ) or die $!;
open( NEWICK, "<", $in )  or die $!;
my $newickString = '';

while (<NEWICK>) {
chomp;
$newickString .= $_;
}

$newickString .= ';';

my $forest = Bio::Phylo::IO->parse(
-format => 'newick',
-string => $newickString
);
my $tree = $forest->first;

my $treedrawer = Bio::Phylo::Treedrawer->new(
-width  => 1920,
-height => 1280,
-shape  => 'CURVY',
-mode   => 'PHYLO',
-format => 'SVG'
);

$treedrawer->set_scale_options(
-width => '100%',
-major => '10%',
-minor => '2%',
-label => 'MYA',
);

$treedrawer->set_tree($tree);
print SVG $treedrawer->draw;
