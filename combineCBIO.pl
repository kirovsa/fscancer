my $cbfs=`find . -name 'data_mutations*'`;
open(CNT,">cnt");
my @cbf=split(/\n/,$cbfs);
foreach my $proj (@cbf) {
my ($dor,$fm,$projname,$file)=split(/\//,$proj);
my ($study,$center,$r)=split("\_",$projname);
my $uid=$study.$center;
#Clean
next if ($seen{$uid}++);
#next if ($proj=~/ccle/i);
#print "$projname\t$study\t$center\t$proj\n";
open(F,$proj);
my $h=<F>;
my @h=split(/\t/,$h);
my $gene;
my $vtype;
for my $i (0..$#h) {
	if ($h[$i] eq "Hugo_Symbol") {
		$gene=$i;
	}
	if ($h[$i] eq "Variant_Type") {
		$vtype=$i;
	}
}
unless ($vtype) { warn $projname; next; }
while (my $b=<F>) {
	my @d=split(/\t/,$b);
	print "$projname\t".$d[$gene]."\t".$d[$vtype]."\n";
	$cnt{$d[$gene]}++;
}
}
foreach my $genek (keys %cnt) {
	print CNT $genek."\t".$cnt{$genek}."\n";
}
