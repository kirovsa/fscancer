my $cbfs=`find . -name 'data_mutations*' |grep -vi ccle|grep -vi pdx|grep -vi cellline|grep -vi cell_line|grep -vi xenograft|grep -vi test`;
#my $cbfs=`find . -name 'data_mutations*' |grep test`;
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
my $h;
do {
	$h=<F>;
} until ($h!~/^#/);
my @h=split(/\t/,$h);
my $gene;
my $vtype;
my $classification;
my $hgvsp;
for my $i (0..$#h) {
	if ($h[$i] eq "Hugo_Symbol") {
		$gene=$i;
	}
	if ($h[$i] eq "Variant_Type") {
		$vtype=$i;
	}
	if ($h[$i] eq "HGVSp") {
		$hgvsp=$i;
	}
	if ($h[$i] eq "Tumor_Sample_Barcode") {
                $sample=$i;
        }
	if ($h[$i] eq "Consequence") {
                $classification=$i;
        }

}
unless ($hgvsp) { warn $projname; next; }
my ($fsstart,$fslen);
while (my $b=<F>) {
	next if ($b=~/^\#/);
	my @d=split(/\t/,$b);
	my ($fsstart,$fslen);
	$d[$vtype]="SNP" if ($d[$classification]=~/inframe/i);
	if ($d[$classification] !~/frameshift/i) {
		$fslen=0;
	}
	else {
		($fsstart,$fslen)=$d[$hgvsp]=~/(\d+)/g;
	}
	unless ($fslen) { $fslen=0; }
	print "$projname\t".$d[$gene]."\t".$d[$sample]."\t".$d[$vtype]."\t".$d[$hgvsp]."\t".$fsstart."\t".$fslen."\n";
	$cnt{$d[$gene]}++;
}
}
foreach my $genek (keys %cnt) {
	print CNT $genek."\t".$cnt{$genek}."\n";
}
