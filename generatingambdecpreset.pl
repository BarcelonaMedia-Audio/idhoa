# header
#
# you need 3 files:
# a header
# the lf matrix coefficients file (passed as parameter of the script)
# the hf matrix coefficients file (passed as parameter of the script)
use strict;
use warnings;



open(INPUTFILE, "<header.studio.ambdec");
open(OUTPUTFILE, ">studio.out");
while(<INPUTFILE>)
{
    my($line) = $_;
    chomp($line);
    print OUTPUTFILE "$line\n";
}
close(INPUTFILE);



# low frequency
print OUTPUTFILE "/lfmatrix/{ \n";
print OUTPUTFILE "order_gain     1.00000  1.00000  1.00000  1.00000 \n";

my $file_name = shift @ARGV;
open(my $file, '<', $file_name) or die $!;
while(<$file>)
{
    my($line) = $_;
    chomp($line);
print OUTPUTFILE "add_row $line\n";
}
close($file);
print OUTPUTFILE "/} \n \n";



## high frequency
print OUTPUTFILE "/hfmatrix/{ \n";
print OUTPUTFILE "order_gain     2.05000  1.80000  1.06250  0.15000 \n";
my $file_name = shift @ARGV;
open(my $file, '<', $file_name) or die $!;
while(<$file>)
{
    my($line) = $_;
    chomp($line);
    print OUTPUTFILE "add_row $line\n";
    }
close($file);
print OUTPUTFILE "/} \n";
print OUTPUTFILE "\n \n";
print OUTPUTFILE "/end";



close(OUTPUTFILE);
