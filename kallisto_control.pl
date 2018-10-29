#!/usr/bin/perl -w

# p.borrill@bham.ac.uk
#
# Aim of script is to run kallisto on RNA-seq for multiple samples to a common reference to calculate expression levels

#### paths and references:
my $path = '/nbi/Research-Groups/NBI/Cristobal-Uauy/';
my $read_path_triticum = '/nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/raw_data/PKG_ENQ-789_60_samples_data_transfer/';

my $ref = "/nbi/group-data/ifs/NBI/Cristobal-Uauy/PB_AFLF/RefSeqv1.1_masked/IWGSC_v1.1_ALL_20170706_transcripts.fasta";
my $index = "/nbi/group-data/ifs/NBI/Cristobal-Uauy/PB_AFLF/RefSeqv1.1_masked/IWGSC_v1.1_ALL_20170706_transcripts_kallisto0.44.0_index";

# NB make index by kallisto index -i IWGSC_v1.1_ALL_20170706_transcripts_kallisto0.44.0_index IWGSC_v1.1_ALL_20170706_transcripts.fasta

my $output_dir = "$path/PB_AFLF/RefSeqv1.1_masked/control_timecourse/kallisto_results";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
my $control_RNA_seq_paired_list = 'input_for_kallisto.txt';


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$output_dir") or die "couldn't move to output directory";

open (INPUT_FILE, "$control_RNA_seq_paired_list") || die "couldn't open the input file $control_RNA_seq_paired_list!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $dir = $array[0];
my $pair_1_R1 = $array[1];
my $pair_1_R2 = $array[2];
my $pair_2_R1 = $array[3];
my $pair_2_R2 = $array[4];
my $output = $array[5];


chdir("$read_path_triticum/$dir") or die "couldn't move to specific read directory";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel kallisto tasks
#
#SBATCH -p nbi-medium,RG-Cristobal-Uauy # partition (queue) e.g. nbi-medium
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks=1 # number of tasks
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --mem 30000 # memory pool for all cores (MB)
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.out # STDOUT
#SBATCH -e $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.err # STDERR
#SBATCH -J JOBNAME_kallisto
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address
SLURM

 my $tmp_file = "$output_dir/tmp/kallisto.$output";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  $SLURM_header =~ s/JOBNAME/$output/g;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum/$dir\n";

  print SLURM "source kallisto-0.44.0\n";
  print SLURM "source samtools-0.1.19\n";

	print SLURM "kallisto quant -i $index -o $output_dir/$output -t 8 --pseudobam $pair_1_R1 $pair_1_R2 $pair_2_R1 $pair_2_R2 | samtools view -Sb - > $output_dir/$output".".bam\n";
	print SLURM "samtools sort $output_dir/$output/pseudoalignments.bam $output_dir/$output".".sorted\n";
  print SLURM "samtools index $output_dir/$output".".sorted.bam $output_dir/$output".".sorted.bai\n";
  print SLURM "rm $output_dir/$output/pseudoalignments.bam\n";

	close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;

}

	    close(INPUT_FILE);
