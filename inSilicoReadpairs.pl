#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Carp;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Getopt::Long 'GetOptions';
use Bio::DB::Sam;
use Data::UUID;
use Pod::Usage;

$Carp::Verbose = 1;

use constant MAX_INSERTS => 2_500;
use constant CACHE_X_READS => 500_000;

{
  my $options = option_builder();
  validateInput($options);
  $options->{'uuid'} = uuid();
  $options->{'read_prefix'} = '@'.$options->{'uuid'}.':';
  $options->{'fai'} = setup_fai($options);
  $options->{'refs'} = get_refs($options);
  $options->{'rand_inserts'} = rand_inserts();
  $options->{'rl'} = 100; # readlength
  generate_reads($options);
}

sub generate_reads {
  my ($options) = @_;
  my $ins_arr_size = @{$options->{'rand_inserts'}} - 1;
  my ($start_1, $end_1, $start_2, $end_2, $insert, $seq_1, $seq_2, $read_name, $nn_s, $tmp_seq);
  my $qual_str = 'C' x $options->{'rl'}; # build a qual_string
  my $read_count = 1;
  my $max_nn_s = $options->{'rl'} * 0.2;
  my (@reads_1, @reads_2);
  foreach my $chr(keys %{$options->{'refs'}}) {
    my $since_success = 0;
    foreach my $segment(@{$options->{'refs'}->{$chr}}) {
      if($since_success > 5_000) {
        warn "Failed to generate a unique readpair in the last 5,000 tries for $chr\n";
        last;
      }
      $since_success++;
      my ($ref_start, $ref_end) = @{$segment};
      my $length = $ref_end - $ref_start + 1;
      my $readpairs_for_depth = int (($length * $options->{'d'} / ($options->{'rl'} * 2)) + 1);
      my $chr_seq = uc $options->{'fai'}->fetch($chr.':'.$ref_start.'-'.$ref_end);
      $chr_seq =~ s/[RYKMSWBDHV]/N/gi;
      for(;$readpairs_for_depth > 0;) {
        $start_1 = int(rand($length)) + 1;
        $insert = $options->{'rand_inserts'}->[$ins_arr_size];

        # no point creating pairs off the end of the reference
        next if($start_1 + $insert + $options->{'rl'} > $length);

        $end_1 = $start_1 + ($options->{'rl'} - 1);
        $seq_1 = substr($chr_seq, $start_1-1, $options->{'rl'}); # zero based

        $tmp_seq = $seq_1;
        $nn_s = ($tmp_seq =~ tr/Nn//);
        next if($nn_s > $max_nn_s);

        $start_2 = $start_1 + $insert;
        $end_2 = $end_1 + $insert;

        $seq_2 = substr($chr_seq, $start_2-1, $options->{'rl'}); # zero based
        $tmp_seq = $seq_2;
        $nn_s = ($tmp_seq =~ tr/Nn//);
        next if($nn_s > $max_nn_s);

        # read 2 needs rev comp
        $seq_2 =~ tr/ACGT/TGCA/;
        $seq_2 = reverse $seq_2;

        # put the source location of the read in the readname
        $read_name = $options->{'read_prefix'}.$chr.':'.($start_1+($ref_start-1)).'_'.($end_1+($ref_start-1)).':'.($start_2+($ref_start-1)).'_'.($end_2+($ref_start-1)).':'.$read_count.'/';

        push @reads_1,  $read_name."1\n"
                        .$seq_1."\n"
                        ."+\n"
                        .$qual_str."\n";
        push @reads_2,  $read_name."2\n"
                        .$seq_2."\n"
                        ."+\n"
                        .$qual_str."\n";

        # reset insert size counter
        $ins_arr_size = @{$options->{'rand_inserts'}} if($ins_arr_size == 0);
        $read_count++;
        $readpairs_for_depth--;
        $ins_arr_size--;

        $since_success = 0;

        if($read_count % CACHE_X_READS == 0) {
          write_out($options, \@reads_1, \@reads_2);
            @reads_1 = ();
            @reads_2 = ();
        }
      }
    }
    #print STDERR 'Finished ref: ',$chr,"\n";
  }
  write_out($options, \@reads_1, \@reads_2) if(@reads_1 > 0);
  return 1;
}

sub write_out {
  my ($options, $reads_1, $reads_2) = @_;
  my $open_type = '>>';
  if(!$options->{'fastq_1'}) {
    # generate a file_path
    $options->{'fastq_1'} = $options->{'o'}.'_1.fastq';
    $options->{'fastq_2'} = $options->{'o'}.'_2.fastq';
    $open_type = '>';
  }
  create_or_append($open_type, $options->{'fastq_1'}, $reads_1);
  create_or_append($open_type, $options->{'fastq_2'}, $reads_2);
  return 1;
}

sub create_or_append {
  my ($open_type, $file, $data) = @_;
  open my $FQ, $open_type, $file or croak 'Failed to create or append to '.$file;
  print $FQ @{$data} or croak 'Failed to write records to '.$file;
  close $FQ or croak 'Failed to close '.$file;
  return 1;
}


sub uuid {
  my $ug= new Data::UUID;
  my $uuid = $ug->create_str();
  $uuid =~ s/\-//g;
  return $uuid;
}

sub get_refs {
  my ($options) = @_;
  my $fai = $options->{'r'}.'.fai';
  croak 'Cannot find fai file for '.$options->{'r'} if(!-e $fai);
  my %refs;
  if(defined $options->{'e'}) {
    open my $GFF, '<', $options->{'e'} or croak 'Failed to open '.$options->{'e'}.' for reading';
      while(my $line=<$GFF>) {
        next if(index($line, '#') == 0);
        chomp $line;
        my ($chr, undef, undef, $start, $stop) = split "\t", $line;
        push @{$refs{$chr}}, [$start, $stop];
      }
    close $GFF;
  }
  else {
    open my $FAI, '<', $fai or croak 'Failed to open '.$fai.' for reading';
    while(my $line=<$FAI>) {
      chomp $line;
      my ($ref, $length) = split /\t/, $line;
      $refs{$ref} = [[1, $length]];
    }
    close $FAI;
  }
  # ref, length, other stuff....
  return \%refs;
}

sub setup_fai {
  my ($options) = @_;
  my $fai = Bio::DB::Sam::Fai->load($options->{'r'});
  return $fai;
}

sub rand_inserts {
  # this is used to approximate a standard curve of insert distributions
  my @insert_ranges = ([100,300],[150,250],[160,240],[170,230],[175,225],[185,215],[195,205]);
  my ($min, $max);
  my @inserts;
  for(1..MAX_INSERTS) {
    ($min, $max) = @{$insert_ranges[int(rand(@insert_ranges))]};
    push @inserts, $min + int(rand($max-$min));
  }
  return \@inserts;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'r|reference=s' => \$opts{'r'},
		'e|exome=s' => \$opts{'e'},
		'd|depth=n' => \$opts{'d'},
		'o|output=s' => \$opts{'o'}
	);

	pod2usage(1) if((keys %opts) == 0);
	pod2usage(0) if($opts{'h'});
	return \%opts;
}

sub validateInput {
  my $opts = shift;
  pod2usage(1) if(! $opts->{'r'} && ! $opts->{'o'});
  return 1;
}

__END__

=head1 NAME

inSilicoReadpairs.pl - Generates a 'perfect' set of fastq inputs from a reference

=head1 SYNOPSIS

inSilicoReadpairs.pl [-h]

  Required inputs:

    --reference  (-r)   Reference file in fasta format with associated fai (fasta index)

    --depth      (-d)   Depth to simulate coverage up to
                         - Depth is based on the reference file
                         - recommend no more that 5x per run, this will aid parallel processing
                           of alignement in downstream analysis
                         - 5x of an exome set is likely to be the equivalent of ~200x in regions
                           of interest

    --output     (-o)   Output location stub (needs to write 2 fastq files so STDOUT will not work)
                         - '_1.fastq' and '_2.fastq' will be automatically appended to stub path

  Optional inputs

    --exome      (-e)   Generate coverage over exome coordinates only based on this gff3 input

    --help       (-h)   Brief documentation

  Examples:

    perl inSilicoReadpairs.pl -r ~/human_37.fa -d 5 -o ~/lane_1
      - will generate wholegenome coverage and output to
        ~/lane_1_1.fastq
        ~/lane_1_2.fastq

    perl inSilicoReadpairs.pl -r ~/human_37.fa -d 5 -o ~/lane_1 -e ~/human_37_exome.gff3
      - will generate exome coverage and output to
        ~/lane_1_1.fastq
        ~/lane_1_2.fastq

=head1 NOTES

  Will use UUID to generate unique prefix for readnames ('-' will be removed)
  Will generate 100 b.p. reads with randomised insert sizes between 250 and 350

=head1 DEPENDENCIES

  Bio::DB::Sam
  Data::UUID

=cut
