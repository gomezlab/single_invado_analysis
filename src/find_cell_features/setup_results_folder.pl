#!/usr/bin/env perl

=head1 Name

setup_results_folder.pl - Moves the experimental data into the proper location in the results folder

=head1 Synopsis

setup_results_folder.pl -cfg adhesion_config_file -debug

=cut

################################################################################
# Global Variables and Modules
################################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use File::Spec::Functions;

use lib "../lib";
use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;

#Perl built-in variable that controls buffering print output, 1 turns off 
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

my $matlab_wrapper;
if (defined $cfg{matlab_executable}) {
    $matlab_wrapper = Math::Matlab::Local->new({ cmd => "$cfg{matlab_executable} -nodisplay -nojvm -nosplash", });
} else {
    $matlab_wrapper = Math::Matlab::Local->new();
}

################################################################################
# Main Program
################################################################################

mkpath($cfg{individual_results_folder});

my @image_sets = (
    [qw(cell_mask_image_prefix cell_mask_file)], 
    [qw(adhesion_image_prefix adhesion_image_file)]
    );

my @matlab_code;

foreach (@image_sets) {
    my $prefix   = $cfg{$_->[0]};
    my $out_file = $cfg{$_->[1]};
    
    my @image_files = <$cfg{exp_data_folder}/$prefix*>;

    if ($opt{debug}) {
        if (scalar(@image_files) > 1) {
            print "Image files found: $image_files[0] - $image_files[$#image_files]\n";
        } else {
            print "Image file found: $image_files[0]\n";
        }
        print "For Config Variable: ",$_->[0],"\n\n";
    }

    push @matlab_code, &create_matlab_code(\@image_files,$prefix,$out_file);
}

my $error_file = catdir($cfg{exp_results_folder},$cfg{matlab_errors_folder},$cfg{setup_errors_file});
&Math::Matlab::Extra::execute_commands($matlab_wrapper,\@matlab_code,$error_file);

################################################################################
#Functions
################################################################################

sub create_matlab_code {
    my @image_files = @{$_[0]};
    my $prefix = $_[1];
    my $out_file = $_[2];

    my @image_stack_count = map { Image::Stack::get_image_stack_number($_) } @image_files;
    
    my @matlab_code;
    if (grep {$_ > 1} @image_stack_count) {
        if (scalar(@image_files) > 1) {
            die "Found more than one image stack in: ", join(", ",@image_files), "\n", 
                "Expected single image stack or multiple non-stacked files\n";
        }
        @matlab_code = &create_matlab_code_stack(\@image_files,$out_file);
    } else {
        @matlab_code = &create_matlab_code_single(\@image_files,$prefix,$out_file);
    }
    return @matlab_code;
}

sub create_matlab_code_stack {
    my @image_files = @{$_[0]};
    my $out_file = $_[1];
    
    my @matlab_code;
    my $min_max_file = catfile(dirname($image_files[0]), $cfg{min_max_file});
    
    $matlab_code[0] .= &create_extr_val_code(\@image_files,$min_max_file);
    
    my $total_stack_images = Image::Stack::get_image_stack_number($image_files[0]);
    foreach my $i_num (1 .. $total_stack_images) {
        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
            
        my $padded_num = sprintf("%0" . length($total_stack_images) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder},$padded_num);
        mkpath($output_path);
        my $final_out_file = catfile($output_path,$out_file);
        $matlab_code[0] .= "write_normalized_image('$image_files[0]','$final_out_file','$min_max_file','I_num',$i_num);\n";
    }
    return @matlab_code;
}

sub create_matlab_code_single {
    my @image_files = @{$_[0]};
    my $prefix = $_[1];
    my $out_file = $_[2];
    
    my @matlab_code;
    my $min_max_file = catfile(dirname($image_files[0]), $cfg{min_max_file});
    
    $matlab_code[0] .= &create_extr_val_code(\@image_files,$min_max_file);
    
    foreach my $file_name (@image_files) {
        my $i_num;
        if ($file_name =~ /$prefix(\d+)\./) {
            $i_num = $1;
        } else {
            warn "Unable to find image number in: $file_name, skipping this image.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
        
        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder},$padded_num);
        mkpath($output_path);
        my $final_out_file = catfile($output_path,$out_file);
        $matlab_code[0] .= "write_normalized_image('$file_name','$final_out_file','$min_max_file');\n";
    }
    return @matlab_code;
}

sub create_extr_val_code {
    my @image_files = @{$_[0]};
    my $min_max_file = $_[1];

    return "find_extr_values('$min_max_file','" . join("','",@image_files) . "');\n";
}
