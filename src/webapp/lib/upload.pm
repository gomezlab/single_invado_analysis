package upload;
use Dancer ':syntax';
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Temp qw/ tempfile /;
use File::Basename qw/ dirname /;
use File::Spec::Functions;
use File::Path qw/ make_path /;
use File::Copy qw/ move /;
use Config::General;

###############################################################################
# Main
###############################################################################

get '/upload' => sub {
	template 'upload';
};

post '/upload' => sub {
	mkdir "uploaded_experiments" if (! -e "uploaded_experiments");

	my $puncta_file = upload('puncta_image') or die $!;
	my ($puncta_fh, $puncta_filename) = tempfile("puncta_XXXXXX",DIR=>"uploaded_experiments");
	while ($puncta_filename =~ /puncta_.*_.*/) {
		unlink $puncta_filename;
		($puncta_fh, $puncta_filename) = tempfile("puncta_XXXXXX",DIR=>"uploaded_experiments");
	}
	$puncta_file->copy_to($puncta_filename);
	
	&is_file_TIFF($puncta_filename);

	my $ecm_filename = $puncta_filename;
	$ecm_filename =~ s/puncta/ecm/;
	my $ecm_file = upload('ecm_image') or die $!;
	$ecm_file->copy_to($ecm_filename);

	if (! &is_file_TIFF($puncta_filename) | ! &is_file_TIFF($ecm_filename)) {
		template 'upload_problem';
	} else {
		#######################################################################
		# Config Processing
		#######################################################################
		my $cfg_filename = $puncta_filename;
		$cfg_filename =~ s/puncta/cfg/;
		$cfg_filename = "$cfg_filename.cfg";

		my %cfg;
		$cfg{submitter_ip} = request->address();
		if (params->{email}) {
			$cfg{email} = params->{email};
		}

		my $out_cfg = new Config::General(\%cfg);
		
		my $header = "<<include ../config/Invado_default.cfg>>";
		open CFG_OUT, ">$cfg_filename";
		print CFG_OUT "$header\n";
		print CFG_OUT $out_cfg->save_string;
		close CFG_OUT;
			
		&organize_uploaded_files($puncta_filename, $ecm_filename, $cfg_filename);
		
		$cfg_filename =~ /cfg_(.*)\.cfg/;

		template 'upload_success', { exp_status_link => "/exp_status/invado_$1" };
	}
};

###############################################################################
# Functions
###############################################################################

sub is_file_TIFF {
	my $file = shift @_;

	my $type_output = `file $file`;
	if ($type_output =~ /TIFF/) {
		return 1;
	} else {
		return 0;
	}
}

sub organize_uploaded_files {
	my $puncta_file = $_[0];
	my $gel_file = $_[1];
	my $cfg_file = $_[2];
	
	$puncta_file =~ /puncta_(.*)/;
	my $out_folder = "invado_$1";
	
	my $base_folder = dirname($puncta_file);
	my $puncta_folder = catdir($base_folder, $out_folder,'Images','puncta');
	my $gel_folder = catdir($base_folder, $out_folder,'Images','gel');
	make_path($puncta_folder, $gel_folder);
	
	move($puncta_file, catfile($puncta_folder,'puncta.tif'));
	move($gel_file, catfile($gel_folder,'gel.tif'));
	move($cfg_file, catdir($base_folder, $out_folder, "analysis.cfg"));
}

true;
