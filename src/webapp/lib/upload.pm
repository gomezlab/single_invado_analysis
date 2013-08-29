package upload;
use Dancer ':syntax';
use Dancer::Session::YAML;
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Temp qw/ tempfile /;
use File::Basename qw/ dirname /;
use File::Spec::Functions;
use File::Path qw/ make_path rmtree/;
use File::Copy qw/ move /;
use File::Find;
use File::Basename;
use Config::General;
use Data::Dumper;
use Storable qw(lock_store lock_nstore lock_retrieve);
use Time::localtime;

my $out_folder = catdir('..','uploaded_experiments');
my $user_exp_info_file = '../user_exp_info.stor';

###############################################################################
# Main
###############################################################################

get '/upload' => sub {
	my $time_spacing = 1;
	if (defined param 'time_spacing') {
		$time_spacing = param 'time_spacing';
	}
	template 'upload', { time_spacing => $time_spacing };
};

post '/upload' => sub {
	if (! -e $out_folder) {
		make_path $out_folder or die $!;
		chmod 0777, $out_folder;
	}

	my $puncta_file = upload('puncta_file') or die $!;
	my $ECM_file = upload('ECM_file') or die $!;

	if (! &is_file_TIFF($puncta_file->tempname) &&
		! &is_file_TIFF($ECM_file->tempname)) {
		template 'upload_problem';
	} else {
		my $exp_out_dir = tempdir("IAS_temp_XXXXXX",DIR=>$out_folder);
		while ($exp_out_dir =~ /IAS_temp_.*_.*/) {
			rmtree($exp_out_dir);
			$exp_out_dir = tempdir("IAS_temp_XXXXXX",DIR=>$out_folder);
		}
		my $exp_ID = basename($exp_out_dir);
		$exp_ID =~ s/_temp_/_/;

		&organize_uploaded_files($exp_out_dir,$puncta_file,$ECM_file);

		#######################################################################
		# Config Processing
		#######################################################################
		my %cfg;
		$cfg{submitter_ip} = request->address();
		
		if (defined session('user_id')) {
			$cfg{session_user_id} = session('user_id');
		}

		my $date_str = `date`;
		chomp($date_str);
		$cfg{sub_date} = $date_str;
		
		$cfg{stdev_thresh} = params->{stdev_thresh_expansion} . " " . params->{stdev_thresh_seed};

		my @copy_if_defined = qw(min_puncta_size max_puncta_size max_ratio email
		note time_spacing exp_note gel_norm_level pixel_size);
		foreach (@copy_if_defined) {
			my $val = params->{$_};
			if (defined $val && $val ne "") {
				$cfg{$_} = params->{$_};
				if ($cfg{$_} eq "on") {
					$cfg{$_} = 1;
				}
				if ($cfg{$_} eq "off") {
					$cfg{$_} = 0;
				}
			}
		}

		my $out_cfg = new Config::General(\%cfg);
		
		my $header = "<<include ../config/webapp_default.cfg>>";
		open CFG_OUT, ">" . catfile($exp_out_dir, "analysis.cfg");
		print CFG_OUT "$header\n";
		print CFG_OUT $out_cfg->save_string;
		close CFG_OUT;
		chmod 0777, catfile($exp_out_dir, "analysis.cfg");
		# &textme_on_upload(%cfg);

		#######################################################################
		# Logged in user processing
		#######################################################################
		if (defined session('user_id')) {
			my %user_exp_data;
			if (-w $user_exp_info_file) {
				%user_exp_data = %{lock_retrieve($user_exp_info_file)};
			}
			push @{$user_exp_data{session('user_id')}}, $exp_ID;
			lock_store \%user_exp_data, $user_exp_info_file;
		}
		
		#######################################################################
		# Get Rid of Temp in Directory Name
		#######################################################################
		my $final_exp_out_dir = $exp_out_dir;
		$final_exp_out_dir =~ s/_temp_/_/;
		move($exp_out_dir, $final_exp_out_dir);

		#######################################################################
		# Return Page
		#######################################################################
		my $exp_status_url = "/exp_status/$exp_ID";
		my $email = param 'email';
		if ($email =~ /gmail/) {
			template 'upload_success', { exp_status_link => $exp_status_url,
				gmail => 1};
		} else {
			template 'upload_success', { exp_status_link => $exp_status_url };
		}
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
	my $out_folder = $_[0];
	my $puncta_file = $_[1];
	my $ECM_file = $_[2];
	
	my $puncta_folder = catdir($out_folder,'Images','puncta');
	make_path($puncta_folder);
	$puncta_file->copy_to(catfile($puncta_folder,'data.tif'));

	my $ECM_folder = catdir($out_folder,'Images','ECM');
	make_path($ECM_folder);
	$ECM_file->copy_to(catfile($ECM_folder,'data.tif'));

	chmod 0777, $out_folder;
	find(\&change_file_perm, $out_folder);
}

sub change_file_perm {
	chmod 0777, $_;
}

sub textme_on_upload {
	my %cfg = @_;
	if (localtime->hour() > 9 && localtime->hour() < 20) {
		my $send_text = "";
		if (defined $cfg{session_user_id}) {
			$send_text = "$cfg{session_user_id}";
		} elsif (defined $cfg{email}) {
			$send_text = "$cfg{email}";
		} else {
			$send_text = "$cfg{submitter_ip}";
		}
		system("curl -u l490C7fX0y:CzO9k6JVFAauVOaN2w81 -d email=airgram\@berginski.com --data-urlencode msg=\"$send_text\" https://api.airgramapp.com/1/send > /dev/null 2>/dev/null");
	}
}	

true;
