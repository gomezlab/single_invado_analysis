#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use POSIX;
use CGI qw/:all/;
use CGI::Carp;

my $prefix = catdir('..','..','Documents');
my $all_exp_dir = catdir($prefix, 'focal_adhesions', 'results', 'latest');
my $col_num = 3;
my $image_size = 200;

my $script_name = basename($0);
my $exp_name;

my $q = new CGI;
if (grep $_ eq 'col_num', $q->param()) {
    $col_num = $q->param('col_num');
    $col_num = 1 if $col_num < 1;
}
if (grep $_ eq 'image_size', $q->param()) {
    $image_size = $q->param('image_size');
    $image_size = 1 if $image_size < 1;
}
if (grep $_ eq 'exp_name', $q->param()) {
    $exp_name = $q->param('exp_name');
}

print $q->header('text/html');

my @exp_list = &get_experiment_list;
my @frame_list = &get_last_movie_frame_list;

###############################################################################
# Main Program
###############################################################################
if ($exp_name) {
    &browse_exp;
} else {
    &browse_all;
}

###############################################################################
# Functions
###############################################################################

#######################################
# General
#######################################

sub convert_to_title {
    my ($text) = @_;

    $text =~ s/_/ /g;
    $text =~ s/\b(\w)/\u$1/g;

    return $text;
}

sub rm_prefix {
    my $text = shift;
    $text =~ s/$prefix(.*)/$1/;
    return $text;
}

sub get_last_movie_frame {
    my $base_dir = $_[0];
    
    my $folder = "edge_track";
    $folder = $_[1] if scalar(@_) > 1;

    my @movie_files = <$base_dir/$folder/*>;
    @movie_files = sort @movie_files;
    
    if (scalar(@movie_files) == 0) {
        $movie_files[0] = "EMPTY! - $base_dir - $folder";
    }

    return $movie_files[-1];
}

sub convert_folder_to_title {
    my ($title) = @_; 
    if ($title =~ /.*\/(.*?)$/) {
        $title = &convert_to_title($1);
    }
    
    return $title;
}

#######################################
# Browse All
#######################################
sub get_experiment_list {
    my @exp_list;
    foreach (<$all_exp_dir/*>) {
        push @exp_list, $_ if (-d $_);
    }
    return sort @exp_list;
}

sub get_last_movie_frame_list {
    my $folder = "edge_track";
    $folder = $_[0] if scalar(@_) > 0;

    my @frame_list;
    foreach (@exp_list) {
        push @frame_list, &get_last_movie_frame(catdir($_,'movies','all'), $folder);
    }
    return @frame_list;
}

sub browse_all {
    print $q->start_html('Focal Adhesion Analysis Results');
    print $q->h1('Focal Adhesion Analysis Results');
    print $q->h2('Display Options:'), "\n";
    print "Column Count: <A HREF=$script_name?col_num=", $col_num + 1, "&image_size=$image_size>+</A> 
                      <A HREF=$script_name?col_num=", $col_num - 1, "&image_size=$image_size>-</A>"; 
    print "<BR>";
    print "Image Size: <A HREF=$script_name?col_num=$col_num&image_size=", $image_size + 50, ">+</A> 
                      <A HREF=$script_name?col_num=$col_num&image_size=", $image_size - 50, ">-</A>"; 
    print "<BR><BR>";
    
    my @rows;
    foreach my $i (0 .. $#exp_list/$col_num) {
        for my $j (0 .. $col_num - 1) {
            my $index = $i * $col_num + $j;
        
            next if ($index > $#exp_list);
    
            my $title = &convert_folder_to_title($exp_list[$index]);
            my $this_name = $exp_list[$index];
            $this_name =~ s/.*\/(.*?)/$1/;

            $rows[$i] .= "\n<td><A HREF=$script_name?exp_name=$this_name>" . $q->h3($title) . "</A><BR>\n";
            
            my $this_file = &rm_prefix($frame_list[$index]);
            $rows[$i] .= $q->a({-href => "$script_name?exp_name=$this_name"},$q->img({src => $this_file, width => $image_size. "px"})) . "</td>"; 
        }
        $rows[$i] .= "\n"
    }
    
    print $q->table({-border => 1},
                    Tr(\@rows));
    print $q->end_html;
}

#######################################
# Browse Experiment
#######################################

sub browse_exp {
    my $base_dir = catdir($all_exp_dir,$exp_name);
    my $base_no_pre = rm_prefix($base_dir);

    if (not( -e $base_dir)) {
        &display_cant_find_exp;
        die;
    }
    
    my $title = &convert_to_title($exp_name);

    print $q->start_html('Focal Adhesion Analysis Results: ' . $title);
    print $q->h1('Focal Adhesion Analysis Results: ' . $title);
    print $q->h2('Movies');
    print $q->h3('Unfiltered');
    
    my $movies_base = catdir($base_dir,'movies');
    my $movies_base_no_pre = rm_prefix(catdir($base_dir,'movies'));

    my @movie_rows = "<td><h3><A HREF=\"" . catdir($movies_base_no_pre,'all','edge_track.mov') . "\">Edge Tracking</A></h3></td>\n" . 
                     "<td><h3><A HREF=\"" . catdir($movies_base_no_pre,'all','time_track.mov') . "\">Time Tracking</A></h3></td>";
    
    my $edge_file = &rm_prefix(&get_last_movie_frame(catdir($movies_base,'all')));
    my $time_file = &rm_prefix(&get_last_movie_frame(catdir($movies_base,'all'),'time_track'));
    push @movie_rows, "<td>" . $q->img({src => $edge_file, width => "80%"}) . "</td>" . 
                      "<td>" . $q->img({src => $time_file, width => "80%"}) . "</td>";
    
    print $q->table(Tr(\@movie_rows));
    
    print $q->h3('Lineage Longevity Filtered (alive for 5 minutes)');
    my @movie_rows = "<td><h3><A HREF=\"" . catdir($movies_base_no_pre,'longev_filtered','edge_track.mov') . "\">Edge Tracking</A></h3></td>\n" . 
                     "<td><h3><A HREF=\"" . catdir($movies_base_no_pre,'longev_filtered','time_track.mov') . "\">Time Tracking</A></h3></td>";
    
    $edge_file = &rm_prefix(&get_last_movie_frame(catdir($movies_base,'longev_filtered')));
    $time_file = &rm_prefix(&get_last_movie_frame(catdir($movies_base,'longev_filtered'),'time_track'));
    
    push @movie_rows, "<td>" . $q->img({src => $edge_file, width => "80%"}) . "</td>" . 
                      "<td>" . $q->img({src => $time_file, width => "80%"}) . "</td>";
    
    print $q->table(Tr(\@movie_rows));
    
    
    print $q->h2('Plots');
    print $q->h3('Individual Adhesion Properties');
    
    my $plots_dir = catdir($base_no_pre,'adhesion_props','plots','png');
    my @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'area_vs_dist.png')), width => "80%" }) . "</td>" .
                  "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'area_vs_pax.png')), width => "80%" }) . "</td>";
    
    push @img_row, "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'sig_vs_dist.png')), width => "80%" }) . "</td>";
    print $q->table(Tr({-valign=>"LEFT"},\@img_row));
    
    print $q->h3('Lineage Properties');
    @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'longev_vs_pax.png')), width => "80%" }) . "</td>" .
               "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'longev_vs_s_dist.png')), width => "80%" }) . "</td>";
    push @img_row, "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'longev_vs_largest_area.png')), width => "80%" }) . "</td>";
    print $q->table(Tr(\@img_row));
    
    print $q->h3('Pixel Value Properties');
    @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'pix_max.png')), width => "80%" }) . "</td>" .
               "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'pix_average.png')), width => "80%" }) . "</td>";
    print $q->table(Tr(\@img_row));
    
    print $q->h3('Misc Files');
    my @list = ($q->a({href => catdir($base_no_pre, 'tracking_seq.csv')}, 'Tracking Matrix') .
        ": This file specifies each of the focal adhesion lineages. Each row of the spreadsheet specifies the focal adhesions that make up a single lineage. The focal adhesions in each lineage are identified using the numbers assigned by the matlab command bwlabel minus one. Negative numbers are used to indicate lineages that have not been born or that have died.");
    print $q->ul(li(\@list));
    print $q->end_html;
}

sub display_cant_find_exp {
    print $q->start_html('Sorry');
    print $q->h1('Can\'t find the experiment named:', $exp_name);
    print $q->end_html;
}
