#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Spec::Functions;
use POSIX;
use CGI qw/:all/;
use CGI::Carp;

my $prefix = catdir('..','..','Documents');
my $dir = catdir($prefix, 'focal_adhesions', 'results', 'latest');
my $col_num = 3;
my $image_size = 200;
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

#######################################
# Browse All
#######################################
sub get_experiment_list {
    my @exp_list;
    foreach (<$dir/*>) {
        push @exp_list, $_ if (-d $_);
    }
    return sort @exp_list;
}

sub get_last_movie_frame_list {
    my $folder = "edge_track";
    $folder = $_[0] if scalar(@_) > 0;

    my @frame_list;
    foreach (@exp_list) {
        push @frame_list, &get_last_movie_frame($_, $folder);
    }
    return @frame_list;
}

sub get_last_movie_frame {
    my $base_dir = $_[0];
    
    my $folder = "edge_track";
    $folder = $_[1] if scalar(@_) > 1;

    my @movie_files = <$base_dir/*/$folder/*>;
    @movie_files = sort @movie_files;
    
    return $movie_files[-1];
}

sub convert_folder_to_title {
    my ($title) = @_; 
    if ($title =~ /.*\/(.*?)$/) {
        $title = &convert_to_title($1);
    }
    
    return $title;
}

sub browse_all {
    print $q->start_html('Focal Adhesion Analysis Results');
    print $q->h1('Focal Adhesion Analysis Results');
    print $q->h2('Display Options:'), "\n";
    print "Column Count: <A HREF=browse.pl?col_num=", $col_num + 1, "&image_size=$image_size>+</A> 
                      <A HREF=browse.pl?col_num=", $col_num - 1, "&image_size=$image_size>-</A>"; 
    print "<BR>";
    print "Image Size: <A HREF=browse.pl?col_num=$col_num&image_size=", $image_size + 50, ">+</A> 
                      <A HREF=browse.pl?col_num=$col_num&image_size=", $image_size - 50, ">-</A>"; 
    print "<BR><BR>";
    
    my @rows;
    foreach my $i (0 .. $#exp_list/$col_num) {
        for my $j (0 .. $col_num - 1) {
            my $index = $i * $col_num + $j;
        
            next if ($index > $#exp_list);
    
            my $title = &convert_folder_to_title($exp_list[$index]);
            my $this_name = $exp_list[$index];
            $this_name =~ s/.*\/(.*?)/$1/;

            $rows[$i] .= "\n<td><A HREF=browse.pl?exp_name=$this_name>" . $q->h3($title) . "</A><BR>\n";
            
            my $this_file = &rm_prefix($frame_list[$index]);
            $rows[$i] .= $q->img({src => $this_file, width => $image_size. "px"}) . "</td>"; 
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
    my $base_dir = catdir($dir,$exp_name);
    my $base_no_pre = rm_prefix($base_dir);

    if (not( -e $base_dir)) {
        &display_cant_find_exp;
        die;
    }
    
    my $title = &convert_to_title($exp_name);

    print $q->start_html('Focal Adhesion Analysis Results: ' . $title);
    print $q->h1('Focal Adhesion Analysis Results: ' . $title);
    print $q->h2('Movies');
    
    my @movie_rows = "<td><h3><A HREF=\"" . catdir($base_no_pre,'edge_track.mov') . "\">Edge Tracking</A></h3></td>\n" . 
                     "<td><h3><A HREF=\"" . catdir($base_no_pre,'time_track.mov') . "\">Time Tracking</A></h3></td>";
    
    my $edge_file = &rm_prefix(&get_last_movie_frame($base_dir));
    my $time_file = &rm_prefix(&get_last_movie_frame($base_dir,'time_track'));

    push @movie_rows, "<td>" . $q->img({src => $edge_file, width => "300px"}) . "</td>" . 
                      "<td>" . $q->img({src => $time_file, width => "300px"}) . "</td>";
    
    print $q->table(Tr(\@movie_rows));
    
    print $q->h2('Plots');
    print $q->h3('Individual Adhesion Properties');
    
    my $plots_dir = catdir($base_no_pre,'adhesion_props','plots','png');
    my @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'area_vs_dist.png')), width => "400px" }) . "</td>" .
                  "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'area_vs_pax.png')), width => "400px" }) . "</td>";
    
    push @img_row, "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'sig_vs_dist.png')), width => "400px" }) . "</td>";
    print $q->table(Tr(\@img_row));
    
    print $q->h3('Lineage Properties');
    @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'longev_vs_pax.png')), width => "400px" }) . "</td>" .
               "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'longev_vs_s_dist.png')), width => "400px" }) . "</td>";
    print $q->table(Tr(\@img_row));
    
    print $q->h3('Pixel Value Properties');
    @img_row = "<td>" . $q->img({src => rm_prefix(catdir($plots_dir,'pix_max.png')), width => "400px" }) . "</td>"; 
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
