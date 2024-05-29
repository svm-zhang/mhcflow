#!/bin/sh

perl -e "use List::Util"
perl -e "use List::MoreUtils"
perl -e "use Bio::DB::Sam"
perl -e "use Data::Dumper"
perl -e "use PolysolverUtils qw(get_hla_frequencies)"
perl -e "use PolysolverUtils qw(population_prior_likelihood_component)"
perl -e "use PolysolverUtils qw(get_err_hash_illumina)"
perl -e "use PolysolverUtils qw(get_score_hash_log)"
perl -e "use PolysolverUtils qw(get_likelihood_log)"
perl -e "use PolysolverUtils qw(md_to_match_string)"