# -*- coding: utf-8 -*-
###########################################################################################
#########      cuRRBS: customised Reduced Representation Bisulfite Sequencing     #########
###########################################################################################
#
# Created by Daniel E. Martin-Herranz, Antonio J.M. Ribeiro and Thomas M. Stubbs.
#
# Copyright (C) 2016,2017 D.E. Martin-Herranz, A.J.M. Ribeiro, T.M. Stubbs.
#
# This file is part of cuRRBS.
#
# cuRRBS is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, either 
# version 3 of the License, or (at your option) any later version.
#
# cuRRBS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with cuRRBS. 
# If not, see <http://www.gnu.org/licenses/>.
#
# USAGE (see help page): python cuRRBS.py -h

__version_info__ = ('1','02')
__version__ = '.'.join(__version_info__)

try:
    import sys
    import argparse
    import os
    import datetime
    import csv
    import numpy
    import bisect
    import itertools
    import copy
    import math
    import multiprocessing
    import random
    import time
    from collections import OrderedDict
except ImportError, e:
    print(
        "\n\nERROR: The python module '{}' is not installed in your system.\n"
            "       Please install it and try again.\n".format(
               e.message.split()[-1])
           )
    sys.exit()

DESCRIPTION = ("""

##########################################################################################
########      cuRRBS: customised Reduced Representation Bisulfite Sequencing     #########
##########################################################################################

Created by Daniel E. Martin-Herranz, Antonio J.M. Ribeiro and Thomas M. Stubbs.

Copyright (C) 2016,2017 D.E. Martin-Herranz, A.J.M. Ribeiro, T.M. Stubbs.

cuRRBS comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to 
redistribute it under the conditions in the GPU GPL-3.0 license.

cuRRBS v{0}

#########################################################################################
"""
).format(__version__)

NO_RESULTS_MESSAGE = ("""
#########################################################################################

No individual enzymes or enzyme combinations were found that satisfy the used settings. 
No output file will be created. 
Please rerun the software with a lower C_Score value (or a higher C_NF/1000 value).

""")

SUCCESSFUL_RUN_MESSAGE = ("""
#########################################################################################

cuRRBS finished its execution successfully!

The output can be found in the following file:
    {}

#########################################################################################

If cuRRBS was useful for your research, please cite our paper and follow us on:

	 https://github.com/demh
	 https://github.com/jose-mr
	 https://github.com/t0mstubbst

#########################################################################################

"""
)

class MyLogger:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = file(filename, 'a')

    def write(self, text):
        if not "Progress:" in text:
            self.logfile.write(text)
        self.stdout.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()

    def flush(self):
        self.stdout.flush()

class ArgumentParserError(Exception): pass

class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        """Default message for all argument errors"""
        print("\nERROR: {}".format(message))
        print("\nRun with the -h option to see all the options and "
              "documentation\n")
        sys.exit()

def file_exists(parser, arg):
    """If file does not exists raise error"""
    if not os.path.exists(arg):
        raise argparse.ArgumentTypeError("File doesn't exist: {}".format(arg))
    else:
        return arg

#def folder_already_exists(parser, arg):
#    """If folder already exists and it is not empty raise error"""
#    if os.path.exists(arg):
#        if os.listdir(arg):
#            raise argparse.ArgumentTypeError(
#                "Folder already exists and it is not empty: {}".format(arg))
#    else:
#        return arg

def check_range(parser, arg, lim_min=None, lim_max=None, ttype=int):
    """If user input is not inside this range, raise error"""
    arg = ttype(arg)
    if lim_min != None and lim_max != None:
        if not lim_min <= arg <= lim_max:
            raise argparse.ArgumentTypeError(
                "should be between {} and {}.".format(lim_min, lim_max))
    elif not lim_max:
        if arg <= lim_min:
            raise argparse.ArgumentTypeError(
                "should be greater than {}.".format(lim_min))

    return arg

parser = ThrowingArgumentParser(
    description=DESCRIPTION,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

#parser = argparse.ArgumentParser(description="")

required = parser.add_argument_group("required arguments")
parser._action_groups = [parser._action_groups[2], parser._action_groups[1]]

parser.add_argument(
    '-V', '--version', 
    action='version', 
    version="{prog}s ({version})".
    format(prog="%(prog)", version=__version__)
)

required.add_argument(
    "-o", metavar="Output folder",
    required=True,
    help="path to the main output folder",
#    type=lambda x: folder_already_exists(parser, x),
)

required.add_argument(
    "-p", metavar="Pre-computed files",
    required=True,
    help="path to the folder that contains the pre-computed files",
    type=lambda x: file_exists(parser, x),
)

required.add_argument(
    "-e", metavar="Enzymes to check", 
    required=True,
    help=("path to the text file which contains the enzymes to check in the "
          "pipeline (e.g. utils/enzymes_to_check_CpG.txt)"),
    type=lambda x: file_exists(parser, x),
)

required.add_argument(
    "-a", metavar="Annotation for sites",
    required=True,
    help=("path to the CSV file which contains the annotation for the "
          "sites of interest (e.g. examples/epigenetic_clock_human_hg38_sites_annotation.csv)"),
    type=lambda x: file_exists(parser, x),
)

required.add_argument(
    "-r", metavar="Read length", 
    type=lambda x: check_range(parser, x, 30, 300),
    required=True,
    help=("read length (in bp) that will be used during the sequencing "
          "for the customised RRBS experiment. This determines whether "
          "a site of interest is 'seen' in a size-selected fragment after the sequencing "
          "(i.e. only if it is close to one of the ends of the fragment). "
          "RANGE: 30-300"),
)

required.add_argument(
    "-s", metavar="adapters Size", 
    type=lambda x: check_range(parser, x, 0, None ),
    required=True,
    help=("total size (in bp) of the adapters used for the customised RRBS "
          "experiment. This will be used to calculate the experimental "
          "size range. e.g. in  the original RRBS protocol (Gu et al., "
          "Nature Methods, 2011) this value is 120 for single-end adapters"),
)

required.add_argument(
    "-c", metavar="C_Score constant", 
    type=lambda x: check_range(parser, x, 0, 1, ttype=float),
    required=True,
    help=("value for the C_Score constant (Score threshold). It must be a number "
          "(integer or float) in the interval (0,1]. Only those enzyme combinations "
          "with a Score > C_Score * max_Score are reported"),
)

required.add_argument(
    "-g", metavar="Genome size", 
    required=True,
    type=float,
    help=("size of the genome used to generate the pre-computed files (in "
	 "Mega-basepairs). The values for the genomes already available "
         "can be found in utils/genome_sizes.txt (e.g. for human hg38: 3088.286401)"),
)


parser.add_argument(
    "-k", metavar="C_NF/1000 constant", 
    default=0.2,
    type=lambda x: check_range(parser, x, 0, 1, ttype=float),
    help="value for the C_NF/1000 constant (NF threshold). It must be a number "
         "(integer or float) in the interval (0,1]. Only those enzyme combinations "
         "with a NF/1000 <= C_NF/1000 * ref_NF/1000 are reported, where ref_NF/1000 "
         "is the NF/1000 that would be generated in a whole-genome "
         "bisulfite-sequencing (WGBS) experiment. DEFAULT: 0.2"
)

parser.add_argument(
    "-d", metavar="experimental error", 
    default=20,
    type=lambda x: check_range(parser, x, 5, 500),
    help="experimental error assumed when performing the size selection step "
         "(in bp). When this value is increased, less size ranges are tested "
         "and the software is faster, but it is also more likely to miss a "
    "more optimal size range than the one reported. DEFAULT: 20. RANGE: 5-500"
)

parser.add_argument(
    "-t", metavar="outpuT size",
    default=30,
    type=lambda x: check_range(parser, x, 0, None ),
    help="maximum number of individual enzymes or enzyme combinations that "
         "will be reported in the final output file. DEFAULT: 30"
)


parser.add_argument(
    "-m", metavar="isoschizomers annotation file",
    default=os.path.dirname(os.path.realpath(__file__)) + "/utils/isoschizomers_CpG_annotation.csv",
    help="path to the file containing the information regarding "
         "the different isoschizomer families and methylation sensitivity of "
         "the restriction enzymes. "
         "DEFAULT: path/to/cuRRBS/utils/isoschizomers_CpG_annotation.csv",
    type=lambda x: file_exists(parser, x),

)

parser.add_argument(
    "-i",
    action="store_true",
    help="include the IDs of the sites of interest that will theoretically "
         "be sequenced in the final output file",
)

if len(sys.argv)==1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# make options
FG_MIN = 20
FG_MAX = 1000
NO_THREADS = 4
MAX_COMBS = 2
REF_NF_1000 = args.g *1000./float(args.r) # We calculate REF_NF_1000 = (genome_size / read_length) / 1000
START_TIME = time.time()
# working directory


def read_isoschizomers(filename):
    """Read isoschizomers file and return a dict with the info
       output_dict = {
          'familiy': ['member1', 'member2']
       }"""
    csv_reader = csv.reader(open(filename))
    next(csv_reader, None)
    output_dict = {}
    for family, enzyme, sensitive in csv_reader:
        if family not in output_dict:
            output_dict[family] = []
        if sensitive=='0':
            output_dict[family].append(enzyme)
    return output_dict


def read_sites(filename):
    """Read sites of interest as a list of tuples
    each tuple: (site_name, chr_name, coord, weight)"""
    csv_reader = csv.reader(open(filename))
    next(csv_reader, None)
    return [(row[0],row[1],int(row[2]),float(row[3])) for row in csv_reader]

def find_fragment(cleaving_sites, coordinate, read_length=args.r, to_print=False):    
    """Returns the length of the restriction fragment that contains this site"""
    index_coord = bisect.bisect(cleaving_sites, coordinate) - 1
    if index_coord == len(cleaving_sites):
        return # not found
    start = cleaving_sites[index_coord]
    length = cleaving_sites[index_coord+1] - start
    if FG_MIN > length or length > FG_MAX:
        return # outside required lengths
    if to_print:
        print(coordinate, start+read_length, start+length-read_length, length)
    if (start+read_length) <= coordinate < (start+length-read_length):
        return # outside read_length
    return length

def print_progress(n, total): 
    progress = float(n)/total * 100.
    time_passed = time.time() - START_TIME
    try:
        time_remaining = int((100.-progress)/(progress/time_passed)/60)
    except ZeroDivisionError:
        time_remaining = 0

    remaining_str = "{:5.0f}m".format(time_remaining) if time_remaining else " <1min"

    sys.stdout.write("\rProgress: {:5.1f}%  Time elapsed:{:5.0f}s "
                     "   Time remaining:{}".format(
        progress, time_passed, remaining_str))
    sys.stdout.flush()
        
def read_precomputed(enzymes, all_data, chr_names, first=False):
    """Reads the precomputed restriction sites for a list of enzymes
    
    For the first enzyme it will also read the chr names
    all_data = {'BsaWI': [arr_chr1, arr_chr2, ...], 'AsuII': ... }
    """
    if first:
        START_TIME = time.time()
    enzymes_len = len(enzymes)
    for no, enzyme in enumerate(enzymes):
        chr_arrays = []
        try:
            with open("{0}/{1}_pre_computed.txt".format(args.p, enzyme), 'r') as f:
                string_by_chr = f.read().replace('\n',' ').split('>')[1:]
                for chr_str in string_by_chr:
                    chr_name, chr_sites_str = chr_str.split(' ', 1)
                    chr_arrays.append(numpy.fromstring(
                        chr_sites_str, dtype=int, sep=' '))
                    if first and no==0:
                        chr_names.append(chr_name)
            all_data[enzyme] = chr_arrays
        except IOError:
            print("\nERROR: There is no pre_computed file for enzyme {}".format(
                enzyme))

def get_list_of_enzymes(enzymes_list_file):
    """Read list of restriction enzymes from file"""
    with open(enzymes_list_file, 'rU') as enzyme_file:
        enzymes = [ line.rstrip() for line in enzyme_file if line.rstrip()]
    return enzymes
        

def find_best_ev(scores_all, fg_counts_all, no_sites, max_score):
    """Find best enrichment value based on array of scores and fg counts"""
    # scores and fragment size count are binned 
    # according to the experimental error
    scores = numpy.add.reduceat(
        scores_all, range(0, len(scores_all), args.d))
    fg_counts = numpy.add.reduceat(
        fg_counts_all, range(0, len(scores_all), args.d))

    # calculate all ev values and save the best one
    ev = 10000
    start = end = score = nf = robustness = 0
    for x in range(0, len(scores)-1):
        for y in range(x+1, len(scores)):
            this_score = sum(scores[x:y])+scores_all[args.d*y]
            this_nf = (sum(fg_counts[x:y])+fg_counts_all[args.d*y])
            if this_score/max_score > args.c and (this_nf/1000.)/REF_NF_1000 <= args.k:
                try:
                    this_ev = -math.log10(this_score/this_nf*no_sites/max_score)
                except ValueError:
                    continue
                if this_ev < ev:
                    ev, start, end, score, nf = this_ev, x, y, this_score, this_nf

    # calculate robustness
    diff_sum = 0
    if score:
        for x, y in itertools.product([-1, 0, 1], repeat=2):
            if (start == 0 and x == -1):
                error_ev = ev
            else:
                try:
                    error_score = sum(scores[start+x:end+y])+scores_all[args.d*(end+y)]
                    error_ev = -math.log10(error_score/(sum(
                        fg_counts[start+x:end+y])+fg_counts_all[args.d*(
                            end+y)])*no_sites/max_score)
                except (ValueError, IndexError):
                    error_ev = ev
                
            diff_sum += abs(ev-error_ev)
        if diff_sum:
            omega = diff_sum/ev
            robustness = math.exp(-omega) 
        else:
            robustness = "NA"
    return (start*args.d, end*args.d, ev, score, nf, robustness)

def find_best_cut(enzyme_mixes, sites, restriction_sites_dict, chr_names, enzyme_mixes_dict):
    max_score = sum([site[3] for site in sites]) # site[3] = weight
    no_sites = len(sites)
    for no_mix, enzyme_mix in enumerate(enzyme_mixes):
        #print(enzyme_mix)
        fragment_sizes_by_chr = []
        scores_count = numpy.zeros(FG_MAX-FG_MIN+1)
        sites_found =  [[] for _ in range(FG_MAX-FG_MIN+1)]
        for no_chr, chromosome in enumerate(chr_names):
            chr_restriction_sites = numpy.unique(numpy.concatenate(
                 [restriction_sites_dict[enzyme][no_chr] for
                  enzyme in enzyme_mix]))
            fragment_sizes_by_chr.append(numpy.diff(chr_restriction_sites))
            for site in sites:
                if site[1] == chromosome:
                    #if site[0]=='S006' and chromosome=='chr10': 
                    #    print('oi')
                    #    fragment_size = find_fragment(chr_restriction_sites,
                    #                    site[2], to_print=True)
                    #    print(fragment_size)
                    fragment_size = find_fragment(chr_restriction_sites,
                                                  site[2])
                    if fragment_size:
                        scores_count[fragment_size-FG_MIN] += site[3]
                        #print(site[0])
                        sites_found[fragment_size-FG_MIN].append(site[0])
        fragment_sizes = numpy.concatenate(fragment_sizes_by_chr)
        fg_size_count = numpy.bincount( # TODO not necessary?
            fragment_sizes)[FG_MIN:FG_MAX+1]

        # this function will give the best possible EV for this enzyme mix
        start, end, ev, score, nf, robustness = find_best_ev(
            scores_count,  
            fg_size_count, 
            no_sites,
            max_score, 
        )
        sites_found = [name for sublist in sites_found[start:end+1] 
                       for name in sublist]
        sites_found.sort()
        enzyme_mixes_dict[enzyme_mix] = (
            start, end, score, score/max_score*100, nf, ev, robustness, sites_found)

def main():
    # create output dir
    try:
        os.mkdir(str(args.o))
    except OSError:
        pass
    output_folder_name = "{}/run_{}".format(
        args.o, datetime.datetime.now().strftime("%Y_%m_%d__%H_%M_%S"))
    try:
        os.mkdir(output_folder_name)
    except OSError:
        print("Cannot create output folder: {}".format(output_folder_name))
    
    log_file_name = "{}/LOG".format(output_folder_name)
    logger = MyLogger(sys.stdout, log_file_name)
    sys.stdout = logger

    print(" ".join(sys.argv))
    
    print(DESCRIPTION)

    #os.chdir(str(args.o))
    # read sites of interest file to a list of tuples
    # [ (site_name, chr_name, coord, weight), ... ]
    sites = read_sites(args.a)
    
    # read enzymes to check
    enzymes = get_list_of_enzymes(args.e)
    enzymes_len = len(enzymes)
    
    print
    print("Step 1/2:  Reading {} enzyme pre-digestions".format(
        len(enzymes)))
    
    # reading predigested files with parallelization
    manager = multiprocessing.Manager()
    restriction_sites = manager.dict()
    chr_names = manager.list()
    jobs = []
    for i in range(NO_THREADS):
        p = multiprocessing.Process(
            target=read_precomputed, 
            args=(enzymes[i::NO_THREADS], restriction_sites, chr_names,i==0),
        )
        jobs.append(p)
        p.start()

    while len(restriction_sites.keys())!= enzymes_len:
        current = len(restriction_sites.keys())
        print_progress(current, enzymes_len)
        time.sleep(0.5)

    print_progress(1,1)
    print("\n")

    for job in jobs:
        job.join()
    

    restriction_sites_dict = dict(restriction_sites)
    ## end of reading parallelization
    if  len(enzymes) != len(restriction_sites_dict.keys()):
        sys.exit()

    #make enzyme combinations
    enzyme_mixes = []
    for n in range(1, MAX_COMBS+1):
        enzyme_mixes += itertools.combinations(enzymes, n)
    random.shuffle(enzyme_mixes)
    enzyme_mixes_len = len(enzyme_mixes)

    print("Step 2/2: Finding the best enzyme mix in {} possibilities".format(
        len(enzyme_mixes)))

    # calculate best ev and score for all enzyme combinations 
    enzyme_mixes_dict = manager.dict()
    
    #find_best_cut(enzyme_mixes, sites, restriction_sites_dict,   
    #                                chr_names, enzyme_mixes_dict)

    jobs = []
    for i in range(NO_THREADS):
        p = multiprocessing.Process(
            target=find_best_cut, 
            args=(enzyme_mixes[i::NO_THREADS], sites, restriction_sites_dict, 
                  chr_names, enzyme_mixes_dict),
        )
        jobs.append(p)
        p.start()

    global START_TIME 
    START_TIME = time.time()
    while all([job.is_alive() for job in jobs]):
        current = len(enzyme_mixes_dict.keys())
        print_progress(current, enzyme_mixes_len)
        time.sleep(1)
    
    print_progress(1,1)
    print("\n")

    for job in jobs:
        job.join()
    # end of parallelization

    iso_dict = read_isoschizomers(args.m)

    enzyme_mixes_dict = dict(enzyme_mixes_dict)
    sorted_results = sorted(enzyme_mixes_dict.items(), key=lambda e: len(e[0]))
    sorted_results = sorted(enzyme_mixes_dict.items(), key=lambda e: e[1][5])
    output_lines = []
    sites_id_str = ""
    if args.i:
        sites_id_str=",Sites_IDs"
    for result in sorted_results[:args.t]:
        output = ""
        enzyme_mix = result[0]
        start, end, score, score_per_cent, nf, ev, robustness, sites_found = result[1]
        if score:
            output=("{0},{1}_{2},{3}_{4},{5},{6},{7},{8},{9:.14f},"
                  "{10:.15f},{11}|{12},{13}".format(
                    " AND ".join(["({})".format(' OR '.join(iso_dict[enzyme])) 
                                  for enzyme in enzyme_mix]),
                    FG_MIN+start + args.s, FG_MIN+end + args.s,
                    FG_MIN+start, FG_MIN+end,
                    score,
                    float("{0:.2f}".format(score_per_cent)),
                    float("{0:.3f}".format(nf/1000.)),
                    float("{0:.2f}".format(REF_NF_1000/(nf/1000.))),
                    ev,
                    robustness,
                    args.c, args.k,
                    len(sites_found),

            ))
        if output:
            if args.i:
                output = "{0},{1}".format(output, ";".join(sites_found))
            output_lines.append(output)


    if output_lines:
        output_csv_name = "{}/final_cuRRBS_output.csv".format(output_folder_name)
        with open(output_csv_name, 'w') as output_csv:
            output_csv.write(
              "Enzyme(s),Experimental_size_range,Theoretical_size_range,"
              "Score,%_max_Score,NF/1000,Cost_Reduction_Factor,"
              "Enrichment_Value,Robustness,C_Score|C_NF/1000,Number_of_sites{}\n"
            .format(sites_id_str))
            for line in output_lines:
                output_csv.write("{}\n".format(line))


        print(SUCCESSFUL_RUN_MESSAGE.format(output_csv_name))
    else:
        print(NO_RESULTS_MESSAGE)


if __name__ == "__main__":
    main()
    
