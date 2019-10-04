#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2019 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of GEDLIB.                                           #
#                                                                          #
#   GEDLIB is free software: you can redistribute it and/or modify it      #
#   under the terms of the GNU Lesser General Public License as published  #
#   by the Free Software Foundation, either version 3 of the License, or   #
#   (at your option) any later version.                                    #
#                                                                          #
#   GEDLIB is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU Lesser General Public License for more details.                    #
#                                                                          #
#   You should have received a copy of the GNU Lesser General Public       #
#   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. #
#                                                                          #
#//////////////////////////////////////////////////////////////////////////#

##
# @file process_results.py
# @brief Processes the results of median tests.
#
# @details 
# Usage: 
# ```sh
# $ python process_results.py <RESULT_FILE> <OUTPUT_PREFIX>
# ```
#
# For more information, execute `$ python process_results.py --help`.

'''Processes the results of median tests.'''

import pandas as pd
import numpy as np
import argparse

def aggregate_results(csv_file):
    rs = pd.read_csv(csv_file)
    rs.itrs = rs.itrs.apply(lambda l: [int(x) for x in l.split(";")])
    grouped = rs.groupby(["percent", "init_type", "num_inits", "algo"], as_index=False)
    concat_lists = lambda ll: [x for l in ll for x in l] 
    as_list = lambda ll: [x for x in ll]
    ratio_greater = lambda l, t: float(len([x for x in l if x > t])) / float(len(l))
    agg_rs = grouped.agg({"time": np.mean, "time_init": np.mean, "time_converged": np.mean, "sod": np.mean, "sod_init": np.mean, "sod_converged": np.mean, "itrs": concat_lists, "state": as_list})
    agg_rs["itrs_max"] = agg_rs.itrs.apply(max)
    agg_rs["itrs_min"] = agg_rs.itrs.apply(min)
    agg_rs["itrs_mean"] = agg_rs.itrs.apply(np.mean)
    agg_rs["itrs_lower_quartile"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .25))
    agg_rs["itrs_median"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .5))
    agg_rs["itrs_upper_quartile"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .75))
    agg_rs["ratio_initialized"] = agg_rs.state.apply(lambda l: ratio_greater(l, 0))
    agg_rs["ratio_converged"] = agg_rs.state.apply(lambda l: ratio_greater(l, 1))
    agg_rs["ratio_refined"] = agg_rs.state.apply(lambda l: ratio_greater(l, 2))
    return agg_rs

def agg_rs_vs_percent(agg_rs):
    row_rs = pd.DataFrame([10, 20],columns = ["percent"])
    algos = ["BRANCH_FAST", "REFINE", "IPFP"]
    init_types = ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]
    nums_inits = [1, 2, 4, 8, 16, 32]
    cols = ['time', 'time_init', 'time_converged', 'sod', 'sod_init', 'sod_converged', 'itrs', 'state', 'itrs_max', 'itrs_min', 'itrs_mean', 'itrs_lower_quartile', 'itrs_median', 'itrs_upper_quartile', 'ratio_initialized', 'ratio_converged', 'ratio_refined']
    for algo in algos:
        for init_type in init_types:
            for num_inits in nums_inits:
                if num_inits > 1 and algo != "RANDOM":
                    continue
                prefix = algo + "_" + init_type + "_I" + str(num_inits) + "_"
                for col in cols:
                    row_rs[prefix + col] = agg_rs[(agg_rs.algo == algo) & (agg_rs.init_type == init_type) & (agg_rs.num_inits == num_inits)][col].reset_index()[col]
    return row_rs

def agg_rs_vs_num_inits(agg_rs):
    row_rs = pd.DataFrame([1, 2, 4, 8, 16, 32], columns = ["num_inits"])
    algos = ["BRANCH_FAST", "REFINE", "IPFP"]
    percents = [10, 20]
    cols = ['percent', 'time', 'time_init', 'time_converged', 'sod', 'sod_init', 'sod_converged', 'itrs', 'state', 'itrs_max', 'itrs_min', 'itrs_mean', 'itrs_lower_quartile', 'itrs_median', 'itrs_upper_quartile', 'ratio_initialized', 'ratio_converged', 'ratio_refined']
    for algo in algos:
        for percent in percents:
            prefix = algo + "_RANDOM_P" + str(percent) + "_"
            for col in cols:
                row_rs[prefix + col] = agg_rs[(agg_rs.algo == algo) & (agg_rs.init_type == "RANDOM") & (agg_rs.percent == percent)][col].reset_index()[col]
    return row_rs

parser = argparse.ArgumentParser(description="Processes results of median tests.")
parser.add_argument("results", help = "path to CSV file containing the results of the experiments")
parser.add_argument("prefix", help = "prefix of the CSV files generated by the script")
args = parser.parse_args()

agg_rs = aggregate_results(args.results)
agg_rs_vs_percent(agg_rs).to_csv(args.prefix + "x=percent.csv", index = False)
agg_rs_vs_num_inits(agg_rs).to_csv(args.prefix + "x=num_inits.csv", index = False)
