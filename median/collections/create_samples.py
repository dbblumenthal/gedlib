#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2020 by David B. Blumenthal                              #
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

from subprocess import call
import argparse

parser = argparse.ArgumentParser(description="Create subsamples of datasets.")
parser.add_argument("dataset", help="name of dataset collection XML file")
args = parser.parse_args()

sample_script = "../../data/collections/sample.py"
size_ratios = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
sample_ids = range(5)
for size_ratio in size_ratios:
    for sample_id in sample_ids:
        sample_name = args.dataset.split('.')[0] + "-" + str(int(size_ratio * 100)) + "-" + str(sample_id) + ".xml"
        command = "python " + sample_script + " " + args.dataset + " " + sample_name + " --size_ratio " + str(size_ratio)
        call(command, shell=True)
