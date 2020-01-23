import xml.etree.ElementTree as ET
import argparse
parser = argparse.ArgumentParser(description="Transforms XML collection to .ds format")
parser.add_argument("dataset", help="Path to existing XML collection.")
args = parser.parse_args()
dataset = ET.parse(args.dataset).getroot()
graphs = [graph.attrib["file"] for graph in dataset]
with open(args.dataset.split(".")[0] + ".ds",mode="wt") as ds_file:
    for graph in graphs:
        ds_file.write(graph + "\n")