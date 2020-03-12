import xml.etree.ElementTree as ET
import argparse
parser = argparse.ArgumentParser(description="Transforms XML collection to .ds format")
parser.add_argument("dataset", help="Path to existing XML collection (without .xml suffix) or prefix thereof.")
parser.add_argument("--prefix", help="Path to existing XML collection.", action="store_true")
args = parser.parse_args()
datasets = ["{}.xml".format(args.dataset)]
if args.prefix:
    datasets = ["{}-{}-{}.xml".format(args.dataset, percent, collection_id) for percent in range(10, 91, 10) for collection_id in range(5)]
for dataset in datasets:
    xml_doc = ET.parse(dataset).getroot()
    graphs = [graph.attrib["file"] for graph in xml_doc]
    ds_file_name = dataset[:-4] + ".ds"
    with open(ds_file_name,mode="w") as ds_file:
        for graph in graphs:
            ds_file.write(graph + "\n")