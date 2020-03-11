import pandas as pd

data = pd.read_csv("clinical_data.csv", index_col=0, dtype=str)
samples = data.index
variables = data.columns
for variable in variables:
    collection = open("../../collections/IBD_{}.xml".format(variable), "w")
    collection.write("<?xml version=\"1.0\"?>\n")
    collection.write("<!DOCTYPE GraphCollection SYSTEM \"http://github.com/dbblumenthal/gedlib/data/collections/GraphCollection.dtd\">\n")
    collection.write("<GraphCollection>\n")
    for sample in samples:
        collection.write("\t<graph file=\"{}.gxl\" class=\"{}\"/>\n".format(sample, data.loc[sample, variable]))
    collection.write("</GraphCollection>\n")
    collection.close()