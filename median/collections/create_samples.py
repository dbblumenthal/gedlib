from subprocess import call

datasets = ["Mutagenicity-Correct.xml", "AIDS.xml", "Letter.xml"]
sample_script = "../../data/collections/sample.py"
size_ratios = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
sample_ids = range(5)
for dataset in datasets:
    for size_ratio in size_ratios:
        for sample_id in sample_ids:
            sample_name = dataset.split('.')[0] + "-" + str(int(size_ratio * 100)) + "-" + str(sample_id) + ".xml"
            command = "python " + sample_script + " " + dataset + " " + sample_name + " --balanced --size_ratio " + str(size_ratio)
            call(command, shell=True)
