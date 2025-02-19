high_cov_sites = open(snakemake.input["high_cov_sites"][0], "r").readlines()
high_cov_ids = ['{0}_{1}'.format(x.split("\t")[0].strip(), x.split("\t")[1].strip()) for x in high_cov_sites]
angsd_sites = open(snakemake.input["angsd_sites"][0], "r").readlines()
angsd_ids = ['{0}_{1}'.format(x.split("\t")[0].strip(), x.split("\t")[1].strip()) for x in angsd_sites]


def find_matching_index(list1, list2):
    """
    Function to quickly find index of element in second list matching items in first list
    """
    inverse_index = {element: index for index, element in enumerate(list1)}

    idx = [index for index, element in enumerate(list2) if element in inverse_index]
    return idx


indices = find_matching_index(high_cov_ids, angsd_ids)
sites = [angsd_sites[i] for i in indices]

with open(snakemake.output[0], "w") as fout:
    for site in sites:
        fout.write(site)
