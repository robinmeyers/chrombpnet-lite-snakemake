import json

with open(snakemake.input.split, "r") as infile:
	split = json.load(infile)

with open(snakemake.input.template, "r") as infile:
	parameters = json.load(infile)

parameters['profile_filename'] = "outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.npz".format(sample = snakemake.wildcards.sample, fold = snakemake.wildcards.fold)
parameters['counts_filename'] = "outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.npz".format(sample = snakemake.wildcards.sample, fold = snakemake.wildcards.fold)
parameters['loci'] = snakemake.config['peaks']
parameters['sequences'] = snakemake.config['genome_fa']
parameters['chroms'] = split['test']
parameters['model'] = "outs/{sample}/fold{fold}/{sample}.fold{fold}.accessibility.torch".format(sample = snakemake.wildcards.sample, fold = snakemake.wildcards.fold)

with open(snakemake.output[0], 'w') as outfile:
	outfile.write(json.dumps(parameters, sort_keys=True, 
		indent=4))