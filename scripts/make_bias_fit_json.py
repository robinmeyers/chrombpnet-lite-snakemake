import json

with open(snakemake.input.split, "r") as infile:
	split = json.load(infile)

with open(snakemake.input.template, "r") as infile:
	parameters = json.load(infile)

parameters['name'] = "outs/{sample}/fold{fold}/{sample}.fold{fold}.bias".format(sample = snakemake.wildcards.sample, fold = snakemake.wildcards.fold)
parameters['signals'] = snakemake.config['sample_bigwigs'][snakemake.wildcards.sample]
parameters['loci'] = snakemake.config['nonpeaks']
parameters['sequences'] = snakemake.config['genome_fa']
parameters['training_chroms'] = split['train']
parameters['validation_chroms'] = split['valid']

with open(snakemake.output[0], 'w') as outfile:
	outfile.write(json.dumps(parameters, sort_keys=True, 
		indent=4))