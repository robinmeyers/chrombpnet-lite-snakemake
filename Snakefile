
config = {}

config['genome_fa'] = "reference/hg38.fa"

config['sample_bigwigs'] = {"D3_DMSO" : ["data/D3_DMSO.unstranded.bw"],
					       "D3_dBrd9" : ["data/D3_dBrd9.unstranded.bw"]}

config['peaks'] = "data/filtered.peaks.bed"
config['nonpeaks'] = "data/filtered.nonpeaks.bed"


rule all:
	input:
		expand("outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.attr.npz",
			sample = config['sample_bigwigs'].keys(), fold = range(0, 5)),
		expand("outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.npz",
			sample = config['sample_bigwigs'].keys(), fold = range(0, 5))
	run:
		print("workflow complete!")


rule all_jsons:
	input:
		expand("outs/{sample}/fold{fold}/bias.fit.json", sample = config['sample_bigwigs'].keys(), fold = range(0, 5)),
		expand("outs/{sample}/fold{fold}/fit.json", sample = config['sample_bigwigs'].keys(), fold = range(0, 5)),
		expand("outs/{sample}/fold{fold}/predict.json", sample = config['sample_bigwigs'].keys(), fold = range(0, 5)),
		expand("outs/{sample}/fold{fold}/interpret.json", sample = config['sample_bigwigs'].keys(), fold = range(0, 5)),
	run:
		print("workflow complete!")

rule get_fold_jsons:
	output: "splits/fold_{fold}.json"
	shell:
		"wget -O {output} https://mitra.stanford.edu/kundaje/oak/anusri/chrombpnet_data/input_files/folds/fold_{wildcards.fold}.json"

rule make_bias_fit_json:
	input:
		template = "json_templates/bias.fit.json",
		split = "splits/fold_{fold}.json"
	output: "outs/{sample}/fold{fold}/bias.fit.json"
	script: "scripts/make_bias_fit_json.py"

rule fit_bias_model:
	input: "outs/{sample}/fold{fold}/bias.fit.json"
	output: "outs/{sample}/fold{fold}/{sample}.fold{fold}.bias.torch"
	shell: "bpnet fit -p {input}"

rule make_fit_json:
	input:
		template = "json_templates/accessibility.fit.json",
		split = "splits/fold_{fold}.json"
	output: "outs/{sample}/fold{fold}/fit.json"
	script: "scripts/make_fit_json.py"

rule fit_accessibility_model:
	input:
		json = "outs/{sample}/fold{fold}/fit.json",
		bias_model = "outs/{sample}/fold{fold}/{sample}.fold{fold}.bias.torch"
	output: "outs/{sample}/fold{fold}/{sample}.fold{fold}.accessibility.torch"
	shell: "chrombpnet fit -p {input.json}"


rule make_predict_json:
	input: 
		template = "json_templates/predict.json",
		split = "splits/fold_{fold}.json"
	output: "outs/{sample}/fold{fold}/predict.json"
	script: "scripts/make_predict_json.py"

rule predict_model:
	input:
		json = "outs/{sample}/fold{fold}/predict.json",
		model = "outs/{sample}/fold{fold}/{sample}.fold{fold}.accessibility.torch"
	output: "outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.npz"
	shell: "chrombpnet predict -p {input.json}"

rule make_interpret_json:
	input:
		template = "json_templates/interpret.json",
		split = "splits/fold_{fold}.json"
	output: "outs/{sample}/fold{fold}/interpret.json"
	script: "scripts/make_interpret_json.py"

rule interpret_model:
	input:
		json = "outs/{sample}/fold{fold}/interpret.json",
		model = "outs/{sample}/fold{fold}/{sample}.fold{fold}.accessibility.torch"
	output: "outs/{sample}/fold{fold}/{sample}.fold{fold}.profile.attr.npz"
	shell: "chrombpnet interpret -p {input.json}"



