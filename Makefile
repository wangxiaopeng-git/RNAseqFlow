.PHONY: test validate dry-run run

test:
	python -m pytest -q
validate:
	rnaseqflow validate -c config/config.yaml
dry-run:
	rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --dry-run
run:
	rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --mem-mb 32000
