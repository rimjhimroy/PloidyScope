.PHONY: env extract-real-fixture compare-real-fixture compare-canonical-fixture compare-serpentine-fixture

BASELINE_PYTHON := .tool-baselines/pixy-venv/bin/python
CANONICAL_VCF := /data/users/rchoudhury/Biscutella_serpentine/results/gatk4/final_vcf/missingness_filtered.vcf.gz
CANONICAL_METADATA := /data/users/rchoudhury/Biscutella_serpentine/configs/metadata_serpentine.tsv
CANONICAL_OUTDIR := tests/data/fixture_comparison/fixture
CANONICAL_PREFIX := comparison_fixture
CANONICAL_REGION := Bv1:1-10000000
CANONICAL_POPS := BL11-S2 DKS TNA TDS

env:
	python -m pip install -r requirements.txt

extract-real-fixture:
	$(BASELINE_PYTHON) scripts/extract_fixture_region.py \
		--vcf $(CANONICAL_VCF) \
		--metadata $(CANONICAL_METADATA) \
		--outdir $(CANONICAL_OUTDIR) \
		--output-prefix $(CANONICAL_PREFIX) \
		--region $(CANONICAL_REGION) \
		--selected-pops $(CANONICAL_POPS)

compare-real-fixture:
	$(BASELINE_PYTHON) scripts/compare_real_fixture.py

compare-canonical-fixture: extract-real-fixture compare-real-fixture

compare-serpentine-fixture: compare-canonical-fixture
