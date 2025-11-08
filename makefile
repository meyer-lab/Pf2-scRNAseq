.PHONY: clean test pyright

flist = $(wildcard pf2scrnaseq/figures/figure*.py)
allOutput = $(patsubst pf2scrnaseq/figures/figure%.py, output/figure%.svg, $(flist))

all: $(allOutput)

output/figure%.svg: pf2scrnaseq/figures/figure%.py
	@ mkdir -p ./output
	uv run fbuild $*

test: .venv
	uv run pytest -s -v -x

.venv:
	uv sync

coverage.xml: .venv
	uv run pytest --junitxml=junit.xml --cov=pf2scrnaseq --cov-report xml:coverage.xml

pyright: .venv
	uv run pyright pf2scrnaseq

clean:
	rm -rf output profile profile.svg
	rm -rf factor_cache