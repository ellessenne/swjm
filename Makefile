.PHONY: style example_swjm example_swjm_individual

style:
	R -e "styler::style_dir(filetype = c('r', 'rmd', 'qmd'))"

example_swjm:
	make style
	R -f "01-example-swjm/01a-data-simulation.R"
	cd 01-example-swjm && stata-mp -e 01b-data-analysis.do

example_swjm_individual:
	make style
	R -f "02-example-swjm-individual/02a-data-simulation.R"
	cd 02-example-swjm-individual && stata-mp -e 02b-data-analysis.do
