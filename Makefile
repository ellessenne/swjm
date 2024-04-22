.PHONY: style example_swjm

style:
	R -e "styler::style_dir(filetype = c('r', 'rmd', 'qmd'))"

example_swjm:
	make style
	R -f "01-example-swjm/01a-data-simulation.R"
	cd 01-example-swjm && stata-mp -e 01b-data-analysis.do
