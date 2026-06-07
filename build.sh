#!/usr/bin/env bash

set -euo pipefail

cd plot
./plot.py
bundle exec ./plot.rb
cd ..

build() {
	cd "$1"
	xelatex "$1"
	if [[ -e "$1.bib" ]]; then
		bibtex "$1"
		xelatex "$1"
	fi
	xelatex "$1"
	cd ..
}

build proposal
build seminar
build thesis
