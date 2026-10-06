#!/bin/zsh
# Run the 3-layer model once with nominal parameters. Run from the repository root.

start=$(date +%s)

ogs models/multilayer_irz_nominal/MULTI_BW_line_IRZ.prj -o data/runs/nominal
