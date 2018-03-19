#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - gnuplot.py
inputs:
  system:
    type: string
    inputBinding:
      position: 1
    default: "linux"
  step:
    type: string
    inputBinding:
      position: 2
    default: "gnuplot"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  gnuplot_input_xvg_path:
    type: File
    inputBinding:
      position: 4
  gnuplot_output_png_path:
    type: string
    inputBinding:
      position: 5
    default: "gnuplot.png"
outputs:
  gnuplot_output_png_file:
    type: File
    outputBinding:
      glob: $(inputs.gnuplot_output_png_path)
