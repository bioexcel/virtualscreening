#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - genion.py
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
    default: "genion"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  gio_input_tpr_path:
    type: File
    inputBinding:
      position: 4
  gio_output_gro_path:
    type: string
    inputBinding:
      position: 5
    default: "gio.gro"
  gio_input_top_zip_path:
    type: File
    inputBinding:
      position: 6
  gio_output_top_zip_path:
    type: string
    inputBinding:
      position: 7
    default: "gio_zip.top"
outputs:
  gio_output_gro_file:
    type: File
    outputBinding:
      glob: $(inputs.gio_output_gro_path)
  gio_output_top_zip_file:
    type: File
    outputBinding:
      glob: $(inputs.gio_output_top_zip_path)
