#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - solvate.py
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
    default: "solvate"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  sol_input_solute_gro_path:
    type: File
    inputBinding:
      position: 4
  sol_output_gro_path:
    type: string
    inputBinding:
      position: 5
    default: "sol.gro"
  sol_input_top_zip_path:
    type: File
    inputBinding:
      position: 6
  sol_output_top_zip_path:
    type: string
    inputBinding:
      position: 7
    default: "sol_top.zip"
outputs:
  sol_output_gro_file:
    type: File
    outputBinding:
      glob: $(inputs.sol_output_gro_path)
  sol_output_top_zip_file:
    type: File
    outputBinding:
      glob: $(inputs.sol_output_top_zip_path)
