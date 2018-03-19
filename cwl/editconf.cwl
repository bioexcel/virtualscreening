#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - editconf.py
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
    default: "editconf"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  ec_input_gro_path:
    type: File
    inputBinding:
      position: 4
  ec_output_gro_path:
    type: string
    inputBinding:
      position: 5
    default: "ec.gro"
outputs:
  ec_output_gro_file:
    type: File
    outputBinding:
      glob: $(inputs.ec_output_gro_path)
