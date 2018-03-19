#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - scwrl.py
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
    default: "scwrl"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  scw_input_pdb_path:
    type: File
    inputBinding:
      position: 4
  scw_output_pdb_path:
    type: string
    inputBinding:
      position: 5
    default: "mutated.pdb"
outputs:
  scw_output_pdb_file:
    type: File
    outputBinding:
      glob: $(inputs.scw_output_pdb_path)
