#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - mdrun.py
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
    default: "mdmin"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  md_input_tpr_path:
    type: File
    inputBinding:
      position: 4
  md_output_trr_path:
    type: string
    inputBinding:
      position: 5
    default: "md.trr"
  md_output_gro_path:
    type: string
    inputBinding:
      position: 6
    default: "md.gro"
  md_output_cpt_path:
    type: string?
    inputBinding:
      position: 7
    default: "md.cpt"
outputs:
  md_output_gro_file:
    type: File
    outputBinding:
      glob: $(inputs.md_output_gro_path)
  md_output_trr_file:
    type: File
    outputBinding:
      glob: $(inputs.md_output_trr_path)
  md_output_cpt_file:
    type: File?
    outputBinding:
      glob: $(inputs.md_output_cpt_path)
