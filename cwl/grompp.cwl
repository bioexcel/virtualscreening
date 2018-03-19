#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - grompp.py
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
    default: "gppions"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  gpp_input_gro_path:
    type: File
    inputBinding:
      position: 4
  gpp_input_top_zip_path:
    type: File
    inputBinding:
      position: 5
  gpp_input_mdp_path:
    type: File
    inputBinding:
      position: 6
  gpp_output_tpr_path:
    type: string
    inputBinding:
      position: 7
    default: "grompp.tpr"
  gpp_input_cpt_path:
    type: File?
    inputBinding:
      position: 8
outputs:
  gpp_output_tpr_file:
    type: File
    outputBinding:
      glob: $(inputs.gpp_output_tpr_path)
