#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - rms.py
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
    default: "rms"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  rms_input_gro_path:
    type: File
    inputBinding:
      position: 4
  rms_input_trr_path:
    type: File
    inputBinding:
      position: 5
  rms_output_xvg_path:
    type: string
    inputBinding:
      position: 6
    default: "rms.xvg"
outputs:
  rms_output_xvg_file:
    type: File
    outputBinding:
      glob: $(inputs.rms_output_xvg_path)
