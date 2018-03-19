#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - pdb2gmx.py
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
    default: "pdb2gmx"
  properties_file:
    type: File
    inputBinding:
      position: 3
    default:
      class: File
      location: test/conf_1ps.yaml
  p2g_input_structure_pdb_path:
    type: File
    format: http://edamontology.org/format_1476 #PDB format from EDAM
    inputBinding:
      position: 4
  p2g_output_gro_path:
    type: string
    inputBinding:
      position: 5
    default: "p2g.gro"
  p2g_output_top_zip_path:
    type: string
    inputBinding:
      position: 6
    default: "p2g_top.zip"

outputs:
  p2g_output_gro_file:
    type: File
    outputBinding:
      glob: $(inputs.p2g_output_gro_path)
  p2g_output_top_zip_file:
    type: File
    outputBinding:
      glob: $(inputs.p2g_output_top_zip_path)
