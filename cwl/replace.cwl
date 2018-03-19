#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - replace.py
inputs:
  pdb_structure:
    type: File
    inputBinding:
      position: 1

outputs:
  fixed_pdb:
    type: File
    outputBinding:
      glob: fixed.pdb
