#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  scw_input_pdb_path: File
  mutation: string

outputs:
  gnuplot_output_png_file:
     type: File
     outputSource: gnuplot/gnuplot_output_png_file
  rms_output_xvg_file:
     type: File
     outputSource: rmsd/rms_output_xvg_file
  md_output_gro_file:
    type: File
    outputSource: equilibration/md_output_gro_file
  md_output_trr_file:
    type: File
    outputSource: equilibration/md_output_trr_file

steps:
  mutate_structure:
    run: scwrl.cwl
    in:
      scw_input_pdb_path: scw_input_pdb_path
      scw_mutation: mutation
    out: [scw_output_pdb_file]
  create_topology:
    run: pdb2gmx.cwl
    in:
      p2g_input_structure_pdb_path: mutate_structure/scw_output_pdb_file
    out: [p2g_output_gro_file, p2g_output_top_zip_file]
  create_water_box:
    run: editconf.cwl
    in:
      ec_input_gro_path: create_topology/p2g_output_gro_file
    out: [ec_output_gro_file]
  solvate:
    run: solvate.cwl
    in:
      sol_input_solute_gro_path: create_water_box/ec_output_gro_file
      sol_input_top_zip_path: create_topology/p2g_output_top_zip_file
    out: [sol_output_gro_file, sol_output_top_zip_file]

  ions_preprocess:
    run: grompp.cwl
    in:
      gpp_input_gro_path: solvate/sol_output_gro_file
      gpp_input_top_zip_path: solvate/sol_output_top_zip_file
    out: [gpp_output_tpr_file]

  add_ions:
    run: genion.cwl
    in:
      gio_input_tpr_path: ions_preprocess/gpp_output_tpr_file
      gio_input_gro_path: solvate/sol_output_gro_file
      gio_input_top_zip_path: solvate/sol_output_top_zip_file
    out: [gio_output_gro_file, gio_output_top_zip_file]

  minimization_preprocess:
    run: grompp.cwl
    in:
      gpp_input_gro_path: add_ions/gio_output_gro_file
      gpp_input_top_zip_path: add_ions/gio_output_top_zip_file
    out: [gpp_output_tpr_file]

  minimization:
    run: mdrun.cwl
    in:
      md_input_tpr_path: minimization_preprocess/gpp_output_tpr_file
    out: [md_output_gro_file, md_output_trr_file]

  nvt_dynamics_preprocess:
    run: grompp.cwl
    in:
      gpp_input_gro_path: minimization/md_output_gro_file
      gpp_input_top_zip_path: add_ions/gio_output_top_zip_file
    out: [gpp_output_tpr_file]

  nvt_dynamics:
    run: mdrun.cwl
    in:
      md_input_tpr_path: nvt_dynamics_preprocess/gpp_output_tpr_file
    out: [md_output_gro_file, md_output_trr_file, md_output_cpt_file]

  npt_dynamics_preprocess:
    run: grompp.cwl
    in:
      gpp_input_gro_path: nvt_dynamics/md_output_gro_file
      gpp_input_top_zip_path: add_ions/gio_output_top_zip_file
      gpp_input_cpt_path: nvt_dynamics/md_output_cpt_file
    out: [gpp_output_tpr_file]

  npt_dynamics:
    run: mdrun.cwl
    in:
      md_input_tpr_path: npt_dynamics_preprocess/gpp_output_tpr_file
    out: [md_output_gro_file, md_output_trr_file, md_output_cpt_file]

  equilibration_preprocess:
    run: grompp.cwl
    in:
      gpp_input_gro_path: npt_dynamics/md_output_gro_file
      gpp_input_top_zip_path: add_ions/gio_output_top_zip_file
      gpp_input_cpt_path: npt_dynamics/md_output_cpt_file
    out: [gpp_output_tpr_file]

  equilibration:
    run: mdrun.cwl
    in:
      md_input_tpr_path: equilibration_preprocess/gpp_output_tpr_file
    out: [md_output_gro_file, md_output_trr_file, md_output_cpt_file]

  rmsd:
    run: rms.cwl
    in:
      rms_input_gro_path: equilibration/md_output_gro_file
      rms_input_trr_path: equilibration/md_output_trr_file
    out: [rms_output_xvg_file]

  gnuplot:
    run: gnuplot.cwl
    in:
      gnuplot_input_xvg_path: rmsd/rms_output_xvg_file
    out: [gnuplot_output_png_file]
