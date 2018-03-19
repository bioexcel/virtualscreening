#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  pdb_structure: File

outputs:
  refined_structure:
     type: File
     outputSource: extract_last_snapshot/refined_structure

steps:
  replace_atom_names:
    run: replace.cwl
    in:
      pdb_structure: pdb_structure
    out: [fixed_pdb]

  create_topology:
    run: pdb2gmx.cwl
    in:
      p2g_input_structure_pdb_path: replace_atom_names/fixed_pdb
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
      gpp_input_mdp_path:
        default:
          class: File
          location: ions.mdp
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
      gpp_input_mdp_path:
        default:
          class: File
          location: min.mdp
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
      gpp_input_mdp_path:
        default:
          class: File
          location: nvt.mdp
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
      gpp_input_mdp_path:
        default:
          class: File
          location: npt.mdp
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
      gpp_input_mdp_path:
        default:
          class: File
          location: eq.mdp
    out: [gpp_output_tpr_file]

  equilibration:
    run: mdrun.cwl
    in:
      md_input_tpr_path: equilibration_preprocess/gpp_output_tpr_file
    out: [md_output_gro_file, md_output_trr_file, md_output_cpt_file]

  extract_last_snapshot:
    run: extract.cwl
    in:
      input_gro_path: equilibration/md_output_gro_file
      input_trr_path: equilibration/md_output_trr_file
    out: [refined_structure]
