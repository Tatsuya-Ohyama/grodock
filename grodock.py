#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
grodock.py - docking (structure merge) tool for gromacs
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import re
import numpy as np

from mods.func_prompt_io import check_exist, check_overwrite
from mods.molecule_objects import Atom, Residue



# =============== variable =============== #
RE_QUOTED_TERM = re.compile(r"^[\"'](.+)[\"']$")
RE_BOX = re.compile(r'^(\s+\d+\.\d+){3}\s*$')
THRESHOLD = 0.001



# =============== function =============== #
def read_map(input_file):
	"""
	Function to read map file of correspondence table for ATOMTYPE

	Args:
		input_file (str): map file

	Returns:
		dict: {pdb_atomtype: gro_atomtype, ...}
	"""
	corr_table = {}
	list_term = []
	with open(input_file, "r") as obj_input:
		for line_val in obj_input:
			line_val = line_val.strip().replace(" ", "")
			if "," in line_val:
				list_term += line_val.split(",")
			elif len(line_val) == 0:
				continue
			else:
				list_term.append(line_val)

	flag_duplicate = False
	for term in [v for v in list_term if len(v) != 0]:
		key, val = term.split(":")
		obj_match = RE_QUOTED_TERM.match(key)
		if obj_match:
			key = obj_match.group(1)
		obj_match = RE_QUOTED_TERM.match(val)
		if obj_match:
			val = obj_match.group(1)

		if key in corr_table.keys():
			sys.stderr.write("ERROR: duplicated atomtypes `{0}` are found.\n".format(key))
			flag_duplicate = True

		corr_table[key] = val

	if flag_duplicate:
		sys.exit(1)

	return corr_table


def read_gro_file(input_file):
	"""
	Function to read gro file and generate Atom and Residue objects

	Args:
		input_file (str): gro file path

	Returns:
		list: molecule: [obj_residue, ...]
		list: box_info: [x(float), y(float), z(float)]
	"""
	molecule = []
	box_info = []
	with open(input_file, "r") as obj_input:
		for line_idx, line_val in enumerate(obj_input, 1):
			if line_idx > 2:
				if RE_BOX.search(line_val):
					box_info = [float(v) for v in line_val.strip().split()]
				else:
					residue_name = line_val[5:10].strip()
					residue_index = int(line_val[0:5].strip())

					if len(molecule) == 0 or "{0}.{1}".format(molecule[-1].name, molecule[-1].number) != "{0}.{1}".format(residue_name, residue_index):
						residue = Residue(residue_name, residue_index)
						molecule.append(residue)

					atom = Atom(line_val[10:15].strip(), line_val[15:20].strip())
					atom.coord = [float(v.strip()) for v in [line_val[20:28], line_val[28:36], line_val[36:44]]]
					molecule[-1].append_atom(atom)
	return molecule, box_info



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="grodock.py - docking (structure merge) tool for gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-g", dest="GRO_RECEPTOR_FILE", metavar="GRO_RECEPTOR_FILE.gro", required=True, help="biomolecules structure obtained from `gmx pdb2gmx`")
	parser.add_argument("-p", dest="PDB_LIGAND_FILE", metavar="LIGAND.pdb", required=True, help="ligand in complex structure")
	parser.add_argument("-a", dest="ACPYPE_LIGAND_FILE", metavar="ACPYPE_LIGAND.gro", required=True, help="ligand obtained from `acpype`")
	parser.add_argument("-m", dest="MAP_FILE", metavar="MAP.txt", help="mapping atomtype file (separated by comma or new line)\nEx. PDB_ATOM_TYPE: GRO_ATOM_TYPE, ...\nNote: if no mapping file, mapped in the order of the atoms.")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.gro", required=True, help="output file for docked complex structure")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()

	# check file
	check_exist(args.GRO_RECEPTOR_FILE, 2)
	check_exist(args.PDB_LIGAND_FILE, 2)
	check_exist(args.ACPYPE_LIGAND_FILE, 2)
	if args.MAP_FILE is not None:
		check_exist(args.MAP_FILE, 2)

	# read gro (biomolecule)
	biomolecule, box_info = read_gro_file(args.GRO_RECEPTOR_FILE)

	# read gro (ligand)
	ligand, _ = read_gro_file(args.ACPYPE_LIGAND_FILE)

	# read coresponding map of ATOMTYPE and check consistency of ATOMTYPE
	corr_table = {}
	corr_list = []
	if args.MAP_FILE is not None:
		# when specify map file
		corr_table = read_map(args.MAP_FILE)
		if len(ligand[0].atoms) != len(corr_table.keys()):
			sys.stderr.write("ERROR: number of atoms for ligand in {0} does not match with ATOMTYPE correspondence table.\n".format(args.ACPYPE_LIGAND_FILE))
			sys.exit(1)
	else:
		# when specify no map file
		corr_list = [obj_atom.name for obj_atom in ligand[0].atoms]

	# read pdb (ligand) and replace coordinates
	list_atomtype_ligand = [atom.name for atom in ligand[0].atoms]
	atom_idx = 0
	with open(args.PDB_LIGAND_FILE, "r") as obj_input:
		for line_val in obj_input:
			if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
				atom_name = line_val[12:17].strip()
				coord = [float(v.strip()) / 10 for i, v in enumerate([line_val[30:38], line_val[38:46], line_val[46:54]])]

				if args.MAP_FILE is None:
					# when specify no map file
					atom_name = corr_list[atom_idx]
					ligand[0].atoms[atom_idx].coord = coord
					atom_idx += 1
					if atom_idx > len(corr_list):
						atom_idx = 0

				else:
					# when specify map file
					if atom_name not in corr_table.keys():
						sys.stderr.write("ERROR: ATOMTYPE `{0}` not found in {1}\n".format(atom_name, args.MAP_FILE))
						sys.exit(1)

					atom_name = corr_table[atom_name]
					atom_idx = list_atomtype_ligand.index(atom_name)
					ligand[0].atoms[atom_idx].coord = coord


	# reassign residue index
	system = biomolecule + ligand
	total_atom = 0
	for residue_idx, obj_residue in enumerate(system, 1):
		obj_residue.index = residue_idx
		total_atom += len(obj_residue.atoms)

	# output file
	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)

	with open(args.OUTPUT_FILE, "w") as obj_output:
		obj_output.write("{0} (coordinates for {1}) docked by grodock.py\n".format(args.ACPYPE_LIGAND_FILE, args.PDB_LIGAND_FILE))
		obj_output.write("{0}\n".format(total_atom))
		num_atom = 0
		for obj_residue in system:
			for obj_atom in obj_residue.atoms:
				num_atom += 1
				obj_output.write("{0:>5}{1:<5}{2:>5}{3:>5}{4[0]:>8.3f}{4[1]:>8.3f}{4[2]:>8.3f}\n".format(obj_residue.index, obj_residue.name, obj_atom.name, num_atom, obj_atom.coord))
		obj_output.write("{0[0]:>10.5f} {0[1]:>10.5f} {0[2]:>10.5f}\n".format(box_info))
