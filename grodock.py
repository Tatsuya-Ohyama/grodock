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

from basic_func import check_exist, check_overwrite
from basic_class import TemplateClass
from molecule_objects import Atom, Residue


# =============== variable =============== #
RE_QUOTED_TERM = re.compile(r"^[\"'](.+)[\"']$")
RE_BOX = re.compile(r'^(\s+\d+\.\d+){3}\s*$')
THRESHOLD = 0.001


# =============== function =============== #
def read_map(input_file):
	"""
	ATOMTYPE の対応マップファイルを読み込む関数
	@param input_file: マップファイル
	@return corr_table(dict): ATOMTYPE の対応マップの dict
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

	for term in [v for v in list_term if len(v) != 0]:
		key, val = term.split(":")
		obj_match = RE_QUOTED_TERM.match(key)
		if obj_match:
			key = obj_match.group(1)
		obj_match = RE_QUOTED_TERM.match(val)
		if obj_match:
			val = obj_match.group(1)

		corr_table[key] = val

	return corr_table



def read_gro_file(input_file):
	"""
	Read gro file and generate Atom and Residue objects
	@param input_file: gro file path
	@return molecule[list]: [obj_residue, ...]
	"""
	molecule = []
	with open(input_file, "r") as obj_input:
		for line_idx, line_val in enumerate(obj_input, 1):
			if line_idx > 2:
				if not RE_BOX.search(line_val):
					residue_name = line_val[5:10].strip()
					residue_index = int(line_val[0:5].strip())

					if len(molecule) == 0 or molecule[-1].id != "{0}.{1}".format(residue_name, residue_index):
						residue = Residue()
						residue.name = residue_name
						residue.index = residue_index
						molecule.append(residue)

					atom = Atom()
					atom.name = line_val[11:16].strip()
					atom.coord = [float(v.strip()) for v in [line_val[20:28], line_val[28:36], line_val[36:44]]]
					molecule[-1].append_atom(atom)
	return molecule



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "grodock.py - docking (structure merge) tool for gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-g", dest = "pdb2gmx_file", metavar = "PDB2GMX_FILE.gro", required = True, help = "biomolecules structure obtained from `gmx pdb2gmx`")
	parser.add_argument("-e", dest = "editconf_file", metavar = "EDITCONF_FILE.gro", required = True, help = "biomolecule structure obtained from `gmx editconf`")
	parser.add_argument("-p", dest = "pdb_ligand_file", metavar = "LIGAND.pdb", required = True, help = "ligand in complex structure")
	parser.add_argument("-a", dest = "acpype_ligand_file", metavar = "ACPYPE_LIGAND.gro", required = True, help = "ligand obtained from `acpype`")
	parser.add_argument("-m", dest = "map_file", metavar = "MAP.txt", required = True, help = "mapping atomtype file (separated by comma or new line)\nEx. PDB_ATOM_TYPE: GRO_ATOM_TYPE, ...")
	parser.add_argument("-o", dest = "output_file", metavar = "OUTPUT.gro", required = True, help = "output file for docked complex structure")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.pdb2gmx_file, 2)
	check_exist(args.editconf_file, 2)
	check_exist(args.pdb_ligand_file, 2)
	check_exist(args.acpype_ligand_file, 2)
	check_exist(args.map_file, 2)

	# ATOMTYPE の対応マップ取得
	corr_table = read_map(args.map_file)

	# 平行移動量の取得
	molecule_pdb2gmx = read_gro_file(args.pdb2gmx_file)
	molecule_editconf = read_gro_file(args.editconf_file)
	shift_vector = np.zeros(3)	# 平行移動量

	coord_pdb2gmx = []
	coord_editconf = []
	for residue_pdb2gmx in molecule_pdb2gmx:
		list_residue_editconf = [residue_editconf for residue_editconf in molecule_editconf if residue_pdb2gmx.id == residue_editconf.id]
		if len(list_residue_editconf) == 0:
			sys.stderr.write("ERROR: {0} not found.\n".format(residue_pdb2gmx.id))
			sys.exit(1)
		if len(list_residue_editconf) > 1:
			sys.stderr.write("ERROR: multiple residues having the same id found: {0}.\n".format(residue_pdb2gmx.id))
			sys.exit(1)
		residue_editconf = list_residue_editconf[0]

		for atom_pdb2gmx in residue_pdb2gmx.atoms:
			list_atom_editconf = [atom_editconf for atom_editconf in residue_editconf.atoms if atom_pdb2gmx.name == atom_editconf.name]
			if len(list_atom_editconf) == 0:
				sys.stderr.write("ERROR: {1} in {0} not found.\n".format(residue_pdb2gmx.id, atom_pdb2gmx.name))
				sys.exit(1)
			if len(list_atom_editconf) > 1:
				sys.stderr.write("ERROR: multiple atoms having the same name found in {0}.\n".format(residue_pdb2gmx.id))
				sys.exit(1)
			atom_editconf = list_atom_editconf[0]

			coord_pdb2gmx.append(atom_pdb2gmx.coord)
			coord_editconf.append(atom_editconf.coord)

	coord_pdb2gmx = np.array(coord_pdb2gmx)
	coord_editconf = np.array(coord_editconf)

	# shift しただけかの確認
	diff_coord = coord_pdb2gmx - coord_editconf
	shift_vector = np.mean(diff_coord, axis = 0)
	if not np.all(diff_coord - shift_vector < THRESHOLD):
		# 平行移動量と平行移動量の平均が同一でない場合、回転が加わっている可能性があるため、警告を表示して終了
		sys.stderr.write("ERROR: `editconf` added rotate operation!!!\nThis program cannot dock.\n")
		sys.exit(1)


	# pdb (リガンド) ファイルの読み込み
	atoms = []
	with open(args.pdb_ligand_file, "r") as obj_input:
		for line_val in obj_input:
			if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
				atom = Atom()
				atom.name = line_val[12:17].strip()
				atom.coord = [float(v.strip()) for v in [line_val[30:38], line_val[38:46], line_val[46:54]]]
				atoms.append(atom)
	if len(atoms) != len(corr_table.keys()):
		sys.stderr.write("ERROR: number of atoms for ligand does not match with ATOMTYPE correspondence table.\n")
		sys.exit(1)


	# gro (生体分子) ファイルの読み込み
	pdb_atom_names = [v.name for v in atoms]
	corr_table_rev = {v: k for k, v in corr_table.items()}

	output_line = []
	num_atom = 0
	with open(args.editconf_file, "r") as obj_input:
		# ボックス情報付きの Gqd 構造
		for line_idx, line_val in enumerate(obj_input, 1):
			output_line.append(line_val)
			if line_idx > 2:
				num_atom += 1

	box_info = output_line.pop()
	num_atom -= 1

	# gro (リガンド) ファイルの読み込み
	with open(args.acpype_ligand_file, "r") as obj_input:
		for line_idx, line_val in enumerate(obj_input, 1):
			if line_idx > 2:
				if not RE_BOX.search(line_val):
					atom_name = line_val[11:16].strip()
					pdb_atom_name = corr_table_rev[atom_name]
					coord = [float(v) for v in atoms[pdb_atom_names.index(pdb_atom_name)].coord]
					coord = np.array(coord)
					coord = coord / 10 - shift_vector
					coord = coord.tolist()
					line_val = line_val[:20] + "".join(["{0:>8.3f}".format(v) for v in coord]) + "\n"
					output_line.append(line_val)
					num_atom += 1

	output_line.pop()
	num_atom -= 1
	output_line[1] = "{0}\n".format(num_atom)
	output_line.append(box_info)
	output_line[0] = "{0} + {1} (coordinates for {2}) docked by grodock.py\n".format(args.editconf_file, args.acpype_ligand_file, args.pdb_ligand_file)

	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		for line_val in output_line:
			obj_output.write(line_val)
