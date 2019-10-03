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
	@return molecule(list), box_info(list): [obj_residue, ...], [x(float), y(float), z(float)]
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

					if len(molecule) == 0 or molecule[-1].id != "{0}.{1}".format(residue_name, residue_index):
						residue = Residue()
						residue.name = residue_name
						residue.index = residue_index
						molecule.append(residue)

					atom = Atom()
					atom.name = line_val[11:16].strip()
					atom.coord = [float(v.strip()) for v in [line_val[20:28], line_val[28:36], line_val[36:44]]]
					molecule[-1].append_atom(atom)
	return molecule, box_info



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


	# ファイルのチェック
	check_exist(args.pdb2gmx_file, 2)
	check_exist(args.editconf_file, 2)
	check_exist(args.pdb_ligand_file, 2)
	check_exist(args.acpype_ligand_file, 2)
	check_exist(args.map_file, 2)


	# ATOMTYPE の対応マップ取得
	corr_table = read_map(args.map_file)


	# gro (生体分子) ファイルの読み込み
	molecule_editconf, box_info = read_gro_file(args.editconf_file)


	# 平行移動量の取得
	shift_vector = np.zeros(3)	# 平行移動量
	molecule_pdb2gmx, _ = read_gro_file(args.pdb2gmx_file)

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
	del(molecule_pdb2gmx)

	# shift しただけかの確認
	diff_coord = coord_pdb2gmx - coord_editconf
	shift_vector = np.mean(diff_coord, axis = 0)
	if not np.all(diff_coord - shift_vector < THRESHOLD):
		# 平行移動量と平行移動量の平均が同一でない場合、回転が加わっている可能性があるため、警告を表示して終了
		sys.stderr.write("ERROR: `editconf` added rotate operation!!!\nThis program cannot dock.\n")
		sys.exit(1)


	# gro (生体分子) ファイルの読み込み
	biomolecule = molecule_editconf


	# gro (リガンド) ファイルの読み込み
	ligand, _ = read_gro_file(args.acpype_ligand_file)
	if len(ligand[0].atoms) != len(corr_table.keys()):
		sys.stderr.write("ERROR: number of atoms for ligand in {0} does not match with ATOMTYPE correspondence table.\n".format(args.acpype_ligand_file))
		sys.exit(1)


	# pdb (リガンド) ファイルの読み込みと座標書き換え
	list_atomtype_ligand = [atom.name for atom in ligand[0].atoms]
	with open(args.pdb_ligand_file, "r") as obj_input:
		for line_val in obj_input:
			if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
				atom_name = line_val[12:17].strip()
				coord = [float(v.strip()) / 10 - shift_vector[i] for i, v in enumerate([line_val[30:38], line_val[38:46], line_val[46:54]])]
				if atom_name not in corr_table.keys():
					sys.stderr.write("ERROR: ATOMTYPE `{0}` not found in {1}\n".format(atom_name, args.map_file))
					sys.exit(1)
				atom_name = corr_table[atom_name]
				atom_idx = list_atomtype_ligand.index(atom_name)
				ligand[0].atoms[atom_idx].coord = coord


	# 系内の残基の番号振り直し
	system = biomolecule + ligand
	total_atom = 0
	for residue_idx, obj_residue in enumerate(system, 1):
		obj_residue.index = residue_idx
		total_atom += len(obj_residue.atoms)


	# ファイル出力
	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		obj_output.write("{0} + {1} (coordinates for {2}) docked by grodock.py\n".format(args.editconf_file, args.acpype_ligand_file, args.pdb_ligand_file))
		obj_output.write("{0}\n".format(total_atom))
		num_atom = 0
		for obj_residue in system:
			for obj_atom in obj_residue.atoms:
				num_atom += 1
				obj_output.write("{0:>5}{1:<5}{2:>5}{3:>5}{4[0]:>8.3f}{4[1]:>8.3f}{4[2]:>8.3f}\n".format(obj_residue.index, obj_residue.name, obj_atom.name, num_atom, obj_atom.coord))
		obj_output.write("{0[0]:>10.5f} {0[1]:>10.5f} {0[2]:>10.5f}\n".format(box_info))
