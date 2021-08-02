#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np



# =============== function =============== #
def split_n(line, length):
	"""
	n 文字毎に分割する関数

	Args:
		line (str): 対象文字列
		length (int): 分割後の個別の文字列の長さ

	Returns:
		list
	"""
	return [line[i:i+length] for i in range(0, len(line), length)]



# =============== Atom class =============== #
class Molecule:
	def __init__(self, pdb_file):
		self._chains = []
		self._coords = None

		self.read_pdb(pdb_file)

	@property
	def atoms(self):
		return [atom for chain in self._chains for residue in chain.residues for atom in residue.atoms]

	@property
	def residues(self):
		return [residue for chain in self._chains for residue in chain.residues]

	@property
	def chains(self):
		return self._chains

	@property
	def coords(self):
		return self._coords

	@property
	def connect_info(self):
		return [atom.connect_info for chain in self._chains for residue in chain.residues for atom in residue.atoms if atom.connect_info is not None]


	def read_pdb(self, pdb_file):
		"""
		PDB ファイルを読み込むメソッド

		Args:
			pdb_file (str): PDB ファイルのパス

		Returns:
			self
		"""
		obj_residue = None
		obj_chain = None
		coords = []
		list_obj_atoms = []
		list_obj_atom_indexes = []
		with open(pdb_file, "r") as obj_input:
			for line_val in obj_input:
				if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
					# ATOM レコードの場合

					# Atom オブジェクトを作成する
					obj_atom = Atom(line_val[12:16], int(line_val[6:11].strip()))
					list_obj_atoms.append(obj_atom)

					if len(line_val) >= 60:
						obj_atom.set_occupancy(float(line_val[54:60].strip()))

					if len(line_val) >= 66:
						obj_atom.set_temp_factor(float(line_val[60:66].strip()))

					if len(line_val) >= 78:
						obj_atom.set_element(line_val[76:78])

					if len(line_val) >= 80:
						obj_atom.set_charge(line_val[78:80])

					if obj_residue is None or obj_residue.unique != "{0}-{1}-{2}".format(line_val[17:20], line_val[21:22], int(line_val[22:26].strip())):
						# 操作中の Residue オブジェクトがない場合か、現在の原子が別の Residue オブジェクトの場合
						# Residue オブジェクトを作成する
						obj_residue = Residue(line_val[17:20], int(line_val[22:26].strip()))
						if obj_chain is None or obj_chain.name != line_val[21:22]:
							# Chain オブジェクトがない場合
							# Chain オブジェクトを作成する
							obj_chain = Chain(line_val[21:22])
							self._chains.append(obj_chain)
						obj_chain.append_residue(obj_residue)
						obj_residue.set_chain(obj_chain)
					obj_residue.append_atom(obj_atom)
					obj_atom.set_residue(obj_residue)
					obj_atom.set_chain(obj_chain)

					# 座標を追加する
					coords.append([
						float(line_val[30:38].strip()),
						float(line_val[38:46].strip()),
						float(line_val[46:54].strip())
					])

				elif line_val.startswith("TER"):
					# TER レコードの場合
					obj_chain = None

				elif line_val.startswith("CONECT"):
					# CONECT レコードの場合
					if len(list_obj_atom_indexes) == 0:
						list_obj_atom_indexes = [atom.number for atom in list_obj_atoms]
					connect_info = [int(v) for v in split_n(line_val[6:].strip(), 5)]
					list_obj_atom_connect = [list_obj_atoms[list_obj_atom_indexes.index(atom_index)] for atom_index in connect_info if atom_index in list_obj_atom_indexes]
					list_obj_atom_connect[0].set_connect(list_obj_atom_connect[1:])

		# 座標を ndarray に変換する
		self._coords = np.array(coords)
		return self


	def translate(self, shift_vector):
		"""
		全体の座標を平行移動するメソッド

		Args:
			shift_vector (ndarray): [x(float), y(float), z(float)]

		Returns:
			self
		"""
		self._coords += shift_vector
		return self


	def rotate(self, phi, theta, psi):
		"""
		全体の座標を回転するメソッド

		Args:
			vector (ndarray): []

		Returns:
			self
		"""
		# r_phi = np.array([
		# 	[1.0, 0.0, 0.0],
		# 	[0.0,  np.cos(phi), np.sin(phi), 0.0],
		# 	[0.0, -np.sin(phi), np.cos(phi), 0.0]
		# ])
		# r_theta = np.array([
		# 	[np.cos(theta), 0.0, -np.sin(theta)],
		# 	[0.0, 1.0, 0.0],
		# 	[np.sin(theta), 0.0,  np.cos(theta)]
		# ])
		# r_psi = np.array([
		# 	[np.cos(psi), -np.sin(psi), 0.0],
		# 	[np.sin(psi),  np.cos(psi), 0.0],
		# 	[0.0, 0.0, 1.0]
		# ])
		# r = r_psi.dot(r_theta)
		# print(1000, r)
		# r = r.dot(r_phi)
		# print(2000, r)
		r = rotM([phi, theta, psi])
		r = r.T
		new_coords = []
		for c in self._coords:
			new_coords.append(np.dot(r, c))
		self._coords = np.array(new_coords)
		return self


	def remove(self, remove_object):
		"""
		所属オブジェクトを削除するメソッド

		Args:
			remove_object (Chain/Residue/Atom object): 削除オブジェクト

		Returns:
			self
		"""
		idx_remove = None
		if isinstance(remove_object, Chain):
			# Chain オブジェクトの場合
			new_chain = []
			for chain in self._chains:
				if chain is remove_object:
					idx_remove = [i for i, atom1 in enumerate(self.atoms) for residue in chain.residues for atom2 in residue.atoms if atom1 is atom2]
				else:
					new_chain.append(chain)

			self._chains = new_chain

		elif isinstance(remove_object, Residue):
			# Residue オブジェクトの場合
			for residue in self.residues:
				if residue is remove_object:
					idx_remove = [i for i, atom1 in enumerate(self.atoms) for atom2 in residue.atoms if atom1 is atom2]
					residue.chain.remove_residue(remove_object)
					break

		elif isinstance(remove_object, Atom):
			# Atom オブジェクトの場合
			for i, atom in enumerate(self.atoms):
				if atom is remove_object:
					atom.residue.remove_atom(remove_object)
					idx_remove = i
					break

		else:
			sys.stderr.write("ERROR: Unknown object.\n")
			sys.exit(1)

		self._coords = np.delete(self._coords, idx_remove, 0)
		return self


	def write_pdb(self, output_file):
		"""
		PDB ファイルとして出力するメソッド

		Args:
			output_file (str): PDB ファイル

		Returns:
			self
		"""
		idx_atom = 0
		with open(output_file, "w") as obj_output:
			for chain in self._chains:
				for residue in chain.residues:
					for atom in residue.atoms:
						obj_output.write("HETATM{0:>5} {1:4} {2:3} {3:>1}{4:>4}    {5[0]:>8.3f}{5[1]:>8.3f}{5[2]:>8.3f}{6}{7}          {8}{9}\n".format(
							atom.number,
							atom.name,
							residue.name,
							chain.name,
							residue.number,
							self._coords[idx_atom],
							atom.occupancy,
							atom.temp_factor,
							atom.element,
							atom.charge
						))
						idx_atom += 1
				obj_output.write("TER\n")
			if len(self.connect_info) != 0:
				for connect in self.connect_info:
					obj_output.write("CONECT" + "".join(["{0:>5}".format(v) for v in connect]) + "\n")
			obj_output.write("END\n")
		return self


class Chain:
	def __init__(self, name):
		self._name = None
		self._residues = []

		self.set_name(name)

	@property
	def name(self):
		return self._name

	@property
	def residues(self):
		return self._residues


	def set_name(self, name):
		"""
		Chain に名前を設定するメソッド

		Args:
			name (str): Chain 名

		Returns:
			self
		"""
		self._name = name
		return self


	def append_residue(self, obj_Residue):
		"""
		残基オブジェクトを追加するメソッド

		Args:
			obj_Residue (obj_Residue): 残基オブジェクト

		Returns:
			self
		"""
		self._residues.append(obj_Residue)
		return self


	def set_residues(self, list_obj_Residues):
		"""
		残基オブジェクトを設定するメソッド

		Args:
			list_obj_Residues (list): [obj_Residue, ...]

		Returns:
			self
		"""
		self._residues = list_obj_Residues
		return self


	def remove_residue(self, remove_Residue):
		"""
		Residue オブジェクトを削除するメソッド

		Args:
			remove_Residue (obj_Residue): Residue オブジェクト

		Returns:
			self
		"""
		self._residues = [residue for residue in self._residues if residue is not remove_Residue]
		return self


class Residue:
	"""
	残基クラス (情報のみ)
	"""
	def __init__(self, name, number):
		self._name = None
		self._number = None
		self._atoms = []
		self._chain = None

		self.set_name(name)
		self.set_number(number)


	@property
	def name(self):
		return self._name

	@property
	def number(self):
		return self._number

	@property
	def unique(self):
		return "{0}-{1}-{2}".format(self._name, self._chain.name, self._number)

	@property
	def atoms(self):
		return self._atoms

	@property
	def chain(self):
		return self._chain


	def set_name(self, name):
		"""
		残基名を設定するメソッド

		Args:
			name (str): 残基名

		Returns:
			self
		"""
		self._name = name
		return self


	def set_number(self, number):
		"""
		残基インデックスを設定するメソッド

		Args:
			number (int): 残基インデックス

		Returns:
			self
		"""
		self._number = number
		return self


	def set_chain(self, obj_Chain):
		"""
		所属する Chain オブジェクトを設定するメソッド

		Args:
			obj_Chain (obj_Chain): 所属する Chain オブジェクト

		Returns:
			self
		"""
		self._chain = obj_Chain
		return self


	def append_atom(self, obj_atom):
		"""
		所属原子を追加するメソッド

		Args:
			obj_atom (obj_Atom): Atom オブジェクト

		Returns:
			self
		"""
		self._atoms.append(obj_atom)
		return self


	def set_atoms(self, list_obj_atoms):
		"""
		所属原子を設定するメソッド

		Args:
			list_obj_atoms (list): [obj_Atom, ...]

		Returns:
			self
		"""
		self._atoms = list_obj_atoms
		return self


	def remove_atom(self, remove_Atom):
		"""
		Atom オブジェクトを削除するメソッド

		Args:
			remove_Atom (obj_Atom): Atom オブジェクト

		Returns:
			self
		"""
		self._atoms = [atom for atom in self._atoms if atom is not remove_Atom]
		return self


class Atom:
	def __init__(self, name, number):
		self._name = None
		self._number = None
		self._chain = None
		self._residue = None
		self._occupancy = None
		self._temp_factor = None
		self._element = None
		self._charge = None
		self._connect = []

		self.set_name(name)
		self.set_number(number)

	@property
	def name(self):
		return self._name

	@property
	def number(self):
		return self._number

	@property
	def chain(self):
		return self._chain

	@property
	def residue(self):
		return self._residue

	@property
	def occupancy(self):
		occupancy = self._occupancy
		if occupancy is None:
			occupancy = 1.0
		return "{0:>6.2f}".format(occupancy)

	@property
	def temp_factor(self):
		temp_factor = self._temp_factor
		if temp_factor is None:
			temp_factor = 0.0
		return "{0:>6.2f}".format(temp_factor)

	@property
	def element(self):
		element = self._element
		if element is None:
			element = ""
		return "{0:>2}".format(element)

	@property
	def charge(self):
		charge = self._charge
		if charge is None:
			charge = ""
		return "{0:>2}".format(charge)

	@property
	def connect(self):
		return self._connect

	@property
	def connect_info(self):
		if len(self._connect) != 0:
			return [self.number] + [v.number for v in self._connect]
		return None


	def set_name(self, name):
		"""
		原子名を設定するメソッド

		Args:
			name (str): 原子名

		Returns:
			self
		"""
		self._name = name
		return self


	def set_number(self, number):
		"""
		原子インデックスを設定するメソッド

		Args:
			number (int): 原子インデックス

		Returns:
			self
		"""
		self._number = number
		return self


	def set_chain(self, obj_Chain):
		"""
		所属する Chain オブジェクトを設定するメソッド

		Args:
			obj_Chain (obj_Chain): 所属する Chain オブジェクト

		Returns:
			self
		"""
		self._chain = obj_Chain
		return self


	def set_residue(self, obj_Residue):
		"""
		所属する Residue オブジェクトを設定するメソッド

		Args:
			obj_Residue (obj_Residue): 所属する Residue オブジェクト

		Returns:
			self
		"""
		self._residue = obj_Residue
		return self


	def set_occupancy(self, occupancy):
		"""
		occupancy を設定するメソッド

		Args:
			occupancy (float): occupancy

		Returns:
			self
		"""
		self._occupancy = occupancy
		return self


	def set_temp_factor(self, temp_factor):
		"""
		temp factor を設定するメソッド

		Args:
			temp_factor (float): temp factor

		Returns:
			self
		"""
		self._temp_factor = temp_factor
		return self


	def set_element(self, element):
		"""
		元素を設定するメソッド

		Args:
			element (str): 元素名

		Returns:
			self
		"""
		self._element = element
		return self


	def set_charge(self, charge):
		"""
		電荷を設定するメソッド

		Args:
			charge (str): 電荷

		Returns:
			self
		"""
		self._charge = charge
		return self


	def set_connect(self, list_connect_obj_Atom):
		"""
		接続情報を設定するメソッド

		Args:
			list_connect_objAtom (list): [obj_Atom, ...]

		Returns:
			self
		"""
		self._connect = list_connect_obj_Atom
		return self
