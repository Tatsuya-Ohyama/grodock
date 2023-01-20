# grodock

## 概要
Gromacs 構造ファイルを結合するプログラム





## 使用方法
```sh
$ grodock.py [-h] -g GRO_RECEPTOR_FILE.gro -p LIGAND.pdb -a ACPYPE_LIGAND.gro [-m MAP.txt] -o OUTPUT.gro [-O]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して、終了する。
* `-g GRO_RECEPTOR_FILE.gro`
	: `gmx pdb2gmx` で出力した生体分子構造ファイル
* `-p LIGAND.pdb`
	: リガンド構造ファイル
* `-a ACPYPE_LIGAND.gro`
	: `acpype` で出力したリガンド構造ファイル
* `-m MAP.txt`
	: ATOMTYPE のマッピングファイル (コンマか、改行で区切る)
	: 例: `PDB_ATOM_TYPE: GRO_ATOM_TYPE, ...`
	: マッピングファイルが指定されない場合、構造ファイルの原子の順番でマッピングされる。
* `-o OUTPUT.gro`
	: ドッキング構造
* `-O`
	: プロンプトを出さずに上書きする。


## 動作要件
* Python3
	* numpy


## License
The MIT License (MIT)

Copyright (c) 2021 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 2.4 (2023-01-20)
* マッピングファイルが指定されていない場合の原子数の不一致に対してエラーメッセージを表示するようにした。

### Ver. 2.3.2 / Ver. 2.3.3 (2022-06-24)
* 使用していないモジュールを削除した。

### Ver. 2.3 (2021-08-31)
* モジュール改変に伴うバグを修正した。

### Ver. 2.2 (2021-08-02)
* 公開した。
