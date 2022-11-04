# pflow
<font color="steelblue" size="4">

A public layer for:
1. "pwmat" / "atom.config"
2. "poscar" / "vasp"
3. "cssr"
4. "json"
5. "xsf"
6. "mcsqs"
7. "prismatic"
8. "yaml"
9. "fleur-inpgen"

</font>

![pic_1](./test_data/pics/pic_1.png)

# 1. Usage
```shell
# 1. 安装
$ pip install -e .
# 2. 更改脚本中的 python 解释器路径
# 3. 添加环境变量
$ export PATH=<your_path>/pflow/click:$PATH
# 4. 格式转换
## 4.1. pwmat -> poscar
$ convertFormat.py  
The path of input file: : ./atom.config
The format of input file: : pwmat
The path of output file: : POSCAR
The format of output file: : poscar

## 4.2. pwmat -> mcsqs
$ convertFormat.py  
The path of input file: : ./atom.config
The format of input file: : pwmat
The path of output file: : rndstr.in
The format of output file: : mcsqs
```