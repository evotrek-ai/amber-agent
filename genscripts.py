import sys
import os
from typing import Dict

leapin_for_biochemMol = """
source leaprc.protein.ff19SB
source leaprc.DNA.bsc1
source leaprc.RNA.OL3
source leaprc.gaff
source leaprc.water.tip3p 
loadamberparams frcmod.ionsjc_tip3p
set default PBRadii mbondi2

rec = loadpdb receptor.pdb
saveamberparm rec rec.gas.prmtop rec.gas.inpcrd

loadamberparams ligand.frcmod
lig = loadmol2 ligand_amber.mol2
saveamberparm lig lig.gas.prmtop lig.gas.inpcrd

com = combine {rec lig}
saveamberparm com com.gas.prmtop com.gas.inpcrd

solvateBox com TIP3PBOX 10
addions com Na+ 0 
addions com Cl- 0 
saveamberparm com wat.prmtop wat.inpcrd
quit
"""

leapin_for_chemMol = """
source leaprc.water.tip3p 
loadamberparams frcmod.ionsjc_tip3p
loadamberparams ligand.frcmod
set default PBRadii mbondi2

obj = loadmol2 ligand_amber.mol2
saveamberparm obj gas.prmtop gas.inpcrd

solvateBox obj TIP3PBOX 10
addions obj Na+ 0 
addions obj Cl- 0 
saveamberparm obj wat.prmtop wat.inpcrd
quit
"""

leapin_for_bioMol = """
source leaprc.protein.ff19SB
source leaprc.DNA.bsc1
source leaprc.RNA.OL3
source leaprc.water.tip3p 
loadamberparams frcmod.ionsjc_tip3p
set default PBRadii mbondi2

obj = loadpdb receptor.pdb
saveamberparm obj gas.prmtop gas.inpcrd

solvateBox obj TIP3PBOX 10
addions obj Na+ 0 
addions obj Cl- 0 
saveamberparm obj wat.prmtop wat.inpcrd
quit
"""

prepare_for_biochemMol = """
pdb4amber -i {receptor}.pdb -o receptor_tmp.pdb -d -y
awk '!($3 == "P" || $3 == "OP1" || $3 == "OP2")' receptor_tmp.pdb > receptor.pdb
pdb4amber -i {ligand}.pdb -o ligand.pdb -d -y
obabel ligand.pdb  -h --partialcharge gasteiger -O ligand.mol2
antechamber -i ligand.mol2 -fi mol2 -o ligand_amber.mol2 -fo mol2 -s 2 -c bcc
parmchk2 -i ligand_amber.mol2 -f mol2 -o ligand.frcmod
"""

prepare_for_bioMol = """
pdb4amber -i {receptor}.pdb -o receptor_tmp.pdb -d -y
awk '!($3 == "P" || $3 == "OP1" || $3 == "OP2")' receptor_tmp.pdb > receptor.pdb
"""

prepare_for_chemMol = """
pdb4amber -i {ligand}.pdb -o ligand.pdb -d -y
obabel ligand.pdb  -h --partialcharge gasteiger -O ligand.mol2
antechamber -i ligand.mol2 -fi mol2 -o ligand_amber.mol2 -fo mol2 -s 2 -c bcc
parmchk2 -i ligand_amber.mol2 -f mol2 -o ligand.frcmod
"""

minimization_template = """
minimization
    &cntrl
        imin=1,                    ! 执行能量最小化
        maxcyc={maxcyc},           ! 最小化最大步数
        ncyc={ncyc},               ! 陡最速下降法步数
        ntb=1,                     ! 常体积周期性边界条件
        ntp=0,                     ! 不进行压力控制
        cut={cut},                 ! 非键相互作用截断距离
        {restraint_mask_line}      ! 可选的原子约束
    /
"""

heating_template = """
heating
    &cntrl
        imin=0,                    ! 执行分子动力学
        ntx=1,                     ! 坐标格式，1=坐标
        irest=0,                   ! 不从之前的模拟继续
        nstlim={nstlim},           ! 模拟步数
        dt={dt},                   ! 时间步长(皮秒)
        ntc=2,                     ! SHAKE约束H-重原子键
        ntf=2,                     ! 不计算H-重原子键力
        ntb=1,                     ! 常体积周期性边界条件
        cut={cut},                 ! 非键相互作用截断距离
        ntt=3,                     ! Langevin恒温器
        gamma_ln=2.0,              ! 碰撞频率(皮秒^-1)
        tempi={tempi},             ! 初始温度(K)
        temp0={temp0},             ! 目标温度(K)
        {restraint_mask_line}      ! 可选的原子约束
        ntpr=1000,                 ! 能量信息输出频率
        ntwx=1000,                 ! 轨迹保存频率
        ntwr=10000,                ! 重启文件保存频率
    /
    &wt
        type='TEMP0',              ! 温度控制
        istep1=0,                  ! 起始步
        istep2={heat_steps},       ! 终止步
        value1={tempi},            ! 起始温度
        value2={temp0},            ! 终止温度
    /
    &wt type='END' /
"""

equilibration_template = """
equilibration
    &cntrl
        imin=0,                    ! 执行分子动力学
        ntx=5,                     ! 坐标和速度格式
        irest=1,                   ! 从之前的模拟继续
        nstlim={nstlim},           ! 模拟步数
        dt={dt},                   ! 时间步长(皮秒)
        ntc=2,                     ! SHAKE约束H-重原子键
        ntf=2,                     ! 不计算H-重原子键力
        ntb=2,                     ! 常压周期性边界条件
        ntp=1,                     ! 各向同性压力控制
        taup=2.0,                  ! 压力弛豫时间(皮秒)
        cut={cut},                 ! 非键相互作用截断距离
        ntt=3,                     ! Langevin恒温器
        gamma_ln=2.0,              ! 碰撞频率(皮秒^-1)
        temp0={temp0},             ! 目标温度(K)
        {restraint_mask_line}      ! 可选的原子约束
        ntpr=1000,                 ! 能量信息输出频率
        ntwx=1000,                 ! 轨迹保存频率
        ntwr=10000,                ! 重启文件保存频率
    /
"""

production_template = """
production
    &cntrl
        imin=0,                    ! 执行分子动力学
        ntx=5,                     ! 坐标和速度格式
        irest=1,                   ! 从之前的模拟继续
        nstlim={nstlim},           ! 模拟步数
        dt={dt},                   ! 时间步长(皮秒)
        ntc=2,                     ! SHAKE约束H-重原子键
        ntf=2,                     ! 不计算H-重原子键力
        ntb=2,                     ! 常压周期性边界条件
        ntp=1,                     ! 各向同性压力控制
        taup=2.0,                  ! 压力弛豫时间(皮秒)
        cut={cut},                 ! 非键相互作用截断距离
        ntt=3,                     ! Langevin恒温器
        gamma_ln=2.0,              ! 碰撞频率(皮秒^-1)
        temp0={temp0},             ! 目标温度(K)
        ntpr={ntpr},               ! 能量信息输出频率
        ntwx={ntwx},               ! 轨迹保存频率
        ntwr={ntwr},               ! 重启文件保存频率
        iwrap=1,                   ! 将分子重新包装入主单元
    /
"""

run_amber_template = """#!/bin/bash

set -e

#======================================================
# prepare inputs
{prepare}

#======================================================
# tleap
tleap -f leap.scr

#======================================================
# run md simulation on GPU

pmemd.cuda -O -i min.in -o min.out -p wat.prmtop -c wat.inpcrd -r min.rst

pmemd.cuda -O -i heat.in -o heat.out -p wat.prmtop -c min.rst -r heat.rst

pmemd.cuda -O -i equi.in -o equi.out -p wat.prmtop -c heat.rst -r equi.rst

pmemd.cuda -O -i prod.in -o prod.out -p wat.prmtop -c equi.rst -x prod.nc -r prod.rst
"""

def gen_tleap(receptor: str | None, ligand: str | None):
    has_ligand = ligand is not None
    if has_ligand and ligand.endswith('.pdb'):
        ligand, _ = os.path.splitext(ligand)
    has_receptor = receptor is not None
    if has_receptor and receptor.endswith('.pdb'):
        receptor, _ = os.path.splitext(receptor)
    params = {
        "receptor": receptor,
        "ligand": ligand,
    }
    prepare = ""
    leapin = ""
    if has_receptor and has_ligand:
        prepare = prepare_for_biochemMol.format(**params)
        leapin = leapin_for_biochemMol
    elif has_receptor:
        prepare = prepare_for_biochemMol.format(**params)
        leapin = leapin_for_bioMol
    elif has_ligand:
        prepare = prepare_for_chemMol.format(**params)
        leapin = leapin_for_chemMol

    # generate leap.src
    with open('leap.scr', 'w') as file:
        file.write(leapin)
        
    return prepare
        

def gen_minimization(parameters: Dict | None):
    if parameters is None:
        parameters = {}
    default_params = {
        "maxcyc": 5000,
        "ncyc": 1000,
        "cut": 10.0,
        "restraint_mask_line": ""
    }
    params = {**default_params, **parameters}
    script = minimization_template.format(**params)
    with open('min.in', 'w') as file:
        file.write(script)
        
def gen_heating(parameters: Dict | None):
    if parameters is None:
        parameters = {}
    default_params = {
        "nstlim": 25000,
        "dt": 0.002,
        "cut": 10.0,
        "tempi": 0.0,
        "temp0": 300.0,
        "heat_steps": 25000,
        "restraint_mask_line": ""
    }
    params = {**default_params, **parameters}
    script = heating_template.format(**params)
    with open('heat.in', 'w') as file:
        file.write(script)
        

        
def gen_equilibration(parameters: Dict | None):
    if parameters is None:
        parameters = {}
    default_params = {
        "nstlim": 250000,
        "dt": 0.002,
        "cut": 10.0,
        "temp0": 300.0,
        "restraint_mask_line": ""
    }
    params = {**default_params, **parameters}
    script = equilibration_template.format(**params)
    with open('equi.in', 'w') as file:
        file.write(script)
        
def gen_production(parameters: Dict | None):
    if parameters is None:
        parameters = {}
    default_params = {
        "nstlim": 5000000,
        "dt": 0.002,
        "cut": 10.0,
        "temp0": 300.0,
        "ntpr": 5000,
        "ntwx": 5000,
        "ntwr": 50000
    }
    params = {**default_params, **parameters}
    script = production_template.format(**params)
    with open('prod.in', 'w') as file:
        file.write(script)
        

def gen_scripts(ligand: str | None, min_params: Dict | None, heat_params: Dict | None, equi_params: Dict | None, prod_params: Dict | None):
    receptor = 'combined.nosmall.pdb'
    if not os.path.isfile(receptor):
        receptor = None

    prepare = gen_tleap(receptor, ligand)
    gen_minimization(min_params)
    gen_heating(heat_params)
    gen_equilibration(equi_params)
    gen_production(prod_params)
    script = run_amber_template.format(prepare=prepare)
    with open('run_amber.sh', 'w') as file:
        file.write(script)


if __name__ == "__main__":
    workdir = sys.argv[1]
    ligand = sys.argv[2]
    os.chdir(workdir)
    gen_scripts(ligand, None, {"nstlim": 5000}, {"nstlim": 5000}, {"nstlim": 5000})