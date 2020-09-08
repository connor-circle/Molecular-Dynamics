# Molecular-Dynamics
This project is about Molecular Dynamics code written with C language.

The force field uesd in the project is [CHARMM Force Field](https://www.charmm.org/) which is composed of bonded terms and non-bonded terms.

The bonded terms include [bond stretch](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Bond_term.c),[angle bend](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Angle_term.c),
[angle torsion](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Torsion_tem.c),[improper term](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Improper_term.c) and
[Urey-Bradley term](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Torsion_tem.c) inclued in angle torsion file

The non-bonded terms are comprised of [electrostatic interaction](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Ewald.c) and [van der Waals interaction](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Vdw_term.c)

The [verlet velocity](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Velocity.c) is applied in the project and The temperature is controlled by [Berendsen](https://github.com/Yuanruikang/Molecular-Dynamics/blob/master/Berendsen.c)
