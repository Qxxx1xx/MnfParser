# 模态中性文件（mnf）数据组织结构

模态中性文件解析功能模块

基于模态分析与综合法的柔性体建模方法，首先需要获得柔性体的模态信息模态中性文件解析功能模块针对有限元分析软件所生成的模态中性文件（MNF)进行解析，获取柔性体质量矩阵、惯量矩阵、模态矩阵等信息。MNF文件将不同类型的数据分块组织起来，共分为29个数据块，但并不是每个MNF文件都需要完全存储所有的数据块。事实上，在这29组数据中，只有一部分对于柔性体模型是必不可少的，其它数据中，一部分数据可以根据对模型的不同需要在生成MNF文件时选择输出，另有一部分数据如惯性时不变矩阵等可以根据已有的数据计算得到。

**MNF文件数据组织结构及各数据块含义如下：**

- Version Code——MNF文件版本号（必含）
- Header——文件头信息，内容包括文件名字、创建日期等（必含）
- Content Summary ——本文件所含内容纲要，即该文件含有哪些数据块（必含）
- Nodal Coordinates ——节点坐标数据，含节点编号（必含）
- Global body Property——柔性体整体属性，包括质量、质心位置、惯性张量三部分数据（必含）
- Eigenvalues——特征值相关数据，包括特征值、频率、广义质量、广义刚度（必含）
- Modes——模态数据（必含）
- Nodal Mass——节点质量数据（必含）
- Nodal Inertias——节点惯量数据，只有对界面点才拥有该部分数据（可选）
- Units——单位信息，一般为四种单位：质量、长度、时间、力（必含）
- Generalized Stiffness ——广义刚度，与Eigenvalues数据块中的相同（必含）
- Generalized Mass ——广义质量，与Eigenvalues数据块中的相同（必含）
- Element Faces——有限元单元数据（必含）
- Generalized Damping ——广义阻尼（可选）
- Mode Shape Transformation——模态变换（可选）
- Interface Nodes——界面点（可选）
- Modal Stress——模态应力（可选）
- Invariants 1-9——九个惯性时不变矩阵（可选）
- Modal Preload——模态预载（可选）
- Modal Loads——模态载荷（可选）
- Modal Strain——模态应力（可选）

**MNF文件被设计为一个跨平台的二进制文件，在数据存储上有固定的编码格式，总结如下：**

- 整形数据使用4个字节，浮点型数据一律采用双精度8个字节，没有单精度数据；
- 对于字符型数据，每个字符用4个字节编码：实际上字符型数据在MNF文件中很少，主要用于存储模型单位制；
- 数据字节顺序采用大端法，即对于多字节数据，最高有效字节在最前面。

**MNF文件编码格式**：

0. 数据块名称——距离文件头的字节数（十六进制）——数据说明

1. Version Code ——0——从文件头开始前4个字节，数值为版本号*10

2. Header——4——Title，从第0+4个字节后开始；Date，从第0x200+4字节后开始；Input file name，从第0×280+4字节后开始；Input file type，从第0x480+4字节后开始；Comments，

3. Content Summary—— ——

4. Nodal Coords——0xa10+4——前4个字节是节点数目，接下来的四个字节是坐标维数，一般为三，即每个点有三个坐标，然后是节点坐标数据的开始，每个节点坐标数据有28个字节，前4个字节（1个int）是节点编号，然后是三个double型的节点坐标，占24个字节

5. 

6. Global body——紧跟Nodal Coords开始——分三部分：mass(1)，Center of mass(3)，Inertia tensor(9),共13个double型数据，占13*8个字节

7. Eigenvalues——紧跟数据块6——数据块6结束后，紧跟一个4个字节的整数，值为模态阶数，然后是Eigenvalues数据，每一行由1个int和5个double组成，int是模态数，其后紧跟的5个double依次是：eigenvalue，Hertz，radians，Gen.Mass，Gen.Stiff

8. 节点编号——紧跟7——4个字节int是代表模态阶数，再4个字节int代表节点数目。随后就是一整块数据，记录了所有的节点编号，每个编号为一个int值，占4个字节

9. Modes——紧跟8——一个int代表模态阶数，然后是该阶模态所有数值，顺序排列，每个点有六个double型的值，但该部分没有标出节点编号，即：01 double double double double double double double double……依次排列下去然后是二阶的02 double double double double double double double double……等等

10. Nodal mass——紧跟9——1个int，代表节点数目，再4个字节，值为1，含 意不详，然后是正式的节点质量数据，格式为：1个int型节点编号，紧跟1个double型节点质量值

11. Nodal inertias——紧跟10——前四个字节是一个int数据，其值为节点数据包个数（节点数据包见下文），接下来的就是节点惯量值，其规律为节点编号（4字节）+ 两个double数据（8x2=16字节），计20个字节，称之为一个节点数据包

12. Units——紧跟11——Units：排列顺序是Mass Units, Length Units, Time Units, Force Units每种单位用4\*16=64个字节编码，共计4\*64个字节

13. Element Faces——紧跟12——在面数据开始之前有16个字节，最前8个字节含义不详，随后的8个字节中，前4个字节（int型）代表面数目，然后4个字节（int型），代表到该块数据结束时的int型数据个数，然后开始面数据，格式：1个int，值为每个面的顶点数目，然后紧跟相应数目的顶点编号数据，编号从0个开始

14. Mode Shape Transformation——紧跟13——第一个4个字节为int表示模态阶数；随后有模态阶数个int，字节数为(模态阶数*4)；最后是一个矩阵，大小是(模态阶数$\times$模态阶数)，数据类型为double，不知道是行优先还是列优先。

15. Interface Nodes——紧跟14——第一个数据表示界面节点个数，类型为int。然后紧跟相应数目的顶点编号数据，编号从0个开始。