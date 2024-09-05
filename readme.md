# MiniSPICE
> 该项目脱胎于[`parser`](https://github.com/rid-sun/practices/tree/main/parser)，历时一月使用C++特性（继承、封装、抽象、多态）重新改写整体架构，给出了一个<u>**简单理解SPICE开发流程**</u>很好的demo。:smile: :smile: :smile:

> 主要参考自：
> * https://www.sfu.ca/~ljilja/cnl/projects/Homotopy/  
> * https://ngspice.sourceforge.io/

---------------------------------------------------------------
## 新特性
- 以**节点为中心**生成KCL方程组，并使用**雅可比矩阵填充形式**而不是直流伴随模型。
- 使用**SUPERLU_DIST**和**KLU**作为底层解法器，分布式并行求解+串行求解。
- setup阶段实现**多种**矩阵格式在**全流程中仅一次转换**，减小引入求解器而带来的格式转换的开销。
- 支持**多种直流分析算法**，传统NR迭代算法、伪瞬态分析算法。
- 支持结果在前端**图形化输出**。

-------------------------------------------------------------

## TODO
### 功能模块
- [ ] 弧长法的同伦分析支持
- [ ] 时间区域的瞬态分析支持
- [ ] 前端界面结果绘制的Python支持
### 求解器
- [ ] CUDA加速的左视**MiniSPICE—LU**研发

-----------------------------------------------------

## 环境依赖
* Intel HPC Toolkit
* Cmake
* Python
* ParMETIS

### 构建
为保证环境一致性，建议使用docker构建。本项目在windows环境下，全面通过测试，大概需要**15G的镜像存储空间 + 410s的项目构建时间**
> `windows + docker + vscode`开发

Linux系统上，使用`OpenMPI + GCC`也可以构建，但需要对项目中的[`src/solver/CMakeLists.txt`](./src/solver/CMakeLists.txt)中的依赖库进行修改，保证一致性。

下面简述`windows + docker`的项目构建流程
1. 安装`Docker`，官网下载，可自定义镜像存储路径
2. 下载本项目到本地
    ```bash
    git clone https://github.com/rid-sun/MiniSPICE.git
    ```
3. 启动`dokcer`，而后按下面指令依序执行
    ```bash
    cd MiniSPICE
    docker build . --file Dockerfile --tag minispice
    ------waitting------
    docker run --mount type=bind,source="%CD%",target=/root/minispice -it minispice
    cd /root/minispice/3rd_lib
    bash install.sh
    ------waitting------
    cd ../
    bash build.sh
    ```
4. 构建完毕

### 运行
```bash
cd bin
mpirun -np num ./minispice -f ../testcase/testcase1/Netlist1.txt -o netlist1
```
其中`num`时进程数；项目构建时，默认使用`KLU`作为解法器，可以通过`-DUSE_SUPERLU=ON`使用`SuperLU_DIST`。
