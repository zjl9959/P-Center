# P-Center

## 仓库文件目录

    deploy/ 项目部署
      instance/ 存放算例
        analyze.py 用于分析算法或配置性能
        benchmark.py 用于批量测试
        config.txt 用于配置批量测试参数
    doc/ 相关文件
    PCenter/ 代码实现
        main.cc
        pcenter_solver.cc
        pcenter_solver.h

## 使用方法

将源代码编译成功后，在**deploy**文件夹下会生成__PCenter.*__可执行程序，在shell中使用命令行传参的方式调用，或者配置好config.txt文件后执行benchmark.py脚本进行批量测试。
注意：_需要提前将算例拷贝到deploy/instance/文件夹下_

### 命令行参数

+ -i: instance name
+ -p: facility number
+ -k: consider range
+ -t: time out seconds
+ -s: random seed value
+ -m: max iteration steps
+ 示例：
    > PCenter -i pr226.tsp -p 40 -t 3
