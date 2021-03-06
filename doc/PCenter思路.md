# PCentre思路

## 问题定义

### 输入

+ 无向图$G(v,(e,w))$
+ 服务节点数目P

### 约束

+ 每个节点由一个离他最近的中心节点服务
+ 所有节点都要有一个服务节点

### 目标

+ 最小化最长服务边

## 算法框架

### 产生初始解

+ 初始解通过逐次挑选中心节点来构造
1. 第一次随便挑选一个节点作为服务节点$P_{1}$
2. 其次从最长服务边中随机挑选一个离用户节点前$k$近的点作为下一个中心节点
3. 重复1，2过程直到选择完所有中心节点

### 局部搜索邻域动作

+ 邻域动作：将最长服务边上的**服务节点**用**代替节点**代替

    > _为减小邻域结构，这里代替节点为最长服务边上离用户节点第k近的用户节点_
+ 邻域动作评估：
1. 找出k个代替节点$N_{k}$
2. 尝试增加一个代替节点做为服务节点
3. 对所有服务节点，尝试删除其中一个f
4. 如果2，3执行，包含以下节点n的服务边将会受到影响
    > 1. 满足条件$D_{kn}<D_{pn}$的节点，_p为n的当前服务节点_
    > 2. 服务边上含有f，失去服务节点的用户会选择其次近的服务节点
    > + 1中的节点只影响状态更新，不影响动作评估
+ 动作更新
1. 新的目标函数值为，老的目标函数值和新增目标函数值中的较大者
2. 对于被删去的服务节点，其服务节点为新增服务节点和离它次近的服务节点中最近节点
3. 对于新增的服务节点，其服务节点为自己，其次近的服务节点为以前最近的服务节点
4. 对于用户节点，其服务点是否更新看它离新增服务节点的距离是否大于当前服务边长

### 禁忌表

+ 对节点交换对动作进行禁忌，禁忌表中记录解禁时间

    > 解禁时间 = 当前时间 + 随机禁忌步长

+ 如果禁忌动作同时比最小的非禁忌动作和历史动作还要好，对该动作进行解禁

----

### 实现细节

+ 如何寻找最长服务边？
    > 1.扫描FDTable, 记录最长的服务边
    > [-]2.在更新目标函数值的时候记录最长服务边
