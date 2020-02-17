## 微机保护算法实现

### 保护动作过程

1. 非故障状态下，采用基于2采样点的**复数**算法，计算电压和电流的复数形式

该算法默认采样得到的值是稳态基频(50Hz)正弦量，因此直接用正弦量的 **相量形式** 进行推导，用两个采样点计算正弦量的两个特征值————幅值和初相。因此，该算法得到的复数的模是最大值，而非有效值。

2. 通过突变量的形式开启故障状态下主程序（中断）

3. 进入故障程序(EeFault)

由于故障时，波形发生畸变，因而采样得到的值已经不再是稳态基频正弦量，而是富含多种频率的量，因此不能使用 `第1步` 中的 **2点复数算法**。这里采样的是全波傅立叶算法。即，进入故障程序后，等待一个周期的采样值，利用傅立叶算法，根据这一个周期 N个点 的采样值，计算各个频率下正弦量的 实部和虚部，进而得到不同频率正弦量幅值和初相角。

4. 利用`第3步`得到的基频(50Hz)电压和电流复数值，进而进行不同保护的判据计算。

对于距离保护，根据故障状态下得到的基频电压和电流的复数形式，计算故障时的阻抗，又称测量阻抗。

对于差动保护，根据本测故障状态下得到的电流复数值 及 对测发送来的电流复数值 进行差分计算。

5. 全波傅立叶公式的前提是：采样值是利用 直流分量及不同频率稳态正弦分量 叠加而成，当采样值中包含衰减分量时，则误差较大。可以对采样值进行FIR一阶差分滤波，过滤衰减直流分量，然后再送给故障程序，计算故障基频复数值。

### 已完成

1. 非故障状态下的电压、电流采样
2. 非故障状态下根据电压、电流复数值计算P、Q、R、X
3. 故障状态下 不同频率的 电压、电流采样（基于全波傅立叶）
4. 计算故障状态下不同频率的阻抗

### 未完成 
1. FIR一阶差分滤波
2. 具体距离保护的实现
3. 相位比较判据
4. 动态展示，电压、电流、保护动作开出

### 


