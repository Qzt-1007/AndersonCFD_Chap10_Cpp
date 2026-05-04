# Anderson's CFD Chap10(C++)

这个仓库是安德森《计算流体力学基础及其应用》第10章算例的C++实现。

主程序为main.cpp，核心的麦考马克算法通过Maccormack.h实现。

对计算粘性应力和热流量的函数作出了优化：在边界处采用线性插值方法（书中给出的算法是一阶差分）

## 计算结果展示

![物面压力分布](P_on_plate.svg)

![平板后缘压力剖面](P_profile_back.svg)

![平板后缘温度剖面](T_profile_back.svg)

![平板后缘水平速度剖面](U_profile_back.svg)

上述图像与作者给出的参考图像吻合。
