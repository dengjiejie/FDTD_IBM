本程序采用IBM+FDTD模拟声场中圆柱散射的问题，进行如下设定：
  1、计算域为(x,y)=(-10,10)，以减少边界处反射波的影响
  2、将欧拉网格的间距设为△x=△y=0.02，欧拉网格点共有1001*1001个。即左下角的(-10,10)为第0行、0列个网格点。
  3、拉格朗日点共取145个。这样拉格朗日点间的弧线间距为0.021，与△x、△y大致匹配。
  4、初始条件下，u、v等值都设为0
  5、初始条件下，p值按照论文中给定进行赋初值，以模拟一个在(4,0)处产生高斯波，随后波源消失的情况。
  6、程序中，设置时间步长为0.005，总共10s，有2000个时间步

下面对每个类的功能进行简要介绍：

  1、主函数main所在的类为：TestProject类
  主要获取时间步信息，之后自动调用类的构造函数TestProject()进行计算
 
  2、ReadFile类
  读取时间步信息

  3、LagPoint类
  对拉格朗日点进行处理：
      writeLagFile()：生成拉格朗日点；在生成拉格朗日点的过程中，为了保持上下对称，将前73个拉格朗日点先生成，之后对称的生成后面的点。
      readLagFile()：读取点信息

   4、EulaPoint类
  对欧拉网格点进行处理：
      eulaGeneration()：生成欧拉网格

   5、SetCondition类
  设置初始条件和边界条件：
      initCondition()
      boundCondition()

   6、TimeOperation类
   进行每个时间步的处理：
      virForce()：计算附加力
      solvePUV()：用FDTD交替求解p、u、v
      saveXYP()：保存坐标点以及相应的p值，生成文件puv.dat
      saveValues_p()
      saveValues_u()
      saveValues_v()
      saveF()：保存附加力

   7、ParDefinition类
    存放整个程序所用到的变量值，在各个方法中需要用到时，直接通过类名.变量名调用即可。


    注意：1、本程序中欧拉点文件gridxyz.dat、拉格朗日点文件lp.dat、以及时间步文件in.dat的路径均为"src/cn/edu/tju/main/"

          2、程序运行出来的结果保存在"src/cn/edu/tju/main/values"文件夹中，在还未运行程序之前，需要将该文件夹中的内容清空，以免结果出错。由于跑出来的结果很大，所以只选取了部分文件放在程序里的values文件夹中。




