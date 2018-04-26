# guass-quadrature
这个程序可以用于计算任意给定代数精度和误差下高斯型求积公式的节点和相应的权重。我们的程序首先运用经典Runge-Kutta法求解Prüfer微分方程，得到了给定次数的Legendre多项式的零点，即为相应求积公式的节点。再通过求权重的公式计算出各个节点对应的权重。由于这个程序中使用了GMP包，因此理论上可以达到任意的精度。程序的运行方法如下：在"data.txt"文件的第一行输入所需的代数精度，第二行输入允许的最大误差，再运行程序，结果会输入到"answer.txt"文件中，第一列是节点，第二列是权重。
