# Implementation-strategy-of-APT-MCMC-algorithm-in-R-environment
Jin 等人提出了一种自动并行回火马尔可夫链蒙特卡洛（Automatic Parallel Tempering Markov Chain Monte Carlo, APT-MCMC）算法 ，可以高效采样复杂的分布，同时将该算法实现为一个名叫“Nii-C”的代码，Nii-C 代码可以从 https://github.com/shengjin/nii-c.git 中获取。这个算法使得在处理如系外行星轨道参数推断这类复杂的高维贝叶斯推断问题时，显著缩短获得可靠后验分布所需的计算时间，提高了分析效率。不过，该代码是 C 语言编写的，无法运用到R环境中。因此想到把这个算法引入到R环境中，因此想了三个实现策略：1.纯R语言版本实现；2.R语言与C++语言混合版本实现；3.R外部接口版本实现。注意，这三种版本的实现都是基于Jin等人实现APT-MCMC算法所使用的Nii-C代码的逻辑。APT-MCMC算法更详细的介绍的参考文献是：Jin S, Jiang W, Wu D H. Automatic parallel tempering markov chain monte carlo with Nii-C[J]. The Astrophysical Journal Supplement Series, 2024, 274(1): 10.

## 各个版本的性能测试配置
基准实验配置也遵循 Jin 等人在其研究里的描述，包括：
1. 模型与数据：包含 15 个轨道参数的双行星系统模型，基于其指定的模拟天体测量数据集进行估计。参数的模拟真值参考其论文 Table 1。
2. Nii-C 运行参数：启动 8 条具有特定温度 β 值序列的并行回火链；设定总运行长度为 300 万次 MCMC 迭代，舍弃前 50 万次作为预烧期（burn-in）。一个不同之处在于，本研究在整个 300 万次迭代期间均激活了 Nii-C 的自动调整提议分布功能，而不是仅在前 100 万次迭代期间激活。
在以下三个版本中，关于模型的设置在user_logll.R和user_prior.R，MCMC迭代参数的设置在input.ini文件中。同时请注意，任何C语言文件都是 https://github.com/shengjin/nii-c.git 中获取的，即Nii-C代码，由Jin等人所创作。

## 纯R语言版本
纯R语言版本的APT-MCMC代码在 pure R language version 文件夹里，main.R是主文件。通过 Rstudio 运行 main.R 文件，就开始了对这个模型15个参数的估计。
需注意的是，这时的版本拟合出的参数有部分参数不收敛，并且计算速度慢，还需要进一步改进。

## R与C++混合版本
R与C++混合版本的APT-MCMC代码在 R and C++ mixed version 文件夹里，main.R是主文件，cppfiles是将C++代码包装成了一个库(在 main.R 中可以看到library(cppfiles))。通过 Rstudio 运行 main.R 文件，就开始了对这个模型15个参数的估计。
需注意的是，这时的版本得到的参数估计结果很差，几乎没有一个拟合的好，因此说明代码仍需要进行修改。

## R外部接口版本
R外部接口版本的APT-MCMC代码在 R External Interface version 文件夹里，首先要对这些C语言文件进行编译，得到main.exe，然后再运行runc.R，就可以实现在R语言环境下运行APT-MCMC算法的代码了，其实主要思路只是在R语言环境中来运行编译好的C语言程序。
