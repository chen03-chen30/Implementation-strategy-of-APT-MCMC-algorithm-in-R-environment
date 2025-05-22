# Implementation-strategy-of-APT-MCMC-algorithm-in-R-environment
Jin 等人提出了一种自动并行回火马尔可夫链蒙特卡洛（Automatic Parallel Tempering Markov Chain Monte Carlo, APT-MCMC）算法 ，可以高效采样复杂的分布，同时将该算法实现为一个名叫“Nii-C”的代码，Nii-C 代码可以从 https://github.com/shengjin/nii-c.git 中获取。这个算法使得在处理如系外行星轨道参数推断这类复杂的高维贝叶斯推断问题时，显著缩短获得可靠后验分布所需的计算时间，提高了分析效率。不过，该代码是 C 语言编写的，无法运用到R环境中。因此想到把这个算法引入到R环境中，因此想了三个实现策略：1.纯R语言版本实现；2.R语言与C++语言混合版本实现；3.R外部接口版本实现。注意，这三种版本的实现都是基于Jin等人实现APT-MCMC算法所使用的Nii-C代码的逻辑。APT-MCMC算法更详细的介绍的参考文献是：Jin S, Jiang W, Wu D H. Automatic parallel tempering markov chain monte carlo with Nii-C[J]. The Astrophysical Journal Supplement Series, 2024, 274(1): 10.

## 纯R语言版本
纯R语言版本的
