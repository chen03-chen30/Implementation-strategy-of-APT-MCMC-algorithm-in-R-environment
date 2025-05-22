#CORRECT

#通过Box-Muller变换方法，将两个独立的均匀分布随机数转换为一个标准正态分布随机数
r8_normal_01 <- function(){
  r8_pi <- 3.141592653589793
  r <- runif(2)
  x <- sqrt(-2 * log(r[1])) * cos(2 * r8_pi * r[2])
  return(x)
}

r8_normal_ab <- function(a, b){
  value <- a + b * r8_normal_01()
  return(value)
}

r8_logunif_ab <- function(a, b){
  r8_unif <- runif(1)
  value <- a * exp(r8_unif * log(b/a))
  return(value)
}

#输出双精度a，b，得到在a～b之间的均匀分布的随机数
r8_unif_ab <- function(a, b){
  r8_unif <- runif(1)
  value <- r8_unif * (b - a) + a
  return(value)
}

#输入整数a，b，得到在a～b之间的均匀分布的随机数
i4_unif_ab <- function(a, b){
  r8_unif <- runif(1)
  value <- round(r8_unif * (b - a)) + a
  return(value)
}

#输入整数a，返回一个在0～a之间的随机整数
i4_unif_0a <- function(a){
  r8_unif <- runif(1)
  value <- floor(r8_unif * a)
  return(value)
}

normal_pdf <- function(x, mean, sd){
  inv_sqrt_2pi <- 0.3989422804014327
  a <- (x - mean)/sd
  return(inv_sqrt_2pi / sd * exp(-0.5 * a * a))
}


