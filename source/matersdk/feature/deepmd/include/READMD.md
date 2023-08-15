# 1. In brief
## 1.1. Note
$$r_{ji} = r_j - r_i$$

$$x = \frac{r-r_s}{r_c-r_s}$$

## 1.2. $s(r)$ -- `smooth function`
$$
s(r) = 
\begin{aligned}
& \frac{1}{r}, \quad r < r_s    \\
& \frac{1}{r} \cdot [x^3(-6x^2+15x-10) + 1], \quad r_s \leq r < r_c \\
& 0, \quad r \geq r_c
\end{aligned}
$$

## 1.3. $\tilde{R}$
$$
\tilde{R} = (s(r), \frac{s(r)x_{ji}}{r}, \frac{s(r)y_{ji}}{r}, \frac{s(r)z_{ji}}{r})
$$



# 2. 便于求导的形式
## 2.1. $s(r)$ -- represented with `switching function`
$$
switch\_func(r) = 
\begin{aligned}
& 1, \quad r<r_s \\
& x^3(-6x^2+15x-10) + 1, \quad r_s \leq r < r_c \\
& 0, \quad r > r_c
\end{aligned}
$$

## 2.2. $\tilde{R}$ -- represented with `switching function`
$$
\tilde{R} = (\frac{switch\_func(r)}{r}, \frac{switch\_func(r)x_{ji}}{r^2_{ji}}, \frac{switch\_func(r)y_{ji}}{r^2_{ji}}, \frac{switch\_func(r)z_{ji}}{r^2_{ji}})
$$