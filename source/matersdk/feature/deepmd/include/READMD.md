# 1. In brief
## 1.1. Note
$$r_{ji} = r_j - r_i$$

$$uu = \frac{r-r_s}{r_c-r_s}$$

## 1.2. $s(r)$ -- `smooth function`
$$
s(r) = 
\begin{aligned}
& \frac{1}{r}, \quad r < r_s    \\
& \frac{1}{r} \cdot [uu^3(-6uu^2+15uu-10) + 1], \quad r_s \leq r < r_c \\
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
switchFunc(r) = 
\begin{aligned}
& 1, \quad r < r_s    \\
& uu^3(-6uu^2+15uu-10) + 1, \quad r_s \leq r < r_c \\
& 0, \quad r \geq r_c
\end{aligned}
$$

## 2.2. $\tilde{R}$ -- represented with `switching function`
$$
\tilde{R} = (\frac{switchFunc(r)}{r}, \frac{switchFunc(r)x_{ji}}{r^2_{ji}}, \frac{switchFunc(r)y_{ji}}{r^2_{ji}}, \frac{switchFunc(r)z_{ji}}{r^2_{ji}})
$$