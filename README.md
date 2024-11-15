
## Language
- [English](#english)
- [中文](#中文)
---
# English

# Twinkle_CPU
The CPU version of [Twinkle](https://github.com/AsterLight0626/Twinkle), producing exactly the same results.

If you use Twinkle (or part of Twinkle) in your research work, we request citing our paper: Wang, S. et al., 2024, in prep.

If you incorporate Twinkle (or part of Twinkle) in your code, we request specifying it in the README file (or equivalent), including the GitHub link ( https://github.com/AsterLight0626/Twinkle ) and asking the user to cite our paper: Wang, S. et al., 2024, in prep.

# Tutorial

The program design is consistent with Twinkle.

It is recommended to use the simple mode of calculating a single source at a time, that is, `/src/calculation/demo_single_point.cpp`

```
Magnification = LittleStar.Mag(s, q, y1, y2, Rs);
```
where $s$ is the binary separation, $q$ is the mass ratio, $y_1$ is the horizontal coordinate of the source center, $y_2$ is the vertical coordinate of the source center, and $R_s$ is the radius of the source (normalized by the Einstein Radius).

---
# 中文

# Twinkle_CPU
[Twinkle](https://github.com/AsterLight0626/Twinkle) 程序的CPU版，给出与 Twinkle 完全一致的结果

如果您在研究工作中使用了Twinkle（或Twinkle的一部分），请引用我们的文章：Wang, S. et al., 2024, in prep.

如果您将Twinkle（或Twinkle的一部分）整合到您的代码中，请您在README文件（或等效文件）中指明这一点。内容包括GitHub链接 (https://github.com/AsterLight0626/Twinkle) ，并要求用户引用我们的论文：Wang, S. et al., 2024, in prep.


# 使用教程 

程序设计与 Twinkle 一致。

推荐使用单次计算一个源的简便模式，即 `/src/calculation/demo_single_point.cpp`

```
Magnification = LittleStar.Mag(s, q, y1, y2, Rs);
```
其中 $s$ 是透镜体之间的距离，$q$ 是透镜体质量比，$y1$ 是源中心的横坐标，$y2$ 是源中心的纵坐标，$Rs$ 是源的半径（由 Einstein Radius 归一化）。
