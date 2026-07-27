# Stochia 指令输入教程

表达式在主界面中央底部输入。左边填写新变量名，右边填写表达式，然后点击“定义变量”。

## 基本规则

- 变量名区分大小写，`X` 与 `x` 是不同变量。
- 必须先创建表达式引用的变量。
- 函数名和逗号使用英文字符。
- 支持括号和 `+ - * / ^`。
- 常量为 `pi` 和 `e`；大写 `E` 可以作为变量名。

## 普通变换

```text
Z = X + Y
W = log(abs(X) + 0.001)
R = sqrt(X^2 + Y^2)
```

## 可变参数批量操作

```text
S = sum(X,Y,Z,W)
P = product(X,Y,Z)
P2 = prod(X,Y,Z)
L = min(X,Y,Z,W)
H = max(X,Y,Z,W)
```

`sum`、`product/prod`、`min`、`max` 至少需要两个参数，并可继续增加参数。

## 独立同分布操作

```text
S = iid_sum(X,10)
P = iid_product(X,5)
L = iid_min(X,20)
H = iid_max(X,20)
```

第一个参数必须是已经存在的随机变量，第二个参数必须是 1 到 1,000,000 之间的整数。

请特别注意：

```text
X + X
```

使用同一个 `X` 样本，因此等于 `2*X`。而：

```text
iid_sum(X,2)
```

会为每一次结果生成两个相互独立、但服从相同分布的 `X` 副本。

IID 的源变量也可以是已有的变换变量：

```text
Z = X + Y
T = iid_sum(Z,10)
```

每个 `Z` 副本会重新独立采样它依赖的整条变量链。

## 期望与方差

```text
M = E(X)
M2 = mean(X)
V = Var(X)
Z = (X-E(X))/sqrt(Var(X))
```

`E(X)` 与 `mean(X)` 等价，返回已经定义变量的理论期望；`Var(X)` 返回理论方差。它们在表达式中是确定常数，不会额外抽取随机样本，并会参与解析化简。例如正态变量的标准化表达式仍会被识别为 `Normal(0,1)`。

如果变量没有可用的理论分布，或者对应的矩不存在/不为有限数，系统会给出明确错误，不会在每次表达式求值时重复运行蒙特卡洛。

## 当前解析规则

以下情况会给出具体的闭式分布：

```text
X ~ Normal(mu,sigma)
iid_sum(X,n) ~ Normal(n*mu, sqrt(n)*sigma)
```

```text
X ~ Exponential(lambda)
iid_sum(X,n) ~ Gamma(shape=n, scale=1/lambda)
```

```text
X ~ Gamma(shape,scale)
iid_sum(X,n) ~ Gamma(n*shape, scale)
```

```text
X ~ Bernoulli(p)
iid_sum(X,n) ~ Binomial(n,p)
```

```text
U ~ Uniform(a,b)
iid_sum(U,n) ~ a*n + (b-a)*IrwinHall(n)
```

均匀分布 IID 和使用精确的分段多项式 PDF/CDF，而不是卷积嵌套或蒙特卡洛。当前交互式精确求值支持 `n <= 128`。

Poisson、具有相同 `p` 的 Binomial 和 Negative Binomial，以及 Chi-square 的 IID 和也会识别对应分布族。多个相互独立的正态变量相加或相减仍会识别为正态分布。同尺度的 Gamma/Exponential 变量相加会合并形状参数。

对于没有常见分布名称、但存在理论积分公式的两个独立连续变量，Stochia 会计算：

```text
Z = X + Y
f_Z(z) = ∫ f_X(x) f_Y(z-x) dx
```

```text
Z = X - Y
f_Z(z) = ∫ f_X(x) f_Y(x-z) dx
```

```text
Z = X * Y
f_Z(z) = ∫ f_X(x) f_Y(z/x) / |x| dx
```

```text
Z = X / Y
f_Z(z) = ∫ |y| f_X(zy) f_Y(y) dy
```

以及独立变量的最值：

```text
F_max(z) = ∏ F_i(z)
F_min(z) = 1 - ∏(1-F_i(z))
```

理论结果从“密度 / PMF”和“累积分布 CDF”查看；“模拟对照”使用蒙特卡洛样本验证理论结果。

## 编辑已有变量

选中左侧变量，点击“编辑所选”：

- 基础变量可以修改分布类型和参数。
- 变换变量可以修改表达式。
- 变量名保持不变，因此引用它的下游变量无需重建。
- 保存修改后，所有下游解析分布会自动重新计算。

## 定义域注意事项

- `log(x)` 要求 `x > 0`。
- `sqrt(x)` 要求 `x >= 0`。
- 除法要求除数不为零。
- 如果某些随机样本违反定义域，蒙特卡洛会给出数学错误。
