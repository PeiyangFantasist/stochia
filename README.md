# Stochia · 随机空间

Stochia 是一个以“让随机变量变得可触摸”为目标的概率论桌面实验室。当前版本是可运行的 Qt 6 MVP，重点验证随机变量建模、分布可视化、变量变换、蒙特卡洛传播与项目保存这条核心工作流。

## 当前能力

- 13 个内置分布：Bernoulli、Binomial、Geometric、Poisson、Negative Binomial、Uniform、Normal、Exponential、Gamma、Beta、Chi-square、Student t、F
- 理论 PDF / PMF / CDF 自绘图形
- 理论分布与经验直方图对照
- 表达式：`+ - * / ^`、括号、`sin cos log exp sqrt abs`、`pi`、`e`
- 可变参数操作：`sum`、`product/prod`、`min`、`max`
- IID 批量构造：`iid_sum`、`iid_product`、`iid_min`、`iid_max`
- 数字特征算子：`E(X)` / `mean(X)` 与 `Var(X)`
- 变换变量，例如 `Z = X + Y`、`W = log(abs(X)) + Y^2`
- 解析分布推断：正态闭包、指数和到 Gamma、Uniform IID 和到缩放 Irwin–Hall、同尺度 Gamma 和、离散分布求和闭包
- 连续变量加减乘除、最值的理论密度/CDF 积分
- 递归依赖采样与循环依赖检查
- 均值、方差及模拟统计
- 已有变量原位编辑及下游理论自动重算
- 内置“帮助”菜单与表达式教程
- `.pvis` JSON 工程保存与加载
- 中英双语数学工作台界面

## 构建

要求：

- CMake 3.21+
- C++17 编译器
- Qt 6.4+，包含 Widgets

Windows（Qt Creator）：

1. 在 Qt Creator 中打开根目录的 `CMakeLists.txt`。
2. 选择安装了 Qt 6 的 Desktop Kit。
3. 构建并运行 `Stochia`。

命令行：

```powershell
cmake -S . -B build -DCMAKE_PREFIX_PATH="C:\Qt\6.8.0\msvc2022_64"
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

如果使用预编译的 Windows 发布包，解压后直接双击 `Stochia.exe`。包内已经包含 Qt 6.8.3 和 MSVC 运行库，不需要单独安装 Qt。

## 使用

启动后自带三个示例变量：

```text
X ~ Normal(mu=0, sigma=1)
Y ~ Exponential(lambda=2)
Z = X + Y
```

使用左下角“添加分布”创建基础随机变量。中下方表达式栏用于构造新变量。顶部视图菜单可切换理论密度、CDF 和模拟对照；右侧可设置采样次数并运行新的随机实验。

批量操作示例：

```text
S = sum(X,Y,Z)
P = product(X,Y,Z)
M = max(X,Y,Z)
T = iid_sum(X,10)
U = iid_max(X,20)
Q = (X-E(X))/sqrt(Var(X))
```

`X+X` 使用同一个随机样本；`iid_sum(X,2)` 会为每次结果生成两个相互独立的 `X` 副本。完整说明可从应用内“帮助”菜单查看。

## 工程结构

```text
src/core/Distribution.*      分布目录、密度/CDF、矩与采样
src/core/Expression.*        递归下降表达式解析器与 AST
src/core/ProbabilityModel.*  变量图、依赖检查与蒙特卡洛传播
src/ui/PlotWidget.*          无第三方依赖的 Qt 自绘图形
src/ui/MainWindow.*          三栏工作台与 .pvis 项目系统
tests/core_tests.cpp         数学核心回归测试
examples/*.pvis              可直接打开的示例项目
```

核心对象和依赖采样流程详见 [`docs/architecture.md`](docs/architecture.md)。
表达式和 IID 用法详见 [`docs/expression-guide.md`](docs/expression-guide.md)。

## 已知边界与下一步

当前版本会为可识别分布给出具体闭式，并为一般独立连续变量的加减乘除和最值建立理论积分。超出解析引擎覆盖范围的表达式仍自动回退到蒙特卡洛。自定义 PDF/CDF、多元联合分布、条件分布、LLN/CLT 动画和更完整的符号化简属于后续里程碑。
