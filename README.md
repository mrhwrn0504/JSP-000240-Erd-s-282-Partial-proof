# 奇数贪心轨道：分子有界情形的部分结果与开放问题 / Odd greedy orbits with bounded numerators

**孙雨晨奖 · 编号 GSP-000240**

## 状态 / Status

本仓库**不证明完整猜想**。它记录的是：

- 已证明（初等，待复核）：A ≤ 4 时不存在无限有界奇数贪心轨道。
- 已证明（计算机辅助）：A ≤ 15 时不存在无限有界轨道，由两套独立实现核对（`code/trim.cpp`、`code/btrim.cpp`）。
- 已证明（计算机辅助）：八周期词 W 在素数 2、11、13 上 p-adic 可容许，并给出 14 步整数种子。
- 归约未完成：A = 16 归约到 W（来源待核对）。
- **未解决**：完整猜想；分子无界的强形式（与完整猜想等价）。

The full conjecture is **open**. See `paper/en/paper.tex` (English) and `paper/cn/` (中文) for the statement of every result and its status.

## 仓库结构 / Layout

```
paper/en/   英文论文源码（LaTeX, amsart）与 PDF
paper/cn/   中文论文 PDF 与 HTML 源文件
code/       修剪证书与辅助脚本
SHA256SUMS.txt  代码与论文文件的校验和
LICENSE     MIT
```

## 复现 / Reproduce

```bash
g++ -O2 -fopenmp code/trim.cpp -o trim && ./trim 15 720720 /dev/null   # 期望：存活集合为空
g++ -O2 -fopenmp code/btrim.cpp -o btrim && ./btrim 15 720720          # 期望：alive=0
```

（上面两条命令以论文附录所述参数为准。）

## 发表前待办 / Before submission

1. 证明或否定"切换点展开两两不同"引理（子贪心路线的缺口）。
2. 由独立人员重跑 A ≤ 15 的证书。
3. 逐条核对参考文献的原始来源。

## 引用 / Citation

请在引用时注明"部分结果，猜想未解决"。Cite as partial results; the conjecture remains open.
