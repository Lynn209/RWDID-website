# RWDID GitHub Pages 网页计算器

这是一个纯静态、浏览器本地运行的 RWDID / RWDID-C 研究计算器。网页部署方式参考 `Lynn209/betaSTEER_web`：`index.html` 放在仓库根目录，使用 GitHub Pages 的 **Deploy from a branch → main → /(root)**，不需要 Shiny、服务器、Node 构建或 GitHub Actions。

## 仓库结构

```text
index.html
input-guide.html
assets/
  app.js
  scoring.js
  styles.css
data/
  评分参数、问卷映射和 friendly13 JSON
analysis_code/
  RWDID_complete_analysis.R
scripts/
  export_friendly13_palette.R
docs/
.nojekyll
README.md
README_CN.md
```

## 上传与部署

1. 新建 GitHub 公共仓库，例如 `RWDID-website`。
2. 将本文件夹内的**全部内容**上传到仓库根目录，确保 `index.html` 直接位于根目录。
3. 提交到 `main` 分支。
4. 进入 **Settings → Pages**。
5. Build and deployment 的 Source 选择 **Deploy from a branch**。
6. Branch 选择 `main`，Folder 选择 `/(root)`，点击 Save。
7. GitHub 会显示网页链接，格式类似 `https://你的用户名.github.io/RWDID-website/`。

## 上传前生成精确 friendly13 配色

本项目不再写入未经核实的近似 HEX 颜色。上传前请在网站文件夹根目录执行一次：

```bash
Rscript scripts/export_friendly13_palette.R
```

该脚本调用 `cols4all::c4a(palette = "friendly13", n = 13, type = "cat", format = "hex")`，将包返回的13个原始 HEX 值和包版本写入 `data/friendly13_palette.json`。将生成后的 JSON 一并上传。网页会校验 palette 名称、13色数量、HEX 格式及唯一性；不会静默替换成其他颜色。未生成时网页会明确显示警告，而不是声称已经使用 friendly13。

## 本地预览

浏览器通过 `fetch()` 读取 JSON 文件，因此不要直接双击 `index.html`。在网站根目录运行：

```bash
python -m http.server 8000
```

然后打开 `http://localhost:8000/`。

## 评分内容

- 15 个固定饮食组分。
- 85 项 Table S1 原始条目审计，其中网页采用 82 项有效问卷输入。
- RWDID：15 个二分类组分 × 训练集 Cox β 系数，按 β 总和归一化，范围 0–1。
- RWDID-C：15 个连续比例组分等权求和，范围 0–15。
- 不进行插补；空白不自动当作零摄入。
- 网页只应用冻结的评分参数，不重新拟合 Cox 或 RCS 模型。

## 单一完整分析代码

`analysis_code/RWDID_complete_analysis.R` 是唯一公开的完整分析代码文件，包含11个按顺序排列的代码模块：饮食数据处理、70/30队列划分、Model 2最低风险参考类别、15组分条目/份量审计、RCS锚点构建与区间协调、最终二分类评分、连续型RWDID-C、训练集β加权及测试集验证。

注意：公开仓库不包含UK Biobank个体级数据及项目预处理数据对象；所需输入已在代码文件开头列明。整合时仅将队列划分脚本的输出路径统一到下游代码使用的 `results/Step_03/ukb_70_30_split`，随机划分算法和seed未改变。当前Step_04_2_1 V1.3与Step_04_2_2 V1.4的输出路径相互匹配，训练代码会导出测试集验证所需的冻结RCS knots。

## 数据保护与用途

公开仓库不包含 UKB 个体级数据、EID 或个体训练评分。网页计算在浏览器本地进行，不主动上传输入。仅供研究使用，不是临床风险预测或诊断工具。
