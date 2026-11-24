# RWDID 网站 · Minimal Calculator Edition

## 页面结构
- `index.html`：首页只保留 RWDID 标题、简短 Instructions、15组分问卷、进度和双评分结果。
- `methods.html`：两种评分公式及冻结的15组分参数表。
- `resources.html`：单一完整 Analysis Code 和必要的映射/参数下载。
- `input-guide.html`：简短填写说明。

首页不再放置三步说明卡片、完整参数表、Friendly13 色板展示、审稿回复或研究流程介绍。

## Friendly13 配色
网页的分类强调色严格使用 `cols4all::friendly13` 的13色。GitHub Actions 在部署时执行 `scripts/export_friendly13_palette.R`，生成并校验 `data/friendly13_palette.json`。没有生成有效色板时，部署检查会失败，不会静默换用其他调色板。

## 部署
1. 解压 ZIP，将其中**全部文件和文件夹**上传到 `Lynn209/RWDID-website` 仓库根目录并覆盖同名文件。
2. 确认根目录有 `index.html`、`assets/`、`data/`、`analysis_code/`、`scripts/` 和 `.github/workflows/pages.yml`。
3. GitHub 仓库进入 **Settings → Pages**，将 Source 设为 **GitHub Actions**。
4. 提交到 `main` 分支，进入 **Actions**，确认 `Build and deploy RWDID GitHub Pages` 的 build 和 deploy 均成功。
5. 打开 `https://lynn209.github.io/RWDID-website/`，使用 Ctrl+F5 强制刷新。

## 评分定义
- RWDID：15个二分类组分按训练集 β 系数加权，范围0–1。
- RWDID-C：15个连续比例组分等权求和，范围0–15。
- 网页只应用冻结参数，不重新拟合 Cox 或 RCS 模型。
- 空白输入不视为零摄入；不进行缺失值插补。

## 隐私与用途
问卷输入仅在浏览器本地处理，不上传至服务器。RWDID 与 RWDID-C 是研究指数，不是临床风险预测或诊断工具。
