# RWDID GitHub Pages — v4 UKB 问卷匹配版

本网站是静态浏览器计算器，用于计算冻结的15组分 RWDID 和 RWDID-C。

## 本次问卷修订
- 按照 Oxford WebQ 原始交互方式：先勾选昨天是否食用某项目，再选择对应数量。
- 未勾选的分类项目按0摄入处理；下拉框不再人为添加“0 / Not consumed”。
- 每个 Field ID 显示对应的UKB原始响应文字，而不是统一的数字范围。
- 82个活动输入正确区分为76个分类 Oxford WebQ 字段和6个连续型 UKB Showcase 总重量字段。
- 26102、26131、26133、26151、26152、26153为连续 g/day 字段，不是分类问卷选项。
- 对 Data-Coding 100001 的PDF与Showcase编码差异作了说明：PDF v1.5显示2和3+，但发布的Showcase编码只有200=2+，且冻结R代码将200转换为2.0；网站遵循已发布编码和冻结分析。
- 增加混合蔬菜/沙拉和寿司的防重复计入提示。
- 已包含经核对的 cols4all friendly13 13色文件。

## 部署
1. 解压ZIP。
2. 将文件夹中的全部内容上传至 `Lynn209/RWDID-website` 仓库根目录并覆盖对应文件。
3. 确认 `index.html`、`assets/`、`data/`、`docs/`、`scripts/` 和 `.github/workflows/pages.yml` 均位于根目录。
4. GitHub → Settings → Pages 保持 Source = GitHub Actions。
5. 提交至 main，在 Actions 中确认 `Build and deploy RWDID GitHub Pages` 成功。
6. 强制刷新网站。

详细的82个字段响应映射见 `docs/ukb_questionnaire_response_mapping.md`。
