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


## V6：面包卷类型联动修复
Field 101160（面包卷摄入量）现在与 UKB Field 20094（面包卷类型，多选，Data-Coding 81）联动，提供 White、Mixed、Wholemeal、Seeded、Other 五种选项。严格复现原始 Step_02_1_1：只有当 Field 20094 恰好仅选择 White（代码1）时，Field 101160 的数量才计入最终白面包卷伪组分；其他单一类型或多类型组合均计为该伪组分的0。Field 20094 只是联动问卷字段，不是新增评分项。82个评分输入、15个组分、评分区间、β系数、RWDID-C参数和公式全部保持冻结。详见 `docs/ukb_questionnaire_response_mapping.md` 和 `CHANGELOG_UKB_QUESTIONNAIRE_FIX.md`。


## V7.0.0 面包卷类型联动修复

V7 将精制谷物组中的面包卷项目显示为“Bread roll intake”，并始终明确显示 UKB Field 20094 的五个多选类型：White、Mixed、Wholemeal、Seeded、Other。选择至少一种类型后，才显示 Field 101160 的数量选项。

评分严格复现冻结 R 代码：仅当 Field 20094 恰好单选 White（code 1）时，Field 101160 的数量才计入白面包衍生项目；其他单选及任何多选组合均不计入。Field 20094 不是新增评分项目。

本版本为 CSS、JavaScript 和 JSON 请求增加 `v=7.0.0` 缓存版本标记，并在网页顶部显示 `V7.0.0`，便于核实 GitHub Pages 实际部署版本。请将 ZIP 内文件直接上传到仓库根目录，不要把整个外层文件夹作为二级目录上传。保留仓库已有的 GitHub Pages workflow。
