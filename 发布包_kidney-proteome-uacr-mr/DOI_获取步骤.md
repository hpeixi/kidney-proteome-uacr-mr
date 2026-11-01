# 获取 GitHub + Zenodo 的仓库地址与 DOI

当前状态：发布包已就绪（`repo/` 目录），但**仓库地址与 DOI 尚未生成**，因为这两项必须由你的 GitHub 与 Zenodo 账号完成，我无法代你注册或授权。

## 一、建 GitHub 仓库（约 10 分钟）

1. 登录 GitHub，新建仓库，建议命名为 `kidney-proteome-uacr-mr`，设为 Public。
2. 在本机 `repo/` 目录下执行：

```bash
git init
git add .
git commit -m "v1.0.0: cis-pQTL MR and colocalization for albuminuria"
git branch -M main
git remote add origin https://github.com/<你的用户名>/kidney-proteome-uacr-mr.git
git push -u origin main
```

3. 把仓库地址替换进 `metadata/CITATION.cff` 的 `repository-code` 字段。

## 二、连 Zenodo 取 DOI（约 10 分钟）

1. 用 GitHub 账号登录 https://zenodo.org ，进入 Settings → GitHub，授权并打开该仓库的开关。
2. 回到 GitHub，创建 Release，例如 tag `v1.0.0`，标题与说明随意。
3. 发布后 Zenodo 会自动归档并生成 DOI（形如 `10.5281/zenodo.1234567`）。
4. 把 DOI 写回 `metadata/CITATION.cff` 的 `doi` 字段与 `metadata/DOI.txt`，并重新提交一次。

## 三、需要同步到稿件的两处

1. Declarations 的“Availability of data and materials”改为：
   `Analysis code and derived result tables are available at https://github.com/<用户名>/kidney-proteome-uacr-mr (DOI: 10.5281/zenodo.XXXXXXX).`
2. 若期刊要求，可在参考文献中加一条数据引用（Zenodo 会给出建议格式）。

## 四、注意事项

- Zenodo 的 DOI 分“版本 DOI”和“全版本 DOI”，引用时建议用全版本 DOI。
- 仓库里不要放受控数据（KPMP、UK Biobank 个体数据等），本项目的输入均为公开汇总统计。
- 建议在仓库根目录保留 `LICENSE`（代码 MIT）与 `metadata/DATA_DICTIONARY.md`（数据字典），这两项已生成。
