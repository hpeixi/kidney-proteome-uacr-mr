# 上传发布包到 GitHub（无需安装 git）

目标仓库：https://github.com/hpeixi/kidney-proteome-uacr-mr
发布包：`kidney_bioinfo_topics/发布包_kidney-proteome-uacr-mr.zip`（35 个文件，约 1 MB）

---

## 方式 A：网页拖拽（推荐，约 3 分钟）

1. 把压缩包解压到任意文件夹，得到：`README.md`、`LICENSE`、`.zenodo.json`、`DOI_获取步骤.md`、`analysis/`、`results/`、`figures/`、`metadata/`。
2. 浏览器打开仓库首页，点右侧的 **uploading an existing file** 链接（或 Add file → Upload files）。
3. 把解压后的**全部内容**拖进上传区——拖文件夹和文件，**不要拖压缩包本身**。
4. 检查文件列表里是否有 `.zenodo.json`。Windows 默认不显示点开头的文件，若没看到，先在同目录地址栏输入 `%cd%` 打开资源管理器并勾选“显示隐藏文件”，再单独把它拖进去一次。
5. 下方 Commit changes，提交说明填 `v1.0.0: analysis code and derived tables`，点 Commit。

上传完成后，仓库首页应能看到 `analysis`、`results`、`figures`、`metadata` 四个目录和 4 个根文件。

---

## 方式 B：GitHub Desktop（以后还要改代码时用）

1. 安装 GitHub Desktop 并登录你的账号。
2. File → Add local repository → 选择解压后的文件夹 → 提示“not a git repository”时选 **create a repository**。
3. 右上角 **Publish repository**，取消勾选 “Keep this code private”。
4. 以后修改后 Commit → Push 即可。

---

## 上传后拿 DOI（Zenodo）

1. 打开 https://zenodo.org ，用 GitHub 账号登录。
2. 右上头像 → **Settings** → **GitHub** → 找到 `kidney-proteome-uacr-mr`，把开关打开。
3. 回 GitHub 仓库 → **Releases** → **Create a new release** → Tag 填 `v1.0.0` → **Publish release**。
4. 回到 Zenodo 该仓库页面，会出现 DOI（形如 `10.5281/zenodo.1234567`）。
5. 把 **仓库地址 + DOI** 发我，我填进 `metadata/CITATION.cff`、`metadata/DOI.txt`，并改稿件 Declarations 的数据可用性声明。

---

## 注意事项

- 不要上传压缩包本身，GitHub 不会自动解压。
- 仓库里不要放 KPMP、UK Biobank 等受控数据；本项目输入全是公开汇总统计。
- 若网页上传中断，重试即可；已上传的文件会保留。
- Release 的 tag 建议固定为 `v1.0.0`，与 `CITATION.cff` 里的 `version` 一致。
