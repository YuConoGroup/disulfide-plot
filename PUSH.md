# 推送到 GitHub

本目录已经是一个独立的 git 仓库，main 分支已有初始 commit。

## 新建空仓库后推送

```bash
cd /data/work/cd/project/disulfide-plot
git remote add origin git@github.com:<你的用户名>/<仓库名>.git
git push -u origin main
```

## 如果仓库已有 README/初始 commit

```bash
git pull --rebase origin main
git push -u origin main
```

## 检查

```bash
git remote -v
git log --oneline
git status
```
