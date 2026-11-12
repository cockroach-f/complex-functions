# 复变函数从零到精通

从零开始的复变函数讲义（华中科技大学《复变函数与积分变换》第 1–6 章：复数与复变函数、解析函数、复积分、级数、留数、共形映射），零基础起步，配 manim 动画可视化、知识点归纳总结与练习题。

**在线阅读：** <https://cockroach-f.github.io/complex-functions/>

姊妹站点：[ISAR 目标识别入门与文献带读](https://cockroach-f.github.io/isar-tutorial/)

## 进度

- [x] 第 1 章 复数与复变函数（正文 + 5 段动画 + 练习）
- [x] 第 2 章 解析函数（正文 + 4 段动画 + 练习）
- [x] 第 3 章 复变函数的积分（正文 + 3 段动画 + 练习）
- [x] 第 4 章 解析函数的级数表示（正文 + 3 段动画 + 练习）
- [x] 第 5 章 留数及其应用（正文 + 4 段动画 + 练习）
- [x] 第 6 章 共形映射（正文 + 4 段动画 + 练习）
- [x] 全课程总结（主线、题型速查、后续课程联系、自测清单）

## 本地预览

```bash
pip install mkdocs-material
mkdocs serve
```

浏览器打开提示的本地地址（默认 `http://127.0.0.1:8000`）。

## 目录结构

```
docs/
  index.md              # 前言与课程知识地图
  ch1.md … ch6.md       # 正文各章
  summary.md            # 全课程总结
  formulas.md           # 公式速查
  references.md         # 参考资料与致谢
  assets/videos/        # manim 渲染的 MP4（720p）
  assets/posters/       # 视频封面
  assets/figs/          # matplotlib 示意图
animations/             # manim 场景源码（common.py + chN_scenes.py）
tools/
  render_videos.sh      # 渲染全部动画并压缩到 docs/assets/videos
  make_figures.py       # 生成静态示意图
```

## 重新渲染动画

依赖 [Manim Community](https://www.manim.community/)（≥0.18）、ffmpeg 和 TeX 发行版：公式用 `MathTex` 排版，中文用 `Tex` + **ctex**（xelatex 编译）。最小安装可用 [TinyTeX](https://yihui.org/tinytex/)，再补装宏包：

```bash
tlmgr install ctex xecjk fandol standalone preview dvisvgm mathtools jknapltx fundus-calligra \
  wasy wasy-type1 tipa microtype psnfss physics doublestroke relsize setspace ragged2e babel-english cm-super
```

`tools/render_videos.sh` 会自动把 `TINYTEX_BIN`（默认 `/workspace/tools/.TinyTeX/bin/x86_64-linux`）加入 PATH。

```bash
pip install manim            # 或 conda install -c conda-forge manim
tools/render_videos.sh                       # 全部场景
tools/render_videos.sh ch1_scenes.py C1_MapSquare   # 单个场景
python tools/make_figures.py                 # 静态图
```

## 部署

推送到 `main` 后，GitHub Actions（`.github/workflows/pages.yml`）自动执行 `mkdocs build --strict` 并发布到 GitHub Pages。

## 资料来源与版权

内容组织参考 [N-E-E/hust-complex_functions](https://github.com/N-E-E/hust-complex_functions) 中整理的华中科技大学课程 PPT、教材与习题全解。**本仓库未包含原 PPT / PDF**；讲解文字、例题解答、示意图与动画均为重新编写制作，仅供学习交流。
