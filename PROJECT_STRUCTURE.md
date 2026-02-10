# 项目结构说明

## 核心文件

```
reference_genome_builder/
├── agent_main.py          # 智能体模式入口（自然语言交互）
├── agent_core.py          # 智能体核心逻辑（LLM集成、任务规划）
├── llm_client.py          # LLM客户端封装（支持OpenAI兼容API）
├── main.py                # 传统模式入口（命令行参数）
├── download.py            # 参考基因组下载模块
├── index_builder.py       # 索引构建模块
├── requirements.txt       # Python依赖
├── LICENSE                # MIT许可证
└── README.md              # 项目说明文档
```

## 文档文件

```
├── ARCHITECTURE.md        # 架构设计文档
├── AGENT_FEATURES.md      # 智能体特性说明
├── LLM_SETUP.md          # LLM配置指南
├── QUICKSTART.md         # 快速开始指南
├── TEST_CASES.md         # 测试用例集合
└── PROJECT_STRUCTURE.md  # 本文件
```

## 配置和示例

```
├── config.yaml.example   # 配置文件示例
├── example.py            # Python API使用示例
└── scripts/
    └── test_all_tools.sh # 批量测试脚本
```

## 输出目录（运行时生成，已加入.gitignore）

```
work_dir/
└── outputs/              # 统一输出目录
    ├── human_GRCh38.p14/
    │   ├── genome/      # 参考基因组文件
    │   ├── indices/     # 索引文件
    │   └── annotation/  # 注释文件（可选）
    └── ...
```

## 忽略的文件（.gitignore）

- `__pycache__/` - Python缓存
- `outputs/` - 输出目录（包含下载的基因组和索引）
- `test_results_*/` - 测试结果
- `*.log`, `*.out` - 日志文件
- `genomes/`, `indices/` - 旧版输出目录
- `.env`, `config.yaml` - 配置文件（可能包含敏感信息）
