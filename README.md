# 参考基因组构建智能体

一个智能化的参考基因组下载和索引构建工具，支持从多个数据源下载参考基因组，并使用多种主流比对工具构建索引。

## 🎯 智能体架构特性

本工具采用**真正的智能体架构**，具备以下核心能力：

- 🧠 **LLM集成**：使用大语言模型理解自然语言请求，自动规划任务
- 🎯 **智能决策**：根据使用场景自动推荐最佳数据源和索引工具
- 📋 **任务规划**：将用户需求自动分解为可执行的任务序列
- 🔄 **错误恢复**：自动重试失败任务，支持降级策略
- 📊 **状态管理**：完整的任务状态跟踪和历史记录
- 🛠️ **工具调用**：动态选择和调用工具，支持任务间结果传递

**两种使用模式**：
1. **传统模式** (`main.py`)：命令行参数，精确控制
2. **智能体模式** (`agent_main.py`)：自然语言交互，智能决策

详见 [ARCHITECTURE.md](ARCHITECTURE.md) 了解完整的架构设计。

## 功能特性

- 🔽 **多数据源支持**：支持从 NCBI RefSeq、Ensembl、UCSC 下载参考基因组
- 🧬 **多物种支持**：支持人类、小鼠等常见模式生物的参考基因组
- 🔧 **多工具支持**：支持 BWA、Bowtie2、STAR、HISAT2、minimap2 等主流比对工具的索引构建
- ⚡ **自动化流程**：一键完成下载、解压缩、索引构建的完整流程
- 📊 **进度显示**：下载和构建过程实时显示进度
- ✅ **智能检测**：自动检测已存在的文件和索引，避免重复下载和构建

## 安装

### 前置要求

1. Python >= 3.7
2. 需要安装相应的比对工具（根据使用需求选择）：
   - **BWA**: `conda install -c bioconda bwa` 或 `brew install bwa`
   - **Bowtie2**: `conda install -c bioconda bowtie2` 或 `brew install bowtie2`
   - **STAR**: `conda install -c bioconda star` 或 `brew install star`
   - **HISAT2**: `conda install -c bioconda hisat2` 或 `brew install hisat2`
   - **minimap2**: `conda install -c bioconda minimap2` 或 `brew install minimap2`

### 安装依赖

```bash
pip install -r requirements.txt
```

## 使用方法

### 🚀 快速开始：智能体模式（推荐）

使用自然语言描述您的需求，智能体会自动规划并执行：

```bash
# RNA-seq分析准备
python agent_main.py "我需要为RNA-seq分析准备人类参考基因组和STAR索引"

# WGS分析准备
python agent_main.py "下载人类参考基因组并构建BWA索引"

# 获取推荐方案
python agent_main.py --recommend "WGS全基因组测序分析"
```

智能体会自动：
1. 理解您的需求和使用场景
2. 推荐最佳的数据源和工具
3. 规划任务步骤
4. 自动执行所有任务

### 📋 传统模式：精确控制

如果您需要精确控制每个步骤，可以使用传统模式：

### 1. 列出可用的参考基因组

```bash
python main.py list-genomes
```

按物种过滤：
```bash
python main.py list-genomes --species human
```

### 2. 列出可用的索引工具

```bash
python main.py list-tools
```

检查特定参考基因组的工具可用性：
```bash
python main.py list-tools --genome genomes/human_GRCh38.p14.fa
```

### 3. 下载参考基因组

从 NCBI RefSeq 下载人类参考基因组：
```bash
python main.py download --source ncbi_refseq --species human
```

从 Ensembl 下载小鼠参考基因组：
```bash
python main.py download --source ensembl --species mouse
```

指定工作目录：
```bash
python main.py download --source ucsc --species human --work-dir /path/to/workdir
```

保留压缩文件（不解压缩）：
```bash
python main.py download --source ncbi_refseq --species human --no-decompress
```

### 4. 构建参考基因组索引

构建 BWA 索引：
```bash
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 8
```

构建多个索引：
```bash
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa bowtie2 star --threads 8
```

构建 STAR 索引（自定义参数）：
```bash
python main.py index --genome genomes/human_GRCh38.p14.fa --tools star --threads 8 --sjdb-overhang 99 --genome-sa-index-n-bases 14
```

### 5. 完整流程（下载 + 构建索引）

一键完成下载和索引构建：
```bash
python main.py pipeline --source ensembl --species human --tools bwa bowtie2 star --threads 8
```

## 目录结构

运行后会在工作目录下创建以下目录结构：

```
work_dir/
├── genomes/              # 下载的参考基因组文件
│   ├── human_GRCh38.p14.fa.gz
│   └── human_GRCh38.p14.fa
├── indices/              # 构建的索引文件
│   ├── human_GRCh38.p14.bwt      # BWA索引
│   ├── human_GRCh38.p14.1.bt2    # Bowtie2索引
│   └── human_GRCh38.p14_star/   # STAR索引目录
└── build_results.json     # 构建结果记录
```

## 支持的数据源和物种

### 数据源

- **NCBI RefSeq**: 官方参考序列数据库
- **Ensembl**: 欧洲生物信息学研究所维护的基因组数据库
- **UCSC**: 加州大学圣克鲁兹分校的基因组浏览器数据库

### 支持的物种

**内置支持（使用 download 工具）**：
- **human**: 人类（Homo sapiens）
- **mouse**: 小鼠（Mus musculus）

**任意物种（使用 download_from_url 工具）**：
- 支持从 NCBI RefSeq、Ensembl 等数据源下载任意物种的参考基因组
- LLM 会自动查找并下载指定物种的参考基因组
- 例如：E.coli、酵母、果蝇、拟南芥等

**示例**：
```bash
# LLM 会自动查找 E.coli 的下载URL并下载
python agent_main.py "我需要为RNA-seq分析准备大肠杆菌参考基因组和STAR索引"

# LLM 会自动查找酵母的下载URL并下载
python agent_main.py "下载酵母参考基因组并构建BWA索引"
```

## 支持的索引工具

| 工具 | 用途 | 适用场景 |
|------|------|----------|
| BWA | 短读长序列比对 | WGS、WES、ChIP-seq |
| Bowtie2 | 超快速短读长序列比对 | WGS、WES、ChIP-seq |
| STAR | RNA-seq 比对 | 转录组分析 |
| HISAT2 | RNA-seq 比对 | 转录组分析（内存占用更小） |
| minimap2 | 长读长序列比对 | PacBio、Nanopore 数据 |

## 配置说明

### STAR 索引参数

- `--sjdb-overhang`: Splice junction database overhang，通常设置为 read 长度 - 1
  - 对于 100bp 的 reads，设置为 99
  - 对于 150bp 的 reads，设置为 149

- `--genome-sa-index-n-bases`: 基因组 SA 索引的碱基数
  - 对于小基因组（< 4.5 Gb），使用默认值 14
  - 对于大基因组（> 4.5 Gb），需要设置为 `log2(基因组大小)/2 - 1`

## 示例工作流

### 示例 1: WGS 数据分析准备

```bash
# 1. 下载人类参考基因组
python main.py download --source ncbi_refseq --species human

# 2. 构建 BWA 索引（用于比对）
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 16
```

### 示例 2: RNA-seq 数据分析准备

```bash
# 1. 下载人类参考基因组
python main.py download --source ensembl --species human

# 2. 构建 STAR 索引（用于 RNA-seq 比对）
python main.py index --genome genomes/human_GRCh38.fa --tools star --threads 16 --sjdb-overhang 99
```

### 示例 3: 多工具索引构建

```bash
# 一次性构建所有常用索引
python main.py pipeline \
  --source ensembl \
  --species human \
  --tools bwa bowtie2 star hisat2 \
  --threads 16 \
  --sjdb-overhang 99
```

## 注意事项

1. **磁盘空间**：参考基因组文件通常很大（人类参考基因组约 3GB），索引文件可能更大（STAR 索引可能超过 30GB），请确保有足够的磁盘空间。

2. **内存需求**：构建索引时，特别是 STAR 索引，需要大量内存。建议至少 32GB RAM。

3. **网络连接**：下载参考基因组需要稳定的网络连接。如果下载中断，程序会重新下载。

4. **工具安装**：确保所需的比对工具已正确安装并在 PATH 中可用。

5. **文件完整性**：程序会自动检查已存在的文件，但不会验证文件完整性。如果怀疑文件损坏，请删除后重新下载。

## 故障排除

### 问题：工具未找到

**解决方案**：
- 确保工具已正确安装
- 检查工具是否在 PATH 中：`which bwa`
- 使用 conda 安装：`conda install -c bioconda <tool_name>`

### 问题：STAR 索引构建失败

**可能原因**：
- 内存不足
- `genomeSAindexNbases` 参数设置不当（对于大基因组）

**解决方案**：
- 增加可用内存或使用更少的线程
- 对于大基因组，调整 `--genome-sa-index-n-bases` 参数

### 问题：下载速度慢

**解决方案**：
- 检查网络连接
- 考虑使用镜像站点
- 使用 `wget` 或 `curl` 手动下载后放入 `genomes/` 目录

## 许可证

本项目遵循 MIT 许可证。

## 贡献

欢迎提交 Issue 和 Pull Request！

## 🤖 智能体架构说明

### 什么是智能体架构？

智能体架构使工具具备**自主决策**和**智能规划**能力，而不仅仅是执行预设命令。

#### 核心能力对比

| 特性 | 传统工具 | 智能体 |
|------|---------|--------|
| **交互方式** | 命令行参数 | 自然语言 |
| **决策能力** | 用户指定 | 自动决策 |
| **任务规划** | 手动步骤 | 自动规划 |
| **错误处理** | 手动重试 | 自动恢复 |
| **推荐功能** | 无 | 智能推荐 |
| **状态管理** | 无 | 完整状态跟踪 |

#### 智能体工作流程

```
用户自然语言请求
    ↓
[LLM理解意图] ← 智能决策
    ↓
[任务规划] ← 自动分解为步骤
    ↓
[工具调用] ← 动态选择工具
    ↓
[状态管理] ← 跟踪执行状态
    ↓
[错误恢复] ← 自动重试/降级
    ↓
结果返回
```

#### 使用示例

**传统模式**（需要明确指定所有参数）：
```bash
python main.py pipeline --source ensembl --species human --tools star --threads 8 --sjdb-overhang 99
```

**智能体模式**（自然语言，自动决策）：
```bash
python agent_main.py "我需要为RNA-seq分析准备人类参考基因组和STAR索引"
```

智能体会自动：
- ✅ 识别这是RNA-seq场景
- ✅ 推荐Ensembl数据源（更适合RNA-seq）
- ✅ 选择STAR工具（RNA-seq最佳选择）
- ✅ 自动设置合适的参数（sjdb_overhang=99）
- ✅ 规划并执行所有步骤

### 架构文档

详细架构说明请参考：[ARCHITECTURE.md](ARCHITECTURE.md)

### LLM配置（必需）

智能体模式**统一使用LLM**进行智能决策，不再支持规则模式。

**配置LLM（必需）**：

1. 安装依赖：
```bash
pip install openai
```

2. 设置环境变量：
```bash
# 必需：API 密钥
export LLM_API_KEY='your-api-key'

# 可选：API 基础 URL（默认使用 OpenAI）
export LLM_BASE_URL='https://dashscope.aliyuncs.com/compatible-mode/v1'  # 阿里云 DashScope
# 或
export LLM_BASE_URL='https://api.openai.com/v1'  # OpenAI（默认）

# 可选：模型名称（默认 gpt-3.5-turbo）
export LLM_MODEL='qwen-plus'  # 阿里云 DashScope
# 或
export LLM_MODEL='gpt-4'  # OpenAI
```

**支持的 LLM 服务**：
- OpenAI API（默认）
- 阿里云 DashScope（兼容模式）
- 其他 OpenAI 兼容的 API 服务

**注意**：如果未配置 LLM，程序将无法运行并提示错误。

## 更新日志

### v1.1.0
- ✨ 新增智能体架构支持
- 🧠 LLM集成，支持自然语言交互
- 🎯 智能推荐功能
- 📋 自动任务规划
- 🔄 错误自动恢复机制
- 📊 完整状态管理

### v1.0.0
- 初始版本
- 支持从 NCBI RefSeq、Ensembl、UCSC 下载参考基因组
- 支持 BWA、Bowtie2、STAR、HISAT2、minimap2 索引构建
- 完整的命令行接口
- 自动化流程支持
