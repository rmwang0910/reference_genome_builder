# 快速开始指南

本指南将帮助您快速上手使用参考基因组构建智能体。

## 安装

### 1. 安装 Python 依赖

```bash
pip install -r requirements.txt
```

### 2. 安装比对工具（根据需要选择）

**使用 conda（推荐）**：
```bash
conda install -c bioconda bwa bowtie2 star hisat2 minimap2
```

**使用 Homebrew（macOS）**：
```bash
brew install bwa bowtie2 star hisat2 minimap2
```

**验证安装**：
```bash
python main.py list-tools
```

## 快速开始

### 场景 1: 下载人类参考基因组

```bash
python main.py download --source ncbi_refseq --species human
```

下载完成后，文件会保存在 `genomes/` 目录下。

### 场景 2: 构建 BWA 索引

假设您已经有了参考基因组文件 `genomes/human_GRCh38.p14.fa`：

```bash
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 8
```

索引文件会保存在 `indices/` 目录下。

### 场景 3: 一键完成下载和索引构建

```bash
python main.py pipeline \
  --source ensembl \
  --species human \
  --tools bwa bowtie2 \
  --threads 8
```

这会自动完成：
1. 下载参考基因组
2. 解压缩
3. 构建指定的索引

## 常用命令

### 查看可用的参考基因组

```bash
# 查看所有
python main.py list-genomes

# 只查看人类
python main.py list-genomes --species human
```

### 查看可用的索引工具

```bash
python main.py list-tools
```

### 下载参考基因组

```bash
# 从 NCBI RefSeq 下载
python main.py download --source ncbi_refseq --species human

# 从 Ensembl 下载
python main.py download --source ensembl --species mouse

# 从 UCSC 下载
python main.py download --source ucsc --species human
```

### 构建索引

```bash
# 构建单个索引
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 8

# 构建多个索引
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa bowtie2 star --threads 8

# 构建 STAR 索引（自定义参数）
python main.py index \
  --genome genomes/human_GRCh38.p14.fa \
  --tools star \
  --threads 8 \
  --sjdb-overhang 99 \
  --genome-sa-index-n-bases 14
```

## 典型工作流

### WGS（全基因组测序）分析准备

```bash
# 1. 下载参考基因组
python main.py download --source ncbi_refseq --species human

# 2. 构建 BWA 索引
python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 16
```

### RNA-seq 分析准备

```bash
# 1. 下载参考基因组
python main.py download --source ensembl --species human

# 2. 构建 STAR 索引（假设 reads 长度为 100bp）
python main.py index \
  --genome genomes/human_GRCh38.fa \
  --tools star \
  --threads 16 \
  --sjdb-overhang 99
```

### 多工具索引构建

```bash
# 一次性构建所有常用索引
python main.py pipeline \
  --source ensembl \
  --species human \
  --tools bwa bowtie2 star hisat2 \
  --threads 16 \
  --sjdb-overhang 99
```

## 输出目录结构

运行后，您的工作目录结构如下：

```
work_dir/
├── genomes/              # 参考基因组文件
│   ├── human_GRCh38.p14.fa.gz    # 压缩文件
│   └── human_GRCh38.p14.fa       # 解压后的FASTA文件
├── indices/              # 索引文件
│   ├── human_GRCh38.p14.bwt      # BWA索引
│   ├── human_GRCh38.p14.1.bt2    # Bowtie2索引
│   └── human_GRCh38.p14_star/    # STAR索引目录
└── build_results.json    # 构建结果记录
```

## 注意事项

1. **磁盘空间**：确保有足够的磁盘空间（参考基因组约 3GB，索引可能更大）
2. **内存**：构建索引，特别是 STAR，需要大量内存（建议 32GB+）
3. **网络**：下载需要稳定的网络连接
4. **时间**：索引构建可能需要较长时间，请耐心等待

## 获取帮助

查看完整帮助信息：

```bash
python main.py --help
python main.py download --help
python main.py index --help
python main.py pipeline --help
```

## 故障排除

### 工具未找到

确保工具已安装并在 PATH 中：
```bash
which bwa
which bowtie2
```

### 下载失败

检查网络连接，或手动下载后放入 `genomes/` 目录。

### 索引构建失败

- 检查内存是否足够
- 检查参考基因组文件是否完整
- 查看错误日志获取详细信息

## 下一步

- 阅读完整的 [README.md](README.md) 了解详细功能
- 查看 [example.py](example.py) 了解 Python API 使用方法
- 根据需要修改 [config.yaml.example](config.yaml.example) 配置文件
