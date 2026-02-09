# 参考基因组构建智能体 - 测试用例集合

本文档包含所有比对工具的测试用例，用于全面测试参考基因组构建智能体的功能。

## 使用说明

在运行测试前，请确保：

1. 已设置 LLM 环境变量：
```bash
export LLM_API_KEY='your-api-key'
export LLM_BASE_URL='https://dashscope.aliyuncs.com/compatible-mode/v1'  # 可选
export LLM_MODEL='qwen-plus'  # 可选
```

2. 已安装所需的比对工具（根据测试用例选择）：
```bash
conda install -c bioconda bwa bowtie2 star hisat2 minimap2
```

3. 确保有足够的磁盘空间（参考基因组和索引文件可能很大）

---

## 1. BWA 索引测试（短读长比对）

### 测试用例 1.1: Human + BWA (WGS场景)
```bash
python agent_main.py "我需要为WGS分析准备人类参考基因组并构建BWA索引" --work-dir test_bwa_human
```

**预期结果**：
- ✓ 下载人类参考基因组（从 NCBI RefSeq 或 Ensembl）
- ✓ 解压缩文件
- ✓ 构建 BWA 索引（单线程，不支持多线程）
- ✓ 索引文件：`test_bwa_human/indices/human_*.bwt`

**检查项**：
- [ ] 文件下载成功
- [ ] 索引文件生成（.bwt, .amb, .ann, .pac, .sa）
- [ ] 日志中显示 "BWA索引构建完成"

---

### 测试用例 1.2: Mouse + BWA
```bash
python agent_main.py "下载小鼠参考基因组并构建BWA索引，用于全基因组测序" --work-dir test_bwa_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 BWA 索引
- ✓ 索引文件：`test_bwa_mouse/indices/mouse_*.bwt`

---

### 测试用例 1.3: E.coli + BWA
```bash
python agent_main.py "我需要为大肠杆菌构建BWA索引用于短读长序列比对" --work-dir test_bwa_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 BWA 索引（测试 URL 下载功能）
- ✓ 索引文件：`test_bwa_ecoli/indices/ecoli_genomic.bwt`

**检查项**：
- [ ] LLM 正确识别 E.coli 并使用 `download_from_url` 工具
- [ ] 下载 URL 正确
- [ ] 路径解析正确（genomes/ 和 indices/ 目录）

---

## 2. Bowtie2 索引测试（短读长比对）

### 测试用例 2.1: Human + Bowtie2
```bash
python agent_main.py "下载人类参考基因组并构建Bowtie2索引" --work-dir test_bowtie2_human
```

**预期结果**：
- ✓ 下载人类参考基因组
- ✓ 解压缩文件
- ✓ 构建 Bowtie2 索引（支持多线程）
- ✓ 索引文件：`test_bowtie2_human/indices/human_*.1.bt2` 等

**检查项**：
- [ ] 索引文件生成（.1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2）
- [ ] 命令中包含 `--threads` 参数（如果指定了线程数）

---

### 测试用例 2.2: Mouse + Bowtie2
```bash
python agent_main.py "为小鼠参考基因组构建Bowtie2索引" --work-dir test_bowtie2_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 Bowtie2 索引
- ✓ 索引文件：`test_bowtie2_mouse/indices/mouse_*.1.bt2`

---

### 测试用例 2.3: E.coli + Bowtie2
```bash
python agent_main.py "我需要为E.coli构建Bowtie2索引" --work-dir test_bowtie2_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 Bowtie2 索引
- ✓ 索引文件：`test_bowtie2_ecoli/indices/ecoli_genomic.1.bt2`

---

## 3. STAR 索引测试（RNA-seq比对）

### 测试用例 3.1: Human + STAR (有注释文件场景)
```bash
python agent_main.py "我需要为RNA-seq分析准备人类参考基因组和STAR索引" --work-dir test_star_human
```

**预期结果**：
- ✓ 下载人类参考基因组
- ✓ 解压缩文件
- ✓ 构建 STAR 索引（可能包含 GTF 文件）
- ✓ 索引目录：`test_star_human/indices/human_*_star/`

**检查项**：
- [ ] 索引目录包含 Genome, SA, SAindex 等文件
- [ ] 如果提供了 GTF 文件，命令中包含 `--sjdbGTFfile` 和 `--sjdbOverhang`
- [ ] 如果没有 GTF 文件，命令中不包含 `--sjdbOverhang`（避免错误）

---

### 测试用例 3.2: Mouse + STAR
```bash
python agent_main.py "为小鼠RNA-seq分析构建STAR索引" --work-dir test_star_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 STAR 索引
- ✓ 索引目录：`test_star_mouse/indices/mouse_*_star/`

---

### 测试用例 3.3: E.coli + STAR (无注释文件场景，测试sjdbOverhang处理)
```bash
python agent_main.py "我需要为RNA-seq分析准备大肠杆菌参考基因组和STAR索引" --work-dir test_star_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 STAR 索引（**不包含 sjdbOverhang 参数**）
- ✓ 索引目录：`test_star_ecoli/indices/ecoli_genomic_star/`

**检查项**：
- [ ] **重要**：命令中**不包含** `--sjdbOverhang` 参数（E.coli 没有注释文件）
- [ ] 日志中显示 "未提供GTF注释文件，构建基础STAR索引"
- [ ] 索引构建成功，无错误

---

## 4. HISAT2 索引测试（RNA-seq比对）

### 测试用例 4.1: Human + HISAT2
```bash
python agent_main.py "下载人类参考基因组并构建HISAT2索引用于转录组分析" --work-dir test_hisat2_human
```

**预期结果**：
- ✓ 下载人类参考基因组
- ✓ 解压缩文件
- ✓ 构建 HISAT2 索引（支持多线程）
- ✓ 索引文件：`test_hisat2_human/indices/human_*.1.ht2` 等

**检查项**：
- [ ] 索引文件生成（.1.ht2 到 .8.ht2）
- [ ] 命令中包含 `-p` 参数（线程数）

---

### 测试用例 4.2: Mouse + HISAT2
```bash
python agent_main.py "为小鼠构建HISAT2索引用于RNA-seq" --work-dir test_hisat2_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 HISAT2 索引
- ✓ 索引文件：`test_hisat2_mouse/indices/mouse_*.1.ht2`

---

### 测试用例 4.3: E.coli + HISAT2
```bash
python agent_main.py "我需要为E.coli构建HISAT2索引" --work-dir test_hisat2_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 HISAT2 索引
- ✓ 索引文件：`test_hisat2_ecoli/indices/ecoli_genomic.1.ht2`

---

## 5. minimap2 索引测试（长读长比对）

### 测试用例 5.1: Human + minimap2
```bash
python agent_main.py "下载人类参考基因组并构建minimap2索引用于长读长序列比对" --work-dir test_minimap2_human
```

**预期结果**：
- ✓ 下载人类参考基因组
- ✓ 解压缩文件
- ✓ 构建 minimap2 索引（单线程，不支持多线程）
- ✓ 索引文件：`test_minimap2_human/indices/human_*.mmi`

**检查项**：
- [ ] 索引文件生成（.mmi）
- [ ] minimap2 索引构建是单线程的

---

### 测试用例 5.2: Mouse + minimap2
```bash
python agent_main.py "为小鼠构建minimap2索引用于PacBio数据比对" --work-dir test_minimap2_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 minimap2 索引
- ✓ 索引文件：`test_minimap2_mouse/indices/mouse_*.mmi`

---

### 测试用例 5.3: E.coli + minimap2
```bash
python agent_main.py "我需要为E.coli构建minimap2索引用于Nanopore数据" --work-dir test_minimap2_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 minimap2 索引
- ✓ 索引文件：`test_minimap2_ecoli/indices/ecoli_genomic.mmi`

---

## 6. 多工具组合测试

### 测试用例 6.1: Human + 多个工具
```bash
python agent_main.py "下载人类参考基因组并构建BWA、Bowtie2和STAR索引" --work-dir test_multi_human
```

**预期结果**：
- ✓ 下载人类参考基因组（只下载一次）
- ✓ 解压缩文件（只解压一次）
- ✓ 构建 BWA 索引
- ✓ 构建 Bowtie2 索引
- ✓ 构建 STAR 索引
- ✓ 所有索引文件都在 `test_multi_human/indices/` 目录

**检查项**：
- [ ] LLM 正确规划任务，避免重复下载和解压缩
- [ ] 使用 `{{task_X.result}}` 引用前一个任务的结果

---

### 测试用例 6.2: Mouse + 多个工具
```bash
python agent_main.py "为小鼠构建BWA和STAR索引" --work-dir test_multi_mouse
```

**预期结果**：
- ✓ 下载小鼠参考基因组
- ✓ 解压缩文件
- ✓ 构建 BWA 索引
- ✓ 构建 STAR 索引

---

### 测试用例 6.3: E.coli + 多个工具
```bash
python agent_main.py "我需要为E.coli构建BWA、Bowtie2和STAR索引" --work-dir test_multi_ecoli
```

**预期结果**：
- ✓ 从 URL 下载 E.coli 参考基因组
- ✓ 解压缩文件
- ✓ 构建 BWA 索引
- ✓ 构建 Bowtie2 索引
- ✓ 构建 STAR 索引（**不包含 sjdbOverhang**）

**检查项**：
- [ ] STAR 索引构建时正确处理 sjdbOverhang 参数

---

## 7. 不同数据源测试

### 测试用例 7.1: NCBI RefSeq
```bash
python agent_main.py "从NCBI RefSeq下载人类参考基因组并构建BWA索引" --work-dir test_ncbi
```

**预期结果**：
- ✓ 从 NCBI RefSeq 下载
- ✓ 构建 BWA 索引

**检查项**：
- [ ] LLM 正确识别数据源要求
- [ ] 下载 URL 来自 NCBI RefSeq

---

### 测试用例 7.2: Ensembl
```bash
python agent_main.py "从Ensembl下载小鼠参考基因组并构建STAR索引" --work-dir test_ensembl
```

**预期结果**：
- ✓ 从 Ensembl 下载
- ✓ 构建 STAR 索引

**检查项**：
- [ ] LLM 正确识别数据源要求
- [ ] 下载 URL 来自 Ensembl

---

### 测试用例 7.3: UCSC
```bash
python agent_main.py "从UCSC下载人类参考基因组并构建Bowtie2索引" --work-dir test_ucsc
```

**预期结果**：
- ✓ 从 UCSC 下载
- ✓ 构建 Bowtie2 索引

**检查项**：
- [ ] LLM 正确识别数据源要求
- [ ] 下载 URL 来自 UCSC

---

## 8. 特殊场景测试

### 测试用例 8.1: 指定线程数
```bash
python agent_main.py "下载人类参考基因组并构建BWA索引，使用8个线程" --work-dir test_threads
```

**预期结果**：
- ✓ 下载人类参考基因组
- ✓ 解压缩文件
- ✓ 构建 BWA 索引（注意：BWA index 不支持多线程，会忽略线程数）

**检查项**：
- [ ] 日志中显示 "BWA index 命令不支持多线程" 的提示
- [ ] 命令中不包含 `-t` 参数

---

### 测试用例 8.2: WES场景
```bash
python agent_main.py "我需要为外显子组测序准备人类参考基因组和BWA索引" --work-dir test_wes
```

**预期结果**：
- ✓ LLM 识别这是 WES 场景
- ✓ 推荐使用 BWA 或 Bowtie2
- ✓ 下载并构建索引

**检查项**：
- [ ] LLM 正确识别使用场景
- [ ] 推荐的工具合适（BWA/Bowtie2）

---

### 测试用例 8.3: ChIP-seq场景
```bash
python agent_main.py "为ChIP-seq分析准备小鼠参考基因组和Bowtie2索引" --work-dir test_chipseq
```

**预期结果**：
- ✓ LLM 识别这是 ChIP-seq 场景
- ✓ 推荐使用 BWA 或 Bowtie2
- ✓ 下载并构建索引

---

## 9. 其他物种测试（测试URL下载功能）

### 测试用例 9.1: 酵母 (Saccharomyces cerevisiae)
```bash
python agent_main.py "下载酵母参考基因组并构建BWA索引" --work-dir test_yeast
```

**预期结果**：
- ✓ LLM 识别酵母物种
- ✓ 使用 `download_from_url` 工具查找并下载
- ✓ 构建 BWA 索引

**检查项**：
- [ ] LLM 正确查找酵母的下载 URL
- [ ] URL 来自 NCBI RefSeq 或 Ensembl
- [ ] 下载成功

---

### 测试用例 9.2: 果蝇 (Drosophila melanogaster)
```bash
python agent_main.py "我需要为果蝇构建STAR索引用于RNA-seq分析" --work-dir test_fly
```

**预期结果**：
- ✓ LLM 识别果蝇物种
- ✓ 使用 `download_from_url` 工具查找并下载
- ✓ 构建 STAR 索引（可能包含 GTF 文件）

**检查项**：
- [ ] LLM 正确查找果蝇的下载 URL
- [ ] STAR 索引构建正确处理注释文件

---

### 测试用例 9.3: 拟南芥 (Arabidopsis thaliana)
```bash
python agent_main.py "下载拟南芥参考基因组并构建Bowtie2索引" --work-dir test_arabidopsis
```

**预期结果**：
- ✓ LLM 识别拟南芥物种
- ✓ 使用 `download_from_url` 工具查找并下载
- ✓ 构建 Bowtie2 索引

---

## 10. 推荐功能测试

### 测试用例 10.1: WGS推荐
```bash
python agent_main.py --recommend "WGS全基因组测序分析"
```

**预期结果**：
- ✓ 推荐数据源：ncbi_refseq（最常用）
- ✓ 推荐工具：bwa 或 bowtie2
- ✓ 提供推荐理由

**检查项**：
- [ ] 推荐方案合理
- [ ] 理由充分

---

### 测试用例 10.2: RNA-seq推荐
```bash
python agent_main.py --recommend "RNA-seq转录组分析"
```

**预期结果**：
- ✓ 推荐数据源：ensembl（适合RNA-seq）或 ncbi_refseq
- ✓ 推荐工具：star 或 hisat2
- ✓ 提供推荐理由和参数建议

**检查项**：
- [ ] 推荐 STAR 或 HISAT2
- [ ] 如果推荐 STAR，包含 sjdb_overhang 参数建议

---

### 测试用例 10.3: 长读长推荐
```bash
python agent_main.py --recommend "PacBio长读长序列分析"
```

**预期结果**：
- ✓ 推荐数据源：ncbi_refseq
- ✓ 推荐工具：minimap2
- ✓ 提供推荐理由

**检查项**：
- [ ] 推荐 minimap2
- [ ] 理由说明适合长读长数据

---

## 11. 错误处理测试

### 测试用例 11.1: 不存在的物种
```bash
python agent_main.py "下载不存在的物种参考基因组并构建BWA索引" --work-dir test_error_species
```

**预期结果**：
- ✓ LLM 尝试查找 URL
- ✓ 如果找不到，给出错误提示或使用本地文件假设

---

### 测试用例 11.2: 无效的工具名称
```bash
python agent_main.py "下载人类参考基因组并构建不存在的工具索引" --work-dir test_error_tool
```

**预期结果**：
- ✓ 给出错误提示：不支持的索引工具

---

## 批量测试脚本

创建一个测试脚本来批量运行测试：

```bash
#!/bin/bash
# test_all_tools.sh

# 设置环境变量
export LLM_API_KEY='your-api-key'
export LLM_BASE_URL='https://dashscope.aliyuncs.com/compatible-mode/v1'
export LLM_MODEL='qwen-plus'

# 测试结果目录
TEST_DIR="test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p $TEST_DIR

# 定义测试用例（格式：请求|工作目录|描述）
declare -a tests=(
    "我需要为WGS分析准备人类参考基因组并构建BWA索引|test_bwa_human|BWA Human"
    "下载小鼠参考基因组并构建Bowtie2索引|test_bowtie2_mouse|Bowtie2 Mouse"
    "我需要为RNA-seq分析准备大肠杆菌参考基因组和STAR索引|test_star_ecoli|STAR E.coli"
    "下载人类参考基因组并构建HISAT2索引用于转录组分析|test_hisat2_human|HISAT2 Human"
    "为小鼠构建minimap2索引用于PacBio数据比对|test_minimap2_mouse|minimap2 Mouse"
)

# 运行测试
for test in "${tests[@]}"; do
    IFS='|' read -r request work_dir description <<< "$test"
    echo "=========================================="
    echo "测试: $description"
    echo "请求: $request"
    echo "工作目录: $work_dir"
    echo "=========================================="
    
    python agent_main.py "$request" --work-dir "$TEST_DIR/$work_dir" 2>&1 | tee "$TEST_DIR/${work_dir}.log"
    
    # 检查退出状态
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "✓ 测试通过: $description"
    else
        echo "✗ 测试失败: $description"
    fi
    
    echo ""
    sleep 2
done

echo "所有测试完成，结果保存在: $TEST_DIR"
```

## 测试检查清单

每个测试用例完成后，检查以下项目：

### 通用检查项
- [ ] 任务规划正确（LLM 正确理解需求）
- [ ] 文件下载成功（如果适用）
- [ ] 文件解压缩成功（如果适用）
- [ ] 索引构建成功
- [ ] 日志中没有错误信息
- [ ] 路径解析正确（genomes/ 和 indices/ 目录）

### 下载检查项
- [ ] 文件存在于 `genomes/` 目录
- [ ] 文件大小合理（不是0字节）
- [ ] 对于非内置物种，LLM 正确使用 `download_from_url`

### 解压缩检查项
- [ ] 解压后的文件存在
- [ ] 文件格式正确（.fa 或 .fasta）

### 索引构建检查项
- [ ] 索引文件存在于 `indices/` 目录
- [ ] 索引文件后缀正确：
  - BWA: `.bwt`, `.amb`, `.ann`, `.pac`, `.sa`
  - Bowtie2: `.1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`
  - STAR: 目录包含 `Genome`, `SA`, `SAindex`
  - HISAT2: `.1.ht2` 到 `.8.ht2`
  - minimap2: `.mmi`

### 特殊检查项
- [ ] **BWA**: 命令中不包含 `-t` 参数（不支持多线程）
- [ ] **STAR (无注释文件)**: 命令中不包含 `--sjdbOverhang` 参数
- [ ] **STAR (有注释文件)**: 命令中包含 `--sjdbGTFfile` 和 `--sjdbOverhang`
- [ ] **Bowtie2/HISAT2**: 命令中包含线程参数（如果指定了线程数）

## 注意事项

1. **磁盘空间**：确保有足够的磁盘空间（人类参考基因组约 3GB，索引可能更大）
2. **内存**：STAR 索引构建需要大量内存（建议 32GB+）
3. **网络**：下载需要稳定的网络连接
4. **时间**：索引构建可能需要较长时间，请耐心等待
5. **LLM API 费用**：使用 LLM 会产生 API 调用费用

## 故障排除

### 问题：LLM 规划错误
- 检查 LLM_API_KEY 是否正确设置
- 检查网络连接
- 查看日志中的 LLM 响应

### 问题：下载失败
- 检查网络连接
- 检查 URL 是否正确
- 查看下载日志

### 问题：索引构建失败
- 检查工具是否已安装：`which bwa`, `which STAR` 等
- 检查内存是否足够
- 检查参考基因组文件是否完整
- 查看错误日志获取详细信息

### 问题：路径错误
- 检查工作目录权限
- 检查文件是否存在
- 查看日志中的路径信息
