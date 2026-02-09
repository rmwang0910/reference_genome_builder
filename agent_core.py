#!/usr/bin/env python
"""
参考基因组构建智能体核心模块

实现真正的智能体架构，包括：
1. LLM集成 - 智能决策和任务规划
2. 工具调用能力 - 动态选择和调用工具
3. 状态管理 - 跟踪任务状态
4. 错误恢复 - 自动重试和替代方案
5. 智能推荐 - 根据用户需求推荐最佳方案
"""

import os
import sys
import logging
import json
import re
from pathlib import Path
from typing import Optional, Dict, List, Any, Callable
from enum import Enum
from dataclasses import dataclass, asdict
from datetime import datetime

logger = logging.getLogger(__name__)

# 导入本地 LLM 客户端
try:
    from llm_client import LLMClient, create_llm_client
    HAS_LLM = True
except ImportError as e:
    HAS_LLM = False
    logger.error(f"无法导入 LLM 客户端: {e}")


class TaskStatus(Enum):
    """任务状态"""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class Task:
    """任务数据类"""
    id: str
    name: str
    description: str
    status: TaskStatus
    result: Optional[Any] = None
    error: Optional[str] = None
    retry_count: int = 0
    max_retries: int = 3
    
    def to_dict(self) -> dict:
        """转换为字典"""
        data = asdict(self)
        data['status'] = self.status.value
        return data


@dataclass
class AgentState:
    """智能体状态"""
    current_task: Optional[str] = None
    tasks: List[Task] = None
    context: Dict[str, Any] = None
    history: List[Dict[str, Any]] = None
    
    def __post_init__(self):
        if self.tasks is None:
            self.tasks = []
        if self.context is None:
            self.context = {}
        if self.history is None:
            self.history = []
    
    def add_task(self, task: Task):
        """添加任务"""
        self.tasks.append(task)
    
    def get_task(self, task_id: str) -> Optional[Task]:
        """获取任务"""
        for task in self.tasks:
            if task.id == task_id:
                return task
        return None
    
    def update_task_status(self, task_id: str, status: TaskStatus, result: Any = None, error: str = None):
        """更新任务状态"""
        task = self.get_task(task_id)
        if task:
            task.status = status
            if result is not None:
                task.result = result
            if error is not None:
                task.error = error
    
    def to_dict(self) -> dict:
        """转换为字典"""
        return {
            'current_task': self.current_task,
            'tasks': [task.to_dict() for task in self.tasks],
            'context': self.context,
            'history': self.history
        }


class ReferenceGenomeAgent:
    """参考基因组构建智能体（增强版）"""
    
    def __init__(self, work_dir: Optional[str] = None, use_llm: bool = True):
        """
        初始化智能体
        
        Args:
            work_dir: 工作目录
            use_llm: 是否使用LLM进行智能决策（必须为True，否则会抛出异常）
        """
        if not use_llm:
            raise ValueError("use_llm 必须为 True，规则模式已被移除")
        
        if not HAS_LLM:
            raise ImportError(
                "LLM 模块不可用。请确保：\n"
                "1. 已安装 openai 库: pip install openai\n"
                "2. 已设置 LLM_API_KEY 环境变量"
            )
        
        if work_dir is None:
            self.work_dir = Path.cwd()
        else:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)
        
        self.genomes_dir = self.work_dir / "genomes"
        self.indices_dir = self.work_dir / "indices"
        
        # 初始化状态
        self.state = AgentState()
        
        # 初始化LLM（必须成功）
        self.use_llm = True
        try:
            self.llm_client = LLMClient()
            logger.info("LLM已初始化，启用智能决策")
        except Exception as e:
            logger.error(f"无法初始化LLM: {e}")
            raise ValueError(
                f"LLM 初始化失败: {e}\n"
                "请确保已设置以下环境变量：\n"
                "  - LLM_API_KEY: API 密钥\n"
                "  - LLM_BASE_URL: API 基础 URL（可选，默认 OpenAI）\n"
                "  - LLM_MODEL: 模型名称（可选，默认 gpt-3.5-turbo）"
            )
        
        # 注册工具
        self.tools = self._register_tools()
        
        logger.info(f"智能体初始化完成，工作目录: {self.work_dir}")
        logger.info(f"LLM模式: 启用（统一使用LLM，规则模式已移除）")
    
    def _register_tools(self) -> Dict[str, Callable]:
        """注册可用工具"""
        from download import ReferenceGenomeDownloader
        from index_builder import IndexBuilder
        
        tools = {}
        
        # 下载工具（从预定义数据源下载，仅支持human和mouse）
        def download_tool(source: str, species: str, **kwargs):
            downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
            result_path = downloader.download_from_source(source, species)
            # 返回绝对路径，确保后续任务能正确找到文件
            return result_path.resolve()
        
        # 从URL下载工具（支持任意物种）
        def download_from_url_tool(url: str, output_filename: Optional[str] = None, **kwargs):
            """从URL下载参考基因组，支持任意物种"""
            downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
            result_path = downloader.download_from_url(url, output_filename)
            # 返回绝对路径，确保后续任务能正确找到文件
            return result_path.resolve()
        
        # 解压缩工具
        def decompress_tool(file_path: str, **kwargs):
            downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
            # 处理路径：如果是相对路径，尝试在 genomes_dir 中查找
            file_path_obj = Path(file_path)
            if not file_path_obj.is_absolute():
                # 相对路径，尝试在 genomes_dir 中查找
                potential_path = self.genomes_dir / file_path_obj
                if potential_path.exists():
                    file_path_obj = potential_path
                else:
                    # 如果不在 genomes_dir，尝试相对于工作目录
                    file_path_obj = (self.work_dir / file_path_obj).resolve()
            else:
                file_path_obj = file_path_obj.resolve()
            
            if not file_path_obj.exists():
                raise FileNotFoundError(
                    f"文件不存在: {file_path_obj}\n"
                    f"请检查文件路径是否正确。已下载的文件在: {self.genomes_dir}"
                )
            
            return downloader.decompress_file(file_path_obj)
        
        # 索引构建工具
        def build_index_tool(genome_fasta: str, tool: str, **kwargs):
            # 处理路径：如果是相对路径，尝试在 genomes_dir 中查找
            genome_path_obj = Path(genome_fasta)
            if not genome_path_obj.is_absolute():
                # 相对路径，尝试在 genomes_dir 中查找
                potential_path = self.genomes_dir / genome_path_obj
                if potential_path.exists():
                    genome_path_obj = potential_path
                else:
                    # 如果不在 genomes_dir，尝试相对于工作目录
                    genome_path_obj = (self.work_dir / genome_path_obj).resolve()
            else:
                genome_path_obj = genome_path_obj.resolve()
            
            if not genome_path_obj.exists():
                raise FileNotFoundError(
                    f"参考基因组文件不存在: {genome_fasta}\n"
                    f"尝试的路径: {genome_path_obj}\n"
                    f"请检查文件路径是否正确。已下载的文件在: {self.genomes_dir}"
                )
            
            builder = IndexBuilder(genome_path_obj, output_dir=self.indices_dir)
            if tool == 'bwa':
                return builder.build_bwa_index(**kwargs)
            elif tool == 'bowtie2':
                return builder.build_bowtie2_index(**kwargs)
            elif tool == 'star':
                return builder.build_star_index(**kwargs)
            elif tool == 'hisat2':
                return builder.build_hisat2_index(**kwargs)
            elif tool == 'minimap2':
                return builder.build_minimap2_index(**kwargs)
            else:
                raise ValueError(f"不支持的索引工具: {tool}")
        
        tools['download'] = download_tool
        tools['download_from_url'] = download_from_url_tool
        tools['decompress'] = decompress_tool
        tools['build_index'] = build_index_tool
        
        return tools
    
    def plan_task(self, user_request: str) -> List[Task]:
        """
        使用LLM规划任务（统一使用LLM，规则模式已移除）
        
        Args:
            user_request: 用户请求（自然语言）
            
        Returns:
            任务列表
        """
        return self._plan_with_llm(user_request)
    
    def _plan_with_llm(self, user_request: str) -> List[Task]:
        """使用LLM规划任务（统一使用LLM，不再回退到规则模式）"""
        system_prompt = """你是一个参考基因组构建专家。你的任务是分析用户需求，规划出需要执行的任务步骤。

可用的工具：
1. download - 从预定义数据源下载参考基因组（仅支持human和mouse）
   - 参数: source (ncbi_refseq/ensembl/ucsc), species (human/mouse)
   - 适用场景：human 和 mouse 这两个内置物种
   
2. download_from_url - 从URL下载参考基因组（支持任意物种）
   - 参数: url (下载URL), output_filename (可选，输出文件名)
   - 适用场景：所有其他物种（如 e.coli、大肠杆菌、ecoli、酵母、果蝇等）
   - 你需要根据物种名称查找合适的下载URL
   - 常见数据源URL模式：
     * NCBI RefSeq: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../..._genomic.fna.gz
     * Ensembl: https://ftp.ensembl.org/pub/release-XXX/fasta/SPECIES/dna/...fa.gz
     * 例如 E.coli (NCBI): https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2/GCF_000005825.2_ASM582v2_genomic.fna.gz
     * 例如 E.coli (Ensembl): https://ftp.ensembl.org/pub/release-110/fasta/escherichia_coli_str_k_12_substr_mg1655_gca_000005825/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005825.ASM584v2.dna.toplevel.fa.gz
   
3. decompress - 解压缩文件
   - 参数: file_path
   - 注意：如果文件已经是 .fa 或 .fasta 格式（未压缩），不需要此步骤
   
4. build_index - 构建索引
   - 参数: 
     * genome_fasta (必需): 参考基因组FASTA文件路径
     * tool (必需): 索引工具名称 (bwa/bowtie2/star/hisat2/minimap2)
     * threads (可选): 线程数，默认4
     * sjdb_overhang (可选): 仅用于STAR，当提供了sjdb_gtf_file时使用，通常为read长度-1
     * sjdb_gtf_file (可选): 仅用于STAR，GTF/GFF注释文件路径
     * genome_sa_index_n_bases (可选): 仅用于STAR，基因组SA索引碱基数，默认14
   - 注意：
     * genome_fasta 可以是本地文件的绝对路径或相对路径
     * 对于STAR索引：
       - 如果提供了 sjdb_gtf_file，可以同时提供 sjdb_overhang（推荐用于RNA-seq）
       - 如果没有提供 sjdb_gtf_file，不能使用 sjdb_overhang 参数
       - 对于没有注释文件的基因组（如细菌），构建STAR索引时不要提供 sjdb_overhang

请仔细分析用户需求，规划出详细的任务步骤。每个任务应该包括：
- id: 任务ID（如 task_1, task_2）
- name: 任务名称
- description: 任务描述
- tool: 使用的工具名称
- parameters: 工具参数（字典格式）

重要规则：
1. 物种判断和工具选择：
   - 如果用户提到 "human"、"人类"、"homo sapiens" → 使用 download 工具，source="ncbi_refseq" 或 "ensembl"，species="human"
   - 如果用户提到 "mouse"、"小鼠"、"mus musculus" → 使用 download 工具，source="ncbi_refseq" 或 "ensembl"，species="mouse"
   - 如果用户提到其他物种（如 "e.coli"、"大肠杆菌"、"ecoli"、"yeast"、"酵母"、"drosophila"、"果蝇" 等）：
     * 优先使用 download_from_url 工具，查找并下载该物种的参考基因组
     * 你需要根据你的知识或搜索能力，找到该物种在 NCBI RefSeq 或 Ensembl 的下载URL
     * 如果无法找到URL，可以假设用户已有本地文件，使用环境变量或默认路径

2. URL查找指南（对于非 human/mouse 物种）：
   - NCBI RefSeq: 访问 https://www.ncbi.nlm.nih.gov/genome 搜索物种，找到参考基因组下载链接
   - Ensembl: 访问 https://www.ensembl.org 搜索物种，在物种页面找到 "Download DNA sequence (FASTA)"
   - 常见物种示例URL（你可以参考这些模式）：
     * E.coli K-12 (NCBI): https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2/GCF_000005825.2_ASM582v2_genomic.fna.gz
     * E.coli (Ensembl): https://ftp.ensembl.org/pub/release-110/fasta/escherichia_coli_str_k_12_substr_mg1655_gca_000005825/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005825.ASM584v2.dna.toplevel.fa.gz

3. 本地文件路径规则（仅当无法下载时使用）：
   - 优先使用环境变量：
     * REFBUILDER_ECOLI_GENOME（用于 e.coli）
     * REFBUILDER_CUSTOM_GENOME（通用自定义基因组）
     * REFBUILDER_TOY_GENOME（用于 toy/测试）
   - 如果没有环境变量，使用默认路径：
     * e.coli: "ecoli.fa"
     * toy/测试: "toy.fa"
     * 其他: 根据物种名称推断，如 "species_name.fa"

4. 任务规划规则：
   - 对于 human/mouse：task_1=download, task_2=decompress（如果需要）, task_3+=build_index
   - 对于其他物种（优先下载）：task_1=download_from_url, task_2=decompress（如果需要）, task_3+=build_index
   - 对于其他物种（无法下载，使用本地文件）：task_1=build_index（跳过 download 和 decompress）

5. 索引工具选择：
   - 如果用户没有指定，根据用途推荐：
     * WGS/WES: bwa 或 bowtie2
     * RNA-seq: star 或 hisat2
     * 长读长: minimap2
   
6. STAR索引特殊注意事项：
   - 对于有注释文件的真核生物（如human、mouse），构建STAR索引时可以：
     * 提供 sjdb_gtf_file 参数（GTF/GFF注释文件路径）
     * 提供 sjdb_overhang 参数（通常为read长度-1，如100表示101bp reads）
   - 对于没有注释文件的物种（如细菌、E.coli），构建STAR索引时：
     * 不要提供 sjdb_gtf_file 参数
     * 不要提供 sjdb_overhang 参数（STAR会报错）
     * 只提供基本参数：genome_fasta, tool="star", threads

6. 参数引用：
   - 使用 {{task_X.result}} 引用前一个任务的结果（用于 download/download_from_url → decompress → build_index 的流程）
   - 对于直接使用本地文件的情况，genome_fasta 应该使用实际文件路径字符串，不要使用 {{task_X.result}}"""
        
        prompt = f"""用户需求: {user_request}

⚠️ 关键提醒：
- 对于 human 或 mouse：使用 download 工具（source + species）
- 对于其他物种（如 e.coli、大肠杆菌、ecoli、yeast、酵母等）：
  * 优先使用 download_from_url 工具，提供该物种的下载URL
  * 你需要根据物种名称，查找 NCBI RefSeq 或 Ensembl 的下载链接
  * 如果无法找到URL，可以假设用户已有本地文件，使用环境变量或默认路径

请以JSON格式输出任务规划，格式如下：
{{
  "tasks": [
    {{
      "id": "task_1",
      "name": "任务名称",
      "description": "任务描述",
      "tool": "工具名称",
      "parameters": {{"参数名": "参数值"}}
    }}
  ]
}}

只返回JSON，不要包含其他文本。"""
        
        try:
            # 使用 generate_json 方法直接获取 JSON
            plan_data = self.llm_client.generate_json(
                prompt=prompt,
                system_prompt=system_prompt,
                max_tokens=2000
            )
            
            tasks = []
            for task_data in plan_data.get('tasks', []):
                task = Task(
                    id=task_data['id'],
                    name=task_data['name'],
                    description=task_data['description'],
                    status=TaskStatus.PENDING,
                    result=None
                )
                task.tool = task_data.get('tool')
                task.parameters = task_data.get('parameters', {})
                tasks.append(task)
            
            logger.info(f"LLM规划了 {len(tasks)} 个任务")
            return tasks
            
        except Exception as e:
            logger.error(f"LLM规划失败: {e}")
            raise ValueError(
                f"无法使用LLM规划任务: {e}\n"
                "请检查：\n"
                "1. LLM_API_KEY 是否正确设置\n"
                "2. 网络连接是否正常\n"
                "3. API 服务是否可用"
            )
    
    def _plan_with_rules(self, user_request: str) -> List[Task]:
        """
        基于规则的任务规划（已废弃，不再使用）
        
        注意：此方法已不再被调用，统一使用 LLM 进行任务规划。
        保留此方法仅用于参考。
        """
        tasks = []
        
        # 简单的关键词匹配
        request_lower = user_request.lower()
        
        # 特殊场景 1：用户明确提到 toy / 小基因组 / 测试 且要求不要下载
        # 这里我们假设用户已经在工作目录准备好了一个小的参考基因组文件，
        # 路径可以通过环境变量 REFBUILDER_TOY_GENOME 配置，默认使用 toy.fa。
        if (
            ('toy' in request_lower or '小基因组' in request_lower or 'mini' in request_lower or '测试' in request_lower)
            and ('不要下载' in request_lower or '不下载' in request_lower or '无需下载' in request_lower)
        ):
            genome_path = os.getenv("REFBUILDER_TOY_GENOME", "toy.fa")
            
            # 检测索引工具（默认还是 bwa）
            tools = []
            if 'bwa' in request_lower:
                tools.append('bwa')
            if 'bowtie2' in request_lower or 'bowtie' in request_lower:
                tools.append('bowtie2')
            if 'star' in request_lower:
                tools.append('star')
            if 'hisat2' in request_lower:
                tools.append('hisat2')
            if 'minimap2' in request_lower or 'minimap' in request_lower:
                tools.append('minimap2')
            if not tools:
                tools = ['bwa']
            
            for i, tool in enumerate(tools, 1):
                task = Task(
                    id=f'task_{i}',
                    name=f'构建{tool}索引（toy基因组）',
                    description=f'使用{tool}为本地 toy 基因组构建索引（{genome_path}）',
                    status=TaskStatus.PENDING
                )
                task.tool = 'build_index'
                task.parameters = {
                    'genome_fasta': genome_path,
                    'tool': tool,
                    'threads': 4
                }
                tasks.append(task)
            
            logger.info(f"基于规则规划了 {len(tasks)} 个 toy 基因组任务（不进行下载）")
            return tasks
        
        # 特殊场景 2：大肠杆菌 / E.coli 等非内置物种
        # 规则：不尝试从远程下载，而是认为用户已经在本地准备好 E.coli 参考基因组，
        # 使用环境变量 REFBUILDER_ECOLI_GENOME 或 REFBUILDER_CUSTOM_GENOME 指定路径，默认 ecoli.fa。
        if ('大肠杆菌' in request_lower) or ('e.coli' in request_lower) or ('ecoli' in request_lower):
            genome_path = (
                os.getenv("REFBUILDER_ECOLI_GENOME")
                or os.getenv("REFBUILDER_CUSTOM_GENOME")
                or "ecoli.fa"
            )
            
            # 工具选择：如果用户提到 star / hisat2，则认为是 RNA-seq；否则默认 bwa
            tools = []
            if 'star' in request_lower:
                tools.append('star')
            if 'hisat2' in request_lower:
                tools.append('hisat2')
            if 'bwa' in request_lower:
                tools.append('bwa')
            if not tools:
                # 根据是否提到 rna/转录组 来默认选择
                if 'rna' in request_lower or '转录组' in request_lower:
                    tools = ['star']
                else:
                    tools = ['bwa']
            
            for i, tool in enumerate(tools, 1):
                task = Task(
                    id=f'task_{i}',
                    name=f'构建{tool}索引（E.coli 基因组）',
                    description=f'使用{tool}为本地 E.coli 参考基因组构建索引（{genome_path}）',
                    status=TaskStatus.PENDING
                )
                task.tool = 'build_index'
                task.parameters = {
                    'genome_fasta': genome_path,
                    'tool': tool,
                    'threads': 4
                }
                tasks.append(task)
            
            logger.info(f"基于规则规划了 {len(tasks)} 个 E.coli 基因组任务（不进行下载）")
            return tasks
        
        # 一般场景：需要下载内置的人/小鼠参考基因组
        # 检测物种
        species = 'human'
        if 'mouse' in request_lower or '小鼠' in request_lower:
            species = 'mouse'
        
        # 检测数据源
        source = 'ncbi_refseq'
        if 'ensembl' in request_lower:
            source = 'ensembl'
        elif 'ucsc' in request_lower:
            source = 'ucsc'
        
        # 检测索引工具
        tools = []
        if 'bwa' in request_lower:
            tools.append('bwa')
        if 'bowtie2' in request_lower or 'bowtie' in request_lower:
            tools.append('bowtie2')
        if 'star' in request_lower:
            tools.append('star')
        if 'hisat2' in request_lower:
            tools.append('hisat2')
        if 'minimap2' in request_lower or 'minimap' in request_lower:
            tools.append('minimap2')
        
        # 如果没有指定工具，根据关键词推荐
        if not tools:
            if 'rna' in request_lower or '转录组' in request_lower:
                tools = ['star']
            elif '长读长' in request_lower or 'pacbio' in request_lower or 'nanopore' in request_lower:
                tools = ['minimap2']
            else:
                tools = ['bwa']  # 默认
        
        # 创建任务
        task1 = Task(
            id='task_1',
            name='下载参考基因组',
            description=f'从{source}下载{species}参考基因组',
            status=TaskStatus.PENDING
        )
        task1.tool = 'download'
        task1.parameters = {'source': source, 'species': species}
        tasks.append(task1)
        
        task2 = Task(
            id='task_2',
            name='解压缩文件',
            description='解压缩下载的参考基因组文件',
            status=TaskStatus.PENDING
        )
        task2.tool = 'decompress'
        task2.parameters = {'file_path': '{{task_1.result}}'}
        tasks.append(task2)
        
        # 为每个索引工具创建任务
        for i, tool in enumerate(tools, 3):
            task = Task(
                id=f'task_{i}',
                name=f'构建{tool}索引',
                description=f'使用{tool}构建参考基因组索引',
                status=TaskStatus.PENDING
            )
            task.tool = 'build_index'
            task.parameters = {
                'genome_fasta': '{{task_2.result}}',
                'tool': tool,
                'threads': 4
            }
            tasks.append(task)
        
        logger.info(f"基于规则规划了 {len(tasks)} 个任务")
        return tasks
    
    def recommend_solution(self, use_case: str) -> Dict[str, Any]:
        """
        根据使用场景推荐最佳方案（统一使用LLM）
        
        Args:
            use_case: 使用场景描述（如 "WGS分析", "RNA-seq分析"）
            
        Returns:
            推荐方案字典
        """
        return self._recommend_with_llm(use_case)
    
    def _recommend_with_llm(self, use_case: str) -> Dict[str, Any]:
        """使用LLM推荐方案（统一使用LLM，不再回退到规则模式）"""
        system_prompt = """你是一个生物信息学专家。你的任务是根据用户的使用场景，推荐最佳的参考基因组数据源和索引工具。

可用的数据源：
- ncbi_refseq: 官方参考序列数据库，最常用
- ensembl: 欧洲生物信息学研究所维护，适合RNA-seq分析
- ucsc: 加州大学圣克鲁兹分校维护

可用的索引工具：
- bwa: 短读长序列比对，适合WGS、WES、ChIP-seq
- bowtie2: 超快速短读长序列比对，适合WGS、WES
- star: RNA-seq比对的最佳工具
- hisat2: RNA-seq比对，内存占用更小
- minimap2: 长读长序列比对，适合PacBio、Nanopore数据

请根据使用场景，推荐最合适的数据源和工具，并提供理由。"""
        
        prompt = f"""使用场景: {use_case}

请以JSON格式输出推荐方案：
{{
  "source": "数据源名称",
  "source_reason": "选择理由",
  "tools": ["工具列表"],
  "tools_reason": "工具选择理由",
  "parameters": {{"threads": 8, "sjdb_overhang": 100}}
}}

只返回JSON，不要包含其他文本。"""
        
        try:
            return self.llm_client.generate_json(
                prompt=prompt,
                system_prompt=system_prompt,
                max_tokens=1000
            )
        except Exception as e:
            logger.error(f"LLM推荐失败: {e}")
            raise ValueError(
                f"无法使用LLM生成推荐方案: {e}\n"
                "请检查：\n"
                "1. LLM_API_KEY 是否正确设置\n"
                "2. 网络连接是否正常\n"
                "3. API 服务是否可用"
            )
    
    def _recommend_with_rules(self, use_case: str) -> Dict[str, Any]:
        """
        基于规则的推荐（已废弃，不再使用）
        
        注意：此方法已不再被调用，统一使用 LLM 进行推荐。
        保留此方法仅用于参考。
        """
        use_case_lower = use_case.lower()
        
        # 默认推荐
        recommendation = {
            'source': 'ncbi_refseq',
            'source_reason': 'NCBI RefSeq是官方参考序列数据库，最常用',
            'tools': ['bwa'],
            'tools_reason': 'BWA是短读长序列比对的标准工具',
            'parameters': {'threads': 8}
        }
        
        # 根据使用场景调整
        if 'rna' in use_case_lower or '转录组' in use_case_lower:
            recommendation['tools'] = ['star']
            recommendation['tools_reason'] = 'STAR是RNA-seq比对的最佳工具'
            recommendation['parameters']['sjdb_overhang'] = 100
        
        elif '长读长' in use_case_lower or 'pacbio' in use_case_lower or 'nanopore' in use_case_lower:
            recommendation['tools'] = ['minimap2']
            recommendation['tools_reason'] = 'minimap2是长读长序列比对的最佳工具'
        
        elif 'wes' in use_case_lower or '外显子' in use_case_lower:
            recommendation['tools'] = ['bwa']
            recommendation['tools_reason'] = 'BWA适合外显子组测序数据'
        
        return recommendation
    
    def execute_task(self, task: Task) -> Any:
        """
        执行单个任务
        
        Args:
            task: 任务对象
            
        Returns:
            任务执行结果
        """
        logger.info(f"执行任务: {task.name} ({task.id})")
        task.status = TaskStatus.IN_PROGRESS
        
        try:
            # 验证任务参数（download 工具只支持 human 和 mouse）
            if task.tool == 'download':
                species = task.parameters.get('species', '').lower()
                if species not in ['human', 'mouse']:
                    raise ValueError(
                        f"download 工具只支持 species='human' 或 'mouse'，不支持 '{species}'。\n"
                        f"对于非内置物种，请使用 download_from_url 工具，提供该物种的下载URL。\n"
                        f"或者如果用户已有本地文件，可以直接使用 build_index 工具。"
                    )
            
            # 解析参数中的任务引用
            parameters = self._resolve_parameters(task.parameters)
            
            # 调用工具
            tool_func = self.tools.get(task.tool)
            if not tool_func:
                raise ValueError(f"未知的工具: {task.tool}")
            
            result = tool_func(**parameters)
            
            task.status = TaskStatus.COMPLETED
            # 确保路径是绝对路径的字符串形式
            if isinstance(result, Path):
                task.result = str(result.resolve())
            else:
                task.result = result
            
            logger.info(f"任务完成: {task.name}")
            return result
            
        except Exception as e:
            task.status = TaskStatus.FAILED
            task.error = str(e)
            task.retry_count += 1
            
            logger.error(f"任务失败: {task.name} - {e}")
            
            # 自动重试
            if task.retry_count < task.max_retries:
                logger.info(f"重试任务: {task.name} (第{task.retry_count}次)")
                return self.execute_task(task)
            else:
                raise
    
    def _resolve_parameters(self, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """解析参数中的任务引用"""
        resolved = {}
        for key, value in parameters.items():
            if isinstance(value, str) and value.startswith('{{') and value.endswith('}}'):
                # 解析任务引用 {{task_X.result}}
                match = re.search(r'task_(\d+)\.result', value)
                if match:
                    task_id = f"task_{match.group(1)}"
                    task = self.state.get_task(task_id)
                    if task and task.result:
                        resolved[key] = task.result
                    else:
                        raise ValueError(f"无法解析任务引用: {value}")
                else:
                    resolved[key] = value
            else:
                resolved[key] = value
        return resolved
    
    def run(self, user_request: str) -> Dict[str, Any]:
        """
        运行智能体（主入口）
        
        Args:
            user_request: 用户请求（自然语言）
            
        Returns:
            执行结果
        """
        logger.info("=" * 80)
        logger.info("参考基因组构建智能体启动")
        logger.info("=" * 80)
        logger.info(f"用户请求: {user_request}")
        
        # 1. 规划任务
        logger.info("\n步骤1: 任务规划")
        tasks = self.plan_task(user_request)
        for task in tasks:
            self.state.add_task(task)
            logger.info(f"  - {task.name}: {task.description}")
        
        # 2. 执行任务
        logger.info("\n步骤2: 执行任务")
        results = {}
        for task in tasks:
            try:
                result = self.execute_task(task)
                results[task.id] = {
                    'status': task.status.value,
                    'result': task.result
                }
            except Exception as e:
                results[task.id] = {
                    'status': task.status.value,
                    'error': task.error
                }
                logger.error(f"任务 {task.id} 执行失败: {e}")
                # 根据错误类型决定是否继续
                if 'download' in task.tool:
                    logger.error("下载失败，无法继续后续任务")
                    break
        
        # 3. 保存状态
        state_file = self.work_dir / "agent_state.json"
        with open(state_file, 'w', encoding='utf-8') as f:
            json.dump(self.state.to_dict(), f, indent=2, ensure_ascii=False)
        
        logger.info(f"\n状态已保存到: {state_file}")
        logger.info("=" * 80)
        logger.info("智能体执行完成")
        logger.info("=" * 80)
        
        return {
            'tasks': results,
            'state': self.state.to_dict()
        }
