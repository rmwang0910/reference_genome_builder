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

# 尝试导入 LLM 相关模块
try:
    # 添加 BioLitKG 路径（如果存在）
    biolitkg_path = Path(__file__).parent.parent.parent / "AI" / "BioLitKG"
    if biolitkg_path.exists() and str(biolitkg_path) not in sys.path:
        sys.path.insert(0, str(biolitkg_path))
    
    from core.llm.openai import OpenAIProvider
    from core.config import get_config
    HAS_LLM = True
except ImportError:
    HAS_LLM = False
    logger.warning("LLM模块不可用，将使用基于规则的决策")


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
            use_llm: 是否使用LLM进行智能决策
        """
        if work_dir is None:
            self.work_dir = Path.cwd()
        else:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)
        
        self.genomes_dir = self.work_dir / "genomes"
        self.indices_dir = self.work_dir / "indices"
        
        # 初始化状态
        self.state = AgentState()
        
        # 初始化LLM（如果可用）
        self.use_llm = use_llm and HAS_LLM
        self.llm_client = None
        
        if self.use_llm:
            try:
                config = get_config()
                if config.llm.api_key:
                    self.llm_client = OpenAIProvider(config.llm)
                    logger.info("LLM已初始化，启用智能决策")
                else:
                    logger.warning("未设置 LLM_API_KEY，将使用基于规则的决策")
                    self.use_llm = False
            except Exception as e:
                logger.warning(f"无法初始化LLM: {e}，将使用基于规则的决策")
                self.use_llm = False
        
        # 注册工具
        self.tools = self._register_tools()
        
        logger.info(f"智能体初始化完成，工作目录: {self.work_dir}")
        logger.info(f"LLM模式: {'启用' if self.use_llm else '禁用'}")
    
    def _register_tools(self) -> Dict[str, Callable]:
        """注册可用工具"""
        from download import ReferenceGenomeDownloader
        from index_builder import IndexBuilder
        
        tools = {}
        
        # 下载工具
        def download_tool(source: str, species: str, **kwargs):
            downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
            return downloader.download_from_source(source, species)
        
        # 解压缩工具
        def decompress_tool(file_path: str, **kwargs):
            downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
            return downloader.decompress_file(Path(file_path))
        
        # 索引构建工具
        def build_index_tool(genome_fasta: str, tool: str, **kwargs):
            builder = IndexBuilder(Path(genome_fasta), output_dir=self.indices_dir)
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
        tools['decompress'] = decompress_tool
        tools['build_index'] = build_index_tool
        
        return tools
    
    def plan_task(self, user_request: str) -> List[Task]:
        """
        使用LLM规划任务
        
        Args:
            user_request: 用户请求（自然语言）
            
        Returns:
            任务列表
        """
        if self.use_llm:
            return self._plan_with_llm(user_request)
        else:
            return self._plan_with_rules(user_request)
    
    def _plan_with_llm(self, user_request: str) -> List[Task]:
        """使用LLM规划任务"""
        prompt = f"""你是一个参考基因组构建专家。根据用户的需求，规划出需要执行的任务步骤。

用户需求: {user_request}

可用的工具：
1. download - 下载参考基因组
   - 参数: source (ncbi_refseq/ensembl/ucsc), species (human/mouse)
2. decompress - 解压缩文件
   - 参数: file_path
3. build_index - 构建索引
   - 参数: genome_fasta, tool (bwa/bowtie2/star/hisat2/minimap2), threads

请分析用户需求，规划出详细的任务步骤。每个任务应该包括：
- id: 任务ID（如 task_1, task_2）
- name: 任务名称
- description: 任务描述
- tool: 使用的工具名称
- parameters: 工具参数（字典格式）

请以JSON格式输出，格式如下：
{{
  "tasks": [
    {{
      "id": "task_1",
      "name": "下载参考基因组",
      "description": "从NCBI RefSeq下载人类参考基因组",
      "tool": "download",
      "parameters": {{"source": "ncbi_refseq", "species": "human"}}
    }},
    {{
      "id": "task_2",
      "name": "解压缩文件",
      "description": "解压缩下载的参考基因组文件",
      "tool": "decompress",
      "parameters": {{"file_path": "{{task_1.result}}"}}
    }}
  ]
}}

注意：
- 如果用户没有指定数据源，推荐使用 ncbi_refseq（最常用）
- 如果用户没有指定索引工具，根据用途推荐：
  - WGS/WES: bwa 或 bowtie2
  - RNA-seq: star 或 hisat2
  - 长读长: minimap2
- 使用 {{task_X.result}} 引用前一个任务的结果
"""
        
        try:
            response = self.llm_client.generate(prompt, max_tokens=2000)
            
            # 解析JSON响应
            json_match = re.search(r'\{.*\}', response, re.DOTALL)
            if json_match:
                plan_data = json.loads(json_match.group())
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
            logger.error(f"LLM规划失败: {e}，使用基于规则的规划")
            return self._plan_with_rules(user_request)
    
    def _plan_with_rules(self, user_request: str) -> List[Task]:
        """基于规则的任务规划"""
        tasks = []
        
        # 简单的关键词匹配
        request_lower = user_request.lower()
        
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
        根据使用场景推荐最佳方案
        
        Args:
            use_case: 使用场景描述（如 "WGS分析", "RNA-seq分析"）
            
        Returns:
            推荐方案字典
        """
        if self.use_llm:
            return self._recommend_with_llm(use_case)
        else:
            return self._recommend_with_rules(use_case)
    
    def _recommend_with_llm(self, use_case: str) -> Dict[str, Any]:
        """使用LLM推荐方案"""
        prompt = f"""你是一个生物信息学专家。根据用户的使用场景，推荐最佳的参考基因组数据源和索引工具。

使用场景: {use_case}

请推荐：
1. 最佳数据源（ncbi_refseq/ensembl/ucsc）及理由
2. 最佳索引工具（bwa/bowtie2/star/hisat2/minimap2）及理由
3. 推荐的参数设置

请以JSON格式输出：
{{
  "source": "ncbi_refseq",
  "source_reason": "理由",
  "tools": ["bwa"],
  "tools_reason": "理由",
  "parameters": {{"threads": 8, "sjdb_overhang": 100}}
}}
"""
        
        try:
            response = self.llm_client.generate(prompt, max_tokens=1000)
            json_match = re.search(r'\{.*\}', response, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
        except Exception as e:
            logger.error(f"LLM推荐失败: {e}")
        
        return self._recommend_with_rules(use_case)
    
    def _recommend_with_rules(self, use_case: str) -> Dict[str, Any]:
        """基于规则的推荐"""
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
            # 解析参数中的任务引用
            parameters = self._resolve_parameters(task.parameters)
            
            # 调用工具
            tool_func = self.tools.get(task.tool)
            if not tool_func:
                raise ValueError(f"未知的工具: {task.tool}")
            
            result = tool_func(**parameters)
            
            task.status = TaskStatus.COMPLETED
            task.result = str(result) if isinstance(result, Path) else result
            
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
