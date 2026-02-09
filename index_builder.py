#!/usr/bin/env python
"""
参考基因组索引构建模块

支持多种索引工具：
- BWA (Burrows-Wheeler Aligner)
- Bowtie2
- STAR (Spliced Transcripts Alignment to a Reference)
- HISAT2
- minimap2
"""

import os
import sys
import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any
import json

logger = logging.getLogger(__name__)


class IndexBuilder:
    """参考基因组索引构建器"""
    
    # 支持的索引工具及其命令
    SUPPORTED_TOOLS = {
        'bwa': {
            'command': 'bwa',
            'index_cmd': 'bwa index',
            'description': 'BWA (Burrows-Wheeler Aligner) - 用于短读长序列比对',
            'index_suffixes': ['.amb', '.ann', '.bwt', '.pac', '.sa']
        },
        'bowtie2': {
            'command': 'bowtie2-build',
            'index_cmd': 'bowtie2-build',
            'description': 'Bowtie2 - 超快速、内存高效的短读长序列比对工具',
            'index_suffixes': ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        },
        'star': {
            'command': 'STAR',
            'index_cmd': 'STAR --runMode genomeGenerate',
            'description': 'STAR - 用于RNA-seq数据的快速比对工具',
            'index_suffixes': ['/Genome', '/SA', '/SAindex']
        },
        'hisat2': {
            'command': 'hisat2-build',
            'index_cmd': 'hisat2-build',
            'description': 'HISAT2 - 用于RNA-seq数据的快速比对工具',
            'index_suffixes': ['.1.ht2', '.2.ht2', '.3.ht2', '.4.ht2', '.5.ht2', '.6.ht2', '.7.ht2', '.8.ht2']
        },
        'minimap2': {
            'command': 'minimap2',
            'index_cmd': 'minimap2 -d',
            'description': 'minimap2 - 用于长读长序列比对',
            'index_suffixes': ['.mmi']
        }
    }
    
    def __init__(self, genome_fasta: Path, output_dir: Optional[Path] = None):
        """
        初始化索引构建器
        
        Args:
            genome_fasta: 参考基因组FASTA文件路径
            output_dir: 索引输出目录，如果为None则使用genome_fasta所在目录
        """
        self.genome_fasta = Path(genome_fasta)
        
        if not self.genome_fasta.exists():
            raise FileNotFoundError(f"参考基因组文件不存在: {self.genome_fasta}")
        
        if output_dir is None:
            self.output_dir = self.genome_fasta.parent
        else:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"参考基因组文件: {self.genome_fasta}")
        logger.info(f"索引输出目录: {self.output_dir}")
    
    def check_tool_available(self, tool: str) -> bool:
        """
        检查索引工具是否可用
        
        Args:
            tool: 工具名称
            
        Returns:
            工具是否可用
        """
        if tool not in self.SUPPORTED_TOOLS:
            logger.error(f"不支持的工具: {tool}")
            return False
        
        tool_info = self.SUPPORTED_TOOLS[tool]
        command = tool_info['command']
        
        try:
            # 检查命令是否可用
            if tool == 'star':
                result = subprocess.run([command, '--version'], 
                                      capture_output=True, text=True, timeout=10)
            elif tool == 'bwa':
                result = subprocess.run([command], 
                                      capture_output=True, text=True, timeout=10)
            else:
                result = subprocess.run([command, '--version'], 
                                      capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0 or 'version' in result.stderr.lower() or 'version' in result.stdout.lower():
                logger.info(f"工具 {tool} 可用")
                return True
            else:
                logger.warning(f"工具 {tool} 可能不可用")
                return False
        except FileNotFoundError:
            logger.warning(f"工具 {tool} 未安装或不在PATH中")
            return False
        except subprocess.TimeoutExpired:
            logger.warning(f"检查工具 {tool} 超时")
            return False
        except Exception as e:
            logger.warning(f"检查工具 {tool} 时出错: {e}")
            return False
    
    def build_bwa_index(self, prefix: Optional[str] = None, threads: int = 1) -> Path:
        """
        构建BWA索引
        
        Args:
            prefix: 索引文件前缀，如果为None则使用genome_fasta的basename
            threads: 线程数
            
        Returns:
            索引文件路径（.bwt文件）
        """
        if not self.check_tool_available('bwa'):
            raise RuntimeError("BWA 工具不可用，请先安装 BWA")
        
        if prefix is None:
            prefix = self.genome_fasta.stem
        
        index_prefix = self.output_dir / prefix
        
        # 检查索引是否已存在
        bwt_file = index_prefix.with_suffix('.bwt')
        if bwt_file.exists():
            logger.info(f"BWA索引已存在: {bwt_file}")
            return bwt_file
        
        logger.info(f"开始构建BWA索引: {self.genome_fasta}")
        logger.info(f"索引前缀: {index_prefix}")
        
        try:
            cmd = ['bwa', 'index', '-p', str(index_prefix), str(self.genome_fasta)]
            if threads > 1:
                cmd.extend(['-t', str(threads)])
            
            logger.info(f"执行命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            logger.info("BWA索引构建完成")
            logger.debug(f"输出: {result.stdout}")
            
            return bwt_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"BWA索引构建失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            raise
    
    def build_bowtie2_index(self, prefix: Optional[str] = None, threads: int = 1) -> Path:
        """
        构建Bowtie2索引
        
        Args:
            prefix: 索引文件前缀，如果为None则使用genome_fasta的basename
            threads: 线程数
            
        Returns:
            索引文件路径（.1.bt2文件）
        """
        if not self.check_tool_available('bowtie2'):
            raise RuntimeError("Bowtie2 工具不可用，请先安装 Bowtie2")
        
        if prefix is None:
            prefix = self.genome_fasta.stem
        
        index_prefix = self.output_dir / prefix
        
        # 检查索引是否已存在
        bt2_file = index_prefix.with_suffix('.1.bt2')
        if bt2_file.exists():
            logger.info(f"Bowtie2索引已存在: {bt2_file}")
            return bt2_file
        
        logger.info(f"开始构建Bowtie2索引: {self.genome_fasta}")
        logger.info(f"索引前缀: {index_prefix}")
        
        try:
            cmd = ['bowtie2-build', '--threads', str(threads), 
                   str(self.genome_fasta), str(index_prefix)]
            
            logger.info(f"执行命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            logger.info("Bowtie2索引构建完成")
            logger.debug(f"输出: {result.stdout}")
            
            return bt2_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Bowtie2索引构建失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            raise
    
    def build_star_index(self, prefix: Optional[str] = None, 
                        sjdb_overhang: Optional[int] = None, 
                        sjdb_gtf_file: Optional[str] = None,
                        threads: int = 1,
                        genome_sa_index_n_bases: int = 14) -> Path:
        """
        构建STAR索引
        
        Args:
            prefix: 索引目录名称，如果为None则使用genome_fasta的basename + "_star"
            sjdb_overhang: splice junction database overhang（通常为read长度-1）
                         注意：只有当提供了 sjdb_gtf_file 时才需要此参数
            sjdb_gtf_file: GTF/GFF注释文件路径（可选）
                         如果提供，STAR会使用注释信息构建splice junction database
            threads: 线程数
            genome_sa_index_n_bases: 基因组SA索引的碱基数（对于大基因组需要调整）
            
        Returns:
            索引目录路径
        """
        if not self.check_tool_available('star'):
            raise RuntimeError("STAR 工具不可用，请先安装 STAR")
        
        if prefix is None:
            prefix = self.genome_fasta.stem + "_star"
        
        index_dir = self.output_dir / prefix
        index_dir.mkdir(parents=True, exist_ok=True)
        
        # 检查索引是否已存在
        genome_file = index_dir / "Genome"
        if genome_file.exists():
            logger.info(f"STAR索引已存在: {index_dir}")
            return index_dir
        
        logger.info(f"开始构建STAR索引: {self.genome_fasta}")
        logger.info(f"索引目录: {index_dir}")
        
        # 构建命令
        cmd = [
            'STAR',
            '--runMode', 'genomeGenerate',
            '--genomeDir', str(index_dir),
            '--genomeFastaFiles', str(self.genome_fasta),
            '--runThreadN', str(threads),
            '--genomeSAindexNbases', str(genome_sa_index_n_bases)
        ]
        
        # 只有当提供了GTF文件时，才添加sjdb相关参数
        if sjdb_gtf_file:
            gtf_path = Path(sjdb_gtf_file)
            if not gtf_path.exists():
                raise FileNotFoundError(f"GTF注释文件不存在: {sjdb_gtf_file}")
            
            cmd.extend(['--sjdbGTFfile', str(gtf_path)])
            
            # 如果提供了sjdb_overhang，添加该参数
            if sjdb_overhang is not None:
                cmd.extend(['--sjdbOverhang', str(sjdb_overhang)])
                logger.info(f"SJDB overhang: {sjdb_overhang}")
            else:
                # 如果没有提供但需要GTF，使用默认值
                cmd.extend(['--sjdbOverhang', '100'])
                logger.info(f"SJDB overhang: 100 (默认值)")
            
            logger.info(f"使用GTF注释文件: {sjdb_gtf_file}")
        else:
            # 没有GTF文件时，不能使用sjdbOverhang
            if sjdb_overhang is not None:
                logger.warning(
                    f"未提供GTF注释文件，忽略 sjdb_overhang={sjdb_overhang} 参数。"
                    "STAR在没有注释文件时不能使用 --sjdbOverhang 参数。"
                )
        
        try:
            logger.info(f"执行命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            logger.info("STAR索引构建完成")
            logger.debug(f"输出: {result.stdout}")
            
            return index_dir
            
        except subprocess.CalledProcessError as e:
            logger.error(f"STAR索引构建失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            # 清理失败的索引目录
            if index_dir.exists():
                shutil.rmtree(index_dir)
            raise
    
    def build_hisat2_index(self, prefix: Optional[str] = None, threads: int = 1) -> Path:
        """
        构建HISAT2索引
        
        Args:
            prefix: 索引文件前缀，如果为None则使用genome_fasta的basename
            threads: 线程数
            
        Returns:
            索引文件路径（.1.ht2文件）
        """
        if not self.check_tool_available('hisat2'):
            raise RuntimeError("HISAT2 工具不可用，请先安装 HISAT2")
        
        if prefix is None:
            prefix = self.genome_fasta.stem
        
        index_prefix = self.output_dir / prefix
        
        # 检查索引是否已存在
        ht2_file = index_prefix.with_suffix('.1.ht2')
        if ht2_file.exists():
            logger.info(f"HISAT2索引已存在: {ht2_file}")
            return ht2_file
        
        logger.info(f"开始构建HISAT2索引: {self.genome_fasta}")
        logger.info(f"索引前缀: {index_prefix}")
        
        try:
            cmd = ['hisat2-build', '-p', str(threads), 
                   str(self.genome_fasta), str(index_prefix)]
            
            logger.info(f"执行命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            logger.info("HISAT2索引构建完成")
            logger.debug(f"输出: {result.stdout}")
            
            return ht2_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"HISAT2索引构建失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            raise
    
    def build_minimap2_index(self, prefix: Optional[str] = None) -> Path:
        """
        构建minimap2索引
        
        Args:
            prefix: 索引文件名，如果为None则使用genome_fasta的basename + ".mmi"
            
        Returns:
            索引文件路径
        """
        if not self.check_tool_available('minimap2'):
            raise RuntimeError("minimap2 工具不可用，请先安装 minimap2")
        
        if prefix is None:
            index_file = self.output_dir / (self.genome_fasta.stem + ".mmi")
        else:
            index_file = self.output_dir / prefix
            if not index_file.suffix == '.mmi':
                index_file = index_file.with_suffix('.mmi')
        
        # 检查索引是否已存在
        if index_file.exists():
            logger.info(f"minimap2索引已存在: {index_file}")
            return index_file
        
        logger.info(f"开始构建minimap2索引: {self.genome_fasta}")
        logger.info(f"索引文件: {index_file}")
        
        try:
            cmd = ['minimap2', '-d', str(index_file), str(self.genome_fasta)]
            
            logger.info(f"执行命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            logger.info("minimap2索引构建完成")
            logger.debug(f"输出: {result.stdout}")
            
            return index_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"minimap2索引构建失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            raise
    
    def build_all_indices(self, tools: Optional[List[str]] = None, 
                         threads: int = 1, **kwargs) -> Dict[str, Any]:
        """
        构建所有支持的索引（或指定的索引工具）
        
        Args:
            tools: 要构建的索引工具列表，如果为None则构建所有可用工具的索引
            threads: 线程数
            **kwargs: 其他参数（如sjdb_overhang用于STAR）
            
        Returns:
            构建结果字典
        """
        if tools is None:
            tools = list(self.SUPPORTED_TOOLS.keys())
        
        results = {}
        
        for tool in tools:
            try:
                if tool == 'bwa':
                    index_path = self.build_bwa_index(threads=threads)
                elif tool == 'bowtie2':
                    index_path = self.build_bowtie2_index(threads=threads)
                elif tool == 'star':
                    sjdb_overhang = kwargs.get('sjdb_overhang', 100)
                    genome_sa_index_n_bases = kwargs.get('genome_sa_index_n_bases', 14)
                    index_path = self.build_star_index(
                        sjdb_overhang=sjdb_overhang,
                        threads=threads,
                        genome_sa_index_n_bases=genome_sa_index_n_bases
                    )
                elif tool == 'hisat2':
                    index_path = self.build_hisat2_index(threads=threads)
                elif tool == 'minimap2':
                    index_path = self.build_minimap2_index()
                else:
                    logger.warning(f"未知的工具: {tool}")
                    continue
                
                results[tool] = {
                    'status': 'success',
                    'index_path': str(index_path)
                }
                
            except Exception as e:
                logger.error(f"构建 {tool} 索引失败: {e}")
                results[tool] = {
                    'status': 'failed',
                    'error': str(e)
                }
        
        return results
    
    def list_available_tools(self) -> List[str]:
        """
        列出所有可用的索引工具
        
        Returns:
            可用工具列表
        """
        available = []
        for tool in self.SUPPORTED_TOOLS.keys():
            if self.check_tool_available(tool):
                available.append(tool)
        return available
