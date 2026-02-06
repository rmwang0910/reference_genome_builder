#!/usr/bin/env python
"""
参考基因组构建智能体主程序

提供命令行接口，支持：
1. 下载参考基因组
2. 构建参考基因组索引
"""

import os
import sys
import logging
import argparse
import json
from pathlib import Path
from typing import Optional, List

# 配置logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from download import ReferenceGenomeDownloader
from index_builder import IndexBuilder


class ReferenceGenomeAgent:
    """参考基因组构建智能体"""
    
    def __init__(self, work_dir: Optional[str] = None):
        """
        初始化智能体
        
        Args:
            work_dir: 工作目录，如果为None则使用当前目录
        """
        if work_dir is None:
            self.work_dir = Path.cwd()
        else:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # 统一结果根目录：work_dir/results
        self.results_root = self.work_dir / "results"
        self.results_root.mkdir(parents=True, exist_ok=True)
        
        # 兼容旧逻辑的回退目录（不推荐直接使用）
        self.genomes_dir = self.work_dir / "genomes"
        self.indices_dir = self.work_dir / "indices"
        self.annotations_dir = self.work_dir / "annotations"
        
        # 当前这次运行对应的结果子目录，如 results/human_GRCh38
        self.current_run_root: Optional[Path] = None
        
        logger.info(f"工作目录: {self.work_dir}")
        logger.info(f"结果根目录: {self.results_root}")
    
    def _get_run_name(self, source: str, species: str) -> str:
        """
        根据数据源和物种生成统一的结果目录名，例如：human_GRCh38
        """
        try:
            genome_info = ReferenceGenomeDownloader.DATA_SOURCES[source.lower()][species.lower()]
            return f"{species}_{genome_info['name']}"
        except Exception:
            # 找不到就退回一个通用名字
            return f"{species}_{source}"
    
    def download_genome(self, source: str, species: str, 
                       decompress: bool = True) -> Path:
        """
        下载参考基因组
        
        Args:
            source: 数据源（'ncbi_refseq', 'ensembl', 'ucsc'）
            species: 物种（'human', 'mouse'）
            decompress: 是否自动解压缩
            
        Returns:
            参考基因组FASTA文件路径
        """
        logger.info("=" * 80)
        logger.info("参考基因组下载")
        logger.info("=" * 80)
        
        # 为本次运行创建统一结果目录：results/human_GRCh38 之类
        run_name = self._get_run_name(source, species)
        run_root = self.results_root / run_name
        genome_dir = run_root / "genome"
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        # 更新当前运行上下文
        self.current_run_root = run_root
        self.genomes_dir = genome_dir
        self.indices_dir = run_root / "indices"
        self.annotations_dir = run_root / "annotation"
        self.indices_dir.mkdir(parents=True, exist_ok=True)
        self.annotations_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"本次运行目录: {run_root}")
        logger.info(f"基因组目录: {self.genomes_dir}")
        logger.info(f"索引目录: {self.indices_dir}")
        logger.info(f"注释目录: {self.annotations_dir}")
        
        downloader = ReferenceGenomeDownloader(output_dir=str(self.genomes_dir))
        
        # 下载压缩文件
        compressed_file = downloader.download_from_source(source, species)
        logger.info(f"下载完成: {compressed_file}")
        
        # 如果需要，解压缩
        if decompress and compressed_file.suffix == '.gz':
            genome_fasta = downloader.decompress_file(compressed_file)
            logger.info(f"解压缩完成: {genome_fasta}")
            return genome_fasta
        else:
            return compressed_file
    
    def download_annotation(self, source: str, species: str,
                           fmt: str = "gtf",
                           decompress: bool = True) -> Path:
        """
        下载注释文件（GTF/GFF）
        
        Args:
            source: 数据源（'ncbi_refseq', 'ensembl', 'ucsc'）
            species: 物种（'human', 'mouse'）
            fmt: 注释格式（'gtf' 或 'gff3'）
            decompress: 是否自动解压缩
            
        Returns:
            注释文件路径（解压后）
        """
        logger.info("=" * 80)
        logger.info("基因注释下载")
        logger.info("=" * 80)
        
        # 如果还没有 current_run_root（例如单独下载注释），也根据 source/species 建一个
        if self.current_run_root is None:
            run_name = self._get_run_name(source, species)
            self.current_run_root = self.results_root / run_name
        
        self.annotations_dir = self.current_run_root / "annotation"
        self.annotations_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"注释目录: {self.annotations_dir}")
        
        downloader = ReferenceGenomeDownloader(output_dir=str(self.annotations_dir))
        
        # 下载压缩文件
        compressed_file = downloader.download_annotation_from_source(source, species, fmt=fmt)
        logger.info(f"下载完成: {compressed_file}")
        
        # 如果需要，解压缩
        if decompress and compressed_file.suffix == '.gz':
            annotation_path = downloader.decompress_file(compressed_file)
            logger.info(f"解压缩完成: {annotation_path}")
            return annotation_path
        else:
            return compressed_file
    
    def build_index(self, genome_fasta: Path, tools: List[str], 
                   threads: int = 1, **kwargs) -> dict:
        """
        构建参考基因组索引
        
        Args:
            genome_fasta: 参考基因组FASTA文件路径
            tools: 要构建的索引工具列表
            threads: 线程数
            **kwargs: 其他参数
            
        Returns:
            构建结果字典
        """
        logger.info("=" * 80)
        logger.info("参考基因组索引构建")
        logger.info("=" * 80)
        
        # 如果 genome_fasta 位于 results/<run_name>/genome 下，则将索引放在对应的 results/<run_name>/indices 下
        try:
            genome_fasta_resolved = genome_fasta.resolve()
        except Exception:
            genome_fasta_resolved = genome_fasta
        
        if self.results_root in genome_fasta_resolved.parents:
            rel = genome_fasta_resolved.relative_to(self.results_root)
            # 期望结构：<run_name>/genome/...
            run_name = rel.parts[0]
            run_root = self.results_root / run_name
            self.current_run_root = run_root
            self.indices_dir = run_root / "indices"
            self.indices_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"检测到结果目录: {run_root}，索引将写入: {self.indices_dir}")
        
        builder = IndexBuilder(genome_fasta, output_dir=self.indices_dir)
        
        # 检查可用工具
        available_tools = builder.list_available_tools()
        logger.info(f"可用的索引工具: {available_tools}")
        
        # 过滤出可用的工具
        tools_to_build = [t for t in tools if t in available_tools]
        unavailable_tools = [t for t in tools if t not in available_tools]
        
        if unavailable_tools:
            logger.warning(f"以下工具不可用，将跳过: {unavailable_tools}")
        
        if not tools_to_build:
            raise RuntimeError("没有可用的索引工具")
        
        # 构建索引
        results = builder.build_all_indices(tools=tools_to_build, threads=threads, **kwargs)
        
        # 打印结果摘要
        logger.info("\n索引构建结果:")
        for tool, result in results.items():
            if result['status'] == 'success':
                logger.info(f"  ✓ {tool}: {result['index_path']}")
            else:
                logger.error(f"  ✗ {tool}: {result.get('error', '未知错误')}")
        
        return results
    
    def run_pipeline(self, source: str, species: str, tools: List[str],
                    decompress: bool = True, threads: int = 1,
                    with_annotation: bool = False,
                    annotation_format: str = "gtf",
                    **kwargs) -> dict:
        """
        运行完整流程：下载 + 构建索引
        
        Args:
            source: 数据源
            species: 物种
            tools: 索引工具列表
            decompress: 是否解压缩
            threads: 线程数
            **kwargs: 其他参数
            
        Returns:
            完整流程结果字典
        """
        logger.info("=" * 80)
        logger.info("参考基因组构建完整流程")
        logger.info("=" * 80)
        
        # 1. 下载参考基因组
        genome_fasta = self.download_genome(source, species, decompress=decompress)
        
        # 2. （可选）下载注释文件
        annotation_path = None
        if with_annotation:
            try:
                annotation_path = self.download_annotation(
                    source,
                    species,
                    fmt=annotation_format,
                    decompress=decompress
                )
            except Exception as e:
                logger.warning(f"注释文件下载失败（不会中断索引构建）: {e}")
                annotation_path = None
        
        # 3. 构建索引
        index_results = self.build_index(genome_fasta, tools, threads=threads, **kwargs)
        
        # 4. 汇总结果
        results = {
            'genome_fasta': str(genome_fasta),
            'index_results': index_results,
            'work_dir': str(self.work_dir)
        }
        if annotation_path:
            results['annotation'] = str(annotation_path)
        
        # 保存结果到JSON文件
        results_file = self.work_dir / "build_results.json"
        with open(results_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        logger.info(f"\n结果已保存到: {results_file}")
        logger.info("=" * 80)
        logger.info("流程完成！")
        logger.info("=" * 80)
        
        return results


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='参考基因组构建智能体',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 下载人类参考基因组
  python main.py download --source ncbi_refseq --species human

  # 构建BWA索引
  python main.py index --genome genomes/human_GRCh38.p14.fa --tools bwa --threads 8

  # 完整流程：下载并构建多个索引
  python main.py pipeline --source ensembl --species human --tools bwa bowtie2 star --threads 8

  # 列出可用的参考基因组
  python main.py list-genomes

  # 列出可用的索引工具
  python main.py list-tools
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='子命令')
    
    # download 命令
    download_parser = subparsers.add_parser('download', help='下载参考基因组')
    download_parser.add_argument('--source', required=True, 
                               choices=['ncbi_refseq', 'ensembl', 'ucsc'],
                               help='数据源')
    download_parser.add_argument('--species', required=True,
                               choices=['human', 'mouse'],
                               help='物种')
    download_parser.add_argument('--work-dir', default='.',
                               help='工作目录（默认: 当前目录）')
    download_parser.add_argument('--no-decompress', action='store_true',
                               help='不解压缩下载的文件')
    
    # download-annotation 命令
    download_ann_parser = subparsers.add_parser('download-annotation', help='下载基因注释文件（GTF/GFF）')
    download_ann_parser.add_argument('--source', required=True,
                                     choices=['ncbi_refseq', 'ensembl', 'ucsc'],
                                     help='数据源')
    download_ann_parser.add_argument('--species', required=True,
                                     choices=['human', 'mouse'],
                                     help='物种')
    download_ann_parser.add_argument('--format', dest='fmt', default='gtf',
                                     choices=['gtf', 'gff3'],
                                     help='注释格式（默认: gtf）')
    download_ann_parser.add_argument('--work-dir', default='.',
                                     help='工作目录（默认: 当前目录）')
    download_ann_parser.add_argument('--no-decompress', action='store_true',
                                     help='不解压缩下载的文件')
    
    # index 命令
    index_parser = subparsers.add_parser('index', help='构建参考基因组索引')
    index_parser.add_argument('--genome', required=True,
                            help='参考基因组FASTA文件路径')
    index_parser.add_argument('--tools', required=True, nargs='+',
                            choices=['bwa', 'bowtie2', 'star', 'hisat2', 'minimap2'],
                            help='要构建的索引工具')
    index_parser.add_argument('--work-dir', default='.',
                            help='工作目录（默认: 当前目录）')
    index_parser.add_argument('--threads', type=int, default=1,
                            help='线程数（默认: 1）')
    index_parser.add_argument('--sjdb-overhang', type=int, default=100,
                            help='STAR索引的SJDB overhang（默认: 100）')
    index_parser.add_argument('--genome-sa-index-n-bases', type=int, default=14,
                            help='STAR索引的genomeSAindexNbases（默认: 14）')
    
    # pipeline 命令
    pipeline_parser = subparsers.add_parser('pipeline', help='完整流程：下载并构建索引')
    pipeline_parser.add_argument('--source', required=True,
                               choices=['ncbi_refseq', 'ensembl', 'ucsc'],
                               help='数据源')
    pipeline_parser.add_argument('--species', required=True,
                               choices=['human', 'mouse'],
                               help='物种')
    pipeline_parser.add_argument('--tools', required=True, nargs='+',
                               choices=['bwa', 'bowtie2', 'star', 'hisat2', 'minimap2'],
                               help='要构建的索引工具')
    pipeline_parser.add_argument('--work-dir', default='.',
                               help='工作目录（默认: 当前目录）')
    pipeline_parser.add_argument('--threads', type=int, default=1,
                               help='线程数（默认: 1）')
    pipeline_parser.add_argument('--no-decompress', action='store_true',
                               help='不解压缩下载的文件')
    pipeline_parser.add_argument('--with-annotation', action='store_true',
                               help='同时下载注释文件（GTF/GFF）')
    pipeline_parser.add_argument('--annotation-format', default='gtf',
                               choices=['gtf', 'gff3'],
                               help='注释格式（默认: gtf）')
    pipeline_parser.add_argument('--sjdb-overhang', type=int, default=100,
                               help='STAR索引的SJDB overhang（默认: 100）')
    pipeline_parser.add_argument('--genome-sa-index-n-bases', type=int, default=14,
                               help='STAR索引的genomeSAindexNbases（默认: 14）')
    
    # list-genomes 命令
    list_genomes_parser = subparsers.add_parser('list-genomes', 
                                                help='列出可用的参考基因组')
    list_genomes_parser.add_argument('--species',
                                    choices=['human', 'mouse'],
                                    help='过滤特定物种')
    
    # list-tools 命令
    list_tools_parser = subparsers.add_parser('list-tools',
                                              help='列出可用的索引工具')
    list_tools_parser.add_argument('--genome',
                                  help='参考基因组FASTA文件路径（用于检查工具）')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.command == 'download':
            agent = ReferenceGenomeAgent(work_dir=args.work_dir)
            genome_fasta = agent.download_genome(
                args.source, 
                args.species,
                decompress=not args.no_decompress
            )
            logger.info(f"\n参考基因组已下载到: {genome_fasta}")
        
        elif args.command == 'download-annotation':
            agent = ReferenceGenomeAgent(work_dir=args.work_dir)
            ann_path = agent.download_annotation(
                args.source,
                args.species,
                fmt=args.fmt,
                decompress=not args.no_decompress
            )
            logger.info(f"\n注释文件已下载到: {ann_path}")
        
        elif args.command == 'index':
            agent = ReferenceGenomeAgent(work_dir=args.work_dir)
            genome_fasta = Path(args.genome)
            results = agent.build_index(
                genome_fasta,
                args.tools,
                threads=args.threads,
                sjdb_overhang=args.sjdb_overhang,
                genome_sa_index_n_bases=args.genome_sa_index_n_bases
            )
            logger.info(f"\n索引构建完成，结果: {results}")
        
        elif args.command == 'pipeline':
            agent = ReferenceGenomeAgent(work_dir=args.work_dir)
            results = agent.run_pipeline(
                args.source,
                args.species,
                args.tools,
                decompress=not args.no_decompress,
                threads=args.threads,
                with_annotation=args.with_annotation,
                annotation_format=args.annotation_format,
                sjdb_overhang=args.sjdb_overhang,
                genome_sa_index_n_bases=args.genome_sa_index_n_bases
            )
        
        elif args.command == 'list-genomes':
            downloader = ReferenceGenomeDownloader()
            genomes = downloader.list_available_genomes(species=args.species)
            
            print("\n可用的参考基因组:")
            print("=" * 80)
            for source, species_dict in genomes.items():
                print(f"\n数据源: {source}")
                for species, info in species_dict.items():
                    print(f"  物种: {species}")
                    print(f"    名称: {info['name']}")
                    print(f"    描述: {info['description']}")
                    print(f"    URL: {info['url']}")
        
        elif args.command == 'list-tools':
            if args.genome:
                genome_fasta = Path(args.genome)
                if not genome_fasta.exists():
                    logger.error(f"参考基因组文件不存在: {genome_fasta}")
                    sys.exit(1)
                builder = IndexBuilder(genome_fasta)
                available = builder.list_available_tools()
            else:
                builder = IndexBuilder(Path('/dev/null'))  # 临时对象
                available = []
                for tool in builder.SUPPORTED_TOOLS.keys():
                    if builder.check_tool_available(tool):
                        available.append(tool)
            
            print("\n可用的索引工具:")
            print("=" * 80)
            if available:
                for tool in available:
                    info = builder.SUPPORTED_TOOLS[tool]
                    print(f"\n{tool}:")
                    print(f"  描述: {info['description']}")
                    print(f"  命令: {info['command']}")
            else:
                print("\n没有检测到可用的索引工具")
                print("\n支持的索引工具:")
                for tool, info in builder.SUPPORTED_TOOLS.items():
                    print(f"  - {tool}: {info['description']}")
    
    except KeyboardInterrupt:
        logger.info("\n用户中断操作")
        sys.exit(1)
    except Exception as e:
        logger.error(f"错误: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
