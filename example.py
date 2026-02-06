#!/usr/bin/env python
"""
参考基因组构建智能体使用示例

演示如何使用 Python API 进行参考基因组下载和索引构建
"""

import logging
from pathlib import Path
from download import ReferenceGenomeDownloader
from index_builder import IndexBuilder
from main import ReferenceGenomeAgent

# 配置logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def example_download_genome():
    """示例：下载参考基因组"""
    print("\n" + "=" * 80)
    print("示例 1: 下载参考基因组")
    print("=" * 80)
    
    # 创建下载器
    downloader = ReferenceGenomeDownloader(output_dir="genomes")
    
    # 列出可用的参考基因组
    print("\n可用的参考基因组:")
    genomes = downloader.list_available_genomes(species='human')
    for source, species_dict in genomes.items():
        for species, info in species_dict.items():
            print(f"  {source} - {species}: {info['name']} - {info['description']}")
    
    # 从 NCBI RefSeq 下载人类参考基因组
    print("\n从 NCBI RefSeq 下载人类参考基因组...")
    compressed_file = downloader.download_from_source('ncbi_refseq', 'human')
    print(f"下载完成: {compressed_file}")
    
    # 解压缩
    print("\n解压缩文件...")
    genome_fasta = downloader.decompress_file(compressed_file)
    print(f"解压缩完成: {genome_fasta}")
    
    return genome_fasta


def example_build_index(genome_fasta: Path):
    """示例：构建参考基因组索引"""
    print("\n" + "=" * 80)
    print("示例 2: 构建参考基因组索引")
    print("=" * 80)
    
    # 创建索引构建器
    builder = IndexBuilder(genome_fasta, output_dir="indices")
    
    # 列出可用的工具
    print("\n检查可用的索引工具...")
    available_tools = builder.list_available_tools()
    print(f"可用的工具: {available_tools}")
    
    if not available_tools:
        print("没有可用的索引工具，请先安装相应的工具")
        return
    
    # 构建 BWA 索引（如果可用）
    if 'bwa' in available_tools:
        print("\n构建 BWA 索引...")
        try:
            index_path = builder.build_bwa_index(threads=4)
            print(f"BWA索引构建完成: {index_path}")
        except Exception as e:
            print(f"BWA索引构建失败: {e}")
    
    # 构建 Bowtie2 索引（如果可用）
    if 'bowtie2' in available_tools:
        print("\n构建 Bowtie2 索引...")
        try:
            index_path = builder.build_bowtie2_index(threads=4)
            print(f"Bowtie2索引构建完成: {index_path}")
        except Exception as e:
            print(f"Bowtie2索引构建失败: {e}")


def example_full_pipeline():
    """示例：完整流程"""
    print("\n" + "=" * 80)
    print("示例 3: 完整流程（下载 + 构建索引）")
    print("=" * 80)
    
    # 创建智能体
    agent = ReferenceGenomeAgent(work_dir="example_output")
    
    # 运行完整流程
    try:
        results = agent.run_pipeline(
            source='ensembl',
            species='human',
            tools=['bwa', 'bowtie2'],  # 只构建可用的工具
            decompress=True,
            threads=4
        )
        
        print("\n完整流程结果:")
        print(f"参考基因组: {results['genome_fasta']}")
        print(f"索引结果: {results['index_results']}")
        
    except Exception as e:
        print(f"流程执行失败: {e}")


def example_custom_url():
    """示例：从自定义URL下载"""
    print("\n" + "=" * 80)
    print("示例 4: 从自定义URL下载")
    print("=" * 80)
    
    downloader = ReferenceGenomeDownloader(output_dir="genomes")
    
    # 从自定义URL下载
    custom_url = "https://example.com/reference_genome.fa.gz"
    print(f"从自定义URL下载: {custom_url}")
    
    try:
        # 注意：这是一个示例URL，实际使用时需要替换为真实的URL
        # genome_file = downloader.download_from_url(custom_url, output_filename="custom_genome.fa.gz")
        # print(f"下载完成: {genome_file}")
        print("（示例：实际使用时取消注释上面的代码）")
    except Exception as e:
        print(f"下载失败: {e}")


if __name__ == "__main__":
    print("参考基因组构建智能体 - 使用示例")
    print("=" * 80)
    
    # 运行示例
    try:
        # 示例1: 下载参考基因组
        # genome_fasta = example_download_genome()
        
        # 示例2: 构建索引（需要先有参考基因组文件）
        # if 'genome_fasta' in locals():
        #     example_build_index(genome_fasta)
        
        # 示例3: 完整流程
        # example_full_pipeline()
        
        # 示例4: 自定义URL下载
        example_custom_url()
        
        print("\n" + "=" * 80)
        print("提示: 取消注释上面的示例代码来运行实际的示例")
        print("=" * 80)
        
    except KeyboardInterrupt:
        print("\n用户中断")
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
