#!/usr/bin/env python
"""
参考基因组下载模块

支持从多个数据源下载参考基因组：
- NCBI RefSeq
- Ensembl
- UCSC
- 自定义URL
"""

import os
import sys
import logging
import subprocess
import hashlib
from pathlib import Path
from typing import Optional, Dict, List
from urllib.parse import urlparse
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)


class ReferenceGenomeDownloader:
    """参考基因组下载器
    
    职责：
    - 下载参考基因组 FASTA
    - 下载基因注释文件（GTF/GFF）
    - 处理重复下载与校验
    """
    
    # 常用参考基因组数据源URL模板（FASTA）
    DATA_SOURCES = {
        'ncbi_refseq': {
            'human': {
                'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz',
                'name': 'GRCh38.p14',
                'description': 'Human reference genome (GRCh38.p14) from NCBI RefSeq'
            },
            'mouse': {
                'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz',
                'name': 'GRCm39',
                'description': 'Mouse reference genome (GRCm39) from NCBI RefSeq'
            }
        },
        'ensembl': {
            'human': {
                'url': 'https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
                'name': 'GRCh38',
                'description': 'Human reference genome (GRCh38) from Ensembl'
            },
            'mouse': {
                'url': 'https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz',
                'name': 'GRCm39',
                'description': 'Mouse reference genome (GRCm39) from Ensembl'
            }
        },
        'ucsc': {
            'human': {
                'url': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
                'name': 'hg38',
                'description': 'Human reference genome (hg38) from UCSC'
            },
            'mouse': {
                'url': 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz',
                'name': 'mm39',
                'description': 'Mouse reference genome (mm39) from UCSC'
            }
        }
    }
    
    # 常用注释文件数据源URL模板（GTF/GFF）
    # 这里选取了当前较稳定的路径，实际使用中可以根据需要扩展/调整 release 版本
    ANNOTATION_SOURCES = {
        'ensembl': {
            'human': {
                'gtf': {
                    'url': 'https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz',
                    'name': 'GRCh38.110',
                    'description': 'Human GTF annotation (GRCh38, Ensembl 110)'
                }
            },
            'mouse': {
                'gtf': {
                    'url': 'https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz',
                    'name': 'GRCm39.110',
                    'description': 'Mouse GTF annotation (GRCm39, Ensembl 110)'
                }
            }
        },
        'ncbi_refseq': {
            'human': {
                'gff3': {
                    # 与 GRCh38.p14 对应的 GFF3
                    'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz',
                    'name': 'GRCh38.p14',
                    'description': 'Human GFF3 annotation (GRCh38.p14, NCBI RefSeq)'
                }
            },
            'mouse': {
                'gff3': {
                    'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz',
                    'name': 'GRCm39',
                    'description': 'Mouse GFF3 annotation (GRCm39, NCBI RefSeq)'
                }
            }
        },
        'ucsc': {
            'human': {
                'gtf': {
                    # UCSC refGene GTF
                    'url': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz',
                    'name': 'hg38.refGene',
                    'description': 'Human GTF annotation (hg38, UCSC refGene)'
                }
            },
            'mouse': {
                'gtf': {
                    'url': 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.refGene.gtf.gz',
                    'name': 'mm39.refGene',
                    'description': 'Mouse GTF annotation (mm39, UCSC refGene)'
                }
            }
        }
    }
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        初始化下载器
        
        Args:
            output_dir: 输出目录，默认为当前目录下的 genomes 文件夹
        """
        if output_dir is None:
            self.output_dir = Path.cwd() / "genomes"
        else:
            self.output_dir = Path(output_dir)
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"参考基因组下载目录: {self.output_dir}")
    
    def list_available_genomes(self, species: Optional[str] = None) -> Dict[str, List[Dict]]:
        """
        列出可用的参考基因组
        
        Args:
            species: 物种名称（如 'human', 'mouse'），如果为None则列出所有
            
        Returns:
            可用的参考基因组字典
        """
        if species:
            available = {}
            for source, genomes in self.DATA_SOURCES.items():
                if species.lower() in genomes:
                    available[source] = {species: genomes[species.lower()]}
            return available
        return self.DATA_SOURCES
    
    def list_available_annotations(self, species: Optional[str] = None) -> Dict[str, Dict]:
        """
        列出可用的注释文件（GTF/GFF）
        
        Args:
            species: 物种名称（如 'human', 'mouse'），如果为None则列出所有
            
        Returns:
            可用的注释文件字典
        """
        if species:
            available = {}
            for source, species_dict in self.ANNOTATION_SOURCES.items():
                if species.lower() in species_dict:
                    available[source] = {species: species_dict[species.lower()]}
            return available
        return self.ANNOTATION_SOURCES
    
    def download_from_url(self, url: str, output_filename: Optional[str] = None, 
                         verify_checksum: bool = False, checksum: Optional[str] = None) -> Path:
        """
        从URL下载参考基因组
        
        Args:
            url: 下载URL
            output_filename: 输出文件名，如果为None则从URL推断
            verify_checksum: 是否验证校验和
            checksum: 预期的MD5校验和（如果verify_checksum为True）
            
        Returns:
            下载文件的路径
        """
        if output_filename is None:
            # 从URL提取文件名
            parsed_url = urlparse(url)
            output_filename = os.path.basename(parsed_url.path)
            if not output_filename:
                output_filename = "reference_genome.fa.gz"
        
        output_path = self.output_dir / output_filename
        
        # 检查文件是否已存在
        if output_path.exists():
            logger.info(f"文件已存在: {output_path}")
            if verify_checksum and checksum:
                if self._verify_file_checksum(output_path, checksum):
                    logger.info("文件校验和验证通过")
                    return output_path
                else:
                    logger.warning("文件校验和不匹配，将重新下载")
                    output_path.unlink()
            else:
                return output_path
        
        logger.info(f"开始下载: {url}")
        logger.info(f"保存到: {output_path}")
        
        try:
            # 使用requests下载，支持断点续传
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(output_path, 'wb') as f:
                if total_size == 0:
                    f.write(response.content)
                else:
                    with tqdm(total=total_size, unit='B', unit_scale=True, desc="下载进度") as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
            
            logger.info(f"下载完成: {output_path}")
            
            # 验证校验和
            if verify_checksum and checksum:
                if not self._verify_file_checksum(output_path, checksum):
                    logger.error("文件校验和验证失败")
                    output_path.unlink()
                    raise ValueError("文件校验和验证失败")
            
            return output_path
            
        except requests.exceptions.RequestException as e:
            logger.error(f"下载失败: {e}")
            if output_path.exists():
                output_path.unlink()
            raise
    
    def download_from_source(self, source: str, species: str, 
                            output_filename: Optional[str] = None) -> Path:
        """
        从预定义的数据源下载参考基因组
        
        Args:
            source: 数据源名称（'ncbi_refseq', 'ensembl', 'ucsc'）
            species: 物种名称（'human', 'mouse'）
            output_filename: 输出文件名，如果为None则使用默认名称
            
        Returns:
            下载文件的路径
        """
        source = source.lower()
        species = species.lower()
        
        if source not in self.DATA_SOURCES:
            raise ValueError(f"不支持的数据源: {source}。支持的数据源: {list(self.DATA_SOURCES.keys())}")
        
        if species not in self.DATA_SOURCES[source]:
            raise ValueError(f"数据源 {source} 不支持物种 {species}。支持的物种: {list(self.DATA_SOURCES[source].keys())}")
        
        genome_info = self.DATA_SOURCES[source][species]
        url = genome_info['url']
        
        if output_filename is None:
            output_filename = f"{species}_{genome_info['name']}.fa.gz"
        
        logger.info(f"从 {source} 下载 {species} 参考基因组 ({genome_info['name']})")
        logger.info(f"描述: {genome_info['description']}")
        
        return self.download_from_url(url, output_filename)
    
    def download_annotation_from_source(
        self,
        source: str,
        species: str,
        fmt: str = "gtf",
        output_filename: Optional[str] = None
    ) -> Path:
        """
        从预定义的数据源下载基因注释文件（GTF/GFF）
        
        Args:
            source: 数据源名称（'ncbi_refseq', 'ensembl', 'ucsc'）
            species: 物种名称（'human', 'mouse'）
            fmt: 注释格式（'gtf' 或 'gff3'）
            output_filename: 输出文件名，如果为None则使用默认名称
            
        Returns:
            下载文件的路径
        """
        source = source.lower()
        species = species.lower()
        fmt = fmt.lower()
        
        if source not in self.ANNOTATION_SOURCES:
            raise ValueError(
                f"不支持的注释数据源: {source}。支持的数据源: {list(self.ANNOTATION_SOURCES.keys())}"
            )
        
        if species not in self.ANNOTATION_SOURCES[source]:
            raise ValueError(
                f"注释数据源 {source} 不支持物种 {species}。支持的物种: {list(self.ANNOTATION_SOURCES[source].keys())}"
            )
        
        species_info = self.ANNOTATION_SOURCES[source][species]
        if fmt not in species_info:
            raise ValueError(
                f"数据源 {source} 物种 {species} 不支持格式 {fmt}。支持的格式: {list(species_info.keys())}"
            )
        
        ann_info = species_info[fmt]
        url = ann_info['url']
        
        if output_filename is None:
            # e.g. human_GRCh38.110.gtf.gz
            ext = "gtf" if fmt == "gtf" else "gff3"
            output_filename = f"{species}_{ann_info['name']}.{ext}.gz"
        
        logger.info(f"从 {source} 下载 {species} 注释文件 ({ann_info['name']}, {fmt})")
        logger.info(f"描述: {ann_info['description']}")
        
        return self.download_from_url(url, output_filename)
    
    def download_using_wget(self, url: str, output_filename: Optional[str] = None) -> Path:
        """
        使用wget下载（如果系统有wget）
        
        Args:
            url: 下载URL
            output_filename: 输出文件名
            
        Returns:
            下载文件的路径
        """
        if output_filename is None:
            parsed_url = urlparse(url)
            output_filename = os.path.basename(parsed_url.path)
            if not output_filename:
                output_filename = "reference_genome.fa.gz"
        
        output_path = self.output_dir / output_filename
        
        # 检查wget是否可用
        try:
            subprocess.run(['wget', '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("wget 不可用，使用 requests 下载")
            return self.download_from_url(url, output_filename)
        
        logger.info(f"使用 wget 下载: {url}")
        
        try:
            cmd = ['wget', '-O', str(output_path), url]
            subprocess.run(cmd, check=True)
            logger.info(f"下载完成: {output_path}")
            return output_path
        except subprocess.CalledProcessError as e:
            logger.error(f"wget 下载失败: {e}")
            raise
    
    def _verify_file_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """
        验证文件的MD5校验和
        
        Args:
            file_path: 文件路径
            expected_checksum: 预期的MD5校验和
            
        Returns:
            校验和是否匹配
        """
        logger.info("计算文件MD5校验和...")
        md5_hash = hashlib.md5()
        
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        
        actual_checksum = md5_hash.hexdigest()
        matches = actual_checksum.lower() == expected_checksum.lower()
        
        if matches:
            logger.info(f"校验和匹配: {actual_checksum}")
        else:
            logger.warning(f"校验和不匹配。预期: {expected_checksum}, 实际: {actual_checksum}")
        
        return matches
    
    def decompress_file(self, compressed_path: Path, output_path: Optional[Path] = None) -> Path:
        """
        解压缩文件（支持.gz格式）
        
        Args:
            compressed_path: 压缩文件路径
            output_path: 输出文件路径，如果为None则自动推断
            
        Returns:
            解压后的文件路径
        """
        if output_path is None:
            if compressed_path.suffix == '.gz':
                output_path = compressed_path.with_suffix('')
            else:
                output_path = compressed_path.with_suffix('.fa')
        
        logger.info(f"解压缩: {compressed_path} -> {output_path}")
        
        import gzip
        
        try:
            with gzip.open(compressed_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    f_out.write(f_in.read())
            
            logger.info(f"解压缩完成: {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"解压缩失败: {e}")
            raise
