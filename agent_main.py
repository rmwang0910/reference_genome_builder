#!/usr/bin/env python
"""
参考基因组构建智能体 - 智能体模式主程序

支持自然语言交互和智能决策
"""

import os
import sys
import logging
import argparse
from pathlib import Path

# 配置logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from agent_core import ReferenceGenomeAgent


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='参考基因组构建智能体（智能体模式）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 自然语言请求
  python agent_main.py "我需要为RNA-seq分析准备人类参考基因组和STAR索引"

  # 指定工作目录
  python agent_main.py "下载人类参考基因组并构建BWA索引" --work-dir /path/to/workdir

  # 禁用LLM（使用基于规则的决策）
  python agent_main.py "下载小鼠参考基因组" --no-llm

  # 获取推荐方案
  python agent_main.py --recommend "WGS全基因组测序分析"
        """
    )
    
    parser.add_argument('request', nargs='?', help='用户请求（自然语言）')
    parser.add_argument('--work-dir', default='.', help='工作目录（默认: 当前目录）')
    parser.add_argument('--no-llm', action='store_true', help='禁用LLM，使用基于规则的决策')
    parser.add_argument('--recommend', help='获取使用场景的推荐方案')
    
    args = parser.parse_args()
    
    # 创建智能体
    agent = ReferenceGenomeAgent(work_dir=args.work_dir, use_llm=not args.no_llm)
    
    try:
        if args.recommend:
            # 推荐模式
            logger.info("=" * 80)
            logger.info("获取推荐方案")
            logger.info("=" * 80)
            
            recommendation = agent.recommend_solution(args.recommend)
            
            print("\n推荐方案:")
            print(f"  数据源: {recommendation['source']}")
            print(f"  理由: {recommendation['source_reason']}")
            print(f"  索引工具: {', '.join(recommendation['tools'])}")
            print(f"  理由: {recommendation['tools_reason']}")
            print(f"  参数: {recommendation['parameters']}")
        
        elif args.request:
            # 执行模式
            results = agent.run(args.request)
            
            print("\n执行结果:")
            for task_id, result in results['tasks'].items():
                status = result.get('status', 'unknown')
                if status == 'completed':
                    print(f"  ✓ {task_id}: {result.get('result', '完成')}")
                else:
                    print(f"  ✗ {task_id}: {result.get('error', '失败')}")
        
        else:
            parser.print_help()
            print("\n提示: 使用自然语言描述您的需求，例如:")
            print('  python agent_main.py "我需要为RNA-seq分析准备人类参考基因组和STAR索引"')
    
    except KeyboardInterrupt:
        logger.info("\n用户中断操作")
        sys.exit(1)
    except Exception as e:
        logger.error(f"错误: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
