#!/usr/bin/env bash
# 参考基因组构建智能体 - 批量测试脚本
#
# 用途：
#   - 批量运行一组典型测试用例（覆盖不同物种和索引工具）
#   - 每个用例单独输出日志
#   - 生成一个汇总结果，方便快速定位失败用例
#
# 使用示例：
#   cd /Users/warm/华大智造/agent/z7zb_tools/reference_genome_builder
#   nohup bash scripts/test_all_tools.sh > test_all_tools.out 2>&1 &

set -u  # 未定义变量视为错误

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT_DIR}"

echo "当前工作目录: ${ROOT_DIR}"

# ===== 配置 LLM 环境变量（如已在环境中配置，可注释掉或修改） =====
# export LLM_API_KEY="your-api-key"
# export LLM_BASE_URL="https://dashscope.aliyuncs.com/compatible-mode/v1"
# export LLM_MODEL="qwen-plus"

# ===== 测试结果根目录（日志 + 统一 work_dir） =====
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
TEST_ROOT="${ROOT_DIR}/test_results_${TIMESTAMP}"
mkdir -p "${TEST_ROOT}"

# 所有测试用例共享一个 work_dir，这样同一物种/版本的参考基因组和索引只下载/构建一次
SHARED_WORK_DIR="${TEST_ROOT}/work"
mkdir -p "${SHARED_WORK_DIR}"

SUMMARY_FILE="${TEST_ROOT}/summary.txt"
echo "测试结果目录: ${TEST_ROOT}" | tee -a "${SUMMARY_FILE}"
echo "开始时间: $(date)" | tee -a "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# ===== 定义测试用例 =====
# 格式：请求|工作子目录名|描述
TESTS=(
  "我需要为WGS分析准备人类参考基因组并构建BWA索引|bwa_human|BWA Human WGS"
  "下载小鼠参考基因组并构建Bowtie2索引|bowtie2_mouse|Bowtie2 Mouse"
  "我需要为RNA-seq分析准备大肠杆菌参考基因组和STAR索引|star_ecoli|STAR E.coli RNA-seq"
  "下载人类参考基因组并构建HISAT2索引用于转录组分析|hisat2_human|HISAT2 Human RNA-seq"
  "为小鼠构建minimap2索引用于PacBio数据比对|minimap2_mouse|minimap2 Mouse Long-read"
)

echo "测试用例总数: ${#TESTS[@]}" | tee -a "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

PASS_COUNT=0
FAIL_COUNT=0

# ===== 逐个运行测试用例 =====
for entry in "${TESTS[@]}"; do
  IFS='|' read -r REQUEST SUBDIR DESC <<< "${entry}"

  echo "==========================================" | tee -a "${SUMMARY_FILE}"
  echo "测试描述: ${DESC}" | tee -a "${SUMMARY_FILE}"
  echo "请求: ${REQUEST}" | tee -a "${SUMMARY_FILE}"
  echo "日志子目录: ${SUBDIR}" | tee -a "${SUMMARY_FILE}"
  echo "------------------------------------------" | tee -a "${SUMMARY_FILE}"

  LOG_DIR="${TEST_ROOT}/${SUBDIR}"
  mkdir -p "${LOG_DIR}"

  LOG_FILE="${LOG_DIR}/run.log"

  # 运行智能体模式（自然语言）
  echo "运行命令:" | tee -a "${SUMMARY_FILE}"
  echo "  python agent_main.py \"${REQUEST}\" --work-dir \"${SHARED_WORK_DIR}\"" | tee -a "${SUMMARY_FILE}"

  # 使用子 shell，避免 set -e 影响整个脚本；即使失败也继续后续测试
  (
    cd "${ROOT_DIR}"
    python agent_main.py "${REQUEST}" --work-dir "${SHARED_WORK_DIR}"
  ) > "${LOG_FILE}" 2>&1
  EXIT_CODE=$?

  if [ ${EXIT_CODE} -eq 0 ]; then
    echo "结果: ✓ 通过" | tee -a "${SUMMARY_FILE}"
    PASS_COUNT=$((PASS_COUNT + 1))
  else
    echo "结果: ✗ 失败 (exit code=${EXIT_CODE})" | tee -a "${SUMMARY_FILE}"
    echo "  日志文件: ${LOG_FILE}" | tee -a "${SUMMARY_FILE}"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi

  echo "" >> "${SUMMARY_FILE}"
  # 稍作间隔，避免日志太密集
  sleep 2
done

echo "==========================================" | tee -a "${SUMMARY_FILE}"
echo "测试完成时间: $(date)" | tee -a "${SUMMARY_FILE}"
echo "通过用例数: ${PASS_COUNT}" | tee -a "${SUMMARY_FILE}"
echo "失败用例数: ${FAIL_COUNT}" | tee -a "${SUMMARY_FILE}"
echo "汇总文件: ${SUMMARY_FILE}"

if [ ${FAIL_COUNT} -ne 0 ]; then
  echo ""
  echo "有失败的测试用例，请查看对应子目录下的 run.log 获取详情。"
fi

