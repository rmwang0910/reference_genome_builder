# LLM 配置指南

本工具**统一使用 LLM** 进行智能决策，不再支持规则模式。

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

这会安装 `openai` 库（支持 OpenAI 兼容的 API）。

### 2. 配置环境变量

**必需的环境变量**：

```bash
export LLM_API_KEY='your-api-key-here'
```

**可选的环境变量**：

```bash
# API 基础 URL（默认：https://api.openai.com/v1）
export LLM_BASE_URL='https://dashscope.aliyuncs.com/compatible-mode/v1'

# 模型名称（默认：gpt-3.5-turbo）
export LLM_MODEL='qwen-plus'
```

### 3. 运行

```bash
python agent_main.py "我需要为RNA-seq分析准备大肠杆菌参考基因组和STAR索引"
```

## 支持的 LLM 服务

### OpenAI API

```bash
export LLM_API_KEY='sk-...'
export LLM_BASE_URL='https://api.openai.com/v1'  # 默认，可省略
export LLM_MODEL='gpt-4'  # 或 gpt-3.5-turbo
```

### 阿里云 DashScope

```bash
export LLM_API_KEY='sk-...'
export LLM_BASE_URL='https://dashscope.aliyuncs.com/compatible-mode/v1'
export LLM_MODEL='qwen-plus'  # 或 qwen-max, qwen-turbo 等
```

### 其他 OpenAI 兼容服务

任何支持 OpenAI API 格式的服务都可以使用，只需设置相应的 `LLM_BASE_URL` 和 `LLM_MODEL`。

## 验证配置

如果配置正确，运行时会看到：

```
2026-02-06 17:02:57,135 - agent_core - INFO - LLM已初始化，启用智能决策
2026-02-06 17:02:57,135 - agent_core - INFO - LLM模式: 启用（统一使用LLM，规则模式已移除）
```

如果配置错误，会看到明确的错误提示：

```
ValueError: LLM 初始化失败: ...
请确保已设置以下环境变量：
  - LLM_API_KEY: API 密钥
  - LLM_BASE_URL: API 基础 URL（可选，默认 OpenAI）
  - LLM_MODEL: 模型名称（可选，默认 gpt-3.5-turbo）
```

## 故障排除

### 问题：ImportError: 需要安装 openai 库

**解决方案**：
```bash
pip install openai
```

### 问题：未设置 LLM_API_KEY

**解决方案**：
```bash
export LLM_API_KEY='your-api-key'
```

### 问题：API 调用失败

**可能原因**：
1. API 密钥无效
2. 网络连接问题
3. API 服务不可用

**解决方案**：
- 检查 API 密钥是否正确
- 检查网络连接
- 验证 API 服务状态
- 检查 `LLM_BASE_URL` 是否正确

### 问题：JSON 解析失败

**可能原因**：LLM 返回的响应格式不正确

**解决方案**：
- 检查 LLM 服务是否正常
- 尝试使用不同的模型
- 查看日志中的详细错误信息

## 注意事项

1. **LLM 是必需的**：不再支持规则模式，必须配置 LLM 才能运行
2. **API 费用**：使用 LLM 会产生 API 调用费用，请注意成本
3. **网络要求**：需要能够访问 LLM API 服务的网络连接
4. **响应时间**：LLM 调用需要一定时间，请耐心等待
