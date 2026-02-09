#!/usr/bin/env python
"""
LLM客户端模块

支持 OpenAI 兼容的 API，包括：
- OpenAI API
- 阿里云 DashScope（兼容模式）
- 其他 OpenAI 兼容的 API 服务
"""

import os
import logging
from typing import Optional, List, Dict, Any

logger = logging.getLogger(__name__)

# 尝试导入 openai 库
try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False
    logger.warning("未安装 openai 库，LLM 功能将不可用。请运行: pip install openai")


class LLMClient:
    """LLM客户端：支持 OpenAI 兼容的 API"""
    
    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 2000
    ):
        """
        初始化 LLM 客户端
        
        Args:
            api_key: API 密钥（从环境变量 LLM_API_KEY 读取，如果未提供）
            base_url: API 基础 URL（从环境变量 LLM_BASE_URL 读取，如果未提供）
            model: 模型名称（从环境变量 LLM_MODEL 读取，如果未提供）
            temperature: 温度参数（默认 0.7）
            max_tokens: 最大 token 数（默认 2000）
        """
        if not HAS_OPENAI:
            raise ImportError("需要安装 openai 库: pip install openai")
        
        # 从环境变量读取配置
        self.api_key = api_key or os.getenv("LLM_API_KEY")
        self.base_url = base_url or os.getenv("LLM_BASE_URL", "https://api.openai.com/v1")
        self.model = model or os.getenv("LLM_MODEL", "gpt-3.5-turbo")
        self.temperature = temperature
        self.max_tokens = max_tokens
        
        if not self.api_key:
            raise ValueError(
                "未设置 LLM_API_KEY。请设置环境变量或传入 api_key 参数。\n"
                "例如: export LLM_API_KEY='your-api-key'"
            )
        
        # 初始化 OpenAI 客户端
        self.client = OpenAI(
            api_key=self.api_key,
            base_url=self.base_url
        )
        
        logger.info(f"LLM 客户端已初始化: model={self.model}, base_url={self.base_url}")
    
    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: Optional[int] = None,
        temperature: Optional[float] = None,
        stream: bool = False
    ) -> str:
        """
        生成文本
        
        Args:
            prompt: 用户提示
            system_prompt: 系统提示（可选）
            max_tokens: 最大 token 数（覆盖初始化时的设置）
            temperature: 温度参数（覆盖初始化时的设置）
            stream: 是否流式输出
            
        Returns:
            生成的文本
        """
        messages = []
        
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        
        messages.append({"role": "user", "content": prompt})
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                max_tokens=max_tokens or self.max_tokens,
                temperature=temperature if temperature is not None else self.temperature,
                stream=stream
            )
            
            if stream:
                # 流式输出处理
                full_response = ""
                for chunk in response:
                    if chunk.choices[0].delta.content:
                        full_response += chunk.choices[0].delta.content
                return full_response
            else:
                return response.choices[0].message.content.strip()
                
        except Exception as e:
            logger.error(f"LLM 调用失败: {e}")
            raise
    
    def generate_json(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: Optional[int] = None,
        temperature: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        生成 JSON 格式的响应
        
        Args:
            prompt: 用户提示
            system_prompt: 系统提示（可选）
            max_tokens: 最大 token 数
            temperature: 温度参数
            
        Returns:
            解析后的 JSON 字典
        """
        import json
        import re
        
        # 在提示中强调需要 JSON 格式
        json_prompt = f"""{prompt}

请确保你的响应是有效的 JSON 格式。只返回 JSON，不要包含其他文本或代码块标记。"""
        
        response = self.generate(
            prompt=json_prompt,
            system_prompt=system_prompt,
            max_tokens=max_tokens,
            temperature=temperature
        )
        
        # 尝试提取 JSON（可能被代码块包裹）
        json_match = re.search(r'\{.*\}', response, re.DOTALL)
        if json_match:
            try:
                return json.loads(json_match.group())
            except json.JSONDecodeError:
                pass
        
        # 如果提取失败，尝试直接解析
        try:
            return json.loads(response)
        except json.JSONDecodeError as e:
            logger.error(f"JSON 解析失败: {e}\n响应内容: {response[:500]}")
            raise ValueError(f"无法解析 LLM 响应为 JSON: {e}")


def create_llm_client() -> Optional[LLMClient]:
    """
    创建 LLM 客户端（如果配置可用）
    
    Returns:
        LLMClient 实例，如果配置不可用则返回 None
    """
    try:
        return LLMClient()
    except (ImportError, ValueError) as e:
        logger.warning(f"无法创建 LLM 客户端: {e}")
        return None
