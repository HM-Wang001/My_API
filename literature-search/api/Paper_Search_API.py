#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因变异文献搜索 API - Vercel Serverless Function
"""
import sys
import os
import json
import re
import html
import time
import requests
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple, Any
from datetime import datetime


# ============ 你的原始代码（完全未修改） ============
@dataclass
class Variant:
    """变异类 - 搜索策略：基因名 AND (hgvs_c OR hgvs_p OR rsid OR chr+坐标)"""

    chrom: Optional[str] = None
    pos: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    hgvs_g: Optional[str] = None
    gene: Optional[str] = None
    rsid: Optional[str] = None

    def __post_init__(self):
        if self.chrom:
            self.chrom = self.chrom.replace('chr', '').replace('Chr', '')

        # 自动解析输入的变异格式
        self._auto_parse()

    def _auto_parse(self):
        """自动解析各种变异格式"""
        # 如果已经有hgvs_c但没解析，尝试提取
        if self.hgvs_c and not self.hgvs_p:
            # 尝试从hgvs_c中提取p.
            p_match = re.search(r'\(p\.([^)]+)\)', self.hgvs_c)
            if p_match:
                self.hgvs_p = f"p.{p_match.group(1)}"
                # 清理hgvs_c，去掉括号部分
                self.hgvs_c = re.sub(r'\(p\.[^)]+\)', '', self.hgvs_c).strip()

    @classmethod
    def from_hgvs(cls, hgvs: str, gene: str = None, hgvs_p: str = None):
        """从HGVS字符串创建变异
        Args:
            hgvs: HGVS字符串，如 "NM_000484.4:c.2059A>C" 或 "NM_000484.4:c.2059A>C(p.K687Q)"
            gene: 基因名
            hgvs_p: 可选的蛋白水平HGVS，如 "p.K687Q" 或 "p.Lys687Gln"
        """
        variant = cls(hgvs_c=hgvs, gene=gene)

        # 尝试从输入字符串中解析c.和p.
        if 'c.' in hgvs:
            # 提取c.部分（在括号前）
            c_match = re.search(r'c\.([^\s;()]+)', hgvs)
            if c_match:
                variant.hgvs_c = f"c.{c_match.group(1)}"

            # 尝试提取括号内的p.
            p_in_paren = re.search(r'\(p\.([^)]+)\)', hgvs)
            if p_in_paren:
                variant.hgvs_p = f"p.{p_in_paren.group(1)}"

        # 如果提供了额外的hgvs_p参数，使用它
        if hgvs_p:
            if 'p.' in hgvs_p:
                p_match = re.search(r'p\.([^\s;]+)', hgvs_p)
                if p_match:
                    variant.hgvs_p = f"p.{p_match.group(1)}"
            else:
                variant.hgvs_p = f"p.{hgvs_p}"

        return variant

    @classmethod
    def from_vcf(cls, chrom: str, pos: int, ref: str, alt: str, gene: str = None):
        return cls(chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene)

    @classmethod
    def from_rsid(cls, rsid: str, gene: str = None):
        if not rsid.startswith('rs'):
            rsid = 'rs' + rsid
        return cls(rsid=rsid, gene=gene)

    @classmethod
    def from_string(cls, gene: str, variant_str: str):
        """从字符串智能解析变异（支持多种格式）"""
        variant_str = variant_str.strip()

        # 1. 检测RS号
        if re.match(r'rs\d+', variant_str, re.IGNORECASE):
            return cls.from_rsid(variant_str)

        # 2. 检测基因组坐标: chr7:140453136A>T
        genomic_match = re.match(r'(?:chr)?([\dXYM]+):(\d+)([A-Z]+)>([A-Z]+)', variant_str, re.IGNORECASE)
        if genomic_match:
            chrom, pos, ref, alt = genomic_match.groups()
            var = cls.from_vcf(chrom, int(pos), ref, alt)
            var.gene = gene
            return var

        # 3. 检测c.变异
        if 'c.' in variant_str:
            return cls.from_hgvs(variant_str, gene=gene)

        # 4. 检测p.变异
        if 'p.' in variant_str:
            var = cls(hgvs_p=variant_str, gene=gene)
            return var

        # 5. 检测简单蛋白质变异: V600E, L858R
        simple_p_match = re.match(r'^([A-Z])(\d+)([A-Z])$', variant_str)
        if simple_p_match:
            ref_aa, pos, alt_aa = simple_p_match.groups()
            var = cls(hgvs_p=f"p.{ref_aa}{pos}{alt_aa}", gene=gene)
            # 也尝试生成三字母形式
            aa_one_to_three = {
                'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
            }
            if ref_aa in aa_one_to_three and alt_aa in aa_one_to_three:
                var.hgvs_p_alt = f"p.{aa_one_to_three[ref_aa]}{pos}{aa_one_to_three[alt_aa]}"
            return var

        # 6. 检测三字母蛋白质变异: Val600Glu
        three_p_match = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', variant_str)
        if three_p_match:
            ref_aa, pos, alt_aa = three_p_match.groups()
            var = cls(hgvs_p=f"p.{ref_aa}{pos}{alt_aa}", gene=gene)
            # 生成单字母形式
            aa_three_to_one = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            if ref_aa in aa_three_to_one and alt_aa in aa_three_to_one:
                var.hgvs_p_alt = f"p.{aa_three_to_one[ref_aa]}{pos}{aa_three_to_one[alt_aa]}"
            return var

        # 7. 如果都匹配不上，当作普通文本
        var = cls(gene=gene)
        var.hgvs_c = variant_str
        return var

    def get_search_terms(self) -> List[str]:
        """
        生成搜索词：基因名 AND (hgvs_c OR hgvs_p OR rsid OR chr+坐标)
        返回用于PubMed搜索的关键词列表
        """
        terms = []

        # 必须有基因名
        if not self.gene:
            return terms

        gene = self.gene

        # 1. HGVS c. 格式 - 只使用核心部分如 "2059A>C"
        if self.hgvs_c:
            c_match = re.search(r'c\.([^\s;]+)', self.hgvs_c)
            if c_match:
                c_core = c_match.group(1)  # "2059A>C"
                # 只使用核心部分，不包含c.
                terms.append(f'{gene} AND ("{c_core}")')

        # 2. HGVS p. 格式 (包括三字母和单字母)
        if self.hgvs_p:
            p_match = re.search(r'p\.([^\s;]+)', self.hgvs_p)
            if p_match:
                p_core = p_match.group(1)  # "Lys687Gln" 或 "K687Q"
                # 搜索p.完整形式
                terms.append(f'{gene} AND ("p.{p_core}")')
                # 搜索不带p.的形式
                terms.append(f'{gene} AND ("{p_core}")')

                # 尝试生成简写形式（三字母转单字母）
                short = self._amino_acid_to_single(p_core)
                if short:
                    terms.append(f'{gene} AND ("p.{short}")')
                    terms.append(f'{gene} AND ("{short}")')

        # 如果有备用的p.形式（从其他格式转换的）
        if hasattr(self, 'hgvs_p_alt') and self.hgvs_p_alt:
            p_match = re.search(r'p\.([^\s;]+)', self.hgvs_p_alt)
            if p_match:
                p_core = p_match.group(1)
                terms.append(f'{gene} AND ("p.{p_core}")')
                terms.append(f'{gene} AND ("{p_core}")')

        # 3. rsID
        if self.rsid:
            terms.append(f'{gene} AND ("{self.rsid}")')
            terms.append(f'{self.rsid}')  # 也单独搜索rs号

        # 4. 基因组坐标搜索 - 只保留数字部分
        if self.chrom and self.pos:
            # 只使用位置数字，不包含变异信息
            terms.append(f'{gene} AND {self.pos}')

        # 去重并保持顺序
        seen = set()
        unique_terms = []
        for t in terms:
            t_lower = t.lower()
            if t_lower not in seen and len(t) > 3:
                seen.add(t_lower)
                unique_terms.append(t)

        return unique_terms

    def _amino_acid_to_single(self, p_core: str) -> Optional[str]:
        """三字母氨基酸转单字母"""
        aa_3to1 = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
            'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
            'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
            'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
        }

        match = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3})', p_core)
        if match:
            wt, pos, mut = match.groups()
            wt_1 = aa_3to1.get(wt, '?')
            mut_1 = aa_3to1.get(mut, '?')
            return f"{wt_1}{pos}{mut_1}"

        return None

    def get_match_keywords(self) -> List[str]:
        """生成用于在文本中匹配的关键词"""
        keywords = []
        if self.hgvs_c:
            match = re.search(r'c\..*', self.hgvs_c)
            if match:
                keywords.append(match.group())
            match = re.search(r'c.(\d+\w+>\w+)', self.hgvs_c)
            if match:
                keywords.append(match.group(1))
        if self.hgvs_p:
            match = re.search(r'p\..*', self.hgvs_p)
            if match:
                keywords.append(match.group())
            # 提取核心部分
            p_match = re.search(r'p\.([^\s;]+)', self.hgvs_p)
            if p_match:
                keywords.append(p_match.group(1))
        if self.rsid:
            keywords.append(self.rsid)
            keywords.append(self.rsid.replace('rs', ''))
        if self.gene:
            keywords.append(self.gene)
        return list(set(keywords))

    def find_variant_sentences(self, text: str) -> List[Dict]:
        """从文本中提取包含特定变异的句子（严格匹配具体变异）"""
        if not text:
            return []

        text = html.unescape(text)

        # 获取严格匹配模式
        strict_patterns = self._get_strict_match_patterns()

        # 分句
        sentences = re.split(r'(?<=[.!?])\s+', text)

        matches = []
        seen = set()

        for sentence in sentences:
            sentence = sentence.strip()
            if len(sentence) < 20:
                continue

            # 严格检查：必须匹配到具体的变异标识
            matched_pattern = None
            matched_keywords = []

            for pattern_name, pattern in strict_patterns:
                found = pattern.findall(sentence)
                if found:
                    matched_pattern = pattern_name
                    # 收集匹配到的关键词
                    for f in found:
                        if isinstance(f, tuple):
                            matched_keywords.extend([m for m in f if m])
                        else:
                            matched_keywords.append(str(f))
                    break

            if not matched_pattern:
                continue

            # 去重
            key = sentence.lower()[:100]
            if key in seen:
                continue
            seen.add(key)

            # 去重关键词
            matched_keywords = list(set(matched_keywords))[:5]

            # 高亮
            highlighted = self._highlight_matches(sentence, matched_keywords or [matched_pattern])

            matches.append({
                'sentence': sentence,
                'highlighted': highlighted,
                'match_type': matched_pattern,
                'keywords': matched_keywords
            })

        return matches

    def _get_strict_match_patterns(self) -> List[tuple]:
        """生成严格的变异匹配模式"""
        patterns = []

        # 1. HGVS c.
        if self.hgvs_c:
            c_match = re.search(r'c\.([^\s;]+)', self.hgvs_c)
            if c_match:
                c_core = c_match.group(1)
                patterns.append(('hgvs_c', re.compile(
                    rf'c\s*\.?\s*{re.escape(c_core)}',
                    re.IGNORECASE
                )))
                # 不带c.前缀
                patterns.append(('hgvs_c_no_prefix', re.compile(
                    rf'\b{re.escape(c_core)}\b',
                    re.IGNORECASE
                )))

        # 2. HGVS p.
        if self.hgvs_p:
            p_match = re.search(r'p\.([^\s;]+)', self.hgvs_p)
            if p_match:
                p_core = p_match.group(1)
                patterns.append(('hgvs_p', re.compile(
                    rf'p\s*\.?\s*{re.escape(p_core)}',
                    re.IGNORECASE
                )))
                # 不带p.前缀
                patterns.append(('hgvs_p_no_prefix', re.compile(
                    rf'\b{re.escape(p_core)}\b',
                    re.IGNORECASE
                )))

                # 简写形式
                short = self._amino_acid_to_single(p_core)
                if short:
                    patterns.append(('hgvs_p_short', re.compile(
                        rf'\b{re.escape(short)}\b',
                        re.IGNORECASE
                    )))
                    patterns.append(('hgvs_p_short_p', re.compile(
                        rf'p\s*\.?\s*{re.escape(short)}',
                        re.IGNORECASE
                    )))

        # 如果有备用的p.形式
        if hasattr(self, 'hgvs_p_alt') and self.hgvs_p_alt:
            p_match = re.search(r'p\.([^\s;]+)', self.hgvs_p_alt)
            if p_match:
                p_core = p_match.group(1)
                patterns.append(('hgvs_p_alt', re.compile(
                    rf'p\s*\.?\s*{re.escape(p_core)}',
                    re.IGNORECASE
                )))

        # 3. rsID
        if self.rsid:
            patterns.append(('rsid', re.compile(
                rf'\b{re.escape(self.rsid)}\b',
                re.IGNORECASE
            )))

        # 4. 基因组坐标
        if all([self.chrom, self.pos]):
            patterns.append(('genomic_coord', re.compile(
                rf'(?:chr)?{re.escape(str(self.chrom))}[:.]\s*{re.escape(str(self.pos))}',
                re.IGNORECASE
            )))
            if self.ref and self.alt:
                patterns.append(('genomic_variant', re.compile(
                    rf'(?:chr)?{re.escape(str(self.chrom))}[:.]\s*{re.escape(str(self.pos))}\s*{re.escape(self.ref)}>\s*{re.escape(self.alt)}',
                    re.IGNORECASE
                )))

        return patterns

    def _highlight_matches(self, sentence: str, keywords: List[str]) -> str:
        """高亮句子中的匹配关键词"""
        highlighted = sentence
        for kw in sorted(set(keywords), key=len, reverse=True):
            if len(kw) < 2:
                continue
            pattern = re.compile(re.escape(kw), re.IGNORECASE)
            highlighted = pattern.sub(lambda m: f'**{m.group()}**', highlighted)
        return highlighted


class PubMedTextFetcher:
    """获取 PubMed 文献的摘要或全文"""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.session = requests.Session()

    def fetch_text(self, pmid: str) -> Optional[str]:
        """获取文献文本（优先全文， fallback 到摘要）"""
        # 尝试获取 PMC 全文
        text = self._fetch_pmc_text(pmid)
        if text:
            return text

        # fallback 到摘要
        return self._fetch_abstract(pmid)

    def _fetch_pmc_text(self, pmid: str) -> Optional[str]:
        """尝试从 PMC 获取全文"""
        try:
            # 先查 PMC ID
            url = f"{self.BASE_URL}/esearch.fcgi"
            params = {
                'db': 'pmc',
                'term': f'{pmid}[uid]',
                'retmode': 'json',
                'retmax': 1
            }
            if self.api_key:
                params['api_key'] = self.api_key

            response = self.session.get(url, params=params, timeout=30)
            data = response.json()
            idlist = data.get('esearchresult', {}).get('idlist', [])

            if not idlist:
                return None

            # 获取全文 XML
            pmcid = idlist[0]
            url = f"{self.BASE_URL}/efetch.fcgi"
            params = {
                'db': 'pmc',
                'id': pmcid,
                'retmode': 'xml'
            }
            if self.api_key:
                params['api_key'] = self.api_key

            response = self.session.get(url, params=params, timeout=30)

            # 简单提取文本
            text = re.sub(r'<[^>]+>', ' ', response.text)
            text = re.sub(r'\s+', ' ', text)
            return text if len(text) > 500 else None

        except Exception:
            return None

    def _fetch_abstract(self, pmid: str) -> Optional[str]:
        """获取摘要"""
        try:
            url = f"{self.BASE_URL}/efetch.fcgi"
            params = {
                'db': 'pubmed',
                'id': pmid,
                'retmode': 'xml'
            }
            if self.api_key:
                params['api_key'] = self.api_key

            response = self.session.get(url, params=params, timeout=30)

            # 提取摘要文本
            text = re.sub(r'<[^>]+>', ' ', response.text)
            text = re.sub(r'\s+', ' ', text)
            return text if len(text) > 50 else None

        except Exception:
            return None


class PubMedSearcher:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'VariantLiteratureSearcher/1.0'
        })
        self.text_fetcher = PubMedTextFetcher(api_key)

    def search(self, variant: Variant, max_results: int = 20) -> List[Dict]:
        terms = variant.get_search_terms()
        all_results = []

        print(f"\n🔍 生成搜索词:")
        for i, term in enumerate(terms[:5], 1):
            print(f"  {i}. {term}")

        for term in terms[:5]:
            try:
                print(f"\n  → 搜索: {term}")
                results = self._search_term(variant, term, max_results)
                all_results.extend(results)
                print(f"    找到 {len(results)} 篇文献")
                time.sleep(0.34)
            except Exception as e:
                print(f"    搜索出错: {e}")

        # 去重
        seen = set()
        unique_results = []
        for r in all_results:
            pmid = r.get('pmid')
            if pmid and pmid not in seen:
                seen.add(pmid)
                unique_results.append(r)

        return unique_results

    def _search_term(self, variant: Variant, term: str, max_results: int) -> List[Dict]:
        """搜索并严格过滤"""

        search_params = {
            'db': 'pubmed',
            'term': term,
            'retmode': 'json',
            'retmax': max_results,
            'sort': 'relevance'
        }

        idlist = self._do_esearch(search_params)

        if not idlist:
            return []

        # 获取详情
        papers = self._fetch_details(variant, idlist)

        # 严格过滤：只保留包含具体变异的
        filtered = []
        for paper in papers:
            sentences = paper.get('variant_sentences', [])
            if sentences:
                filtered.append(paper)

        return filtered

    def _do_esearch(self, params: dict) -> List[str]:
        """执行 esearch 并返回 ID 列表"""
        url = f"{self.BASE_URL}/esearch.fcgi"
        if self.api_key:
            params['api_key'] = self.api_key

        try:
            response = self.session.get(url, params=params, timeout=30)
            data = response.json()
            return data.get('esearchresult', {}).get('idlist', [])
        except Exception as e:
            print(f"搜索出错: {e}")
            return []

    def _fetch_details(self, variant: Variant, idlist: List[str]) -> List[Dict]:
        ids = ','.join(idlist)

        summary_params = {
            'db': 'pubmed',
            'id': ids,
            'retmode': 'json'
        }
        if self.api_key:
            summary_params['api_key'] = self.api_key

        summary_url = f"{self.BASE_URL}/esummary.fcgi"
        response = self.session.get(summary_url, params=summary_params, timeout=30)
        data = response.json()

        results = data.get('result', {})
        papers = []

        for pmid in idlist:
            if pmid == 'uids':
                continue

            info = results.get(pmid, {})
            if not info:
                continue

            paper = {
                'pmid': pmid,
                'title': info.get('title', 'N/A'),
                'authors': [a.get('name', '') for a in info.get('authors', [])[:3]],
                'journal': info.get('fulljournalname', info.get('source', 'N/A')),
                'year': info.get('pubdate', 'N/A')[:4] if info.get('pubdate') else 'N/A',
                'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                'source': 'PubMed'
            }

            # 提取变异相关句子
            print(f"  📖 分析 PMID:{pmid}...")
            paper['variant_sentences'] = self._extract_sentences(variant, pmid)

            papers.append(paper)

        return papers

    def _extract_sentences(self, variant: Variant, pmid: str) -> List[Dict]:
        """从文献中提取包含变异的句子"""
        text = self.text_fetcher.fetch_text(pmid)
        if not text:
            return []

        sentences = variant.find_variant_sentences(text)

        for s in sentences:
            s['source'] = 'fulltext' if len(text) > 3000 else 'abstract'

        return sentences


class VariantLiteratureSearch:
    def __init__(self, ncbi_api_key: Optional[str] = None):
        self.pubmed = PubMedSearcher(api_key=ncbi_api_key)

    def search(self, variant: Variant, databases: Optional[List[str]] = None) -> Dict:
        if databases is None:
            databases = ['pubmed', 'clinvar']

        results = {
            'query': {
                'terms': variant.get_search_terms(),
                'timestamp': datetime.now().isoformat()
            },
            'results': {}
        }

        if 'pubmed' in databases:
            print("\n🔍 正在搜索 PubMed...")
            results['results']['pubmed'] = self.pubmed.search(variant)

        return results

    def search_and_format(self, variant: Variant,
                          databases: Optional[List[str]] = None) -> str:
        """格式化输出，严格显示变异匹配"""
        results = self.search(variant, databases)

        lines = []
        lines.append("=" * 80)
        lines.append("🔬 基因变异文献搜索报告")
        lines.append("=" * 80)
        lines.append(f"搜索时间: {results['query']['timestamp']}")
        lines.append(f"基因: {variant.gene}")
        lines.append(f"变异: {variant.hgvs_c or variant.hgvs_p or variant.rsid or 'N/A'}")
        lines.append(f"搜索词: {', '.join(results['query']['terms'][:5])}")
        lines.append(f"匹配规则: 必须包含具体变异标识 (c./p./rsID/坐标)")
        lines.append("=" * 80)

        # PubMed 结果
        if 'pubmed' in results['results']:
            lines.append("\n📚 PubMed 文献搜索结果")
            lines.append("-" * 80)
            papers = results['results']['pubmed']

            if papers:
                lines.append(f"✅ 找到 {len(papers)} 篇包含该变异的文献\n")

                for i, paper in enumerate(papers[:20], 1):
                    lines.append(f"\n{'=' * 70}")
                    lines.append(f"[{i}] 📄 {paper['title']}")
                    lines.append(f"    👤 作者: {', '.join(paper['authors'][:3])}")
                    lines.append(f"    📖 期刊: {paper['journal']} ({paper['year']})")
                    lines.append(f"    🔗 PMID: {paper['pmid']}")
                    lines.append(f"    🌐 链接: {paper['url']}")

                    sentences = paper.get('variant_sentences', [])
                    if sentences:
                        lines.append(f"\n    ✅ 包含该变异 (找到 {len(sentences)} 处):")

                        for j, sent in enumerate(sentences[:5], 1):
                            lines.append(f"\n    证据 {j} [匹配类型: {sent['match_type']}]:")
                            lines.append(f"    📝 {sent['highlighted']}")

                            if sent.get('keywords'):
                                lines.append(f"    🔑 匹配关键词: {', '.join(sent['keywords'][:3])}")

                        if len(sentences) > 5:
                            lines.append(f"\n    ... 还有 {len(sentences) - 5} 处匹配 ...")
                    else:
                        lines.append("\n    ⚠️ 未找到具体变异信息")

            else:
                lines.append("\n❌ 未找到包含该具体变异的文献")
                lines.append("\n   可能原因:")
                lines.append("   1. 该变异尚未在文献中报道")
                lines.append("   2. 文献中使用不同的命名方式")
                lines.append("   3. 变异信息只在全文正文中，不在摘要里")

        return "\n".join(lines)


# ============ Vercel API Handler ============
def handler(request: Dict[str, Any]) -> Dict[str, Any]:
    """Vercel Serverless Function 入口"""
    headers = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
        'Access-Control-Allow-Headers': 'Content-Type'
    }

    if request.get('method') == 'OPTIONS':
        return {'statusCode': 200, 'headers': headers, 'body': ''}

    try:
        # 解析参数
        query = request.get('queryStringParameters', {}) or {}
        body = {}
        if request.get('method') == 'POST':
            try:
                body = json.loads(request.get('body', '{}') or '{}')
            except:
                body = {}

        # 获取公共参数
        gene = query.get('gene') or body.get('gene')
        input_type = query.get('type') or body.get('type') or 'auto'  # hgvs, vcf, rsid, auto
        format_type = query.get('format') or body.get('format') or 'json'
        api_key = query.get('api_key') or body.get('api_key') or os.environ.get('NCBI_API_KEY')

        if not gene:
            return {
                'statusCode': 400,
                'headers': headers,
                'body': json.dumps({'error': 'Missing parameter: gene'})
            }

        # ============================================================
        # 根据 type 参数选择不同的方法
        # ============================================================
        variant = None

        if input_type == 'hgvs':
            # 使用 from_hgvs 方法
            # 参数: hgvs (必需), gene (必需), hgvs_p (可选)
            hgvs = query.get('hgvs') or body.get('hgvs')
            hgvs_p = query.get('hgvs_p') or body.get('hgvs_p')

            if not hgvs:
                return {
                    'statusCode': 400,
                    'headers': headers,
                    'body': json.dumps({'error': 'Missing parameter: hgvs (required for type=hgvs)'})
                }
            variant = Variant.from_hgvs(hgvs, gene=gene, hgvs_p=hgvs_p)

        elif input_type == 'vcf':
            # 使用 from_vcf 方法
            # 参数: chrom (必需), pos (必需), ref (必需), alt (必需), gene (必需)
            chrom = query.get('chrom') or body.get('chrom')
            pos = query.get('pos') or body.get('pos')
            ref = query.get('ref') or body.get('ref')
            alt = query.get('alt') or body.get('alt')

            if not all([chrom, pos, ref, alt]):
                return {
                    'statusCode': 400,
                    'headers': headers,
                    'body': json.dumps({
                        'error': 'Missing parameters: chrom, pos, ref, alt (required for type=vcf)'
                    })
                }
            variant = Variant.from_vcf(chrom, int(pos), ref, alt, gene=gene)

        elif input_type == 'rsid':
            # 使用 from_rsid 方法
            # 参数: rsid (必需), gene (必需)
            rsid = query.get('rsid') or body.get('rsid')

            if not rsid:
                return {
                    'statusCode': 400,
                    'headers': headers,
                    'body': json.dumps({'error': 'Missing parameter: rsid (required for type=rsid)'})
                }
            variant = Variant.from_rsid(rsid, gene=gene)

        else:  # 'auto' 或默认
            # 使用 from_string 智能解析
            # 参数: gene (必需), variant_str (必需)
            variant_str = query.get('variant') or body.get('variant')

            if not variant_str:
                return {
                    'statusCode': 400,
                    'headers': headers,
                    'body': json.dumps({'error': 'Missing parameter: variant (required for type=auto)'})
                }
            variant = Variant.from_string(gene, variant_str)

        # 执行搜索
        searcher = VariantLiteratureSearch(ncbi_api_key=api_key)
        results = searcher.search(variant, databases=['pubmed', 'clinvar'])

        # 返回结果
        if format_type == 'text':
            text_output = searcher.search_and_format(variant, databases=['pubmed', 'clinvar'])
            return {
                'statusCode': 200,
                'headers': {**headers, 'Content-Type': 'text/plain; charset=utf-8'},
                'body': text_output
            }

        return {
            'statusCode': 200,
            'headers': {**headers, 'Content-Type': 'application/json; charset=utf-8'},
            'body': json.dumps({
                'success': True,
                'input_type': input_type,
                'query': results['query'],
                'total': len(results['results'].get('pubmed', [])),
                'results': results['results'].get('pubmed', [])
            }, ensure_ascii=False, indent=2)
        }

    except Exception as e:
        return {
            'statusCode': 500,
            'headers': headers,
            'body': json.dumps({'success': False, 'error': str(e)})
        }
