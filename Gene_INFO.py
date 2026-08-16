# app.py - 基因信息查询API服务（支持基因符号和NM号查询）

from flask import Flask, request, jsonify
import requests
import xml.etree.ElementTree as ET
import re
from datetime import datetime

app = Flask(__name__)

# 缓存
cache = {}


def get_gene_id_from_symbol(symbol):
    """
    通过基因符号搜索Gene ID
    """
    try:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            "db": "gene",
            "term": f"{symbol}[Gene Name] AND human[Organism]",
            "retmode": "json",
            "retmax": "1"
        }

        r = requests.get(search_url, params=search_params, timeout=10)
        r.raise_for_status()
        data = r.json()

        if data["esearchresult"]["idlist"]:
            return data["esearchresult"]["idlist"][0]
        return None

    except Exception as e:
        print(f"搜索基因符号 {symbol} 出错: {e}")
        return None


def get_gene_id_from_nm(nm_id):
    """
    通过NM号搜索Gene ID
    """
    try:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            "db": "gene",
            "term": nm_id,
            "retmode": "json",
            "retmax": "1"
        }

        r = requests.get(search_url, params=search_params, timeout=10)
        r.raise_for_status()
        data = r.json()

        if data["esearchresult"]["idlist"]:
            return data["esearchresult"]["idlist"][0]
        return None

    except Exception as e:
        print(f"搜索NM号 {nm_id} 出错: {e}")
        return None


def get_gene_info_by_id(gene_id):
    """
    通过Gene ID获取基因完整信息
    """
    # 检查缓存
    if gene_id in cache:
        return cache[gene_id]

    try:
        # 获取Gene XML
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "xml"
        }

        r = requests.get(fetch_url, params=fetch_params, timeout=10)
        r.raise_for_status()

        # 解析XML提取所有信息
        root = ET.fromstring(r.text)
        info = parse_gene_xml(root)

        # 添加Gene ID
        info['gene_id'] = gene_id

        # 缓存结果
        cache[gene_id] = info
        return info

    except Exception as e:
        print(f"获取基因信息 {gene_id} 出错: {e}")
        return None


def parse_gene_xml(root):
    """解析Gene XML，提取所有信息"""
    info = {
        'gene_symbol': None,
        'gene_name': None,
        'chromosome': None,
        'strand': None,
        'location': None,
        'map_location': None,
        'summary': None,
        'other_names': [],
        'aliases': [],
        'transcripts': [],
        'mane_transcript': None,
        'genomic_context': None,
        'search_query': None  # 存储原始查询
    }

    try:
        # 查找Gene元素
        for gene in root.findall(".//Entrezgene"):
            # 基因符号
            for gene_symbol in gene.findall(".//Gene-track_geneid"):
                info['gene_symbol'] = gene_symbol.text

            # 基因名称
            for gene_name in gene.findall(".//Gene-ref_locus"):
                info['gene_name'] = gene_name.text

            # 其他名称
            for name in gene.findall(".//Gene-ref_syn"):
                if name.text and name.text not in info['other_names']:
                    info['other_names'].append(name.text)

            # 染色体位置
            for location in gene.findall(".//Gene-commentary_comment"):
                for comment in location.findall(".//Gene-commentary"):
                    for source in comment.findall(".//Gene-commentary_source"):
                        for organism in source.findall(".//Org-ref_taxname"):
                            if "chromosome" in organism.text.lower() or "chr" in organism.text.lower():
                                info['chromosome'] = organism.text

            # 染色体定位
            for map_loc in gene.findall(".//Gene-track_map-location"):
                info['map_location'] = map_loc.text

            # 链信息
            for seq_loc in gene.findall(".//Seq-loc"):
                for seq_interval in seq_loc.findall(".//Seq-interval"):
                    for seq_id in seq_interval.findall(".//Seq-id"):
                        for seq_id_other in seq_id.findall(".//Other-source"):
                            for seq_id_name in seq_id_other.findall(".//Other-source_anchor"):
                                if seq_id_name.text and (
                                        'NC_' in seq_id_name.text or 'chr' in seq_id_name.text.lower()):
                                    info['chromosome'] = seq_id_name.text

                    for start in seq_interval.findall(".//Seq-interval_from"):
                        info['location'] = start.text
                    for end in seq_interval.findall(".//Seq-interval_to"):
                        info['location'] = f"{info['location']}-{end.text}" if info['location'] else end.text

                    for strand in seq_interval.findall(".//Na-strand"):
                        value = strand.get("value")
                        if value == "plus":
                            info['strand'] = "+"
                        elif value == "minus":
                            info['strand'] = "-"

            # 基因摘要
            for summary in gene.findall(".//Entrezgene_summary"):
                info['summary'] = summary.text

        # 查找转录本信息
        for gb_seq in root.findall(".//GBSeq"):
            for gb_feature in gb_seq.findall(".//GBFeature"):
                for gb_key in gb_feature.findall(".//GBFeature_key"):
                    if gb_key.text in ["mRNA", "transcript", "cds"]:
                        for product in gb_feature.findall(".//GBFeature_quals/GBQualifier"):
                            for gb_qual_name in product.findall(".//GBQualifier_name"):
                                if gb_qual_name.text == "transcript_id":
                                    for gb_qual_value in product.findall(".//GBQualifier_value"):
                                        if gb_qual_value.text and gb_qual_value.text not in info['transcripts']:
                                            info['transcripts'].append(gb_qual_value.text)

                                        # 检查是否为MANE转录本
                                        if gb_qual_value.text:
                                            mane_patterns = [
                                                r'NM_\d+\.\d+',
                                                r'XR_\d+\.\d+',
                                                r'NR_\d+\.\d+'
                                            ]
                                            for pattern in mane_patterns:
                                                if re.match(pattern, gb_qual_value.text):
                                                    info['mane_transcript'] = gb_qual_value.text

        # 构建基因组上下文
        if info['chromosome'] and info['location']:
            strand_symbol = info['strand'] if info['strand'] else '?'
            info['genomic_context'] = f"{info['chromosome']}({strand_symbol}){info['location']}"

        # 清理
        if info['summary']:
            info['summary'] = info['summary'].strip()
        if info['map_location']:
            info['map_location'] = info['map_location'].strip()

    except Exception as e:
        print(f"解析XML出错: {e}")

    return info


def detect_query_type(query):
    """
    检测查询类型
    返回: 'gene_symbol' 或 'transcript_id'
    """
    query = query.strip()

    # 检查是否为转录本ID (NM_, NR_, XM_, XR_开头)
    transcript_patterns = [
        r'^NM_\d+',
        r'^NR_\d+',
        r'^XM_\d+',
        r'^XR_\d+'
    ]

    for pattern in transcript_patterns:
        if re.match(pattern, query, re.IGNORECASE):
            return 'transcript_id'

    # 默认为基因符号
    return 'gene_symbol'


def get_gene_info(query):
    """
    统一的基因查询接口
    支持基因符号或转录本ID
    """
    query = query.strip()
    query_type = detect_query_type(query)

    # 检查缓存
    cache_key = f"{query_type}:{query}"
    if cache_key in cache:
        return cache[cache_key]

    try:
        gene_id = None

        if query_type == 'transcript_id':
            # 通过转录本ID查询
            gene_id = get_gene_id_from_nm(query)
        else:
            # 通过基因符号查询
            gene_id = get_gene_id_from_symbol(query)

        if not gene_id:
            return None

        # 获取基因信息
        info = get_gene_info_by_id(gene_id)

        if info:
            info['search_query'] = query
            info['query_type'] = query_type
            # 缓存
            cache[cache_key] = info

        return info

    except Exception as e:
        print(f"查询基因出错: {e}")
        return None


@app.route('/')
def index():
    """API首页"""
    return jsonify({
        "service": "基因信息查询API",
        "version": "3.0",
        "description": "支持通过基因符号(如BRCA1)或转录本ID(如NM_001142864)查询",
        "endpoints": {
            "/": "API文档",
            "/search": "查询基因信息（支持基因符号或NM号）",
            "/batch": "批量查询多个基因",
            "/strand": "仅查询正负链",
            "/health": "健康检查",
            "/cache": "查看缓存状态"
        },
        "usage": {
            "/search?q=BRCA1": "通过基因符号查询",
            "/search?q=NM_001142864": "通过转录本ID查询",
            "/strand?q=BRCA1": "仅查询正负链",
            "/batch": "POST请求批量查询"
        },
        "examples": {
            "gene_symbol": "BRCA1, TP53, EGFR, KRAS, BRAF",
            "transcript_id": "NM_001142864, NM_000546, NM_005228"
        }
    })


@app.route('/search')
def search():
    """
    查询基因信息（支持基因符号或转录本ID）
    参数: q (查询字符串)
    示例: /search?q=BRCA1 或 /search?q=NM_001142864
    """
    query = request.args.get('q')
    if not query:
        return jsonify({
            "error": "Missing parameter 'q'",
            "example": "/search?q=BRCA1 或 /search?q=NM_001142864"
        }), 400

    info = get_gene_info(query)

    if not info:
        return jsonify({
            "error": f"未找到基因: {query}",
            "query": query,
            "suggestion": "请检查基因符号或转录本ID是否正确"
        }), 404

    return jsonify({
        "query": query,
        "query_type": info.get('query_type'),
        "gene_id": info.get('gene_id'),
        "gene_symbol": info.get('gene_symbol'),
        "gene_name": info.get('gene_name'),
        "chromosome": info.get('chromosome'),
        "strand": info.get('strand'),
        "location": info.get('location'),
        "map_location": info.get('map_location'),
        "genomic_context": info.get('genomic_context'),
        "other_names": info.get('other_names', []),
        "transcripts": info.get('transcripts', [])[:10],  # 只返回前10个
        "total_transcripts": len(info.get('transcripts', [])),
        "mane_transcript": info.get('mane_transcript'),
        "summary": info.get('summary', '')[:500] + ('...' if len(info.get('summary', '')) > 500 else ''),
        "cached": f"{info.get('query_type')}:{query}" in cache
    })


@app.route('/strand')
def strand_api():
    """
    仅查询正负链
    参数: q (查询字符串)
    示例: /strand?q=BRCA1 或 /strand?q=NM_001142864
    """
    query = request.args.get('q')
    if not query:
        return jsonify({
            "error": "Missing parameter 'q'",
            "example": "/strand?q=BRCA1"
        }), 400

    info = get_gene_info(query)

    if not info:
        return jsonify({
            "error": f"未找到基因: {query}",
            "query": query
        }), 404

    return jsonify({
        "query": query,
        "gene_symbol": info.get('gene_symbol'),
        "strand": info.get('strand'),
        "cached": f"{info.get('query_type')}:{query}" in cache
    })


@app.route('/batch', methods=['POST'])
def batch_api():
    """
    批量查询多个基因
    参数: JSON数组，可以是基因符号或转录本ID
    示例: ["BRCA1", "TP53", "NM_001142864"]
    """
    data = request.get_json()

    if not data:
        return jsonify({
            "error": "Invalid request body",
            "example": '["BRCA1", "TP53", "NM_001142864"]'
        }), 400

    if not isinstance(data, list):
        return jsonify({
            "error": "Request body must be a JSON array"
        }), 400

    results = []
    for query in data[:50]:  # 限制最多50个
        info = get_gene_info(query)

        if info:
            results.append({
                "query": query,
                "query_type": info.get('query_type'),
                "gene_symbol": info.get('gene_symbol'),
                "gene_name": info.get('gene_name'),
                "chromosome": info.get('chromosome'),
                "strand": info.get('strand'),
                "location": info.get('location'),
                "mane_transcript": info.get('mane_transcript'),
                "cached": f"{info.get('query_type')}:{query}" in cache
            })
        else:
            results.append({
                "query": query,
                "error": "Not found"
            })

    return jsonify({
        "total": len(results),
        "success_count": sum(1 for r in results if 'error' not in r),
        "results": results,
        "cached_count": sum(1 for r in results if r.get('cached', False))
    })


@app.route('/health')
def health():
    """健康检查"""
    return jsonify({
        "status": "healthy",
        "cache_size": len(cache),
        "timestamp": datetime.now().isoformat()
    })


@app.route('/cache')
def cache_status():
    """查看缓存状态"""
    cache_items = list(cache.keys())
    return jsonify({
        "cache_size": len(cache),
        "cached_items": cache_items[:20],  # 只显示前20个
        "total_cached": len(cache_items),
        "timestamp": datetime.now().isoformat()
    })


@app.route('/clear_cache', methods=['POST'])
def clear_cache():
    """清空缓存"""
    global cache
    cache = {}
    return jsonify({
        "message": "缓存已清空",
        "timestamp": datetime.now().isoformat()
    })


if __name__ == '__main__':
    print("=" * 60)
    print("基因信息查询API V3.0")
    print("=" * 60)
    print("本地访问: http://localhost:9876")
    print("")
    print("可用接口:")
    print("  1. GET  /search?q=BRCA1          - 通过基因符号查询")
    print("  2. GET  /search?q=NM_001142864   - 通过转录本ID查询")
    print("  3. GET  /strand?q=BRCA1          - 仅查询正负链")
    print("  4. POST /batch                   - 批量查询")
    print("  5. GET  /health                  - 健康检查")
    print("  6. GET  /cache                   - 查看缓存")
    print("  7. POST /clear_cache             - 清空缓存")
    print("=" * 60)
    print("")
    print("使用示例:")
    print("  http://localhost:9876/search?q=BRCA1")
    print("  http://localhost:9876/search?q=NM_001142864")
    print("  http://localhost:9876/strand?q=TP53")
    print("=" * 60)
    app.run(host='0.0.0.0', port=9876, debug=True)