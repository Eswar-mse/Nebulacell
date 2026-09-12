"""
Rule-based (regex) natural language query parser.

No LLM API call — zero cost. Handles a constrained set of query patterns
and returns a structured filter dict, or None if nothing matched.

This intentionally does NOT execute arbitrary code from user input;
it only ever produces a filter dict from a fixed small grammar.
"""

import re

GENE_LEVEL_RE = re.compile(
    r"(?:show|highlight)\s+cells?\s+(?:with|expressing)\s+(high|low)\s+([A-Za-z0-9\-\.]+)",
    re.IGNORECASE,
)

CLUSTER_RE = re.compile(
    r"(?:show|highlight)\s+(?:only\s+)?(?:cluster|domain)\s+(-?\d+)",
    re.IGNORECASE,
)

RESET_RE = re.compile(r"(reset|clear|show all|show everything)", re.IGNORECASE)


def parse_query(query: str) -> dict | None:
    """
    Parse a natural language query into a structured filter dict.

    Returns one of:
      {"type": "gene_level", "gene": str, "level": "high"|"low"}
      {"type": "cluster", "cluster_id": str}
      {"type": "reset"}
      None  (if nothing matched)
    """
    query = query.strip()

    if RESET_RE.search(query):
        return {"type": "reset"}

    m = GENE_LEVEL_RE.search(query)
    if m:
        level, gene = m.group(1).lower(), m.group(2)
        return {"type": "gene_level", "gene": gene, "level": level}

    m = CLUSTER_RE.search(query)
    if m:
        return {"type": "cluster", "cluster_id": m.group(1)}

    return None


def validate_filter(filter_dict: dict, adata) -> tuple[bool, str]:
    """
    Validate a parsed filter against the actual loaded dataset schema.
    Returns (is_valid, error_message).
    """
    if filter_dict is None:
        return False, "Couldn't understand that query. Try: 'show cells with high GeneX' or 'show cluster 3'."

    ftype = filter_dict["type"]

    if ftype == "reset":
        return True, ""

    if ftype == "gene_level":
        gene = filter_dict["gene"]
        if gene not in adata.var_names:
            return False, f"Gene '{gene}' not found in this dataset."
        return True, ""

    if ftype == "cluster":
        cid = filter_dict["cluster_id"]
        valid_ids = set(adata.obs["spatial_domain"].astype(str).unique())
        if cid not in valid_ids:
            return False, f"Cluster '{cid}' not found. Valid clusters: {sorted(valid_ids)}"
        return True, ""

    return False, "Unrecognized filter type."
