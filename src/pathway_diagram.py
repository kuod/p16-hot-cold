"""
p16/CDKN2A cell-cycle pathway diagram.

Provides:
  P16_PATHWAY_MERMAID  — Mermaid flowchart source string
  render_to_png()      — saves a PNG via the mermaid.ink API (no CLI required)
  streamlit_html()     — returns self-contained HTML for st.components.v1.html
"""

import base64
import json
from pathlib import Path

# ── Mermaid source ─────────────────────────────────────────────────────────────
# Nodes use ["label"] syntax; <br/> gives line breaks inside labels.
# Cancer escape routes are red dashed edges converging on their target node.

P16_PATHWAY_MERMAID = """
flowchart TD

    %% ── Upstream inducers ────────────────────────────────────────────────────
    ONC["Oncogenic activation<br/>RAS · MYC · BRAF"]
    DNADMG["DNA damage<br/>Replication stress · IR"]
    TGFB["TGF-β signalling"]

    %% ── Sensor / adaptor layer ───────────────────────────────────────────────
    ARF["p14ARF<br/>(CDKN2A alt. reading frame)"]
    ATM["ATM / ATR<br/>DNA-damage kinases"]
    CHK["CHK1 / CHK2"]

    %% ── INK4 family — direct CDK4/6 inhibitors ──────────────────────────────
    P16["p16^INK4a (CDKN2A)<br/>CDK4/6 inhibitor"]
    P15["p15^INK4b (CDKN2B)<br/>CDK4/6 inhibitor"]
    P18["p18^INK4c (CDKN2C)<br/>CDK4/6 inhibitor"]

    %% ── p53 branch ───────────────────────────────────────────────────────────
    MDM2["MDM2<br/>p53 E3 ubiquitin ligase"]
    P53["p53 (TP53)<br/>transcription factor"]
    P21["p21^CIP1 (CDKN1A)<br/>CDK2 inhibitor"]

    %% ── Cell-cycle machinery ─────────────────────────────────────────────────
    CDK46["CDK4 / CDK6<br/>+ Cyclin D"]
    CDK2["CDK2 / CDK1<br/>+ Cyclin E / A"]
    RBi["Rb (RB1)<br/>hypophosphorylated"]
    RBp["Rb-P<br/>hyperphosphorylated"]
    E2Fi["E2F sequestered<br/>→ G1 arrest"]
    E2Ff["E2F released<br/>→ S-phase entry"]

    %% ── Cell fate outputs ────────────────────────────────────────────────────
    ARREST["Cell cycle arrest<br/>(senescence / quiescence)"]
    PROLIF["Proliferation"]
    SASP["SASP secretion<br/>IL-6 · CXCL8 · MMPs<br/>→ TME remodelling"]

    %% ── Cancer escape nodes (coloured red) ──────────────────────────────────
    C1["CDKN2A del / mut<br/>~30 % pan-cancer<br/>→ p16 loss"]
    C2["CDK4 / CDK6 amplification<br/>→ kinase overactivation"]
    C3["CCND1 amplification<br/>→ Cyclin D excess"]
    C4["RB1 del / mut<br/>~15 % pan-cancer<br/>→ Rb loss"]
    C5["TP53 mutation<br/>~50 % pan-cancer<br/>→ p53 absent"]
    C6["MDM2 amplification<br/>→ p53 constitutively degraded"]
    C7["ATM del / mut<br/>→ DDR impaired"]
    C8["PTEN del / mut<br/>→ PI3K / AKT elevated<br/>→ CDK activation"]

    %% ── Normal signalling flow ───────────────────────────────────────────────
    ONC -->|"ARF induction"| ARF
    ONC -->|"drives Cyclin D"| CDK46
    DNADMG --> ATM
    ATM --> CHK
    CHK -->|"phosphorylates"| P53
    TGFB --> P15
    ARF -->|"sequesters"| MDM2
    MDM2 -.->|"ubiquitinates / degrades"| P53
    P53 -->|"transcribes"| P21
    P21 -->|"inhibits"| CDK2
    P16 -->|"inhibits"| CDK46
    P15 -->|"inhibits"| CDK46
    P18 -->|"inhibits"| CDK46
    CDK46 -->|"phosphorylates Rb"| RBp
    CDK2 -->|"phosphorylates Rb"| RBp
    RBi --> E2Fi
    RBp --> E2Ff
    E2Fi --> ARREST
    E2Ff --> PROLIF
    ARREST --> SASP

    %% ── Cancer escape bypass routes (dashed) ────────────────────────────────
    C1 -. "p16 absent" .-> CDK46
    C2 -. "overactive kinase" .-> CDK46
    C3 -. "excess Cyclin D" .-> CDK46
    C4 -. "Rb absent" .-> E2Ff
    C5 -. "p53 absent" .-> P53
    C6 -. "p53 degraded" .-> MDM2
    C7 -. "DDR blunted" .-> ATM
    C8 -. "AKT promotes CDK" .-> CDK46

    %% ── Feedback loops ───────────────────────────────────────────────────────
    %% SASP downstream effects (separate node to avoid ARREST→SASP cycle)
    SASP_OUT["SASP downstream<br/>NK/CTL clearance · paracrine p21"]

    SASP --> SASP_OUT
    SASP_OUT -->|"paracrine p21 induction"| P21
    P53 -->|"transcribes CDKN2A<br/>(positive feedback)"| P16
    TGFB -->|"TGF-β paradox:<br/>pro-tumourigenic in TME"| SASP

    %% ── Style classes ────────────────────────────────────────────────────────
    classDef inducer    fill:#fde68a,stroke:#d97706,color:#78350f
    classDef sensor     fill:#e0e7ff,stroke:#6366f1,color:#312e81
    classDef suppressor fill:#dcfce7,stroke:#16a34a,color:#14532d
    classDef kinase     fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    classDef rb         fill:#e2e8f0,stroke:#64748b,color:#1e293b
    classDef arrest     fill:#fef9c3,stroke:#ca8a04,color:#78350f
    classDef escape     fill:#fee2e2,stroke:#dc2626,color:#7f1d1d
    classDef prolif     fill:#fca5a5,stroke:#ef4444,color:#7f1d1d

    class ONC,DNADMG,TGFB inducer
    class ARF,ATM,CHK sensor
    class P16,P15,P18,P21,P53,MDM2 suppressor
    class CDK46,CDK2 kinase
    class RBi,RBp,E2Fi,E2Ff rb
    class ARREST,SASP,SASP_OUT arrest
    class PROLIF prolif
    class C1,C2,C3,C4,C5,C6,C7,C8 escape
""".strip()


# ── Rendering helpers ──────────────────────────────────────────────────────────

def render_to_png(output_path: Path, timeout: int = 30) -> None:
    """
    Render P16_PATHWAY_MERMAID to a PNG file via the mermaid.ink public API.
    Requires network access; no local CLI installation needed.
    """
    import requests

    # mermaid.ink accepts base64-encoded JSON {code, mermaid}
    payload = json.dumps({"code": P16_PATHWAY_MERMAID, "mermaid": {"theme": "default"}})
    encoded = base64.urlsafe_b64encode(payload.encode()).decode()
    url = f"https://mermaid.ink/img/{encoded}?type=png"

    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(resp.content)
    print(f"Saved pathway diagram → {output_path}")


def annotate_with_data(alteration_rates: dict) -> str:
    """
    Return a copy of P16_PATHWAY_MERMAID with cancer-escape node labels
    updated to include observed alteration rates from the dataset.

    Parameters
    ----------
    alteration_rates : dict
        Mapping from gene symbol to pan-cancer alteration fraction (0–1).
        Recognised keys: CDKN2A, RB1, TP53, ATM, PTEN.
        Missing keys are silently ignored (label left unchanged).

    Returns
    -------
    str
        Modified Mermaid source with data-driven % annotations.

    Example
    -------
    >>> rates = {"CDKN2A": 0.31, "TP53": 0.48, "RB1": 0.14, "ATM": 0.05, "PTEN": 0.09}
    >>> diagram = annotate_with_data(rates)
    """
    import re

    diagram = P16_PATHWAY_MERMAID

    # Nodes that already contain a "~N % pan-cancer" placeholder
    pct_nodes = {
        "CDKN2A": r'(C1\["CDKN2A del / mut<br/>)~\d+ % pan-cancer',
        "RB1":    r'(C4\["RB1 del / mut<br/>)~\d+ % pan-cancer',
        "TP53":   r'(C5\["TP53 mutation<br/>)~\d+ % pan-cancer',
    }
    for gene, pattern in pct_nodes.items():
        if gene in alteration_rates:
            pct = round(100 * alteration_rates[gene])
            diagram = re.sub(pattern, rf'\g<1>{pct} % pan-cancer', diagram)

    # Nodes that do NOT yet have a rate — append one before the closing quote
    append_nodes = {
        "ATM":  r'(C7\["ATM del / mut<br/>→ DDR impaired)(")',
        "PTEN": r'(C8\["PTEN del / mut<br/>→ PI3K / AKT elevated<br/>→ CDK activation)(")',
    }
    for gene, pattern in append_nodes.items():
        if gene in alteration_rates:
            pct = round(100 * alteration_rates[gene])
            diagram = re.sub(pattern, rf'\g<1><br/>{pct} % pan-cancer\g<2>', diagram)

    return diagram


def streamlit_html(height: int = 900) -> str:
    """
    Return a self-contained HTML string that renders the pathway diagram
    interactively using Mermaid.js (CDN).  Pass to st.components.v1.html().
    """
    escaped = P16_PATHWAY_MERMAID.replace("</", "<\\/")
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
  <script>
    mermaid.initialize({{
      startOnLoad: true,
      theme: 'default',
      flowchart: {{ useMaxWidth: true, htmlLabels: true }}
    }});
  </script>
  <style>
    body {{ margin: 0; padding: 12px; background: #ffffff; font-family: sans-serif; }}
    .mermaid svg {{ max-width: 100% !important; height: auto; }}
    .legend {{ margin-top: 16px; font-size: 12px; color: #555; }}
    .legend span {{ display: inline-block; width: 14px; height: 14px;
                    border-radius: 2px; margin-right: 4px; vertical-align: middle; }}
  </style>
</head>
<body>
  <div class="mermaid">
{P16_PATHWAY_MERMAID}
  </div>
  <div class="legend">
    <b>Legend:</b>&nbsp;
    <span style="background:#fde68a;border:1px solid #d97706"></span>Upstream inducers&nbsp;&nbsp;
    <span style="background:#e0e7ff;border:1px solid #6366f1"></span>DNA-damage sensors&nbsp;&nbsp;
    <span style="background:#dcfce7;border:1px solid #16a34a"></span>Tumour suppressors&nbsp;&nbsp;
    <span style="background:#dbeafe;border:1px solid #3b82f6"></span>Cyclin-dependent kinases&nbsp;&nbsp;
    <span style="background:#e2e8f0;border:1px solid #64748b"></span>Rb / E2F&nbsp;&nbsp;
    <span style="background:#fef9c3;border:1px solid #ca8a04"></span>Arrest / SASP&nbsp;&nbsp;
    <span style="background:#fee2e2;border:1px solid #dc2626"></span>Cancer escape routes (dashed edges)
  </div>
</body>
</html>"""
