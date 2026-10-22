#!/usr/bin/env python3
from pathlib import Path
import argparse
import warnings
import json
import pandas as pd
import plotly.graph_objects as go


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot expression of a selected gene as an interactive violin plot, "
                    "with a dropdown to split samples by any metadata feature."
    )
    parser.add_argument(
        "--expression_matrix",
        type=Path,
        required=True,
        help="Path to the gene expression matrix (genes on rows, samples on columns, tab-separated)."
    )
    parser.add_argument(
        "--gene",
        type=str,
        required=True,
        help="Name of the gene to plot (must match a row index in the expression matrix)."
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Path to the sample metadata file (CSV, samples on rows)."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=False,
        default=None,
        help="Path for the output HTML plot file. Defaults to 'violin_<gene>.html'."
    )
    return parser.parse_args()


def load_and_validate(expression_matrix_path, gene, metadata_path):
    # Load expression matrix
    expr = pd.read_csv(expression_matrix_path, index_col=0, sep="\t")

    # Validate gene
    if gene not in expr.index:
        raise ValueError(
            f"Gene '{gene}' not found in the expression matrix. "
            f"Available genes (first 10): {list(expr.index[:10])}"
        )

    # Load metadata
    meta = pd.read_csv(metadata_path, index_col=0, sep=",")

    # Warn about numeric features (informational, non-fatal)
    numeric_features = [c for c in meta.columns if pd.api.types.is_numeric_dtype(meta[c])]
    if numeric_features:
        warnings.warn(
            f"The following metadata features are numeric and will be included in the "
            f"dropdown but are not ideal for grouping: {numeric_features}. "
            "Consider using categorical features for meaningful violin plots.",
            UserWarning
        )

    # Check sample overlap
    meta_samples = set(meta.index)
    expr_samples = set(expr.columns)
    missing = meta_samples - expr_samples
    if missing:
        raise ValueError(
            f"The following samples from the metadata are not present in the expression matrix: {missing}"
        )

    # Restrict to metadata samples (ignore extra expression matrix columns)
    common_samples = [s for s in meta.index if s in expr_samples]
    gene_expr = expr.loc[gene, common_samples]

    # Build tidy dataframe: one row per sample, expression + all metadata features
    df = meta.loc[common_samples].copy()
    df.insert(0, "expression", gene_expr)
    df.index.name = "sample"
    df = df.reset_index()

    return df


def build_html(df, gene, output_path):
    """
    Produce a self-contained HTML file embedding all data as JSON.
    The page renders a Plotly violin plot and a <select> dropdown that
    lets the user pick any metadata feature; the plot updates client-side.
    """
    features = [c for c in df.columns if c not in ("sample", "expression")]
    samples = df["sample"].tolist()
    expression = df["expression"].tolist()

    # Build a dict: feature -> list of group labels (one per sample, in order)
    feature_groups = {f: df[f].astype(str).tolist() for f in features}

    data_json = json.dumps({
        "gene": gene,
        "samples": samples,
        "expression": expression,
        "features": features,
        "feature_groups": feature_groups,
    })

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Violin plot – {gene}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{
    font-family: Arial, sans-serif;
    margin: 30px;
    background: #fff;
  }}
  #controls {{
    margin-bottom: 18px;
    display: flex;
    align-items: center;
    gap: 12px;
  }}
  #controls label {{
    font-size: 15px;
    font-weight: bold;
    color: #333;
  }}
  #feature-select {{
    font-size: 14px;
    padding: 5px 10px;
    border-radius: 5px;
    border: 1px solid #bbb;
    background: #f9f9f9;
    cursor: pointer;
  }}
  #warn-box {{
    display: none;
    margin-top: 8px;
    padding: 6px 12px;
    background: #fff8e1;
    border-left: 4px solid #f9a825;
    font-size: 13px;
    color: #795548;
    border-radius: 3px;
  }}
  #plot {{
    width: 100%;
    height: 560px;
  }}
</style>
</head>
<body>

<div id="controls">
  <label for="feature-select">Group samples by:</label>
  <select id="feature-select">
    <option value="__none__">— All samples —</option>
  </select>
</div>
<div id="warn-box">⚠ The selected feature is numeric. Violin plots work best with categorical groupings.</div>
<div id="plot"></div>

<script>
const DATA = {data_json};

// Detect which features are numeric (all values parse as numbers)
const numericFeatures = new Set(
  DATA.features.filter(f =>
    DATA.feature_groups[f].every(v => v !== "" && !isNaN(Number(v)))
  )
);

// Populate dropdown
const sel = document.getElementById("feature-select");
DATA.features.forEach(f => {{
  const opt = document.createElement("option");
  opt.value = f;
  opt.textContent = f;
  sel.appendChild(opt);
}});

function getTraces(feature) {{
  if (feature === "__none__") {{
    // Single violin for all samples
    return [{{
      type: "violin",
      y: DATA.expression,
      name: DATA.gene,
      text: DATA.samples,
      hovertemplate: "<b>%{{text}}</b><br>Expression: %{{y}}<extra></extra>",
      box: {{ visible: true }},
      meanline: {{ visible: true }},
      points: "all",
      jitter: 0.3,
      pointpos: 0,
      marker: {{ size: 6 }},
    }}];
  }}

  const groups = DATA.feature_groups[feature];
  const uniqueGroups = [...new Set(groups)].sort();

  return uniqueGroups.map(g => {{
    const indices = groups.reduce((acc, val, i) => {{ if (val === g) acc.push(i); return acc; }}, []);
    return {{
      type: "violin",
      y: indices.map(i => DATA.expression[i]),
      name: g,
      text: indices.map(i => DATA.samples[i]),
      hovertemplate: "<b>%{{text}}</b><br>Expression: %{{y}}<extra></extra>",
      box: {{ visible: true }},
      meanline: {{ visible: true }},
      points: "all",
      jitter: 0.3,
      pointpos: 0,
      marker: {{ size: 6 }},
    }};
  }});
}}

function getLayout(feature) {{
  const groupLabel = feature === "__none__" ? "All samples" : feature;
  return {{
    title: {{ text: `Expression of <i>${{DATA.gene}}</i> — ${{groupLabel}}`, x: 0.5 }},
    yaxis: {{ title: `${{DATA.gene}} expression` }},
    xaxis: {{ title: feature === "__none__" ? "" : feature }},
    legend: {{ title: {{ text: feature === "__none__" ? "" : feature }} }},
    template: "plotly_white",
    violingap: 0.3,
    violinmode: "group",
  }};
}}

function renderPlot(feature) {{
  const warnBox = document.getElementById("warn-box");
  warnBox.style.display = (feature !== "__none__" && numericFeatures.has(feature)) ? "block" : "none";
  Plotly.react("plot", getTraces(feature), getLayout(feature));
}}

sel.addEventListener("change", () => renderPlot(sel.value));

// Initial render: single violin, all samples
renderPlot("__none__");
</script>
</body>
</html>
"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


def main():
    args = parse_args()

    output = args.output if args.output is not None else Path(f"violin_{args.gene}.html")

    df = load_and_validate(args.expression_matrix, args.gene, args.metadata)

    build_html(df, args.gene, output)
    print(f"Plot saved to: {output}")


if __name__ == "__main__":
    main()
