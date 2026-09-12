import numpy as np
import scanpy as sc
import dash
from dash import html, dcc, Input, Output, State
import plotly.graph_objects as go

from app.data.loader import load_demo_dataset, load_h5ad, get_spatial_coords
from app.data.image_overlay import get_histology_image, image_array_to_data_uri
from app.clustering.hdbscan_runner import run_spatial_clustering
from app.callbacks.hover import compute_top_genes, get_gene_expression, get_gene_options
from app.nlquery.parser import parse_query, validate_filter

app = dash.Dash(__name__)
app.title = "NebulaCell"


class AppState:
    def __init__(self):
        self.adata = None
        self.coords = None          # scaled to match the histology image's pixel space (or raw if no image)
        self.top_genes = []
        self.gene_options = []
        self.default_gene = None
        self.dataset_label = ""
        self.image_uri = None
        self.image_width = None
        self.image_height = None

    def load(self, adata, label: str):
        if getattr(adata, "isbacked", False):
            # NOTE: true out-of-core computation isn't implemented for the
            # clustering/hover/gene-coloring pipeline yet — everything downstream
            # needs a dense in-memory array. Backed mode currently only saves
            # time/memory at the initial file peek; we materialize here before
            # running any real computation. Real streaming support is a v2 item.
            print(f"[loader] materializing backed dataset into memory ({adata.n_obs} cells)...")
            adata = adata.to_memory()
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)

        adata = run_spatial_clustering(adata)
        self.adata = adata

        raw_coords = get_spatial_coords(adata)
        img, scale_factor, sample_key = get_histology_image(adata)

        if img is not None:
            self.coords = raw_coords * scale_factor
            self.image_uri = image_array_to_data_uri(img)
            self.image_height, self.image_width = img.shape[0], img.shape[1]
        else:
            self.coords = raw_coords
            self.image_uri = None
            self.image_width = None
            self.image_height = None

        self.top_genes = compute_top_genes(adata, n_top=3)
        self.gene_options = get_gene_options(adata)
        self.default_gene = self.gene_options[0] if self.gene_options else None
        self.dataset_label = label


state = AppState()
state.load(load_demo_dataset(), "demo dataset (Visium fluo, 2800 cells)")


def compute_filter_mask(filter_dict: dict) -> np.ndarray:
    n = state.coords.shape[0]

    if filter_dict is None or filter_dict["type"] == "reset":
        return np.ones(n, dtype=bool)

    if filter_dict["type"] == "gene_level":
        expr = get_gene_expression(state.adata, filter_dict["gene"])
        threshold = np.percentile(expr, 75 if filter_dict["level"] == "high" else 25)
        return expr >= threshold if filter_dict["level"] == "high" else expr <= threshold

    if filter_dict["type"] == "cluster":
        return (state.adata.obs["spatial_domain"].astype(str) == filter_dict["cluster_id"]).to_numpy()

    return np.ones(n, dtype=bool)


def build_figure(color_mode: str, selected_gene: str, mask: np.ndarray, show_image: bool) -> go.Figure:
    if color_mode == "spatial_domain":
        color_vals = state.adata.obs["spatial_domain"].cat.codes
        colorscale = "Turbo"
    elif color_mode == "gene" and selected_gene:
        color_vals = get_gene_expression(state.adata, selected_gene)
        colorscale = "Viridis"
    elif "leiden" in state.adata.obs:
        color_vals = state.adata.obs["leiden"].cat.codes
        colorscale = "Plasma"
    else:
        color_vals = state.adata.obs["spatial_domain"].cat.codes
        colorscale = "Turbo"

    gene_expr = get_gene_expression(state.adata, selected_gene) if selected_gene else None
    hover_text = []
    for i in range(state.coords.shape[0]):
        lines = [f"cell {i}", f"top genes: {state.top_genes[i]}"]
        if gene_expr is not None:
            lines.append(f"{selected_gene}: {gene_expr[i]:.3f}")
        hover_text.append("<br>".join(lines))

    opacities = np.where(mask, 1.0, 0.08)
    marker_size = 5 if (show_image and state.image_uri) else 6

    fig = go.Figure(
        data=go.Scattergl(
            x=state.coords[:, 0],
            y=state.coords[:, 1],
            mode="markers",
            marker=dict(
                size=marker_size,
                color=color_vals,
                colorscale=colorscale,
                opacity=opacities,
                line=dict(width=0),
            ),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
        )
    )

    layout_images = []
    if show_image and state.image_uri:
        layout_images.append(
            dict(
                source=state.image_uri,
                xref="x",
                yref="y",
                x=0,
                y=0,
                sizex=state.image_width,
                sizey=state.image_height,
                sizing="stretch",
                layer="below",
                opacity=1.0,
            )
        )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0a0a0f",
        plot_bgcolor="#0a0a0f",
        yaxis=dict(scaleanchor="x", scaleratio=1, autorange="reversed"),
        margin=dict(l=10, r=10, t=10, b=10),
        height=700,
        images=layout_images,
    )
    return fig


# --- Layout ---
app.layout = html.Div(
    style={
        "backgroundColor": "#0a0a0f",
        "minHeight": "100vh",
        "color": "#e0e0e0",
        "fontFamily": "monospace",
        "padding": "40px",
    },
    children=[
        html.H1("NebulaCell", style={"color": "#00fff2", "textShadow": "0 0 10px #00fff2"}),
        html.P(id="dataset-label", style={"color": "#ff00e6"}),

        html.Div(
            [
                dcc.Input(
                    id="h5ad-path",
                    type="text",
                    placeholder="/path/to/your/data.h5ad",
                    style={
                        "width": "420px", "padding": "8px", "marginRight": "10px",
                        "backgroundColor": "#1a1a24", "color": "#e0e0e0", "border": "1px solid #00fff2",
                    },
                ),
                html.Button("Load dataset", id="load-btn", n_clicks=0),
            ],
            style={"marginBottom": "10px"},
        ),
        html.Div(id="load-feedback", style={"color": "#ff5555", "marginBottom": "20px"}),

        dcc.RadioItems(
            id="color-mode",
            options=[
                {"label": "  Leiden (expression clusters)", "value": "leiden"},
                {"label": "  HDBSCAN (spatial domains)", "value": "spatial_domain"},
                {"label": "  Gene expression", "value": "gene"},
            ],
            value="leiden",
            inline=True,
            style={"marginBottom": "15px"},
            labelStyle={"marginRight": "30px", "color": "#00fff2"},
        ),
        html.Div(
            [
                html.Label("Gene: ", style={"color": "#ff00e6", "marginRight": "10px"}),
                dcc.Dropdown(id="gene-select", style={"width": "300px", "color": "#0a0a0f"}),
            ],
            style={"display": "flex", "alignItems": "center", "marginBottom": "15px"},
        ),
        dcc.Checklist(
            id="show-image",
            options=[{"label": "  Show histology image", "value": "on"}],
            value=["on"],
            style={"marginBottom": "20px", "color": "#00fff2"},
        ),
        html.Div(
            [
                dcc.Input(
                    id="nl-query",
                    type="text",
                    placeholder="Try: show cells with high Xkr4 / show cluster 3",
                    style={
                        "width": "420px", "padding": "8px", "marginRight": "10px",
                        "backgroundColor": "#1a1a24", "color": "#e0e0e0", "border": "1px solid #00fff2",
                    },
                ),
                html.Button("Run query", id="nl-submit", n_clicks=0, style={"marginRight": "10px"}),
                html.Button("Reset", id="nl-reset", n_clicks=0),
            ],
            style={"marginBottom": "10px"},
        ),
        html.Div(id="nl-feedback", style={"color": "#ff5555", "marginBottom": "20px"}),
        dcc.Graph(id="main-scatter"),
    ],
)


@app.callback(
    Output("dataset-label", "children"),
    Output("gene-select", "options"),
    Output("gene-select", "value"),
    Output("load-feedback", "children"),
    Input("load-btn", "n_clicks"),
    State("h5ad-path", "value"),
    prevent_initial_call=False,
)
def load_dataset(n_clicks, path_value):
    if n_clicks == 0:
        return (
            f"Spatial transcriptomics visualization engine — {state.dataset_label}",
            [{"label": g, "value": g} for g in state.gene_options],
            state.default_gene,
            "",
        )

    if not path_value:
        return dash.no_update, dash.no_update, dash.no_update, "Enter a path to a .h5ad file first."

    try:
        adata = load_h5ad(path_value)
        if getattr(adata, "isbacked", False):
            label = f"{path_value} ({adata.n_obs} cells, backed mode — large dataset)"
        else:
            label = f"{path_value} ({adata.n_obs} cells)"
        state.load(adata, label)
        return (
            f"Spatial transcriptomics visualization engine — {state.dataset_label}",
            [{"label": g, "value": g} for g in state.gene_options],
            state.default_gene,
            "",
        )
    except FileNotFoundError:
        return dash.no_update, dash.no_update, dash.no_update, f"File not found: {path_value}"
    except ValueError as e:
        return dash.no_update, dash.no_update, dash.no_update, str(e)
    except Exception as e:
        return dash.no_update, dash.no_update, dash.no_update, f"Failed to load dataset: {e}"


@app.callback(
    Output("main-scatter", "figure"),
    Output("nl-feedback", "children"),
    Output("nl-query", "value"),
    Input("color-mode", "value"),
    Input("gene-select", "value"),
    Input("nl-submit", "n_clicks"),
    Input("nl-reset", "n_clicks"),
    Input("load-btn", "n_clicks"),
    Input("show-image", "value"),
    State("nl-query", "value"),
)
def update_figure(color_mode, selected_gene, run_clicks, reset_clicks, load_clicks, show_image_val, query_text):
    triggered_id = dash.callback_context.triggered_id
    show_image = "on" in (show_image_val or [])

    if triggered_id in ("nl-reset", "load-btn"):
        mask = np.ones(state.coords.shape[0], dtype=bool)
        fig = build_figure(color_mode, selected_gene, mask, show_image)
        return fig, "", "" if triggered_id == "nl-reset" else dash.no_update

    filter_dict = parse_query(query_text) if query_text else None
    is_valid, error_msg = validate_filter(filter_dict, state.adata) if query_text else (True, "")

    mask = compute_filter_mask(filter_dict) if is_valid else np.ones(state.coords.shape[0], dtype=bool)
    fig = build_figure(color_mode, selected_gene, mask, show_image)

    feedback = "" if is_valid else error_msg
    return fig, feedback, dash.no_update


if __name__ == "__main__":
    app.run(debug=True)
