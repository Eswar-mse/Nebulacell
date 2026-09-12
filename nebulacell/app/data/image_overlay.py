"""
Extracts the histology image and scale factor from a Visium-style
AnnData object (adata.uns['spatial']), for overlaying under the
spatial scatter plot.
"""

import base64
from io import BytesIO

import numpy as np
import anndata as ad
from PIL import Image


def get_histology_image(adata: ad.AnnData, resolution: str = "hires"):
    """
    Extract the histology image + matching scale factor from a Visium
    AnnData object.

    Returns (image_array, scale_factor, sample_key) or (None, None, None)
    if no spatial image metadata is present (e.g. non-Visium datasets).

    resolution: "hires" or "lowres" — hires is sharper but larger;
    lowres loads faster for big dashboards.
    """
    if "spatial" not in adata.uns:
        return None, None, None

    sample_keys = list(adata.uns["spatial"].keys())
    if not sample_keys:
        return None, None, None

    sample_key = sample_keys[0]  # single-sample datasets; multi-sample is a later concern
    spatial_data = adata.uns["spatial"][sample_key]

    images = spatial_data.get("images", {})
    scalefactors = spatial_data.get("scalefactors", {})

    if resolution not in images:
        # Fall back to whatever resolution IS available
        available = list(images.keys())
        if not available:
            return None, None, None
        resolution = available[0]

    image_array = images[resolution]
    scale_key = f"tissue_{resolution}_scalef"
    scale_factor = scalefactors.get(scale_key)

    if scale_factor is None:
        return None, None, None

    return image_array, scale_factor, sample_key


def image_array_to_data_uri(image_array: np.ndarray) -> str:
    """Convert a numpy image array (0-1 float or 0-255 uint8) to a base64 PNG data URI for Plotly."""
    arr = image_array
    if arr.dtype != np.uint8:
        arr = (arr * 255).clip(0, 255).astype(np.uint8)

    img = Image.fromarray(arr)
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"
