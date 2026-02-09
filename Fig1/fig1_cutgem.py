#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cutgem (Fig.1 part-1): Cut GEM by QuPath ROI polygons (GeoJSON), after binning.

Input:
  - GEM table (whitespace-delimited; optionally .gz)
  - QuPath ROI GeoJSON (Polygon / MultiPolygon)

Output:
  - <outdir>/<prefix>.QuPath.cut.gem.gz    (tab-delimited gz)
  - <outdir>/<prefix>.QuPath.roi.geojson   (GeoJSON with properties.name set to 1..N)

Example:
  python cutgem.py -g input.gem.gz -l rois.geojson -b 50 -C x y -o sample1
"""

import argparse
import os
import json
import skimage as ski
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame


def _extract_exterior_rings(features):
    """Return list of exterior rings, each ring is a (N,2) ndarray."""
    rings = []
    for feature in features:
        geometry = feature.get("geometry", {})
        gtype = geometry.get("type", None)
        coords = geometry.get("coordinates", None)
        if gtype == "Polygon" and coords and len(coords) > 0:
            rings.append(np.array(coords[0]))
        elif gtype == "MultiPolygon" and coords:
            for poly in coords:
                if poly and len(poly) > 0:
                    rings.append(np.array(poly[0]))
    return rings


def main(args):
    print("[cutgem] script start")

    # Read GEM
    gem_data = pd.read_csv(
        args.gem,
        sep=r"\s+",
        header=0,
        compression="infer",
        comment="#",
        low_memory=True,
        na_filter=False,
        engine="c",
    )
    print(f"[cutgem] GEM loaded: rows={len(gem_data):,} cols={len(gem_data.columns)}")

    # Bin coordinates (same behavior as original)
    gem_data["bx"] = round(gem_data[args.columns[0]] / args.bin, 0).astype(int)
    gem_data["by"] = round(gem_data[args.columns[1]] / args.bin, 0).astype(int)
    gem_data["bcoor"] = gem_data["bx"].map(str) + "_" + gem_data["by"].map(str)

    # Load GeoJSON
    with open(args.label, "r", encoding="utf-8") as f:
        data = json.load(f)

    features = data.get("features", [])
    if len(features) == 0:
        raise ValueError("GeoJSON has no features. Please check your ROI file.")

    # Extract polygons and rasterize to label table
    rings = _extract_exterior_rings(features)
    if len(rings) == 0:
        raise ValueError("No Polygon/MultiPolygon geometry found in the ROI GeoJSON.")

    label_frames = []
    for i, ts in enumerate(rings, start=1):
        # Keep original convention: polygon(ts[:,0], ts[:,1])
        rr_cc = ski.draw.polygon(ts[:, 0], ts[:, 1])
        df1 = DataFrame({"x": rr_cc[0], "y": rr_cc[1]})
        df1["Label"] = i
        label_frames.append(df1)

    label = pd.concat(label_frames, ignore_index=True)
    label["coor"] = label["x"].map(str) + "_" + label["y"].map(str)

    # Update GeoJSON feature names (1..N) and write to a new file (do not overwrite input by default)
    for i in range(1, len(features) + 1):
        data["features"][i - 1].setdefault("properties", {})
        data["features"][i - 1]["properties"]["name"] = i

    # Filter + assign Label
    gem_data = gem_data[gem_data["bcoor"].isin(label["coor"])].copy()
    gem_data["Label"] = gem_data["bcoor"].map(dict(zip(label["coor"], label["Label"])))

    # Drop helper columns
    gem_data.drop(columns=["bcoor", "bx", "by"], inplace=True)

    # Outputs
    outdir = args.outdir
    if not os.path.exists(outdir):
        print(f"[cutgem] create {outdir}")
        os.makedirs(outdir, exist_ok=True)

    out_gem = os.path.join(outdir, f"{args.output}.QuPath.cut.gem.gz")
    out_geojson = os.path.join(outdir, f"{args.output}.QuPath.roi.geojson")

    gem_data.to_csv(out_gem, sep="\t", header=True, index=False, compression="gzip")

    with open(out_geojson, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)

    print(f"[cutgem] GEM saved: {out_gem} (rows={len(gem_data):,})")
    print(f"[cutgem] ROI saved: {out_geojson}")
    print("[cutgem] mission completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cut GEM by QuPath ROI GeoJSON after binning.")
    parser.add_argument("-g", "--gem", help="input the gem file path", type=str, required=True)
    parser.add_argument("-C", "--columns", help="-C x y", type=str, nargs="+", default=["x", "y"])
    parser.add_argument("-l", "--label", help="ROI GeoJSON filepath", type=str, required=True)
    parser.add_argument("-b", "--bin", help="bin size (default 50)", type=int, default=50)
    parser.add_argument("-o", "--output", help="output prefix (default output)", type=str, default="output", nargs="?")
    parser.add_argument("--outdir", help="output directory (default gem_cut)", type=str, default="gem_cut")
    args = parser.parse_args()
    main(args)
