# Fig4 updated files for the new GC workflow

This bundle updates the files that are directly affected by the new GC processing strategy.

## Updated points

- `fig4_gc_processing.R`
  - adds pre-subset relabeling for `seurat_clusters == 9 -> GC_Column`
  - subsets `GC_Antral_Mural / GC_Preantral / GC_Atretic / GC_Column`
  - keeps the user's newer multi-stage GC workflow
  - writes final refined labels into `meta.data$celltype`

- `configs/fig4_combined.yaml`
  - updates GC config to the new workflow
  - updates ssGSEA order to:
    - `GC_Progenitor`
    - `GC_Preantral`
    - `GC_Mitotic`
    - `GC_Cumulus`
    - `GC_Antral_Mural`
    - `GC_Atretic`
  - disables `c2l_sc` cluster-based remapping for `obj_gr`
  - changes Monocle3 rooting from a hard-coded cluster to `celltype == GC_Progenitor`
  - keeps the notebook-aligned `c2l_spatial` two-input setting

- `fig4_monocle3.R`
  - now supports rooting by a metadata column/value pair (`root_by_col`, `root_value`)

## Recommended replacement

Copy these files into your project:

```bash
cp fig4_gc_update_bundle/Fig4/fig4_gc_processing.R Stereo-cell-oocyte/Fig4/
cp fig4_gc_update_bundle/Fig4/fig4_monocle3.R Stereo-cell-oocyte/Fig4/
cp fig4_gc_update_bundle/Fig4/configs/fig4_combined.yaml Stereo-cell-oocyte/Fig4/configs/
cp fig4_gc_update_bundle/Fig4/README.md Stereo-cell-oocyte/Fig4/
```

## Notes

Because the new final GC object preserves `celltype` labels after the post-label exclusion/reclustering step, downstream modules should use the saved `celltype` column directly rather than trying to re-infer labels from the final `seurat_clusters`.
