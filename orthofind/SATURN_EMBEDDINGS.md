# SATURN ESM2 `.pt` embedding notebook

请在 Jupyter 中按顺序运行 [`Generate_SATURN_ESM2_Embeddings.ipynb`](Generate_SATURN_ESM2_Embeddings.ipynb)。该 notebook 使用已有的 `./models/esm2_t36_3B_UR50D.pt`，取 **ESM2 第 33 层**的残基 mean pooling，并为九个物种生成 SATURN 可读取的 PyTorch `.pt` 文件。

- 输入：`orthofind/proteomes/{cro,emu,esc,gga,hsa,mmu,psi,pvi,tse}.fa`。
- 最终输出：`orthofind/proteomes/saturn_esm2_t36_3B_layer33/<species>/<species>_gene_embeddings_ESM2_t36_3B_layer33.pt`。
- `.pt` 内容为 `{raw_gene_id: CPU float32 Tensor}`，每个 tensor 的 shape 为 `(2560,)`。
- 本流程只创建或删除 `saturn_esm2_t36_3B_layer33/` 内的文件，绝不会写入或覆盖现有 TranscriptFormer 的 `*_gene.h5`、`*.esm_input.fa`、`*.gene_id_mapping.tsv`。

默认每个 gene ID 保留原 FASTA 中第一条非空 protein，因而和已有 TranscriptFormer 预处理规则一致。如果 FASTA 的一条记录就是一个 gene，这正是所需行为；若需要对 isoform 按 gene symbol 聚合，应先以 gene-level FASTA 作为输入，或修改 notebook 的 gene-ID 提取规则。

运行前确保 notebook 环境已经安装 `torch`、`fair-esm`、`biopython`、`numpy`、`tqdm` 和（如需 GPU）可用 CUDA PyTorch。Cell 2 中 `OVERWRITE_EXISTING` 默认是 `False`，已有有效 SATURN `.pt` 时会跳过；它不影响任何 TranscriptFormer 文件。
