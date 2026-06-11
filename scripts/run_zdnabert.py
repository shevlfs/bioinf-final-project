import argparse, os, sys
import numpy as np

def seq2kmer(seq, k=6):
    return [seq[x:x + k] for x in range(len(seq) + 1 - k)]

def split_seq(seq, length=512, pad=16):
    res = []
    for st in range(0, len(seq), length - pad):
        end = min(st + 512, len(seq))
        res.append((st, seq[st:end]))
        if end == len(seq):
            break
    return res

def stitch_np_seq(np_seqs, pad=16):
    res = np.array([], dtype=np.float32)
    for seq in np_seqs:
        res = res[:-pad] if res.size else res
        res = np.concatenate([res, seq])
    return res

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="models/zdnabert/6-new-12w-0")
    ap.add_argument("--genome", default="data/genome.fna")
    ap.add_argument("--out", default="results/zdnabert.bed")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--fp16", action="store_true", help="half precision (faster on MPS/GPU)")
    ap.add_argument("--device", default=None, help="mps|cuda|cpu (auto if unset)")
    ap.add_argument("--scaffold", default=None)
    a = ap.parse_args()

    import torch
    from transformers import BertTokenizer, BertForTokenClassification
    from scipy import ndimage

    if a.device:
        device = a.device
    elif torch.backends.mps.is_available():
        device = "mps"
    elif torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"
    sys.stderr.write(f"device: {device}\n")

    tok = BertTokenizer.from_pretrained(a.model_dir)
    model = BertForTokenClassification.from_pretrained(a.model_dir).to(device).eval()
    if a.fp16:
        model = model.half()
    vocab = tok.get_vocab()
    unk = tok.unk_token_id
    pad_id = tok.pad_token_id if tok.pad_token_id is not None else 0

    only = set(a.scaffold.split(",")) if a.scaffold else None
    n_regions = 0
    fout = open(a.out, "w")
    for name, seq in read_fasta(a.genome):
        if only and name not in only:
            continue
        seq = seq.upper()
        kmers = seq2kmer(seq, 6)
        windows = split_seq(kmers, 512, 16)
        win_ids = [[vocab.get(km, unk) for km in w] for _, w in windows]
        preds = [None] * len(win_ids)
        with torch.no_grad():
            for b0 in range(0, len(win_ids), a.batch):
                batch = win_ids[b0:b0 + a.batch]
                lens = [len(x) for x in batch]
                mx = max(lens)
                ids = torch.full((len(batch), mx), pad_id, dtype=torch.long)
                mask = torch.zeros((len(batch), mx), dtype=torch.long)
                for i, x in enumerate(batch):
                    ids[i, :len(x)] = torch.tensor(x, dtype=torch.long)
                    mask[i, :len(x)] = 1
                logits = model(ids.to(device), attention_mask=mask.to(device)).logits
                p1 = torch.softmax(logits, dim=-1)[:, :, 1].float().cpu().numpy()
                for i, L in enumerate(lens):
                    preds[b0 + i] = p1[i, :L]
                sys.stderr.write(f"  {name}: {min(b0+a.batch,len(win_ids))}/{len(win_ids)} windows\r")
        stitched = stitch_np_seq(preds, pad=16)
        sys.stderr.write("\n")
        labeled, maxlab = ndimage.label(stitched > a.threshold)
        for lab in range(1, maxlab + 1):
            idx = np.where(labeled == lab)[0]
            if idx.shape[0] > a.min_len:
                start = int(idx[0])
                end = int(idx[-1]) + 6
                score = float(stitched[idx].max())
                fout.write(f"{name}\t{start}\t{end}\tZDNABERT\t{score:.3f}\t.\n")
                n_regions += 1
    fout.close()
    print(f"Z-DNABERT regions (p>{a.threshold}, len>{a.min_len}): {n_regions} -> {a.out}")

if __name__ == "__main__":
    main()
