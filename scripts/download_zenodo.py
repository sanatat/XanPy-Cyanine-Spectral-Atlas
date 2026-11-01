#!/usr/bin/env python3
"""Download the original archived compound ZIPs by DOI record and verify MD5.

Archive files are NOT included in the GitHub repository. Run this only when
reproduction from the original source is desired. Requires an internet connection.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]

def md5(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dest", type=Path, default=ROOT/"data/raw")
    p.add_argument("--verify-only", action="store_true",help="No downloads; only verify existing archives")
    args=p.parse_args()
    args.dest.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((ROOT/"data/source_manifest.json").read_text(encoding="utf-8"))
    for item in manifest["zip_archives"]:
        dest=args.dest/item["filename"]
        if dest.exists() and md5(dest)==item["md5"]:
            print(f"PASS {dest.name}: original archive MD5 verified")
            continue
        if args.verify_only:
            print(f"MISSING_OR_BAD_HASH {dest.name}")
            raise SystemExit(1)
        if dest.exists():
            raise RuntimeError(f"File exists but MD5 differs: {dest}. Remove it manually to redownload.")
        tmp=dest.with_suffix('.partial')
        req=Request(item["url"],headers={"User-Agent":"XanPy-reproducibility/1.0"})
        print(f"Downloading {item['filename']} from Zenodo...")
        try:
            with urlopen(req,timeout=100) as response, tmp.open('wb') as output:
                for block in iter(lambda: response.read(1024*1024), b''):
                    output.write(block)
            check=md5(tmp)
            if check!=item['md5']:
                raise ValueError(f"MD5 mismatch for {item['filename']}: got {check}")
            tmp.replace(dest)
            print(f"PASS {dest.name}: original archive MD5 verified")
        finally:
            tmp.unlink(missing_ok=True)
    print('All four archive checksums verified.')

if __name__=='__main__': main()
