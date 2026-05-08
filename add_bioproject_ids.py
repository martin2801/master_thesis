"""
add_bioproject_ids.py
---------------------
Reads a sample metadata CSV containing a BioSample accession column
(BioSampleAccession), fetches the linked BioProject ID for each sample
from the NCBI Entrez API, and writes the result to a new CSV with an
extra column: BioprojectAccession.

Requirements:
    pip install biopython pandas

Usage:
    python add_bioproject_ids.py \
        --input  infantis_clean_genome_collection.csv \
        --output infantis_clean_genome_collection_with_bioproject.csv \
        --email  your@email.com          # required by NCBI
        [--api-key YOUR_NCBI_API_KEY]    # optional but raises rate limit 3x->10 req/s
        [--batch-size 100]               # accessions fetched per API call (default 100)
        [--delay 0.4]                    # seconds between batches (default 0.4)
"""

import argparse
import time
import sys
import pandas as pd
from Bio import Entrez
from xml.etree import ElementTree as ET


# ---------------------------------------------------------------------------
# Core lookup
# ---------------------------------------------------------------------------

def fetch_bioproject_for_biosamples(biosample_ids: list[str]) -> dict[str, str]:
    """
    Given a list of BioSample accession strings, returns a dict mapping
    each accession -> BioProject accession (e.g. 'PRJNA123456', or '' if not found).

    Strategy: fetch the BioSample XML records directly. Each BioSample XML
    embeds its linked BioProject accession in a <Link> element — no elink needed.

    Example XML structure:
        <BioSample accession="SAMN45870197" ...>
          ...
          <Links>
            <Link type="entrez" target="bioproject" label="PRJNA123456">123456</Link>
          </Links>
        </BioSample>
    """
    result = {acc: "" for acc in biosample_ids}

    # Step 1: accessions -> UIDs via esearch
    term = " OR ".join(f"{acc}[Accession]" for acc in biosample_ids)
    search_handle = Entrez.esearch(db="biosample", term=term, retmax=len(biosample_ids))
    search_record = Entrez.read(search_handle)
    search_handle.close()

    uid_list = search_record["IdList"]
    if not uid_list:
        return result

    # Step 2: fetch full BioSample XML for all UIDs in one call
    fetch_handle = Entrez.efetch(
        db="biosample",
        id=",".join(uid_list),
        rettype="xml",
        retmode="xml",
    )
    raw_xml = fetch_handle.read()
    fetch_handle.close()

    # Step 3: parse XML - extract accession and BioProject link per sample
    root = ET.fromstring(raw_xml)
    for sample in root.findall(".//BioSample"):
        # Get the SAMN... accession for this sample
        bs_acc = sample.attrib.get("accession", "")
        if not bs_acc or bs_acc not in result:
            continue

        # BioProject is in <Links><Link target="bioproject" label="PRJNAXXXXXX">
        for link in sample.findall(".//Links/Link[@target='bioproject']"):
            label = link.attrib.get("label", "")
            if label.startswith("PRJ"):
                result[bs_acc] = label
                break  # take first BioProject link

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Add a BioprojectAccession column to a BioSample metadata CSV."
    )
    parser.add_argument("--input",         required=True,  help="Path to input CSV")
    parser.add_argument("--output",        required=True,  help="Path for output CSV")
    parser.add_argument("--email",         required=True,  help="Your e-mail (required by NCBI)")
    parser.add_argument("--api-key",       default=None,   help="NCBI API key (optional)")
    parser.add_argument("--batch-size",    type=int, default=100,
                        help="Accessions per API batch (default 100)")
    parser.add_argument("--delay",         type=float, default=0.4,
                        help="Seconds to wait between batches (default 0.4)")
    parser.add_argument("--biosample-col", default="BioSampleAccession",
                        help="Name of the BioSample accession column (default: BioSampleAccession)")
    args = parser.parse_args()

    # Configure Entrez
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # Load CSV
    print(f"Loading {args.input} ...")
    df = pd.read_csv(args.input)

    if args.biosample_col not in df.columns:
        sys.exit(
            f"ERROR: Column '{args.biosample_col}' not found.\n"
            f"Available columns: {list(df.columns)}"
        )

    accessions = df[args.biosample_col].dropna().unique().tolist()
    print(f"Found {len(accessions)} unique BioSample accessions.")

    # Fetch in batches
    bioproject_map: dict[str, str] = {}
    batches = [
        accessions[i : i + args.batch_size]
        for i in range(0, len(accessions), args.batch_size)
    ]

    for idx, batch in enumerate(batches, 1):
        print(f"  Batch {idx}/{len(batches)} ({len(batch)} accessions) ...", end=" ", flush=True)
        try:
            mapping = fetch_bioproject_for_biosamples(batch)
            bioproject_map.update(mapping)
            found = sum(1 for v in mapping.values() if v)
            print(f"resolved {found}/{len(batch)}")
        except Exception as exc:
            print(f"WARNING - batch failed: {exc}. Retrying one-by-one ...")
            for acc in batch:
                try:
                    single = fetch_bioproject_for_biosamples([acc])
                    bioproject_map.update(single)
                    time.sleep(0.2)
                except Exception as exc2:
                    print(f"    Skipping {acc}: {exc2}")
                    bioproject_map.setdefault(acc, "")

        if idx < len(batches):
            time.sleep(args.delay)

    # Add new column and save
    df["BioprojectAccession"] = df[args.biosample_col].map(bioproject_map).fillna("")

    resolved = (df["BioprojectAccession"] != "").sum()
    df.to_csv(args.output, index=False)
    print(f"\nDone. {resolved}/{len(df)} rows have a BioProject accession.")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()