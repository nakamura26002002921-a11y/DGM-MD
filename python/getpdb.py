import sys
import requests


def download_pdb(pdbid: str, out_path: str):
    pdbid = pdbid.upper()
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    print(f"Downloading {pdbid} from:")
    print(url)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(r.content)
    print(f"Saved to: {out_path}")


def main():
    if len(sys.argv) != 3:
        print("Usage:")
        print("  python3 getpdb.py PDBID outputfile")
        print("")
        print("Example:")
        print("  python3 getpdb.py 1ALC 1ALC.pdb")
        sys.exit(1)
    pdbid = sys.argv[1]
    out_path = sys.argv[2]
    try:
        download_pdb(pdbid, out_path)
    except requests.HTTPError as e:
        print(f"HTTP error: {e}")
        sys.exit(1)
    except requests.RequestException as e:
        print(f"Request failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
