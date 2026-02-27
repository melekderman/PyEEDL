import os, glob
from electron import create_mcdc_file

def main():
    base_dir = os.path.dirname(__file__)
    raw_dir  = os.path.join(base_dir, "../raw_data")
    out_dir  = os.path.join(base_dir, "../mcdc_data")
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))
    for p in files:
        create_mcdc_file(p, out_dir)

if __name__ == "__main__":
    main()