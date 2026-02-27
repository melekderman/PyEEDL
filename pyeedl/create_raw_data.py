import os
from electron import save_element_h5

def main():
    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../eedl")
    out_dir = os.path.join(base_dir, "../raw_data")
    os.makedirs(out_dir, exist_ok=True)

    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EEDL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_element_h5(fn, out_dir)

if __name__ == "__main__":
    main()