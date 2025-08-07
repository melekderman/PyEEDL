import os
from electron import save_element_h5

def main():
    IN_DIR  = "../eedl"
    OUT_DIR = "../mcdc-electron"
    os.makedirs(OUT_DIR, exist_ok=True)

    for Z in range(1, 101):
        fn = os.path.join(IN_DIR, f"EEDL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_element_h5(fn, OUT_DIR)
        else:
            print(f"[!] File not found: {fn}")

if __name__ == "__main__":
    main()