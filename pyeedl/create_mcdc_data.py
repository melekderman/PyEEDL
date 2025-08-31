import os, glob
import numpy as np
import h5py


def linear_interpolation(xs_energy_grid, energy_ref, xs_ref):
    xs_new = np.interp(energy_ref, xs_energy_grid, xs_ref, left=0.0).astype("f8")
    if np.any (xs_energy_grid > energy_ref[-1]):
        print("Warning: xs_energy_grid has values larger than energy grid max.")
    return xs_new


base_dir = os.path.dirname(__file__)
raw_dir = os.path.join(base_dir, "raw_data")
files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))

in_path = files[0]
print(f"Reading: {in_path}")

with h5py.File(in_path, "r") as src:
    meta = src.get("metadata")
    Z = int(meta.attrs.get("Z", 0))
    AWR = float(meta.attrs.get("AWR", np.nan))
    Sym = meta.attrs.get("Sym", "")

    # total xs
    total_energy = src.get("/total_xs/cross_section/energy")
    total_xs = src.get("/total_xs/cross_section/cross_section")

    # scattering
    total_scat_xs_energy = src.get("/elastic_scatter/cross_section/total/energy")
    total_scat_xs = src.get("/elastic_scatter/cross_section/total/cross_section")
    la_scat_xs_energy = src.get("/elastic_scatter/cross_section/large_angle/energy")
    la_scat_xs = src.get("/elastic_scatter/cross_section/large_angle/cross_section")

    la_scat_dist_mu = src.get("/elastic_scatter/distributions/large_angle/mu")
    la_scat_dist_energy = src.get("/elastic_scatter/distributions/large_angle/inc_energy")
    la_scat_dist_prob = src.get("/elastic_scatter/distributions/large_angle/probability")

    # bremsstrahlung
    brem_xs_energy = src.get("/bremsstrahlung/cross_section/energy")
    brem_xs = src.get("/bremsstrahlung/cross_section/cross_section")

    brem_avg_loss_energy = src.get("/bremsstrahlung/distributions/loss_inc_energy")
    brem_avg_loss = src.get("/bremsstrahlung/distributions/avg_loss")

    # excitation
    exc_xs_energy = src.get("/excitation/cross_section/energy")
    exc_xs = src.get("/excitation/cross_section/cross_section")

    exc_avg_loss_energy = src.get("/excitation/distributions/loss_inc_energy")
    exc_avg_loss = src.get("/excitation/distributions/avg_loss")

    # ionization
    ion_total_xs_energy = src.get("/ionization/cross_section/energy")
    ion_total_xs = src.get("/ionization/cross_section/cross_section")

    subshell_names = []
    xs_root = src.get("/ionization/cross_section")
    if isinstance(xs_root, h5py.Group):
        for name in xs_root.keys():
            if name == "total":
                continue
            subshell_names.append(name)

    subshell_xs_data = {}
    for name in subshell_names:
        xs = src.get(f"/ionization/cross_section/{name}/cross_section")
        xs_energy = src.get(f"/ionization/cross_section/{name}/energy")
        subshell_xs_data[name] = {
            f"xs_{name}": xs,
            f"xs_energy_{name}": xs_energy
        }

    subshell_dist_data = {}
    for name in subshell_names:
        inc_energy = src.get(f"/ionization/cross_section/{name}/distributions/inc_energy")
        out_energy = src.get(f"/ionization/cross_section/{name}/distributions/out_energy")
        prob = src.get(f"/ionization/cross_section/{name}/distributions/b")
        be = src.get(f"/ionization/cross_section/{name}/binding_energy")
        subshell_dist_data[name] = {
            f"dist_inc_energy_{name}": inc_energy,
            f"dist_out_energy_{name}": out_energy,
            f"dist_prob_{name}": prob,
            f"binding_energy_{name}": be
        }

out_dir = os.path.join(base_dir, "mcdc_data")
os.makedirs(out_dir, exist_ok=True)
sym_str = str(Sym)
out_path = os.path.join(out_dir, f"{sym_str}.h5")

xs_energy_grid = total_scat_xs_energy
xs_sc_total = linear_interpolation(xs_energy_grid, total_scat_xs_energy, total_scat_xs)
xs_sc_la = linear_interpolation(xs_energy_grid, la_scat_xs_energy, la_scat_xs)
xs_brem = linear_interpolation(xs_energy_grid, brem_xs_energy, brem_xs)
xs_exc = linear_interpolation(xs_energy_grid, exc_xs_energy, exc_xs)
xs_ion_total = linear_interpolation(xs_energy_grid, ion_total_xs_energy, ion_total_xs)


with h5py.File(out_path, "w") as h5f:
    h5f.create_dataset("atomic_weight_ratio", data=AWR)
    h5f.create_dataset("atomic_number", data=Z)
    h5f.create_group("electron_reactions")
    h5f.create_dataset("electron_reactions/xs_energy_grid", data=xs_energy_grid)
    h5f.create_group("electron_reactions/total")
    h5f.create_dataset("electron_reactions/total/xs", data=total_xs)
