import os, glob
import numpy as np
import h5py


from function import linear_interpolation, build_pdf

base_dir = os.getcwd()
#base_dir = os.path.dirname(__file__)
raw_dir = os.path.join(base_dir, "../raw_data")
files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))

BARN_TO_CM2 = 1e-24

in_path = files[0]
print(f"Reading: {in_path}")

with h5py.File(in_path, "r") as src:
    Z = src["metadata/Z"]
    AWR = src["metadata/AWR"]
    Sym = src["metadata/Sym"].asstr()[()]

    # total xs
    total_energy = src["/total_xs/cross_section/energy"][:]
    total_xs = src["/total_xs/cross_section/cross_section"][:]


    # scattering
    total_scat_xs_energy = src["/elastic_scatter/cross_section/total/energy"][:]
    total_scat_xs = src["/elastic_scatter/cross_section/total/cross_section"][:]
    la_scat_xs_energy = src["/elastic_scatter/cross_section/large_angle/energy"][:]
    la_scat_xs = src["/elastic_scatter/cross_section/large_angle/cross_section"][:]

    la_scat_dist_mu = src["/elastic_scatter/distributions/large_angle/mu"][:]
    la_scat_dist_energy = src["/elastic_scatter/distributions/large_angle/inc_energy"][:]
    la_scat_dist_prob = src["/elastic_scatter/distributions/large_angle/probability"][:]

    # bremsstrahlung
    brem_xs_energy = src["/bremsstrahlung/cross_section/energy"][:]
    brem_xs = src["/bremsstrahlung/cross_section/cross_section"][:]

    brem_avg_loss_energy = src["/bremsstrahlung/distributions/loss_inc_energy"][:]
    brem_avg_loss = src["/bremsstrahlung/distributions/avg_loss"][:]

    # excitation
    exc_xs_energy = src["/excitation/cross_section/energy"][:]
    exc_xs = src["/excitation/cross_section/cross_section"][:]

    exc_avg_loss_energy = src["/excitation/distributions/loss_inc_energy"][:]
    exc_avg_loss = src["/excitation/distributions/avg_loss"][:]

    # ionization
    ion_total_xs_energy = src["/ionization/cross_section/total/energy"][:]
    ion_total_xs = src["/ionization/cross_section/total/cross_section"][:]

    subshell_names = []
    xs_root = src["/ionization/cross_section"]
    if isinstance(xs_root, h5py.Group):
        for name in xs_root.keys():
            if name == "total":
                continue
            subshell_names.append(name)

    subshell_xs_data = {}
    for name in subshell_names:
        xs = src[f"/ionization/cross_section/{name}/cross_section"][:]
        xs_energy = src[f"/ionization/cross_section/{name}/energy"][:]
        subshell_xs_data[name] = {
            f"xs_{name}": xs,
            f"xs_energy_{name}": xs_energy
        }

    subshell_dist_data = {}
    for name in subshell_names:
        inc_energy = src[f"/ionization/distributions/{name}/inc_energy"][:]
        out_energy = src[f"/ionization/distributions/{name}/out_energy"][:]
        prob = src[f"/ionization/distributions/{name}/b"][:]
        be = src[f"/ionization/distributions/{name}/binding_energy"]
        subshell_dist_data[name] = {
            f"dist_inc_energy_{name}": inc_energy,
            f"dist_out_energy_{name}": out_energy,
            f"dist_prob_{name}": prob,
            f"binding_energy_{name}": be
        }

    out_dir = os.path.join(base_dir, "../mcdc_data")
    os.makedirs(out_dir, exist_ok=True)
    sym_str = Sym
    out_path = os.path.join(out_dir, f"{sym_str}.h5")
    
    xs_energy_grid = total_scat_xs_energy
    xs_sc_total = linear_interpolation(xs_energy_grid, total_scat_xs_energy, total_scat_xs)
    xs_sc_la = linear_interpolation(xs_energy_grid, la_scat_xs_energy, la_scat_xs)
    # print(f" DEBUG: size brem energy: {brem_xs_energy.size}, size brem xs: {brem_xs.size}")
    xs_brem = linear_interpolation(xs_energy_grid, brem_xs_energy, brem_xs)
    xs_exc = linear_interpolation(xs_energy_grid, exc_xs_energy, exc_xs)
    xs_ion_total = linear_interpolation(xs_energy_grid, ion_total_xs_energy, ion_total_xs)
    
    xs_sc_sa = xs_sc_total - xs_sc_la

    # barn -> cm^2
    total_xs_cm2 = np.asarray(total_xs, dtype="f8") * BARN_TO_CM2
    xs_sc_total  = np.asarray(xs_sc_total, dtype="f8") * BARN_TO_CM2
    xs_sc_la     = np.asarray(xs_sc_la,    dtype="f8") * BARN_TO_CM2
    xs_sc_sa     = np.asarray(xs_sc_sa,    dtype="f8") * BARN_TO_CM2
    xs_brem      = np.asarray(xs_brem,     dtype="f8") * BARN_TO_CM2
    xs_exc       = np.asarray(xs_exc,      dtype="f8") * BARN_TO_CM2
    xs_ion_total = np.asarray(xs_ion_total,dtype="f8") * BARN_TO_CM2
    
    with h5py.File(out_path, "w") as h5f:
        h5f.create_dataset("atomic_weight_ratio", data=AWR)
        h5f.create_dataset("atomic_number", data=Z)
        h5f.create_group("electron_reactions")
        h5f.create_dataset("electron_reactions/xs_energy_grid", data=xs_energy_grid)
        # total
        h5f.create_group("electron_reactions/total")
        h5f.create_dataset("electron_reactions/total/xs", data=total_xs)
        # elastic scattering
        # total
        h5f.create_group("electron_reactions/elastic_scattering")
        h5f.create_dataset("electron_reactions/elastic_scattering/xs", data=xs_sc_total)
        # Large Angle
        h5f.create_group("electron_reactions/elastic_scattering/large_angle")
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/xs", data=xs_sc_la)
        h5f.create_group("electron_reactions/elastic_scattering/large_angle/scattering_cosine")
        energy_grid, energy_offset, val, PDF = build_pdf(la_scat_dist_energy, la_scat_dist_mu, la_scat_dist_prob)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/energy_grid", data=energy_grid)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/energy_offset", data=energy_offset)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/value", data=val)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/PDF", data=PDF)
        # Small Angle
        h5f.create_group("electron_reactions/elastic_scattering/small_angle")
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/xs", data=xs_sc_sa)
        # bremsstrahlung
        h5f.create_group("electron_reactions/bremsstrahlung")
        h5f.create_dataset("electron_reactions/bremsstrahlung/xs", data=xs_brem)
        h5f.create_group("electron_reactions/bremsstrahlung/energy_loss")
        h5f.create_dataset("electron_reactions/bremsstrahlung/energy_loss/energy", data=brem_avg_loss_energy)
        h5f.create_dataset("electron_reactions/bremsstrahlung/energy_loss/value", data=brem_avg_loss)
        # excitation
        h5f.create_group("electron_reactions/excitation")
        h5f.create_dataset("electron_reactions/excitation/xs", data=xs_exc)
        h5f.create_group("electron_reactions/excitation/energy_loss")
        h5f.create_dataset("electron_reactions/excitation/energy_loss/energy", data=exc_avg_loss_energy)
        h5f.create_dataset("electron_reactions/excitation/energy_loss/value", data=exc_avg_loss)
        # ionization (total + subshells)
        h5f.create_group("electron_reactions/ionization")
        h5f.create_dataset("electron_reactions/ionization/xs", data=xs_ion_total)
        subs = h5f.create_group("electron_reactions/ionization/subshells")
        for shell in subshell_names:
            g = subs.create_group(shell)
            xs_e = subshell_xs_data[shell][f"xs_energy_{shell}"]
            xs_v = subshell_xs_data[shell][f"xs_{shell}"]
            xs_on_grid = linear_interpolation(xs_energy_grid, xs_e, xs_v)
            g.create_dataset("xs", data=xs_on_grid)
            inc = subshell_dist_data[shell][f"dist_inc_energy_{shell}"]
            out = subshell_dist_data[shell][f"dist_out_energy_{shell}"]
            prob = subshell_dist_data[shell][f"dist_prob_{shell}"]
            eg, off, val, PDF = build_pdf(inc, out, prob)
            dg = g.create_group("product")
            dg.create_dataset("energy_grid", data=eg)
            dg.create_dataset("energy_offset", data=off)
            dg.create_dataset("value", data=val)
            dg.create_dataset("PDF", data=PDF)
            # binding energy
            be = subshell_dist_data[shell][f"binding_energy_{shell}"][()]
            g.create_dataset("binding_energy", data=be)
