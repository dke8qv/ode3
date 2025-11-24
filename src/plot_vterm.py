import numpy as np
import matplotlib.pyplot as plt


M = 10.0     # kg
g = 9.81     # m/s^2
k_drag = 0.1 # same k used in vterm with air
v_terminal = np.sqrt(M * g / k_drag)

# ------------------------------------------------------------
# Helper function to load 6-column text files
# ------------------------------------------------------------
def load_txt(fname):
    data = np.loadtxt(fname)
    # Columns = t, x, y, vx, vy, last
    return {
        "t":  data[:,0],
        "x":  data[:,1],
        "y":  data[:,2],
        "vx": data[:,3],
        "vy": data[:,4],
        "last": data[:,5],
    }


noair = load_txt("vterm_noair.txt")
air   = load_txt("vterm_air.txt")


def compute_energies(d):
    vx, vy, y = d["vx"], d["vy"], d["y"]
    KE = 0.5 * M * (vx**2 + vy**2)
    U  = M * g * y
    Et = KE + U
    return KE, Et

KE_noair, Et_noair = compute_energies(noair)
KE_air,   Et_air   = compute_energies(air)


def plot_energy(t, KE, Et, fname, title):
    plt.figure(figsize=(8,5))
    plt.plot(t, KE, label="Kinetic Energy K(t)")
    plt.plot(t, Et, label="Total Energy E(t)")
    plt.xlabel("Time t [s]")
    plt.ylabel("Energy [J]")
    plt.title(title)
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    plt.close()
    print("Saved", fname)

def plot_trajectory(x, y, fname, title):
    plt.figure(figsize=(8,5))
    plt.plot(x, y)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(title)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    plt.close()
    print("Saved", fname)

def plot_speed(t, vx, vy, fname, title, v_terminal=None):
    v = np.sqrt(vx**2 + vy**2)
    plt.figure(figsize=(8,5))

    # Numerical speed curve
    plt.plot(t, v, label="Speed v(t)")

    # Terminal velocity horizontal line
    if v_terminal is not None:
        plt.axhline(v_terminal, color='red', linestyle='--',
                    label=f"Terminal Velocity = {v_terminal:.2f} m/s")

    plt.xlabel("Time [s]")
    plt.ylabel("Speed [m/s]")
    plt.title(title)
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    plt.close()
    print("Saved", fname)



# Vacuum case
plot_energy(noair["t"], KE_noair, Et_noair,
            "KE_total_noair.png",
            "Vacuum Case: Kinetic and Total Energy vs Time")

plot_trajectory(noair["x"], noair["y"],
                "trajectory_noair.png",
                "Vacuum Trajectory")

# Air case
plot_energy(air["t"], KE_air, Et_air,
            "KE_total_air.png",
            "Air Drag: Kinetic and Total Energy vs Time")

plot_trajectory(air["x"], air["y"],
                "trajectory_air.png",
                "Trajectory With Air Drag")

plot_speed(air["t"], air["vx"], air["vy"],
           "speed_air.png",
           "Speed vs Time With Air Drag",
           v_terminal=v_terminal)

print("All plots complete.")

