import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

import index


def sample_trajectories(max_trajectories=40, max_steps=300):
    """Collect a manageable subset of trajectories for plotting."""
    trajectories = []
    for R in index.rings:
        for peri in index.periaps:
            for theta in index.thetas:
                r0, v0 = index.inbound_state(R, theta, peri, index.epoch0)
                steps = min(int((index.T_MAX / index.DT).decompose()), max_steps)
                path = []
                for _, r_vec in index.advance_two_body_sun(r0, v0, index.epoch0, steps):
                    path.append(r_vec.to(u.AU).value)
                trajectories.append(np.array(path))
                if len(trajectories) >= max_trajectories:
                    return trajectories
    return trajectories


def plot_trajectories(output_path: Path):
    trajectories = sample_trajectories()

    fig, ax = plt.subplots(figsize=(6, 6))
    for path in trajectories:
        ax.plot(path[:, 0], path[:, 1], color="tab:blue", alpha=0.35, linewidth=0.8)

    # Planet positions at epoch0
    for body in index.BODIES:
        r_planet = index.planet_r(body, index.epoch0).value
        ax.scatter(r_planet[0], r_planet[1], s=30, label=body.name)

    # Gate circle for reference
    gate_radius = index.R_OUTER_GATE.to_value(u.AU)
    gate = plt.Circle((0.0, 0.0), gate_radius, color="black", linestyle="--", fill=False, alpha=0.3)
    ax.add_artist(gate)

    belt_radius = index.R_BELT_GATE.to_value(u.AU)
    belt = plt.Circle((0.0, 0.0), belt_radius, color="tab:orange", linestyle="-.", fill=False, alpha=0.5)
    ax.add_artist(belt)

    ax.set_xlabel("x (AU)")
    ax.set_ylabel("y (AU)")
    ax.set_title("Inbound Trajectory Samples (Sun-Centered)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-gate_radius * 1.1, gate_radius * 1.1)
    ax.set_ylim(-gate_radius * 1.1, gate_radius * 1.1)
    ax.legend(loc="upper right", fontsize="small", ncol=2)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    output = Path("plots") / "trajectory_overview.png"
    plot_trajectories(output)
    print(f"Wrote {output.as_posix()}")
