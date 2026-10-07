"""REBOUND checks of the patched-conic model.

`check_flyby` sets up one encounter at a chosen planetocentric impact parameter and compares the heliocentric (E, h) after the flyby with the analytic prediction. `check_chain` places two planets so that a predicted i -> j chain happens, integrates it with real gravity, then repeats with planet i's mass set to zero (removal test) to see whether j is still met. `fly_plan` flies any planned sequence of flybys and returns the outbound asymptote, to check the patched-conic assist.
"""
import math

import rebound

from .constants import PLANETS
from .flyby import Crossing, local_velocity, state_after
from .scores import _t_from_peri, polar, true_anomaly

R_START = 60.0


def _crossing_time(E, h, r, sign, t_peri):
	"""Time the conic (perihelion at t_peri) reaches radius r moving `sign`."""
	return t_peri + _t_from_peri(E, h, r, sign)


def _polar_angle(E, h, r, sign, omega):
	"""Polar angle at radius r on the conic with perihelion direction omega."""
	return polar(h, omega, true_anomaly(E, h, r, sign))


def _helio_Eh(sim, k):
	s, p = sim.particles[0], sim.particles[k]
	x, y = p.x - s.x, p.y - s.y
	vx, vy = p.vx - s.vx, p.vy - s.vy
	r = math.hypot(x, y)
	return 0.5 * (vx * vx + vy * vy) - 1 / r, x * vy - y * vx, r, (x * vx + y * vy) / r


def _new_sim(phases, masses):
	sim = rebound.Simulation()
	sim.G = 1.0
	sim.integrator = "ias15"
	sim.add(m=1.0)
	for p, ph, m in zip(PLANETS, phases, masses):
		sim.add(m=m, a=p.a, e=0.0, f=ph, primary=sim.particles[0])
	sim.N_active = 1 + len(PLANETS)
	sim.testparticle_type = 0
	return sim


def _visitor(v_inf, b):
	"""Incoming in-plane hyperbola at R_START, perihelion along +x.

	Returns (E, h, t_peri) with t_peri the time from start to perihelion.
	"""
	E, h = 0.5 * v_inf ** 2, b * v_inf
	e = math.sqrt(1 + 2 * E * h * h)
	A = 1 / (2 * E)
	F0 = math.acosh((1 + R_START / A) / e)
	return E, h, (e * math.sinh(F0) - F0) * A ** 1.5


def _add_visitor(sim, E, h):
	e = math.sqrt(1 + 2 * E * h * h)
	f0 = -math.acos((h * h / R_START - 1) / e)
	# inc = pi flips the in-plane sense, so the polar angle is -f (retrograde)
	inc = 0.0 if h > 0 else math.pi
	# far out the visitor orbits the barycenter, not the wobbling Sun
	sim.add(primary=sim.com(), a=-1 / (2 * E), e=e, inc=inc, f=f0)


def _placement(E, h, sign, t_peri, omega, planet, b_p):
	"""Phase at t=0 putting `planet` at signed impact parameter b_p."""
	p = PLANETS[planet]
	cr = Crossing(planet, *local_velocity(E, h, p.a, sign))
	t_c = _crossing_time(E, h, p.a, sign, t_peri)
	theta = _polar_angle(E, h, p.a, sign, omega)
	s = b_p * cr.w_mag / cr.v_r          # along-track offset of the planet
	return theta + s / p.a - t_c * p.v_circ / p.a, t_c, cr


def _run(sim, t_end, watch, track=None):
	"""Integrate to t_end, recording each watched planet's encounter.

	For planet j: rmin, the measured planetocentric impact parameter b_p (h_rel / v_inf_rel, conserved along the planetocentric hyperbola, so the sample time need not be exact), and the visitor's heliocentric (E, h) just before and after it is inside 3 r_SOI. If `track` is a list, (t, visitor xy, planet xys) samples relative to the Sun are appended.
	"""
	k = len(sim.particles) - 1
	enc = {j: {"rmin": math.inf, "b_meas": math.nan, "before": None, "after": None}
		   for j in watch}
	inside = {j: False for j in watch}
	for j in watch:
		enc[j]["crossings"] = []
	r_prev = None
	dt = 0.02
	t = sim.t
	while t < t_end:
		ps = sim.particles
		if track is not None:
			track.append((t, (ps[k].x - ps[0].x, ps[k].y - ps[0].y),
						  [(ps[1 + n].x - ps[0].x, ps[1 + n].y - ps[0].y)
						   for n in range(len(PLANETS))]))
		E_now, h_now, r_now, vr_now = _helio_Eh(sim, k)
		Eh = (E_now, h_now)
		if r_prev is not None:
			for j in watch:
				a = PLANETS[j].a
				if (r_prev - a) * (r_now - a) < 0:
					vis = ps[k]
					theta = math.atan2(vis.y - ps[0].y, vis.x - ps[0].x)
					enc[j]["crossings"].append((t, theta, 1 if vr_now > 0 else -1, Eh))
		r_prev = r_now
		near = math.inf
		for j in watch:
			pj, v = ps[1 + j], ps[k]
			x, y = v.x - pj.x, v.y - pj.y
			d = math.hypot(x, y)
			r_soi = PLANETS[j].soi
			now_in = d < 3 * r_soi
			if now_in and not inside[j]:
				enc[j]["before"] = Eh
			if inside[j] and not now_in:
				enc[j]["after"] = Eh
			inside[j] = now_in
			if d < enc[j]["rmin"]:
				vx, vy = v.vx - pj.vx, v.vy - pj.vy
				vinf2 = vx * vx + vy * vy - 2 * pj.m / d
				enc[j]["rmin"] = d
				enc[j]["t_min"] = t
				enc[j]["b_meas"] = (x * vy - y * vx) / math.sqrt(max(vinf2, 1e-30))
			near = min(near, d / r_soi)
		dt = 0.0005 if near < 4 else 0.02
		t = min(t + dt, t_end)
		sim.integrate(t, exact_finish_time=1)
	for j in watch:
		if enc[j]["after"] is None and inside[j] is False and enc[j]["before"] is not None:
			enc[j]["after"] = _helio_Eh(sim, k)[:2]
	return enc


def _shoot(v_inf, b, targets, masses, t_end, iters=None, tol=1e-3, track=None):
	"""Adjust planet phases until each target encounter hits its b_p.

	targets: ordered list of dicts {planet, phase, b, cr, sign}. The analytic placement is only a first guess: the Sun's reflex motion and distant pulls shift the real encounter, and a later encounter is very sensitive to an earlier one. So each pass corrects only the earliest unconverged target: by re-aiming at the flown path's actual ring crossing when it is far, else by a secant step on the measured b_p(phase) once two tries exist, else linearly (along-track offset db * |w| / v_r).
	"""
	E, h, _ = _visitor(v_inf, b)
	phases = [0.0] * len(PLANETS)
	for tg in targets:
		phases[tg["planet"]] = tg["phase"]
	watch = [tg["planet"] for tg in targets]
	last = {}                      # planet -> (phase, b_meas) of previous try
	for _ in range(iters or 12 * len(targets)):
		sim = _new_sim(phases, masses)
		_add_visitor(sim, E, h)
		enc = _run(sim, t_end, watch)
		after = 0.0
		done = True
		for tg in targets:
			j, e = tg["planet"], enc[tg["planet"]]
			err = tg["b"] - e["b_meas"]
			if abs(err) > tol * abs(tg["b"]):
				done = False
				if e["rmin"] < 0.05 * PLANETS[j].a and tg["cr"] is not None:
					prev = last.get(j)
					now = (phases[j], e["b_meas"])
					slope = None
					if prev is not None and prev[0] != now[0]:
						slope = (now[1] - prev[1]) / (now[0] - prev[0])
					lin = tg["cr"].v_r * PLANETS[j].a / tg["cr"].w_mag   # db/dphase
					if slope is None or slope * lin <= 0 or not 0.2 < slope / lin < 5:
						slope = lin
					last[j] = now
					phases[j] += err / slope
				else:
					ph, cr = _aim_at_crossing(e, j, tg["sign"], after, tg["b"])
					if ph is None:
						raise RuntimeError(f"path never crosses {PLANETS[j].name}'s ring")
					phases[j], tg["cr"] = ph, cr
					last.pop(j, None)
				break
			after = e["t_min"]
		if done:
			break
	sim = _new_sim(phases, masses)
	_add_visitor(sim, E, h)
	return phases, _run(sim, t_end, watch, track=track)


def _replay(v_inf, b, phases, masses, t_end, watch, track=None):
	E, h, _ = _visitor(v_inf, b)
	sim = _new_sim(phases, masses)
	_add_visitor(sim, E, h)
	return _run(sim, t_end, watch, track=track)


def check_flyby(v_inf, b, planet, b_p, sign=-1):
	"""Compare analytic post-flyby (E, h) with REBOUND for one encounter.

	The prediction uses the N-body (E, h) on entry, so this isolates the patched-conic deflection itself from approach-phase perturbations.
	"""
	masses = [0.0] * len(PLANETS)
	masses[planet] = PLANETS[planet].mu
	E, h, t_peri = _visitor(v_inf, b)
	ph, t_c, cr = _placement(E, h, sign, t_peri, 0.0, planet, b_p)
	t_end = t_c + 10 * PLANETS[planet].soi / cr.w_mag + 0.5
	_, enc = _shoot(v_inf, b, [dict(planet=planet, phase=ph, b=b_p, cr=cr, sign=sign)],
					masses, t_end)
	e = enc[planet]
	E0, h0 = e["before"]
	cr0 = Crossing(planet, *local_velocity(E0, h0, PLANETS[planet].a, sign))
	E_an, h_an, _, dv = state_after(cr0, e["b_meas"])
	return {"b_target": b_p, "b_meas": e["b_meas"], "rmin": e["rmin"],
			"E_pred": E_an, "E_nbody": e["after"][0],
			"h_pred": h_an, "h_nbody": e["after"][1], "dv_over_v": dv / cr0.v_mag}


def plan_chain(v_inf, b, i, sign_i, bp_i, j, bp_j):
	"""Analytic first guess for an i -> j chain.

	Returns ({planet: (phase, b_p, crossing)}, time at j, sign of j's leg).
	"""
	E, h, t_peri = _visitor(v_inf, b)
	ph_i, t_i, cr_i = _placement(E, h, sign_i, t_peri, 0.0, i, bp_i)
	E2, h2, sg2, _ = state_after(cr_i, bp_i)
	a_i = PLANETS[i].a
	theta_i = _polar_angle(E, h, a_i, sign_i, 0.0)
	e2 = math.sqrt(max(0.0, 1 + 2 * E2 * h2 * h2))
	f_i = math.acos(max(-1.0, min(1.0, (h2 * h2 / a_i - 1) / e2))) * sg2
	omega2 = theta_i - (f_i if h2 > 0 else -f_i)
	t_peri2 = t_i - _crossing_time(E2, h2, a_i, sg2, 0.0)
	sign_j = -1 if (PLANETS[j].a < a_i and sg2 < 0) else +1
	ph_j, t_j, cr_j = _placement(E2, h2, sign_j, t_peri2, omega2, j, bp_j)
	return {i: (ph_i, bp_i, cr_i), j: (ph_j, bp_j, cr_j)}, t_j, sign_j


def _aim_at_crossing(enc_j, j, sign, after_t, b_p):
	"""Phase putting planet j at b_p where the flown path crosses its ring."""
	for t, theta, sg, (E, h) in enc_j["crossings"]:
		if t > after_t and sg == sign:
			cr = Crossing(j, *local_velocity(E, h, PLANETS[j].a, sg))
			s = b_p * cr.w_mag / cr.v_r
			return theta + s / PLANETS[j].a - t * PLANETS[j].v_circ / PLANETS[j].a, cr
	return None, None


def check_chain(v_inf, b, i, sign_i, bp_i, j, bp_j, crit, t_extra=3.0, tracks=None):
	"""Fly a targeted i -> j chain, then repeat with planet i massless.

	The removal run replays the same planet phases with i's mass set to zero. `crit` gives the qualifying |b_p| at j; `tracks`, if a dict, receives "chain" and "removed" position samples.
	"""
	targets, t_j, sign_j = plan_chain(v_inf, b, i, sign_i, bp_i, j, bp_j)
	masses = [0.0] * len(PLANETS)
	masses[i], masses[j] = PLANETS[i].mu, PLANETS[j].mu
	tgt = [dict(planet=i, phase=targets[i][0], b=bp_i, cr=targets[i][2], sign=sign_i),
		   dict(planet=j, phase=targets[j][0], b=bp_j, cr=targets[j][2], sign=sign_j)]
	t_end = t_j + t_extra
	phases, enc = _shoot(v_inf, b, tgt, masses, t_end)
	tr_c = [] if tracks is not None else None
	tr_r = [] if tracks is not None else None
	enc = _replay(v_inf, b, phases, masses, t_end, [i, j], track=tr_c)
	removed = list(masses)
	removed[i] = 0.0
	enc_r = _replay(v_inf, b, phases, removed, t_end, [i, j], track=tr_r)
	if tracks is not None:
		tracks["chain"], tracks["removed"] = tr_c, tr_r
	b_crit_j = crit(tgt[1]["cr"])
	return {"t_i": enc[i]["t_min"], "t_j": enc[j]["t_min"], "phases": phases,
			"chain": enc, "removed": enc_r, "b_crit_j": b_crit_j,
			"dv_i": _dv_frac(enc[i], i), "dv_j": _dv_frac(enc[j], j),
			"j_survives": abs(enc_r[j]["b_meas"]) < b_crit_j}


def _dv_frac(e, j):
	"""|dv|/v at ring j reconstructed from the measured (E, h) jump."""
	if e["before"] is None or e["after"] is None:
		return math.nan
	a = PLANETS[j].a
	(E0, h0), (E1, h1) = e["before"], e["after"]
	out = []
	for E, h in ((E0, h0), (E1, h1)):
		vr2 = max(2 * E + 2 / a - (h / a) ** 2, 0.0)
		out.append((math.sqrt(vr2), h / a))
	# radial sign is ambiguous from (E, h); take the smaller |dv| option
	v0 = out[0]
	best = min(math.hypot(sx * out[1][0] - v0[0], out[1][1] - v0[1]) for sx in (1, -1))
	return best / math.hypot(*v0)


def _asymptote_nbody(sim, k):
	"""Outbound asymptotic velocity vector of particle k (barycentric)."""
	s, p = sim.com(last=k), sim.particles[k]
	x, y, vx, vy = p.x - s.x, p.y - s.y, p.vx - s.vx, p.vy - s.vy
	r = math.hypot(x, y)
	v2 = vx * vx + vy * vy
	rv = x * vx + y * vy
	ex, ey = (v2 - 1 / r) * x - rv * vx, (v2 - 1 / r) * y - rv * vy
	E, h = 0.5 * v2 - 1 / r, x * vy - y * vx
	e = math.hypot(ex, ey)
	if E <= 0:
		return None
	f_inf = math.acos(-1 / e)
	th = math.atan2(ey, ex) + (f_inf if h > 0 else -f_inf)
	v = math.sqrt(2 * E)
	return v * math.cos(th), v * math.sin(th)


def fly_plan(v_inf, b, plan, t_margin=5.0, track=None):
	"""Fly a planned sequence of flybys in REBOUND.

	plan: [(planet, sign, b_p)] in time order. Only planned planets have mass. Phases are found by `_shoot`, starting from the Sun-only path's ring crossings. Returns (encounters, outbound asymptote vector).
	"""
	watch = [i for i, _, _ in plan]
	E, h, _ = _visitor(v_inf, b)
	sim = _new_sim([0.0] * len(PLANETS), [0.0] * len(PLANETS))
	_add_visitor(sim, E, h)
	t_last = 0.0
	probe = _run(sim, 400.0, watch)
	for i, sg, _ in plan:
		ts = [t for t, _, s, _ in probe[i]["crossings"] if s == sg]
		t_last = max(t_last, max(ts) if ts else 0.0)
	targets = []
	after = 0.0
	for i, sg, bp in plan:
		ph, cr = _aim_at_crossing(probe[i], i, sg, after, bp)
		targets.append(dict(planet=i, phase=ph, b=bp, cr=cr, sign=sg))
		after = next(t for t, _, s, _ in probe[i]["crossings"] if s == sg and t > after)
	masses = [0.0] * len(PLANETS)
	for i in watch:
		masses[i] = PLANETS[i].mu
	t_end = 1.5 * t_last + t_margin
	phases, enc = _shoot(v_inf, b, targets, masses, t_end)
	sim = _new_sim(phases, masses)
	_add_visitor(sim, E, h)
	enc = _run(sim, t_end, watch, track=track)
	sim.integrate(t_end + 200.0)
	return enc, _asymptote_nbody(sim, len(sim.particles) - 1)
