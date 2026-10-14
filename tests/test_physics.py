"""Regression tests for the checks done while building the model.

	.venv/bin/python -m pytest -q
"""
import math

import numpy as np
import pytest

from solar_visitor_metric import planets, route


def circular_planets():
	"""Eight coplanar circular planets at the real semimajor axes, all at longitude zero at t = 0."""
	a = np.array([0.387, 0.723, 1.0, 1.524, 5.203, 9.537, 19.19, 30.07])
	return planets.Orbits(a, np.zeros(8), np.broadcast_to(np.eye(3), (8, 3, 3)), np.zeros(8), 0.0)


def test_planet_pass_through_a_planet_is_found():
	"""A path planted through Earth's position at its outbound crossing of 1 AU passes Earth at ~0; Earth a day late misses it."""
	orbits = circular_planets()
	q, v = 0.5, 0.02
	e = 1 + q * v * v / planets.MU
	nu = math.acos((q * (1 + e) - 1) / e)
	F = 2 * math.atanh(math.sqrt((e - 1) / (e + 1)) * math.tan(nu / 2))
	t1 = (e * math.sinh(F) - F) * planets.MU / v ** 3
	n_earth = math.sqrt(planets.MU)
	R = np.eye(3)[None]
	orbits.M0[2] = nu - n_earth * t1
	d, t = planets.passes(q, v, R, np.array([0.0]), orbits)
	assert d[0, 2] < 1e-6 and t[0, 2] == pytest.approx(t1, abs=1e-3)
	orbits.M0[2] -= n_earth
	assert 0.005 < planets.passes(q, v, R, np.array([0.0]), orbits)[0][0, 2] < 0.03


def test_planet_conic_conserves_energy():
	q, v = 0.3, 0.015
	X, V = planets.conic(q, v, planets.random_rotations(np.random.default_rng(1), 4), np.zeros(4), np.linspace(-400, 400, 8)[None].repeat(4, 0))
	energy = 0.5 * np.sum(V * V, -1) - planets.MU / np.linalg.norm(X, axis=-1)
	assert np.allclose(energy, 0.5 * v * v, rtol=1e-9)


def test_planet_rotations_are_uniform():
	R = planets.random_rotations(np.random.default_rng(2), 200_000)
	assert np.allclose(np.einsum("nab,nac->nbc", R, R), np.eye(3), atol=1e-12)
	z = R[:, :, 2]
	assert np.all(np.abs(z.mean(0)) < 0.01) and np.abs((z[:, 2] ** 2).mean() - 1 / 3) < 0.01


def test_planet_statistic_counts_several_passes():
	one = planets.statistic(np.array([1e-3] + [0.5] * 7))[0]
	three = planets.statistic(np.array([0.01] * 3 + [0.5] * 5))[0]
	none = planets.statistic(np.full(8, 0.9))[0]
	assert three < one < 0.05 < none


def test_planet_null_calibration(monkeypatch):
	"""P(bits >= x) <= 2^-x on fresh null visitors, and not far below it."""
	monkeypatch.setattr(planets, "N_CURVE", 40_000)
	monkeypatch.setattr(planets, "N_STAT", 40_000)
	orbits = circular_planets()
	q, v = 0.5, 0.02
	null = planets.Null(q, v, orbits)
	d = planets.null_passes(q, v, orbits, 20_000, (7,))
	bits = np.array([null.bits(x) for x in d])
	for x in (1, 2, 4):
		assert 0.75 < (bits >= x).mean() * 2 ** x < 1.15


def test_route_p_values_match_brute_force():
	rng = np.random.default_rng(1)
	field = route.Field.synthetic(rng)
	idx = np.argsort(field.D)[[0, 3]]
	sub = route.Field(field.pos[idx], field.vel[idx], field.D[idx])
	v, M = 26.3, 200000
	Vs = rng.normal(size=(M, 3))
	Vs = v * Vs / np.linalg.norm(Vs, axis=1)[:, None]
	for frac in (0.1, 0.5):
		p = route.p_values(sub, v, frac * sub.D, -1)
		for k in range(2):
			u = sub.vel[k] - Vs
			c = (u @ sub.pos[k]) / np.linalg.norm(u, axis=1) / sub.D[k]
			mc = (c >= math.cos(math.asin(frac))).mean()
			assert p[k] == pytest.approx(mc, abs=4 * math.sqrt(mc / M) + 1e-5)


def test_route_hipparcos_meets_gaia_at_its_epoch():
	"""Barnard's star, the fastest across the sky, moves 4 arcmin between the catalogs' epochs; propagated, Hipparcos lands within MATCH of Gaia with the same space velocity."""
	cols = "ra dec plx s_plx pmra s_pmra pmdec s_pmdec rv s_rv"
	hip = dict(zip(cols.split(), "269.45402263 4.66828781 548.31 1.51 -798.58 1.72 10328.12 1.22 -110.6 0.2".split()))
	gaia = dict(zip(cols.split(), "269.44850253 4.73942005 546.97594 0.04012 -801.55098 0.03182 10362.39421 0.03607 -110.46822 0.13126".split()))
	h = route._kinematics([hip], cols, route.GAIA_YEAR - route.HIPPARCOS_YEAR)
	g = route._kinematics([gaia], cols, 0.0)
	assert route._near(h.unit, g.unit)[0]
	assert not route._near(route._kinematics([hip], cols, 0.0).unit, g.unit)[0]
	assert np.linalg.norm(h.vel - g.vel) < 1.0


def test_route_aim_hits_star():
	rng = np.random.default_rng(2)
	field = route.Field.synthetic(rng)
	for k in np.argsort(field.D)[:40]:
		V, t = route.aim_at(field.pos[k], field.vel[k], 58.0, -1)
		if V is None:
			continue
		assert np.linalg.norm(V) == pytest.approx(58.0, rel=1e-9)
		d, _ = route.encounters(field, V, -1)
		assert d[k] < 1e-6
		return
	pytest.fail("no reachable star")


def test_route_hops_match_aim():
	"""The hop is where an exactly aimed visitor meets the star; candidates reach well beyond HOP today."""
	rng = np.random.default_rng(3)
	field = route.Field.synthetic(rng)
	v = 26.3
	L = route.hops(field, v, -1)
	ok = np.flatnonzero(L <= route.HOP)
	assert (field.D[ok] > route.HOP).mean() > 0.5
	n = 0
	for k in ok[:100]:
		V, t = route.aim_at(field.pos[k], field.vel[k], v, -1)
		if V is None:
			continue
		n += 1
		assert v * abs(t) * route.KMS_PC_MYR == pytest.approx(L[k], rel=1e-9)
		assert np.linalg.norm(field.pos[k] + (field.vel[k] - V) * t * route.KMS_PC_MYR) < 1e-6
	assert n > 50


def test_route_turn_score_matches_brute_force():
	"""The departure p is the Rice-convolved fraction of the turn circle passing the star this closely."""
	from scipy.stats import rice
	rng = np.random.default_rng(4)
	field = route.Field.synthetic(rng)
	err = route.Errors()
	V_in = 26.3 * route.random_unit(rng)
	V_out = route.sun_turn(V_in, 0.5, 1.0)
	_, _, info = route.turn_score(field, V_in, V_out, err, rng)
	turn = route.turn_circle(V_in, V_out)
	xs = np.linspace(0, 2 * math.pi, 100000, endpoint=False)
	cand = np.flatnonzero(info["weights"])
	for j in cand[np.argsort(info["p"][cand] / info["weights"][cand])[:3]]:
		r, _ = route.miss(field.pos[j], field.vel[j], turn(xs), +1)
		s = info["sigma"][j]
		assert info["p"][j] == pytest.approx(rice.cdf(info["d"][j] / s, r / s).mean(), rel=0.02)
