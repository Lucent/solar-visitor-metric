"""Regression tests for the checks done while building the model.

	.venv/bin/python -m pytest -q
"""
import math

import numpy as np
import pytest

from solar_visitor_metric import nbody, route, scores
from solar_visitor_metric.chain import dv_criterion, from_infinity
from solar_visitor_metric.constants import PLANETS, kms
from solar_visitor_metric.flyby import Crossing, local_velocity, state_after

JUPITER = 4


@pytest.mark.parametrize("b,b_p", [(3.0, 0.01), (3.0, -0.003), (-3.0, 0.01), (-3.0, -0.03)])
def test_patched_conic_matches_nbody(b, b_p):
	"""Post-flyby heliocentric (E, h) agree with REBOUND to < 2%."""
	r = nbody.check_flyby(kms(26.3), b, JUPITER, b_p)
	assert r["b_meas"] == pytest.approx(b_p, rel=2e-3)
	assert r["E_nbody"] == pytest.approx(r["E_pred"], rel=0.02)
	assert r["h_nbody"] == pytest.approx(r["h_pred"], rel=0.02)


def test_deflection_conserves_planetocentric_speed():
	v = kms(26.3)
	E, h = 0.5 * v * v, 0.8 * v
	cr = Crossing(JUPITER, *local_velocity(E, h, PLANETS[JUPITER].a, -1))
	for b_p in (0.002, -0.01, 0.1):
		vr, vt, _ = cr.deflect(b_p)
		w2 = math.hypot(vr, vt - PLANETS[JUPITER].v_circ)
		assert w2 == pytest.approx(cr.w_mag, rel=1e-12)


def test_dv_window_edge_gives_threshold():
	"""At b_p = b_for_dv(f), the flyby changes v by exactly f |v|."""
	v = kms(32.3)
	E, h = 0.5 * v * v, 2.0 * v
	cr = Crossing(JUPITER, *local_velocity(E, h, PLANETS[JUPITER].a, -1))
	b = cr.b_for_dv(0.1)
	_, _, _, dv = state_after(cr, b)
	assert dv / cr.v_mag == pytest.approx(0.1, rel=1e-9)


def test_chain_probabilities_are_a_distribution():
	o = from_infinity(kms(26.3), 1.0, dv_criterion(0.1))
	assert o.dist.sum() + o.impact == pytest.approx(1.0, abs=1e-6)
	assert 0 < o.at_least(2) < o.at_least(1) < 0.01


def test_visit_score_null_calibration():
	"""P(visit bits >= x) <= 2^-x (conservative) and not far below it."""
	s = scores.null_sample(kms(26.3), 40000, seed=11, visits=True)[:, 2]
	for x in (2, 4, 6):
		ratio = (s >= x).mean() * 2 ** x
		assert 0.5 < ratio < 1.15


def test_visit_score_timed_null_calibration():
	"""Same check with real planet phases on a common clock."""
	s = scores.null_sample(kms(26.3), 40000, seed=12, visits=True, timed=True)[:, 2]
	for x in (2, 4, 6):
		ratio = (s >= x).mean() * 2 ** x
		assert 0.5 < ratio < 1.15


def test_timed_walk_places_planet_at_its_phase():
	"""Jupiter set 0.01 rad ahead of the inbound crossing shows up there."""
	v, b, heading, t_peri, ds0 = kms(26.3), 0.8, 0.3, 50.0, 0.01
	E, h = 0.5 * v * v, b * v
	a = PLANETS[JUPITER].a
	f = scores.true_anomaly(E, h, a, -1)
	theta = scores.polar(h, heading, f)
	t_c = t_peri + scores._t_from_peri(E, h, a, -1)
	lam = list(scores.LAMBDA_J2000)
	lam[JUPITER] = math.degrees(theta + ds0 - t_c * a ** -1.5)
	path = scores.fly_timed(v, b, heading, t_peri, lam)
	before = [r for r in path.records if PLANETS[r[0]].a > a]
	assert all(r[4] == 0.0 for r in before)          # nothing deflected earlier
	rec = next(r for r in path.records if r[0] == JUPITER and r[1] == -1)
	cr = Crossing(JUPITER, *local_velocity(E, h, a, -1))
	assert rec[2] == pytest.approx(ds0 / math.pi, rel=1e-9)
	assert rec[3] == pytest.approx(ds0 * a * cr.v_r / cr.w_mag, rel=1e-9)


def test_order_bits_basics():
	assert scores.order_bits([1.0] * 10) == 0.0
	one = scores.order_bits([1e-4] + [0.5] * 9)
	two = scores.order_bits([1e-4, 1e-3] + [0.5] * 8)
	assert two > one > 5


def test_assist_zero_without_flybys():
	path = scores.fly(kms(26.3), 0.8, default_p=0.99)
	assert path.assist() < 1e-6


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


def test_targeted_chain_and_removal():
	"""A Jupiter -> Saturn chain flies as planned; without Jupiter it breaks."""
	import run
	v, crit = kms(26.3), dv_criterion(0.1)
	sg, bp_i, bp_j = run.find_chain(v, 2.0, JUPITER, 5, crit)
	res = nbody.check_chain(v, 2.0, JUPITER, sg, bp_i, 5, bp_j, crit)
	assert res["chain"][JUPITER]["b_meas"] == pytest.approx(bp_i, rel=2e-3)
	assert res["chain"][5]["b_meas"] == pytest.approx(bp_j, rel=2e-3)
	assert res["dv_i"] >= 0.1 and res["dv_j"] >= 0.1
	assert not res["j_survives"]
