"""Plain-SVG sanity-check figures (no plotting dependency).

Colors are CSS custom properties with a dark-mode override, so the files read correctly in either browser theme.
"""
import math
from html import escape

import numpy as np

from .constants import PLANETS, V_UNIT_KMS

LIGHT = {"bg": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#8a8984",
		 "grid": "#e4e3df", "s1": "#2a78d6", "s2": "#eb6834", "s3": "#1baf7a", "sun": "#eda100"}

# Presentation attributes carry the light colors (renders anywhere, incl. viewers without CSS variables); CSS rules override them in dark mode.
STYLE = """<style>
text { font: 12px system-ui, sans-serif; fill: #52514e; }
.title { font-size: 15px; font-weight: 600; fill: #0b0b0b; }
.ink { fill: #0b0b0b; }
.small { font-size: 11px; fill: #8a8984; }
.ring { fill: none; stroke: #e4e3df; stroke-width: 1; }
.grid { fill: none; stroke: #e4e3df; stroke-width: 1; }
.axis { stroke: #8a8984; stroke-width: 1; }
@media (prefers-color-scheme: dark) {
  .bg { fill: #1a1a19; }
  text { fill: #c3c2b7; }
  .title, .ink { fill: #ffffff; }
  .small { fill: #8f8e86; }
  .ring, .grid { stroke: #33322f; }
  .s1 { stroke: #3987e5; } .s2 { stroke: #d95926; } .s3 { stroke: #199e70; }
  .halo { stroke: #1a1a19; }
}
</style>"""

SERIES = ("var(--s1)", "var(--s2)", "var(--s3)")


def _svg(w, h, body):
	return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
			f'viewBox="0 0 {w} {h}">{STYLE}<rect class="bg" width="100%" height="100%" '
			f'fill="{LIGHT["bg"]}"/>{_resolve(body)}</svg>\n')


def _resolve(body):
	"""Inline var(--x) as the light hex; tag stroked series for dark mode."""
	for k in ("s1", "s2", "s3"):
		body = body.replace(f'stroke="var(--{k})"', f'class="{k}" stroke="{LIGHT[k]}"')
	body = body.replace('stroke="var(--bg)"', 'class="halo" stroke="#fcfcfb"')
	for k, v in LIGHT.items():
		body = body.replace(f"var(--{k})", v)
	return body


def _poly(pts, color, width=2, dash=None, clip=None):
	if len(pts) < 2:
		return ""
	d = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
	extra = f' stroke-dasharray="{dash}"' if dash else ""
	extra += f' clip-path="url(#{clip})"' if clip else ""
	return (f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="{width}" '
			f'stroke-linejoin="round" stroke-linecap="round"{extra}/>')


# --- chain figure ---------------------------------------------------------

def chain_figure(path, res, tracks, i, j, v_inf, b, label):
	"""Overview of a targeted i -> j chain plus planet-centered close-ups.

	Solid: real run. Dashed: same initial conditions with planet i massless.
	"""
	W, H = 1180, 620
	ov, R = (300, 340), 250                       # overview center, radius px
	extent = 1.2 * max(PLANETS[i].a, PLANETS[j].a)
	sc = R / extent

	def to_ov(x, y):
		return ov[0] + x * sc, ov[1] - y * sc

	body = [f'<text x="24" y="32" class="title">{escape(label)}</text>',
			f'<text x="24" y="52" class="small">v∞ = {v_inf * V_UNIT_KMS:.1f} km/s, '
			f'b = {b:+.2f} AU (sign = prograde/retrograde), planar, circular orbits. '
			f'Solid: real gravity. Dashed: {PLANETS[i].name} massless (removal test).</text>',
			f'<clipPath id="ovclip"><circle cx="{ov[0]}" cy="{ov[1]}" r="{R}"/></clipPath>']
	for n, p in enumerate(PLANETS):
		if p.a <= extent:
			body.append(f'<circle cx="{ov[0]}" cy="{ov[1]}" r="{p.a * sc:.1f}" class="ring"/>')
	body.append(f'<circle cx="{ov[0]}" cy="{ov[1]}" r="5" fill="var(--sun)"/>')

	for key, color, dash in (("removed", "var(--muted)", "6 5"), ("chain", "var(--s1)", None)):
		pts = [to_ov(*v) for _, v, _ in tracks[key]]
		body.append(_poly(pts, color, 2, dash, "ovclip"))

	# planets at their encounter epochs, with their arcs over the run
	chain = tracks["chain"]
	for k, color in ((i, "var(--s2)"), (j, "var(--s3)")):
		t_k = res["t_i"] if k == i else res["t_j"]
		arc = [to_ov(*pl[k]) for _, _, pl in chain]
		body.append(_poly(arc, color, 1, "2 4"))
		near = min(chain, key=lambda s: abs(s[0] - t_k))
		x, y = to_ov(*near[2][k])
		body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="{color}" '
					f'stroke="var(--bg)" stroke-width="2"/>')
		body.append(f'<text x="{x + 10:.1f}" y="{y - 8:.1f}" class="ink">'
					f'{PLANETS[k].name}</text>')

	# close-ups
	for row, k in enumerate((i, j)):
		body.append(_closeup(res, tracks, k, row == 1, x0=620, y0=90 + row * 260, size=220,
							 tag="ci" + str(row)))
	return open(path, "w").write(_svg(W, H, "".join(body)))


def _closeup(res, tracks, k, is_target, x0, y0, size, tag):
	"""Planet-centered view of one encounter; is_target marks the second
	planet, whose qualifying radius is drawn."""
	enc_c, enc_r = res["chain"][k], res["removed"][k]
	b_meas = enc_c["b_meas"]
	half = max(2.5 * abs(b_meas), 1.5 * res["b_crit_j"] if is_target else 0,
			   4 * PLANETS[k].radius)
	sc = size / (2 * half)
	cx, cy = x0 + size / 2, y0 + size / 2
	out = [f'<clipPath id="{tag}"><rect x="{x0}" y="{y0}" width="{size}" height="{size}"/></clipPath>',
		   f'<rect x="{x0}" y="{y0}" width="{size}" height="{size}" fill="none" class="grid"/>']
	for key, color, dash in (("removed", "var(--muted)", "6 5"), ("chain", "var(--s1)", None)):
		pts = []
		for _, v, pl in tracks[key]:
			dx, dy = v[0] - pl[k][0], v[1] - pl[k][1]
			if abs(dx) < 3 * half and abs(dy) < 3 * half:
				pts.append((cx + dx * sc, cy - dy * sc))
		out.append(_poly(pts, color, 2, dash, tag))
	if is_target:
		rc = res["b_crit_j"] * sc
		out.append(f'<circle cx="{cx}" cy="{cy}" r="{rc:.1f}" fill="none" stroke="var(--s3)" '
				   f'stroke-dasharray="3 3" clip-path="url(#{tag})"/>')
	out.append(f'<circle cx="{cx}" cy="{cy}" r="{max(PLANETS[k].radius * sc, 2):.1f}" '
			   f'fill="{"var(--s3)" if is_target else "var(--s2)"}"/>')
	km = 1.495978707e8
	lines = [f"{PLANETS[k].name}, planet-centered, ±{half * km / 1e6:.2f} M km",
			 f"real: b_p = {b_meas * km / 1e6:+.3f} M km, r_min = {enc_c['rmin'] * km / 1e6:.3f} M km",
			 f"removed: r_min = {enc_r['rmin'] * km / 1e6:.2f} M km"]
	if is_target:
		lines.append(f"dashed ring: qualifying |b_p| = {res['b_crit_j'] * km / 1e6:.3f} M km")
	for n, s in enumerate(lines):
		out.append(f'<text x="{x0 + size + 14}" y="{y0 + 16 + 17 * n}" '
				   f'class="{"small" if n else ""}">{escape(s)}</text>')
	return "".join(out)


# --- probability profile --------------------------------------------------

def profile_figure(path, profiles, crit_label):
	"""P(>=1) and P(>=2) qualifying encounters vs perihelion, per v_inf.

	profiles: list of (name, v_kms, qs, p1_pro, p2_pro, p1_retro, p2_retro). Two stacked panels sharing the x axis (log perihelion), log y each.
	"""
	W, H = 960, 720
	L, Rm, T = 80, 230, 120
	pw, ph = W - L - Rm, 230
	xlo, xhi = math.log10(0.05), math.log10(40)

	def X(q):
		return L + (math.log10(q) - xlo) / (xhi - xlo) * pw

	body = ['<text x="24" y="30" class="title">Chance of qualifying planet encounters '
			'vs perihelion</text>',
			f'<text x="24" y="52" class="small">Null model: random planet phases, planar arrival. '
			f'Criterion: {escape(crit_label)}. Solid prograde, dashed retrograde.</text>']
	panels = [("P(≥1 encounter)", 3, -1, 3), ("P(≥2 encounters: chain)", 4, -8, -4)]
	for row, (title, col, ylo, yhi) in enumerate(panels):
		top = T + row * (ph + 70)
		# auto range from data
		vals = [v for pr in profiles for v in (pr[col] + pr[col + 2]) if v > 0]
		if vals:
			yhi = math.ceil(math.log10(max(vals)))
			ylo = max(yhi - 6, math.floor(math.log10(min(vals))))

		def Y(p, top=top, ylo=ylo, yhi=yhi):
			lp = math.log10(max(p, 10 ** ylo))
			return top + ph - (lp - ylo) / (yhi - ylo) * ph

		body.append(f'<text x="{L}" y="{top - 8}" class="ink">{title}</text>')
		for e in range(ylo, yhi + 1):
			y = Y(10 ** e)
			body.append(f'<line x1="{L}" x2="{L + pw}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>')
			body.append(f'<text x="{L - 8}" y="{y + 4:.1f}" text-anchor="end">1e{e}</text>')
		for p in PLANETS:
			x = X(p.a)
			body.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{top}" y2="{top + ph}" class="grid"/>')
			if row == 0:
				body.append(f'<text x="{x:.1f}" y="{top - 26}" text-anchor="middle" '
							f'class="small">{p.name[:3]}</text>')
		body.append(f'<line x1="{L}" x2="{L + pw}" y1="{top + ph}" y2="{top + ph}" class="axis"/>')
		clip = f"pc{row}"
		body.append(f'<clipPath id="{clip}"><rect x="{L}" y="{top}" width="{pw}" height="{ph}"/></clipPath>')
		for s, (name, v, qs, *ps) in enumerate(profiles):
			pro, retro = ps[col - 3], ps[col - 1]
			for series, dash in ((pro, None), (retro, "5 4")):
				pts = [(X(q), Y(p)) for q, p in zip(qs, series) if p > 0]
				body.append(_poly(pts, SERIES[s], 2, dash, clip))
			if row == 0:
				ly = T + 18 * s
				body.append(f'<line x1="{L + pw + 20}" x2="{L + pw + 44}" y1="{ly}" y2="{ly}" '
							f'stroke="{SERIES[s]}" stroke-width="2"/>')
				body.append(f'<text x="{L + pw + 50}" y="{ly + 4}">{escape(name)} {v:g} km/s</text>')
	for q in (0.1, 0.3, 1, 3, 10, 30):
		x = X(q)
		body.append(f'<text x="{x:.1f}" y="{T + 2 * ph + 70 + 20}" text-anchor="middle">{q:g}</text>')
	body.append(f'<text x="{L + pw / 2}" y="{T + 2 * ph + 70 + 42}" text-anchor="middle">'
				f'perihelion q (AU, log)</text>')
	open(path, "w").write(_svg(W, H, "".join(body)))


# --- scored example paths ---------------------------------------------------

def _conic_points(E, h, omega, f0, f1, n=240):
	e = math.sqrt(max(0.0, 1 + 2 * E * h * h))
	pts = []
	for k in range(n + 1):
		f = f0 + (f1 - f0) * k / n
		den = 1 + e * math.cos(f)
		if den <= 1e-6:
			continue
		r = h * h / den
		th = omega + (f if h > 0 else -f)
		pts.append((r * math.cos(th), r * math.sin(th)))
	return pts


def paths_figure(path_file, examples, extent=6.5):
	"""One panel per scored example: actual path vs the Sun-only path.

	examples: list of (title, Path, lines) with lines = text rows.
	"""
	cols = 2
	pw, ph = 600, 380
	rows = (len(examples) + cols - 1) // cols
	W, H = cols * pw + 20, rows * ph + 70
	body = ['<text x="24" y="30" class="title">Scored example paths (patched conic, planar)</text>',
			'<text x="24" y="50" class="small">Solid: path with planets. Dashed: Sun only. '
			'Dots: planets at their planned flybys.</text>']
	for n, (title, path, lines) in enumerate(examples):
		x0, y0 = 10 + (n % cols) * pw, 70 + (n // cols) * ph
		cx, cy, R = x0 + 150, y0 + 190, 140
		sc = R / extent
		clip = f"pp{n}"
		body.append(f'<clipPath id="{clip}"><rect x="{cx - R}" y="{cy - R}" width="{2 * R}" '
					f'height="{2 * R}"/></clipPath>')
		for p in PLANETS:
			if p.a <= extent * 1.5:
				body.append(f'<circle cx="{cx}" cy="{cy}" r="{p.a * sc:.1f}" class="ring" '
							f'clip-path="url(#{clip})"/>')
		body.append(f'<circle cx="{cx}" cy="{cy}" r="4" fill="var(--sun)"/>')
		E0, h0 = 0.5 * path.v_inf ** 2, path.b * path.v_inf
		e0 = math.sqrt(1 + 2 * E0 * h0 * h0)
		fi = math.acos(-1 / e0) - 1e-3
		ref = [(cx + x * sc, cy - y * sc) for x, y in _conic_points(E0, h0, 0.0, -fi, fi)]
		body.append(_poly(ref, "var(--muted)", 1.5, "5 4", clip))
		for seg in path.segments:
			f1 = seg[4] - 1e-3 if seg is path.segments[-1] else seg[4]
			pts = [(cx + x * sc, cy - y * sc) for x, y in _conic_points(*seg[:4], f1)]
			body.append(_poly(pts, "var(--s1)", 2, None, clip))
		for (i, sg), (x, y) in path.marks.items():
			body.append(f'<circle cx="{cx + x * sc:.1f}" cy="{cy - y * sc:.1f}" r="5" '
						f'fill="var(--s2)" stroke="var(--bg)" stroke-width="1.5"/>')
			body.append(f'<text x="{cx + x * sc + 8:.1f}" y="{cy - y * sc - 6:.1f}" '
						f'class="small">{PLANETS[i].name}</text>')
		body.append(f'<text x="{x0 + 310}" y="{y0 + 70}" class="ink">{escape(title)}</text>')
		for k, s in enumerate(lines):
			body.append(f'<text x="{x0 + 310}" y="{y0 + 92 + 18 * k}" class="small">'
						f'{escape(s)}</text>')
	open(path_file, "w").write(_svg(W, H, "".join(body)))


def tails_figure(path_file, assist_tails, visit_cal, marks):
	"""Left: null tail of the assist measure per v_inf, with examples. Right: visit-score calibration, observed tail x 2^bits (1 = exact).

	assist_tails: list of (name, v_kms, sorted assist sample). visit_cal: (bits grid, ratio list) for one speed. marks: list of (label, assist value).
	"""
	W, H = 980, 470
	body = ['<text x="24" y="30" class="title">Null distributions (sanity checks)</text>',
			'<text x="24" y="50" class="small">Monte Carlo of random planet phases, '
			'arrivals with perihelion inside Neptune.</text>']
	# left panel
	L, T, pw, ph = 80, 90, 480, 300
	xlo, xhi, ylo, yhi = -4, 1, -7, 0

	def X(a):
		return L + (math.log10(a) - xlo) / (xhi - xlo) * pw

	def Y(p):
		return T + ph - (math.log10(p) - ylo) / (yhi - ylo) * ph

	body.append(f'<text x="{L}" y="{T - 12}" class="ink">P(assist ≥ a): change in outbound '
				f'velocity / v∞</text>')
	for e in range(ylo, yhi + 1):
		body.append(f'<line x1="{L}" x2="{L + pw}" y1="{Y(10 ** e):.1f}" y2="{Y(10 ** e):.1f}" class="grid"/>')
		body.append(f'<text x="{L - 8}" y="{Y(10 ** e) + 4:.1f}" text-anchor="end">1e{e}</text>')
	for e in range(xlo, xhi + 1):
		body.append(f'<line x1="{X(10 ** e):.1f}" x2="{X(10 ** e):.1f}" y1="{T}" y2="{T + ph}" class="grid"/>')
		body.append(f'<text x="{X(10 ** e):.1f}" y="{T + ph + 18}" text-anchor="middle">{10 ** e:g}</text>')
	body.append(f'<text x="{L + pw / 2}" y="{T + ph + 38}" text-anchor="middle">a (log)</text>')
	body.append(f'<clipPath id="tl"><rect x="{L}" y="{T}" width="{pw}" height="{ph}"/></clipPath>')
	for s, (name, v, arr) in enumerate(assist_tails):
		n = len(arr)
		grid = np.geomspace(10 ** xlo, 10 ** xhi, 120)
		idx = np.searchsorted(arr, grid, side="left")
		pts = [(X(a), Y((n - k) / n)) for a, k in zip(grid, idx) if n - k > 0]
		body.append(_poly(pts, SERIES[s], 2, None, "tl"))
		body.append(f'<line x1="{L + pw - 150}" x2="{L + pw - 128}" y1="{T + 16 + 18 * s}" '
					f'y2="{T + 16 + 18 * s}" stroke="{SERIES[s]}" stroke-width="2"/>')
		body.append(f'<text x="{L + pw - 122}" y="{T + 20 + 18 * s}">{escape(name)} {v:g}</text>')
	for k, (label, a) in enumerate(marks):
		if a <= 0 or not math.isfinite(a):
			continue
		x = X(a)
		body.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T}" y2="{T + ph}" stroke="var(--ink2)" '
					f'stroke-dasharray="2 3"/>')
		body.append(f'<text x="{x + 4:.1f}" y="{T + ph - 8 - 14 * k}" class="small">{escape(label)}</text>')
	# right panel
	L2, pw2 = 650, 290
	bits, ratio = visit_cal

	def X2(b):
		return L2 + b / max(bits) * pw2

	def Y2(r):
		return T + ph - (math.log2(r) + 2) / 4 * ph   # ratio 1/4 .. 4

	body.append(f'<text x="{L2}" y="{T - 12}" class="ink">Visit score calibration</text>')
	for r in (0.25, 0.5, 1, 2, 4):
		cls = "axis" if r == 1 else "grid"
		body.append(f'<line x1="{L2}" x2="{L2 + pw2}" y1="{Y2(r):.1f}" y2="{Y2(r):.1f}" class="{cls}"/>')
		body.append(f'<text x="{L2 - 8}" y="{Y2(r) + 4:.1f}" text-anchor="end">{r:g}×</text>')
	for b in range(0, int(max(bits)) + 1, 2):
		body.append(f'<text x="{X2(b):.1f}" y="{T + ph + 18}" text-anchor="middle">{b}</text>')
	body.append(f'<text x="{L2 + pw2 / 2}" y="{T + ph + 38}" text-anchor="middle">visit bits</text>')
	body.append(_poly([(X2(b), Y2(max(r, 0.25))) for b, r in zip(bits, ratio) if r > 0],
					  "var(--s1)", 2))
	body.append(f'<text x="{L2}" y="{T + ph + 60}" class="small">observed P(score ≥ x) / 2^-x; '
				f'below 1× = conservative</text>')
	open(path_file, "w").write(_svg(W, H, "".join(body)))


# --- route score ------------------------------------------------------------

def calibration_figure(path_file, title, panels, n_samples=None):
	"""Small multiples: observed null tail x 2^bits per score (1 = exact).

	panels: list of (panel title, [(label, bits array, ratio array)]). Curves stop where fewer than 10 null events are expected (noise).
	"""
	pw, ph, L, T = 300, 230, 70, 100
	W, H = L + len(panels) * (pw + 40) + 20, T + ph + 100
	body = [f'<text x="24" y="30" class="title">{escape(title)}</text>',
			'<text x="24" y="50" class="small">Observed P(score ≥ x) / 2^-x under the null. '
			'1× is exact; below 1× is conservative; above 1× would overstate evidence.</text>']
	for n, (ptitle, series) in enumerate(panels):
		x0 = L + n * (pw + 40)
		bmax = max(float(np.max(b)) for _, b, _ in series)

		def X(b, x0=x0, bmax=bmax):
			return x0 + b / bmax * pw

		def Y(r):
			r = min(max(r, 1 / 16), 4)
			return T + ph - (math.log2(r) + 4) / 6 * ph

		body.append(f'<text x="{x0}" y="{T - 14}" class="ink">{escape(ptitle)}</text>')
		for r in (1 / 16, 1 / 4, 1, 4):
			cls = "axis" if r == 1 else "grid"
			body.append(f'<line x1="{x0}" x2="{x0 + pw}" y1="{Y(r):.1f}" y2="{Y(r):.1f}" class="{cls}"/>')
			if n == 0:
				body.append(f'<text x="{x0 - 8}" y="{Y(r) + 4:.1f}" text-anchor="end">{r:g}×</text>')
		for b in range(0, int(bmax) + 1, 2):
			body.append(f'<text x="{X(b):.1f}" y="{T + ph + 18}" text-anchor="middle">{b}</text>')
		body.append(f'<text x="{x0 + pw / 2}" y="{T + ph + 38}" text-anchor="middle">bits</text>')
		for s, (label, bits, ratio) in enumerate(series):
			pts = [(X(b), Y(r)) for b, r in zip(bits, ratio)
				   if r > 0 and (n_samples is None or n_samples * 2.0 ** -b >= 10)]
			body.append(_poly(pts, SERIES[s], 2))
			ly = T + ph + 58 + 16 * s
			body.append(f'<line x1="{x0}" x2="{x0 + 22}" y1="{ly}" y2="{ly}" stroke="{SERIES[s]}" stroke-width="2"/>')
			body.append(f'<text x="{x0 + 28}" y="{ly + 4}">{escape(label)}</text>')
	open(path_file, "w").write(_svg(W, H + 20, "".join(body)))


def route_map(path_file, field, V_in, V_out, ends, title, lines, slab=3.0):
	"""Route in the plane of the inbound and outbound asymptotes.

	ends: list of (star index, encounter time Myr, label). Stars within `slab` pc of the plane are drawn at their present positions; the origin/destination are drawn moving to where the visitor meets them.
	"""
	from .route import KMS_PC_MYR
	e1 = -V_in / np.linalg.norm(V_in)
	e3 = np.cross(V_in, V_out)
	e3 /= np.linalg.norm(e3)
	e2 = np.cross(e3, e1)
	met = [field.pos[k] + field.vel[k] * t * KMS_PC_MYR for k, t, _ in ends]
	ext = 1.15 * max(max(field.D[k] for k, _, _ in ends),
					 max(max(abs(m @ e1), abs(m @ e2)) for m in met))
	W, H = 980, 640
	cx, cy, R = 330, 350, 270
	sc = R / ext

	def P(v):
		return cx + (v @ e1) * sc, cy - (v @ e2) * sc

	body = [f'<text x="24" y="30" class="title">{escape(title)}</text>',
			f'<text x="24" y="50" class="small">Plane of the inbound and outbound asymptotes; '
			f'gray: stars within ±{slab:g} pc of it, at their present positions.</text>',
			f'<clipPath id="rm"><rect x="{cx - R}" y="{cy - R}" width="{2 * R}" height="{2 * R}"/></clipPath>',
			f'<rect x="{cx - R}" y="{cy - R}" width="{2 * R}" height="{2 * R}" class="grid"/>']
	inslab = (np.abs(field.pos @ e3) < slab) & (field.D < 1.5 * ext)
	for p in field.pos[inslab]:
		x, y = P(p)
		body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2" fill="var(--muted)" clip-path="url(#rm)"/>')
	for V, (k, t, label), color in ((V_in, ends[0], "var(--s1)"), (V_out, ends[1], "var(--s1)")):
		ts = np.linspace(0, t, 50)
		body.append(_poly([P(V * tt * KMS_PC_MYR) for tt in ts], color, 2, None, "rm"))
		track = [P(field.pos[k] + field.vel[k] * tt * KMS_PC_MYR) for tt in ts]
		body.append(_poly(track, "var(--s2)", 1.5, "4 3", "rm"))
		x0, y0 = track[0]
		x1, y1 = track[-1]
		body.append(f'<circle cx="{x0:.1f}" cy="{y0:.1f}" r="4" fill="none" stroke="var(--s2)" stroke-width="1.5"/>')
		body.append(f'<circle cx="{x1:.1f}" cy="{y1:.1f}" r="5" fill="var(--s2)"/>')
		body.append(f'<text x="{x1 + 8:.1f}" y="{y1 - 6:.1f}" class="ink">{escape(label)}</text>')
	body.append(f'<circle cx="{cx}" cy="{cy}" r="5" fill="var(--sun)"/>')
	body.append(f'<text x="{cx + 8}" y="{cy + 16}" class="small">Sun</text>')
	for k, s in enumerate(lines):
		body.append(f'<text x="640" y="{110 + 18 * k}" class="{"ink" if k == 0 else "small"}">{escape(s)}</text>')
	body.append(f'<text x="640" y="{110 + 18 * (len(lines) + 1)}" class="small">Hollow: star now. '
				f'Filled: star when the visitor passes.</text>')
	open(path_file, "w").write(_svg(W, H, "".join(body)))


# --- real visitors ----------------------------------------------------------

def real_figure(path_file, results, extent=6.0):
	"""Top-down ecliptic view of each real visitor's path near the Sun.

	results: list of dicts from real.score plus 'lines' (text rows). Planets are drawn at the visitor's closest-approach time to each of them.
	"""
	pw, ph = 520, 470
	W, H = len(results) * pw + 20, ph + 90
	body = ['<text x="24" y="30" class="title">The three known interstellar objects</text>',
			'<text x="24" y="50" class="small">Ecliptic plane, top-down, ±6 AU. Path: N-body with '
			'JPL Horizons planets. Dots: each planet when the visitor passed closest.</text>']
	for n, r in enumerate(results):
		x0, y0 = 10 + n * pw, 70
		cx, cy, R = x0 + 200, y0 + 210, 190
		sc = R / extent
		clip = f"rv{n}"
		body.append(f'<clipPath id="{clip}"><rect x="{cx - R}" y="{cy - R}" width="{2 * R}" '
					f'height="{2 * R}"/></clipPath>')
		body.append(f'<rect x="{cx - R}" y="{cy - R}" width="{2 * R}" height="{2 * R}" class="grid"/>')
		for p in PLANETS:
			if p.a < extent * 1.5:
				body.append(f'<circle cx="{cx}" cy="{cy}" r="{p.a * sc:.1f}" class="ring" '
							f'clip-path="url(#{clip})"/>')
		pos = r["pos"]
		iso = pos[::10, -1]
		body.append(_poly([(cx + x * sc, cy - y * sc) for x, y, _ in iso], "var(--s1)", 2, None, clip))
		body.append(f'<circle cx="{cx}" cy="{cy}" r="4" fill="var(--sun)"/>')
		for j, (d, t, _) in enumerate(r["close"], start=1):
			k = int(np.argmin(np.abs(r["times"] - t)))
			x, y = pos[k, j, 0], pos[k, j, 1]
			if abs(x) > extent or abs(y) > extent:
				continue
			body.append(f'<circle cx="{cx + x * sc:.1f}" cy="{cy - y * sc:.1f}" r="4" '
						f'fill="var(--s2)" stroke="var(--bg)" stroke-width="1.5"/>')
			body.append(f'<text x="{cx + x * sc + 7:.1f}" y="{cy - y * sc - 5:.1f}" '
						f'class="small">{PLANETS[j - 1].name}</text>')
		body.append(f'<text x="{x0 + 10}" y="{y0 + 425}" class="ink">{escape(r["name"])}</text>')
		for k, s in enumerate(r["lines"]):
			body.append(f'<text x="{x0 + 10}" y="{y0 + 445 + 16 * k}" class="small">{escape(s)}</text>')
	open(path_file, "w").write(_svg(W, H + 60, "".join(body)))
