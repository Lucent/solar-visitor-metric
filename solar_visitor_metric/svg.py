"""Plain-SVG sanity-check figures (no plotting dependency).

Colors are CSS custom properties with a dark-mode override, so the files read correctly in either browser theme.
"""
import math
from html import escape

import numpy as np


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
