import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { CSS2DRenderer, CSS2DObject } from "three/addons/renderers/CSS2DRenderer.js";
import { Line2 } from "three/addons/lines/Line2.js";
import { LineGeometry } from "three/addons/lines/LineGeometry.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";
import * as Astronomy from "astronomy-engine";

const id = (id) => document.getElementById(id);
const data = await (await fetch("data.json")).json();

const J2000 = 2451545;
const PC = 206264.806;
const NOW = Date.now() / 864e5 + 2440587.5;
const ECL = Astronomy.Rotation_EQJ_ECL();
const SPAN = Math.asinh(data.time.span / data.time.tau);
const form = document.forms[0];
const stars = { pos: Float64Array.from(data.stars.pos), vel: Float64Array.from(data.stars.vel) };
const visitors = Object.fromEntries(data.visitors.map((v) => [v.tag, { ...v, path: Float64Array.from(v.path) }]));
const planetNames = data.visitors[0].planets.map((p) => p.name);

// One law for space and time. A distance r is drawn at asinh(r / r0): to scale inside the zoom radius r0, one unit per factor of e beyond it; or, with Linear checked, at r / r0 everywhere. The slider u runs time the same way, t - t_peri = tau sinh(u asinh(span / tau)).
const r0 = () => 10 ** form.zoom.value;
const scale = (r) => form.linear.checked ? r / r0() : Math.asinh(r / r0());
const visitor = () => visitors[form.visitor.value];
const time = () => visitor().tp + data.time.tau * Math.sinh(form.u.value * SPAN);
const slider = (jd) => Math.asinh((jd - visitor().tp) / data.time.tau) / SPAN;

function mapped(x, y, z, out, i) {
	const r = Math.hypot(x, y, z);
	const k = scale(r) / r;
	out[i] = x * k;
	out[i + 1] = y * k;
	out[i + 2] = z * k;
}

// Write true positions (AU) into a geometry's mapped position attribute.
function place(geometry, xyz) {
	const attr = geometry.getAttribute("position");
	for (let i = 0; i < xyz.length; i += 3)
		mapped(xyz[i], xyz[i + 1], xyz[i + 2], attr.array, i);
	attr.needsUpdate = true;
	geometry.computeBoundingSphere();
}

function buffer(n) {
	const g = new THREE.BufferGeometry();
	g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * n), 3));
	return g;
}

// A label sits beside its mark, not on it: anchored by its lower-left corner for a point, and centered above the outline for a sphere.
function label(text, kind) {
	const el = document.createElement("div");
	el.className = `label ${kind}`;
	el.textContent = text;
	const tag = new CSS2DObject(el);
	tag.center.set(kind == "shell" ? 0.5 : 0, 1);
	return tag;
}

const fmt = (x, digits = 2) => x.toLocaleString("en", { maximumSignificantDigits: digits });
const pct = (p) => `${fmt(100 * p)}%`;

function when(jd) {
	const years = (jd - NOW) / 365.25;
	const side = years < 0 ? "ago" : "from now";
	if (Math.abs(years) < 2000)
		return new Date((jd - 2440587.5) * 864e5).toISOString().slice(0, 10);
	if (Math.abs(years) < 1e6)
		return `${Math.round(Math.abs(years)).toLocaleString()} years ${side}`;
	return `${fmt(Math.abs(years) / 1e6)} million years ${side}`;
}

// Round, softly edged points; WebGL's own are hard squares. A ring marks where the visitor meets a star.
function sprite(draw) {
	const c = document.createElement("canvas");
	c.width = c.height = 64;
	const g = c.getContext("2d");
	g.beginPath();
	g.arc(32, 32, 28, 0, 2 * Math.PI);
	draw(g);
	return new THREE.CanvasTexture(c);
}

const DISC = sprite((g) => {
	g.fillStyle = "#fff";
	g.fill();
});
const RING = sprite((g) => {
	g.lineWidth = 7;
	g.strokeStyle = "#fff";
	g.stroke();
});
const point = (size, extra) => new THREE.PointsMaterial({ size, sizeAttenuation: false, map: DISC, transparent: true, depthWrite: false, ...extra });

// --- scene ---------------------------------------------------------------

const scene = new THREE.Scene();
// Depth for the star field alone: stars beyond the Sun fade toward the page, so turning the view reads as near and far.
scene.fog = new THREE.Fog(0, 30, 70);
const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 500);
camera.up.set(0, 0, 1);
camera.position.set(0, -34, 22);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setPixelRatio(devicePixelRatio);
const labels = new CSS2DRenderer();
id("scene").append(renderer.domElement, labels.domElement);
const controls = new OrbitControls(camera, labels.domElement);
controls.enableZoom = false;
// A drag coasts to a stop, so the parallax that shows depth plays on after the hand lets go.
controls.enableDamping = true;
renderer.setAnimationLoop(() => controls.update());
controls.enablePan = false;
controls.addEventListener("change", render);

const materials = {
	// A distance from the Sun in every direction is a sphere: drawn as a faint limb that brightens toward its edge, it is round from every angle, while orbits tilt into ellipses.
	shell: new THREE.ShaderMaterial({
		uniforms: { color: { value: new THREE.Color() } },
		vertexShader: `varying vec3 n; varying vec3 view; void main() { vec4 p = modelViewMatrix * vec4(position, 1.0); n = normalMatrix * normal; view = -p.xyz; gl_Position = projectionMatrix * p; }`,
		fragmentShader: `uniform vec3 color; varying vec3 n; varying vec3 view; void main() { float rim = 1.0 - abs(dot(normalize(n), normalize(view))); gl_FragColor = vec4(color, 0.22 * pow(rim, 7.0)); \n#include <colorspace_fragment>\n }`,
		transparent: true,
		depthWrite: false,
	}),
	orbit: new THREE.LineBasicMaterial(),
	planet: point(7),
	sun: point(12),
	star: point(1.5, { opacity: 0.5 }),
	origin: point(8),
	destination: point(8),
	originMeet: point(18, { map: RING }),
	destinationMeet: point(18, { map: RING }),
	originTrack: new THREE.LineDashedMaterial({ dashSize: 0.15, gapSize: 0.1 }),
	destinationTrack: new THREE.LineDashedMaterial({ dashSize: 0.15, gapSize: 0.1 }),
	originMiss: new THREE.LineBasicMaterial(),
	destinationMiss: new THREE.LineBasicMaterial(),
	visitor: point(11),
	// Stalks drop each mark to the ecliptic, its foot a dot where it would land: their lengths say how far above or below the plane things are, their feet where across it.
	...Object.fromEntries(["planet", "visitor", "origin", "destination"].flatMap((kind) => [
		[`${kind}Stalk`, new THREE.LineBasicMaterial({ transparent: true, opacity: 0.45 })],
		[`${kind}Foot`, point(4, { opacity: 0.7 })],
	])),
	// A curtain of drop lines from the whole path, evenly spaced as drawn: the path's height above the ecliptic everywhere, not just where the visitor is.
	curtain: new LineMaterial({ linewidth: 2, transparent: true, opacity: 0.22 }),
	past: new LineMaterial({ linewidth: 2.5 }),
	ahead: new LineMaterial({ linewidth: 1.5, transparent: true, opacity: 0.35 }),
};

Object.values(materials).forEach((m) => m.fog = m == materials.star);

// Colors live in the stylesheet; the scene reads them, and again when the scheme flips.
function paint() {
	const css = getComputedStyle(document.documentElement);
	const color = (name) => new THREE.Color(css.getPropertyValue(`--${name}`).trim());
	materials.shell.uniforms.color.value = color("muted");
	Object.entries({ orbit: "muted", planet: "ink2", sun: "sun", star: "muted", origin: "s2", destination: "s3", originMeet: "s2", destinationMeet: "s3", originTrack: "s2", destinationTrack: "s3", originMiss: "s2", destinationMiss: "s3", visitor: "s1", past: "s1", ahead: "s1", curtain: "s1",
		planetStalk: "ink2", planetFoot: "ink2", visitorStalk: "s1", visitorFoot: "s1", originStalk: "s2", originFoot: "s2", destinationStalk: "s3", destinationFoot: "s3" })
		.forEach(([m, c]) => materials[m].color = color(c));
	scene.fog.color = color("plane");
	render();
}

// Decades of distance as spheres about the Sun, from 100 AU out (within the planets, the orbits are the scale). The map is radial, so a sphere stays a sphere: a unit sphere scaled like any distance, labeled at the top of its outline.
const sphere = new THREE.SphereGeometry(1, 96, 48);
const shells = [[100, "100 AU"], [1e3, "1,000 AU"], [1e4, "10,000 AU"], [PC, "1 parsec"], [10 * PC, "10 parsecs"], [100 * PC, "100 parsecs"]]
	.map(([R, text]) => {
		const ball = new THREE.Mesh(sphere, materials.shell);
		ball.add(label(text, "shell"));
		ball.userData.R = R;
		scene.add(ball);
		return ball;
	});

// The compass: the ecliptic, its north, and the way to the Galactic center, turned with the camera, in a corner of the view rather than across it. The near half of the circle is solid and the far half dashed, so it never reads flipped.
const GAL = Astronomy.CombineRotation(Astronomy.Rotation_GAL_EQJ(), ECL);
const center = (() => {
	const v = Astronomy.RotateVector(GAL, new Astronomy.Vector(1, 0, 0, Astronomy.MakeTime(0)));
	return new THREE.Vector3(v.x, v.y, v.z);
})();
const RIM = Array.from({ length: 65 }, (_, k) => new THREE.Vector3(Math.cos(k * Math.PI / 32), Math.sin(k * Math.PI / 32), 0));
function compass() {
	const seen = (w) => w.clone().transformDirection(camera.matrixWorldInverse);
	const xy = (v, r = 0.8) => `${(r * v.x).toFixed(3)} ${(-r * v.y).toFixed(3)}`;
	const rim = RIM.map(seen);
	const half = (front) => rim.map((v, k) => `${k && (rim[k - 1].z > 0) == front && (v.z > 0) == front ? "L" : "M"}${xy(v)}`).join("");
	const box = id("compass");
	box.querySelector(".near").setAttribute("d", half(true));
	box.querySelector(".far").setAttribute("d", half(false));
	const north = seen(new THREE.Vector3(0, 0, 1));
	box.querySelector(".north").setAttribute("d", `M0 0L${xy(north)}`);
	box.querySelector(".north + text").setAttribute("transform", `translate(${xy(north, 0.95)})`);
	const gc = seen(center);
	box.querySelector(".center").setAttribute("transform", `translate(${xy(gc)})`);
}

// Planet orbits from one revolution of Astronomy Engine's ephemeris.
const helio = (name, jd) => Astronomy.RotateVector(ECL, Astronomy.HelioVector(name, jd - J2000));
const orbits = planetNames.map((name) => {
	const period = Astronomy.PlanetOrbitalPeriod(name);
	const xyz = Float64Array.from(Array.from({ length: 256 }, (_, k) => helio(name, NOW + k * period / 256))
		.flatMap((v) => [v.x, v.y, v.z]));
	const loop = new THREE.LineLoop(buffer(256), materials.orbit);
	loop.userData.xyz = xyz;
	scene.add(loop);
	return loop;
});
const planets = new THREE.Points(buffer(planetNames.length), materials.planet);
const planetLabels = planetNames.map((name) => label(name, "planet"));
scene.add(planets, ...planetLabels);
scene.add(new THREE.Points(buffer(1), materials.sun));

const field = new THREE.Points(buffer(stars.pos.length / 3), materials.star);
scene.add(field);

// The path bright where the visitor has been, faint where it is going.
const past = new Line2(new LineGeometry(), materials.past);
const ahead = new Line2(new LineGeometry(), materials.ahead);
let route;
const traveler = new THREE.Points(buffer(1), materials.visitor);
const travelerLabel = label("", "visitor");
scene.add(past, ahead, traveler, travelerLabel);

function stalks(kind, n) {
	const lines = new THREE.LineSegments(buffer(2 * n), materials[`${kind}Stalk`]);
	const feet = new THREE.Points(buffer(n), materials[`${kind}Foot`]);
	scene.add(lines, feet);
	return { lines, feet };
}
// Drop the first `count` drawn points of `from` to the plane z = 0.
function drop({ lines, feet }, from, count, visible = true) {
	const p = from.getAttribute("position").array;
	const l = lines.geometry.getAttribute("position").array;
	const f = feet.geometry.getAttribute("position").array;
	for (let i = 0; i < count; i++) {
		l.set([p[3 * i], p[3 * i + 1], p[3 * i + 2], p[3 * i], p[3 * i + 1], 0], 6 * i);
		f.set([p[3 * i], p[3 * i + 1], 0], 3 * i);
	}
	[[lines, 2 * count], [feet, count]].forEach(([o, k]) => {
		o.geometry.getAttribute("position").needsUpdate = true;
		o.geometry.setDrawRange(0, k);
		o.geometry.computeBoundingSphere();
		o.visible = visible;
	});
}
const planetStalks = stalks("planet", planetNames.length);
const travelerStalk = stalks("visitor", 1);
const CURTAIN = 0.8;   // scene units between drop lines along the drawn path
const curtain = new LineSegments2(new LineSegmentsGeometry(), materials.curtain);
scene.add(curtain);

// Each end's best stars, each with its track from now to where the visitor meets it; the best one's meeting is a ring on the path, joined to the star by the miss.
const TRACK = 32;
const ends = Object.fromEntries(["origin", "destination"].map((end) => {
	const best = new THREE.Points(buffer(5), materials[end]);
	const tracks = new THREE.LineSegments(buffer(5 * 2 * (TRACK - 1)), materials[`${end}Track`]);
	const tag = label("", end);
	const meet = new THREE.Points(buffer(1), materials[`${end}Meet`]);
	const gap = new THREE.LineSegments(buffer(2), materials[`${end}Miss`]);
	// The stars' tracks to their meetings read as clutter; kept out of the scene for now.
	scene.add(best, /* tracks, */ tag, meet, gap);
	return [end, { best, tracks, tag, meet, gap, stalks: stalks(end, 5) }];
}));

// --- state -> picture ------------------------------------------------------

function starAt(i, jd, out, k) {
	const dt = jd - data.epoch;
	out[k] = stars.pos[3 * i] + stars.vel[3 * i] * dt;
	out[k + 1] = stars.pos[3 * i + 1] + stars.vel[3 * i + 1] * dt;
	out[k + 2] = stars.pos[3 * i + 2] + stars.vel[3 * i + 2] * dt;
}

function visitorAt(v, u) {
	const n = v.path.length / 3;
	const f = (u + 1) / 2 * (n - 1);
	const i = Math.min(Math.floor(f), n - 2);
	const a = f - i;
	return [0, 1, 2].map((c) => v.path[3 * i + c] * (1 - a) + v.path[3 * i + 3 + c] * a);
}

// How far the picture is to scale: everywhere, or out to the zoom radius.
function extent(au) {
	return au < 0.1 * PC ? `${fmt(au)} AU` : `${fmt(au / PC)} pc`;
}

function zoomed() {
	const v = visitor();
	form.radius.value = form.linear.checked ? "to scale throughout" : `to scale within ${extent(r0())}, logarithmic beyond`;
	shells.forEach((ball) => ball.scale.setScalar(scale(ball.userData.R)));
	orbits.forEach((loop) => place(loop.geometry, loop.userData.xyz));
	route = new Float32Array(v.path.length);
	for (let i = 0; i < route.length; i += 3)
		mapped(v.path[i], v.path[i + 1], v.path[i + 2], route, i);
	const drops = [];
	for (let i = 3, along = 0, next = 0; i < route.length; i += 3) {
		const step = Math.hypot(route[i] - route[i - 3], route[i + 1] - route[i - 2], route[i + 2] - route[i - 1]);
		if (!step)
			continue;
		for (; next <= along + step; next += CURTAIN) {
			const f = (next - along) / step;
			const [x, y, z] = [0, 1, 2].map((c) => route[i - 3 + c] + (route[i + c] - route[i - 3 + c]) * f);
			drops.push(x, y, z, x, y, 0);
		}
		along += step;
	}
	// a new geometry each time: a line geometry draws the segment count it was first rendered with
	curtain.geometry.dispose();
	curtain.geometry = new LineSegmentsGeometry().setPositions(drops);
	Object.entries(ends).forEach(([end, { tracks, meet, gap }]) => {
		const s = v.route[end].stars[0];
		const at = visitorAt(v, slider(s.jd));
		place(meet.geometry, at);
		starAt(s.i, s.jd, at, 3);
		place(gap.geometry, at);
		const xyz = new Float64Array(5 * 2 * (TRACK - 1) * 3);
		v.route[end].stars.forEach((s, j) => {
			const at = Array.from({ length: TRACK }, (_, k) => data.epoch + (s.jd - data.epoch) * k / (TRACK - 1));
			at.slice(1).forEach((jd, k) => {
				starAt(s.i, at[k], xyz, 6 * ((TRACK - 1) * j + k));
				starAt(s.i, jd, xyz, 6 * ((TRACK - 1) * j + k) + 3);
			});
		});
		place(tracks.geometry, xyz);
		tracks.computeLineDistances();
	});
	moved();
}

// A planet label hides while time sweeps its planet across the screen faster than READABLE pixels a second; the camera's own motion does not count, as both ends are projected through the same camera.
const READABLE = 120;
let lastMoved = performance.now();
let movedAt = lastMoved;
function blur(from, to) {
	const box = labels.domElement.getBoundingClientRect();
	const a = from.clone().project(camera);
	const b = to.clone().project(camera);
	const px = Math.hypot((b.x - a.x) * box.width, (b.y - a.y) * box.height) / 2;
	return px / Math.max(movedAt - lastMoved, 1) * 1000 > READABLE;
}

let settle;
function moved() {
	lastMoved = movedAt;
	movedAt = performance.now();
	// once time stops, every label returns
	clearTimeout(settle);
	settle = setTimeout(() => planetLabels.forEach((tag) => tag.element.classList.remove("fast")), 250);
	const v = visitor();
	const jd = time();
	form.date.value = when(jd);
	const p = visitorAt(v, +form.u.value);
	place(traveler.geometry, p);
	const here = traveler.geometry.getAttribute("position").array;
	const split = 3 * (Math.min(Math.floor((+form.u.value + 1) / 2 * (route.length / 3 - 1)), route.length / 3 - 2) + 1);
	// A line geometry draws the segment count it was first rendered with, so each new split gets a new geometry.
	[[past, Float32Array.of(...route.subarray(0, split), ...here)], [ahead, Float32Array.of(...here, ...route.subarray(split))]].forEach(([line, xyz]) => {
		line.geometry.dispose();
		line.geometry = new LineGeometry().setPositions(xyz);
	});
	drop(travelerStalk, traveler.geometry, 1);
	travelerLabel.element.textContent = v.tag;
	travelerLabel.position.fromArray(traveler.geometry.getAttribute("position").array);
	// The ephemeris holds for millennia, not for the million years the stars need.
	const near = Math.abs(jd - NOW) < 3000 * 365.25;
	planets.visible = near;
	planetLabels.forEach((tag) => tag.visible = near);
	if (near) {
		place(planets.geometry, planetNames.flatMap((name) => {
			const h = helio(name, jd);
			return [h.x, h.y, h.z];
		}));
		planetLabels.forEach((tag, k) => {
			const was = tag.position.clone();
			tag.position.fromArray(planets.geometry.getAttribute("position").array, 3 * k);
			tag.element.classList.toggle("fast", blur(was, tag.position));
		});
	}
	drop(planetStalks, planets.geometry, planetNames.length, near);
	const xyz = new Float64Array(stars.pos.length);
	for (let i = 0; i < xyz.length / 3; i++)
		starAt(i, jd, xyz, 3 * i);
	place(field.geometry, xyz);
	Object.entries(ends).forEach(([end, { best, tag, stalks }]) => {
		const shown = v.route[end].stars;
		const at = new Float64Array(15);
		shown.forEach((s, j) => starAt(s.i, jd, at, 3 * j));
		place(best.geometry, at);
		best.geometry.setDrawRange(0, shown.length);
		drop(stalks, best.geometry, shown.length);
		tag.element.textContent = `${end} ${fmt(shown[0].bits)} bits`;
		tag.position.fromArray(best.geometry.getAttribute("position").array);
	});
	render();
}

function render() {
	const up = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 1);
	shells.forEach((ball) => ball.children[0].position.copy(up));
	renderer.render(scene, camera);
	labels.render(scene, camera);
	compass();
}

new ResizeObserver(([{ contentRect: { width, height } }]) => {
	renderer.setSize(width, height, false);
	labels.setSize(width, height);
	materials.past.resolution.set(width, height);
	materials.ahead.resolution.set(width, height);
	materials.curtain.resolution.set(width, height);
	camera.aspect = width / height;
	camera.updateProjectionMatrix();
	render();
}).observe(id("scene"));

// --- the panels: built once; from then on the radios are the state ------------

// How many random arrivals at this speed do as well: out of 100 while that is at least one, as 1 in N beyond.
const chance = (p) => p >= 0.01 ? `${Math.round(100 * p)} in 100` : p > 1e-6 ? `1 in ${fmt(1 / p)}` : `1 in ${fmt(1e-6 / p)} million`;
const REPO = "https://github.com/lucent/solar-visitor-metric";
// A score's method is a reference, as a citation is: a superscript after its name.
const method = (part) => `<sup><a href="${REPO}#:~:text=${encodeURIComponent(part)}">how</a></sup>`;
// A star goes by its SIMBAD name where it has one, its catalog designation otherwise; the link opens it in SIMBAD.
const named = (i) => data.stars.names[i] || data.stars.ids[i];
const simbad = (i) => `https://simbad.cds.unistra.fr/simbad/sim-id?Ident=${encodeURIComponent(data.stars.ids[i])}`;
const star = (i) => `<a href="${simbad(i)}" title="${data.stars.ids[i]}">${named(i)}</a>`;
const best = (v) => Math.max(v.visit, v.route.bits);

function bits(b) {
	return `<data class="bits" value="${b}" style="--bits: ${b}"><b>${fmt(b)}</b> bits</data>`;
}

function stop(s) {
	return `<tr style="--bits: ${s.bits}"><th>${star(s.i)}</th><td>${fmt(s.miss)} ± ${fmt(s.sigma)} pc</td><td>${fmt(s.hop)} pc</td><td><button type="button" value="${s.jd}">${when(s.jd)}</button></td></tr>`;
}

// The trace that came closest to looking planned, in words.
function closest(v) {
	if (best(v) == 0)
		return "Neither trace: every random arrival at its speed does as well.";
	const [kind, b] = v.visit >= v.route.bits ? ["flybys", v.visit] : ["route", v.route.bits];
	return `Closest to planned: its ${kind}, ${fmt(b)} bits. ${chance(2 ** -b)} random arrivals do as well.`;
}

function panel(v) {
	const a = v.assist;
	const plan = v.route.planned;
	return `<input type="radio" name="visitor" form="controls" id="v${v.tag}" value="${v.tag}"><label for="v${v.tag}" class="visitor"><b>${v.tag}</b>${v.short}${bits(best(v))}</label>
<article class="visitor">
<h2><a href="https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=${v.tag}">${v.name}</a></h2>
<p class="lead">${closest(v)}</p>
<p>${fmt(v.vinf, 3)} km/s far out · perihelion ${fmt(v.q, 3)} AU · <button type="button" value="${v.tp}">${when(v.tp)}</button></p>
<h2><span>Flybys${method("Visit. For every planet")}</span>${bits(v.visit)}</h2>
<p class="lead">${chance(2 ** -v.visit)} random arrivals pass the planets as closely.</p>
<p>Each planet at a random point on its orbit; the score covers all eight together.</p>
<h3>Planet passes</h3>
<table class="planets">
<thead><tr><th></th><th>closest</th><th>planet radii</th><th>date</th><th>chance</th></tr></thead>
${v.planets.map((p) => `<tr style="--bits: ${-Math.log2(p.p)}"><th>${p.name}</th><td>${fmt(p.d)} AU</td><td>${Math.round(p.radii).toLocaleString()}</td><td><button type="button" value="${p.jd}">${when(p.jd)}</button></td><td>${pct(p.p)}</td></tr>`).join("\n")}
</table>
<h3><span>Gravity assist${method("Assist. The visitor is flown twice")}</span><data value="${a.direct}">${pct(a.direct)}</data></h3>
<p class="lead">The planets' pull changed its outbound velocity by ${pct(a.direct)}. A spacecraft's gravity assist changes it by 10% or more.</p>
<p>The Sun's reflex wobble, which shifts every visitor's perihelion a little and is nothing a trajectory could exploit, makes the total ${pct(a.total)}.</p>
<h2><span>Route${method("Route. Inbound and outbound")}</span>${bits(v.route.bits)}</h2>
<p class="lead">${chance(v.route.beaten)} random arrivals pass nearby stars as closely before and after the Sun.</p>
<p>A route planned at this speed from the nearest star it can reach, ${star(plan.origin)}, to the next, ${star(plan.destination)}, scores ${fmt(plan.bits)} bits on the same catalog.</p>
<p>Reachable stars are those it could meet within 20 pc of the Sun; the destination is scored given the origin and the observed ${fmt(v.route.turn)}° turn at the Sun. Gaia radial velocities count only when reliable: <a href="${REPO}#:~:text=${encodeURIComponent("A mismeasured star")}">a mismeasured one</a> once gave ʻOumuamua a false origin.</p>
${["origin", "destination"].map((end) => `<h3><span>${end} <small title="reachable stars">${v.route[end].candidates.toLocaleString()}</small></span>${bits(v.route[end].stars[0].bits)}</h3>
<table class="${end}">
<thead><tr><th>best of the reachable</th><th>missed by</th><th>from Sun</th><th>date</th></tr></thead>
${v.route[end].stars.map(stop).join("")}
</table>`).join("\n")}
</article>`;
}

document.querySelector("aside").innerHTML = data.visitors.map(panel).join("\n");

// The verdict: flybys, route and both together for each visitor, then for three missions flown at one speed. The visitor names select them.
// The visitors' scores are the result, badged; the missions are the yardstick, plain.
const cell = (b, badge) => `<td style="--bits: ${b}">${badge ? bits(b) : `${fmt(b)} bits`}<small>${chance(2 ** -b)}</small></td>`;
const cells = (r, badge) => cell(r.flybys, badge) + cell(r.route, badge) + cell(r.together, badge);
const m = data.missions;
const mission = { flybys: m.flyby, route: `${named(m.origin)} to ${named(m.destination)}, turning at the Sun`, both: "Both at once" };
const by = visitors[data.visitors.find((v) => Math.abs(v.vinf - m.speed) < 1).tag];
id("verdict").innerHTML = `<thead><tr><th></th><th>Flybys</th><th>Route</th><th>Together</th></tr></thead>
<tbody class="visitor">${data.visitors.map((v) => `<tr><th><label for="v${v.tag}">${v.tag} ${v.short}</label></th>${cells({ flybys: v.visit, route: v.route.bits, together: v.together }, true)}</tr>`).join("")}</tbody>
<tbody class="planned"><tr><th colspan="4">Simulated missions at ${by.short}'s speed, ${fmt(m.speed)} km/s</th></tr>
${m.rows.map((r) => `<tr><th>${mission[r.kind]}</th>${cells(r)}</tr>`).join("")}</tbody>
<tfoot><tr><td></td>${`<td>${[[1, "coin flip"], [10, "1 in 1,000"], [19.3, "royal flush"]].map(([at, text]) => `<span style="--at: ${at}">${text}</span>`).join("")}</td>`.repeat(3)}</tr></tfoot>`;
document.querySelector(`#v${data.visitors.at(-1).tag}`).checked = true;
id("now").value = NOW;

// The time ruler: decades marked where they fall on the slider, under a band as thick as the log of the time each step covers, so the ends visibly race and the middle crawls.
const tick = (days) => (Math.asinh(days / data.time.tau) / SPAN + 1) / 2;
const band = (u) => Math.log(data.time.tau * SPAN * Math.cosh(u * SPAN));
const thin = band(0);
const thick = band(1);
const profile = Array.from({ length: 41 }, (_, k) => k / 20 - 1).map((u) => [(u + 1) * 50, 50 * (band(u) - thin * 0.6) / (thick - thin * 0.6)]);
id("ruler").style.setProperty("--band", `polygon(${[...profile.map(([x, h]) => `${x}% ${50 - h}%`), ...profile.toReversed().map(([x, h]) => `${x}% ${50 + h}%`)].join(", ")})`);
id("ruler").insertAdjacentHTML("beforeend", [[-1e6, "−1 Myr"], [-1e4, "−10 kyr"], [-100, "−100 yr"], [-1, "−1 yr"], [0, "perihelion"], [1, "+1 yr"], [100, "+100 yr"], [1e4, "+10 kyr"], [1e6, "+1 Myr"]]
	.map(([years, text]) => `<span style="--at: ${tick(years * 365.25)}">${text}</span>`).join(""));

// --- input -----------------------------------------------------------------

// One way time moves: a journey of the slider to u, through its own log time, so crossing a million years passes through every scale on the way. A date link is a short eased journey; Play is a steady one to the end. Starting any journey, or a hand on the slider, ends the one before.
const GLIDE = 900;      // ms for a date link
const PACE = 27500;     // ms per unit of u while playing: the whole slider in under a minute
let journey = 0;

function travel(to, ms, ease, then = () => form.play.checked = false) {
	const from = +form.u.value;
	const start = performance.now();
	const mine = ++journey;
	requestAnimationFrame(function step(now) {
		if (mine != journey)
			return;
		const k = Math.min(1, Math.max(0, now - start) / ms);
		form.u.value = from + (to - from) * ease(k);
		timed();
		if (k < 1)
			requestAnimationFrame(step);
		else
			then();
	});
}

// Play: a steady journey to the end; at the end, from the beginning, as a player does.
function play() {
	if (+form.u.value == 1)
		form.u.value = -1;
	form.play.checked = true;
	travel(1, (1 - form.u.value) * PACE, steady);
}

// Until the viewer zooms, the zoom follows the visitor: r0 is set so the visitor is always drawn at SEEN from the Sun, a third of the way to the view's edge, inverting the one law (asinh, or linear when checked). Space then flows past it: the scale closes as it falls in, to the slider's floor at perihelion, and opens to parsecs as it leaves for the stars.
const SEEN = 5;
// Until the viewer turns the view, it turns with the visitor: its bearing swings by PAN times the visitor's angle from perihelion, measured in its orbit's plane (always within ±180°, so it never wraps), so a pass is seen from changing sides.
const PAN = 0.3;
let following = true;
let steering = true;
const HOME = camera.position.clone();
const home = Math.atan2(HOME.y, HOME.x);
controls.addEventListener("start", () => steering = false);
// Fit, as the thread browser has it: a glide back to the view that follows the visitor, re-engaging both the zoom and the pan.
let fitting = null;
function fit() {
	fitting = { start: performance.now(), from: camera.position.clone(), zoom: +form.zoom.value };
	following = steering = true;
	requestAnimationFrame(function step() {
		if (!form.play.checked)
			timed();
		if (fitting)
			requestAnimationFrame(step);
	});
}
const vec = (a) => new THREE.Vector3(...a);

function timed() {
	const v = visitor();
	const p = vec(visitorAt(v, +form.u.value));
	const k = fitting ? smooth(Math.min(1, (performance.now() - fitting.start) / GLIDE)) : 1;
	if (following) {
		const zoom = Math.log10(p.length() / (form.linear.checked ? SEEN : Math.sinh(SEEN)));
		form.zoom.value = fitting ? fitting.zoom + (zoom - fitting.zoom) * k : zoom;
	}
	if (steering) {
		const q = vec(visitorAt(v, 0));
		const normal = q.clone().cross(vec(visitorAt(v, 0.01))).normalize();
		const anomaly = Math.atan2(normal.dot(q.clone().cross(p)), q.dot(p));
		const bearing = home + PAN * anomaly;
		const flat = Math.hypot(HOME.x, HOME.y);
		const aim = new THREE.Vector3(flat * Math.cos(bearing), flat * Math.sin(bearing), HOME.z);
		camera.position.copy(fitting ? fitting.from.clone().lerp(aim, k).setLength(fitting.from.length() + (aim.length() - fitting.from.length()) * k) : aim);
		camera.lookAt(0, 0, 0);
	}
	if (k == 1)
		fitting = null;
	following ? zoomed() : moved();
}

const smooth = (k) => k < 0.5 ? 4 * k ** 3 : 1 - (2 - 2 * k) ** 3 / 2;
const steady = (k) => k;

function jump(jd) {
	form.play.checked = false;
	travel(slider(jd), GLIDE, smooth);
}

// A visit starts where the inbound path crosses ENTRY, just outside Neptune's orbit and the Kuiper belt: the slider position of the last sample before it.
const ENTRY = 50;
function arrival(v) {
	const n = v.path.length / 3;
	const k = Array.from({ length: n }, (_, k) => k).find((k) => Math.hypot(...v.path.subarray(3 * k, 3 * k + 3)) < ENTRY);
	return -1 + 2 * (k - 1) / (n - 1);
}

id("fit").addEventListener("click", fit);
document.addEventListener("click", (e) => {
	const button = e.target.closest("button[value]");
	if (button)
		jump(+button.value);
});
// A visitor is shown as a visit: a glide to its arrival, then on through the pass.
document.querySelector("aside").addEventListener("change", () => {
	zoomed();
	form.play.checked = true;
	travel(arrival(visitor()), GLIDE, smooth, play);
});
form.u.addEventListener("input", () => {
	journey++;
	form.play.checked = false;
	timed();
});
form.zoom.addEventListener("input", () => {
	following = false;
	zoomed();
});
form.linear.addEventListener("change", () => following ? timed() : zoomed());
form.play.addEventListener("change", () => {
	journey++;
	if (form.play.checked)
		play();
});
labels.domElement.addEventListener("wheel", (e) => {
	e.preventDefault();
	following = false;
	form.zoom.value = +form.zoom.value + e.deltaY * 0.002;
	zoomed();
}, { passive: false });
document.addEventListener("themechange", paint);

// Hover: the nearest star under the pointer, and how far it is from the Sun at the slider's time. The view is for turning and zooming; its links live in the panel.
const ray = new THREE.Raycaster();
ray.params.Points.threshold = 0.08;
labels.domElement.addEventListener("pointermove", (e) => {
	const box = labels.domElement.getBoundingClientRect();
	const x = e.clientX - box.left;
	const y = e.clientY - box.top;
	ray.setFromCamera(new THREE.Vector2(2 * x / box.width - 1, 1 - 2 * y / box.height), camera);
	const [hit] = ray.intersectObject(field);
	const tip = id("tip");
	tip.style.setProperty("--x", x);
	tip.style.setProperty("--y", y);
	if (!hit) {
		tip.textContent = "";
		return;
	}
	const at = new Float64Array(3);
	starAt(hit.index, time(), at, 0);
	const now = Math.hypot(...stars.pos.subarray(3 * hit.index, 3 * hit.index + 3)) / PC;
	tip.textContent = `${fmt(Math.hypot(...at) / PC)} pc from the Sun at this time, ${fmt(now)} pc now; ${named(hit.index)}`;
});

form.u.value = arrival(visitor());
timed();
paint();
play();
