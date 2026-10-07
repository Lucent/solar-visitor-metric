import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { CSS2DRenderer, CSS2DObject } from "three/addons/renderers/CSS2DRenderer.js";
import { Line2 } from "three/addons/lines/Line2.js";
import { LineGeometry } from "three/addons/lines/LineGeometry.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
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

// One law for space and time. A distance r is drawn at asinh(r / r0): to scale inside the zoom radius r0, one unit per factor of e beyond it. The slider u runs time the same way, t - t_peri = tau sinh(u asinh(span / tau)).
const r0 = () => 10 ** form.zoom.value;
const visitor = () => visitors[form.visitor.value];
const time = () => visitor().tp + data.time.tau * Math.sinh(form.u.value * SPAN);
const slider = (jd) => Math.asinh((jd - visitor().tp) / data.time.tau) / SPAN;

function mapped(x, y, z, out, i) {
	const r = Math.hypot(x, y, z);
	const k = Math.asinh(r / r0()) / r;
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

function label(text, kind) {
	const el = document.createElement("div");
	el.className = `label ${kind}`;
	el.textContent = text;
	return new CSS2DObject(el);
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

// --- scene ---------------------------------------------------------------

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 500);
camera.up.set(0, 0, 1);
camera.position.set(0, -34, 22);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setPixelRatio(devicePixelRatio);
const labels = new CSS2DRenderer();
id("scene").append(renderer.domElement, labels.domElement);
const controls = new OrbitControls(camera, labels.domElement);
controls.enableZoom = false;
controls.enablePan = false;
controls.addEventListener("change", render);

const materials = {
	grid: new THREE.LineBasicMaterial(),
	orbit: new THREE.LineBasicMaterial(),
	planet: new THREE.PointsMaterial({ size: 6, sizeAttenuation: false }),
	sun: new THREE.PointsMaterial({ size: 10, sizeAttenuation: false }),
	star: new THREE.PointsMaterial({ size: 2, sizeAttenuation: false }),
	from: new THREE.PointsMaterial({ size: 7, sizeAttenuation: false }),
	to: new THREE.PointsMaterial({ size: 7, sizeAttenuation: false }),
	fromTrack: new THREE.LineDashedMaterial({ dashSize: 0.15, gapSize: 0.1 }),
	toTrack: new THREE.LineDashedMaterial({ dashSize: 0.15, gapSize: 0.1 }),
	visitor: new THREE.PointsMaterial({ size: 9, sizeAttenuation: false }),
	path: new LineMaterial({ linewidth: 2 }),
};

// Colors live in the stylesheet; the scene reads them, and again when the scheme flips.
function paint() {
	const css = getComputedStyle(document.documentElement);
	const color = (name) => new THREE.Color(css.getPropertyValue(`--${name}`).trim());
	Object.entries({ grid: "grid", orbit: "muted", planet: "ink2", sun: "sun", star: "muted", from: "s2", to: "s3", fromTrack: "s2", toTrack: "s3", visitor: "s1", path: "s1" })
		.forEach(([m, c]) => materials[m].color = color(c));
	render();
}

// Decades of distance as rings in the ecliptic. A sphere about the Sun maps to a sphere, so a ring is a unit circle scaled by asinh(R / r0).
const circle = new THREE.BufferGeometry().setFromPoints(
	Array.from({ length: 128 }, (_, k) => new THREE.Vector3(Math.cos(k * Math.PI / 64), Math.sin(k * Math.PI / 64), 0)));
const shells = [[1, "1 AU"], [10, "10 AU"], [100, "100 AU"], [1e3, "1,000 AU"], [1e4, "10,000 AU"], [PC, "1 parsec"], [10 * PC, "10 parsecs"], [100 * PC, "100 parsecs"]]
	.map(([R, text]) => {
		const ring = new THREE.LineLoop(circle, materials.grid);
		const tag = label(text, "shell");
		tag.position.set(Math.SQRT1_2, -Math.SQRT1_2, 0);
		ring.add(tag);
		ring.userData.R = R;
		scene.add(ring);
		return ring;
	});

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

const path = new Line2(new LineGeometry(), materials.path);
const traveler = new THREE.Points(buffer(1), materials.visitor);
const travelerLabel = label("", "visitor");
scene.add(path, traveler, travelerLabel);

// Each end's best stars, each with its track from now to where the visitor meets it.
const TRACK = 32;
const ends = Object.fromEntries(["from", "to"].map((end) => {
	const best = new THREE.Points(buffer(5), materials[end]);
	const tracks = new THREE.LineSegments(buffer(5 * 2 * (TRACK - 1)), materials[`${end}Track`]);
	const tag = label("", end);
	scene.add(best, tracks, tag);
	return [end, { best, tracks, tag }];
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

function zoomed() {
	const v = visitor();
	shells.forEach((ring) => ring.scale.setScalar(Math.asinh(ring.userData.R / r0())));
	orbits.forEach((loop) => place(loop.geometry, loop.userData.xyz));
	const xyz = new Float32Array(v.path.length);
	for (let i = 0; i < xyz.length; i += 3)
		mapped(v.path[i], v.path[i + 1], v.path[i + 2], xyz, i);
	path.geometry.setPositions(xyz);
	Object.entries(ends).forEach(([end, { tracks }]) => {
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

function moved() {
	const v = visitor();
	const jd = time();
	form.date.value = when(jd);
	const p = visitorAt(v, +form.u.value);
	place(traveler.geometry, p);
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
		planetLabels.forEach((tag, k) => tag.position.fromArray(planets.geometry.getAttribute("position").array, 3 * k));
	}
	const xyz = new Float64Array(stars.pos.length);
	for (let i = 0; i < xyz.length / 3; i++)
		starAt(i, jd, xyz, 3 * i);
	place(field.geometry, xyz);
	Object.entries(ends).forEach(([end, { best, tag }]) => {
		const shown = v.route[end].stars;
		const at = new Float64Array(15);
		shown.forEach((s, j) => starAt(s.i, jd, at, 3 * j));
		place(best.geometry, at);
		best.geometry.setDrawRange(0, shown.length);
		tag.element.textContent = `${end} ${fmt(shown[0].bits)} bits`;
		tag.position.fromArray(best.geometry.getAttribute("position").array);
	});
	render();
}

function render() {
	renderer.render(scene, camera);
	labels.render(scene, camera);
}

new ResizeObserver(([{ contentRect: { width, height } }]) => {
	renderer.setSize(width, height, false);
	labels.setSize(width, height);
	materials.path.resolution.set(width, height);
	camera.aspect = width / height;
	camera.updateProjectionMatrix();
	render();
}).observe(id("scene"));

// --- the panels: built once; from then on the radios are the state ------------

function bits(b) {
	return `<data value="${b}">${fmt(b)} bits</data>`;
}

function stop(s) {
	return `<li style="--bits: ${s.bits}"><span>Gaia DR3 ${data.stars.ids[s.i]}: missed by ${fmt(s.miss)} ± ${fmt(s.sigma)} pc, ${fmt(s.hop)} pc from the Sun, <button type="button" value="${s.jd}">${when(s.jd)}</button></span></li>`;
}

function panel(v) {
	const a = v.assist;
	return `<input type="radio" name="visitor" form="controls" id="v${v.tag}" value="${v.tag}"><label for="v${v.tag}">${v.tag}</label>
<article>
<h2>${v.name}</h2>
<p>${fmt(v.vinf, 3)} km/s before it fell toward the Sun; perihelion <button type="button" value="${v.tp}">${when(v.tp)}</button> at ${fmt(v.q, 3)} AU.</p>
<h2>Planet passes ${bits(v.visit)}</h2>
<p>How often a randomly placed planet would have come this close, and the surprise in bits. With eight planets to be near, the combined score corrects for having looked at all of them.</p>
<table>
${v.planets.map((p) => `<tr style="--bits: ${-Math.log2(p.p)}"><th>${p.name}</th><td>${fmt(p.d)} AU</td><td>${fmt(p.hill)} Hill</td><td><button type="button" value="${p.jd}">${when(p.jd)}</button></td><td>${pct(p.p)}</td></tr>`).join("\n")}
</table>
<h2>Gravity assist <data value="${a.total}">${pct(a.total)}</data></h2>
<p>The planets changed its outbound velocity by ${pct(a.total)} of its speed (their direct pull ${pct(a.direct)}, the Sun's wobble ${pct(a.reflex)}). A flyby that steers a spacecraft changes it by 10% or more.</p>
<h2>Star-to-star route ${bits(v.route.bits)}</h2>
<p>${pct(v.route.beaten)} of random visitors at this speed pass the nearby stars at least as tellingly. The route is scored as one itinerary: where it came from, then where the Sun's turn sent it, given that arrival. Its best candidates, each corrected for the stars it could have tried:</p>
${["from", "to"].map((end) => `<p>${end == "from" ? "Came from" : `Heading to, given the arrival and the Sun's ${fmt(v.route.turn)}° turn`}: ${bits(v.route[end].stars[0].bits)} among ${v.route[end].candidates.toLocaleString()} stars it could meet within 20 pc of the Sun.</p>
<ol class="${end}">${v.route[end].stars.map(stop).join("")}</ol>`).join("\n")}
</article>`;
}

document.querySelector("aside").innerHTML = data.visitors.map(panel).join("\n");
document.querySelector(`#v${data.visitors.at(-1).tag}`).checked = true;
id("now").value = NOW;

// --- input -----------------------------------------------------------------

function jump(jd) {
	form.u.value = slider(jd);
	moved();
}

function play() {
	if (!form.play.checked)
		return;
	form.u.value = +form.u.value + 0.0006;
	if (+form.u.value >= 1)
		form.play.checked = false;
	moved();
	requestAnimationFrame(play);
}

document.addEventListener("click", (e) => {
	const button = e.target.closest("button[value]");
	if (button)
		jump(+button.value);
});
document.querySelector("aside").addEventListener("change", zoomed);
form.u.addEventListener("input", moved);
form.zoom.addEventListener("input", zoomed);
form.play.addEventListener("change", () => requestAnimationFrame(play));
labels.domElement.addEventListener("wheel", (e) => {
	e.preventDefault();
	form.zoom.value = +form.zoom.value + e.deltaY * 0.002;
	zoomed();
}, { passive: false });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", paint);

// Hover: the nearest star under the pointer, and how far it is from the Sun at the slider's time.
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
	tip.textContent = `Gaia DR3 ${data.stars.ids[hit.index]}: ${fmt(Math.hypot(...at) / PC)} pc from the Sun then, ${fmt(now)} pc now`;
});

jump(NOW);
zoomed();
paint();
