// Ptree's theme control. The checkbox reflects the displayed theme. Clicks save only a choice that differs from the system; only "Dark" and "Light" override live system changes. A "themechange" event tells the scene to read its colors again.
const toggle = document.getElementById("dark");
const system = matchMedia("(prefers-color-scheme: dark)");

function applyTheme(dark) {
	toggle.checked = dark;
	document.documentElement.classList.toggle("Dark", dark);
	document.dispatchEvent(new Event("themechange"));
}

toggle.checked = document.documentElement.classList.contains("Dark");
toggle.addEventListener("change", () => {
	const dark = toggle.checked;
	if (dark === system.matches)
		localStorage.removeItem("Theme");
	else
		localStorage.setItem("Theme", dark ? "Dark" : "Light");
	applyTheme(dark);
});
system.addEventListener("change", () => {
	if (!["Dark", "Light"].includes(localStorage.getItem("Theme")))
		applyTheme(system.matches);
});
