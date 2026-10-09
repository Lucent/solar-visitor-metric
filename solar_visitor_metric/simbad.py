"""Readable star names from SIMBAD, keyed by catalog designation ("Gaia DR3 …", "HIP …"), cached in data/."""
import csv
import os
import re
from urllib.parse import urlencode
from urllib.request import urlopen

TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
CHUNK = 1000
GREEK = dict(zip("alf bet gam del eps zet eta tet iot kap lam mu. nu. ksi omi pi. rho sig tau ups phi chi psi ome".split(),
				 "αβγδεζηθικλμνξοπρστυφχψω"))


def names(ids, path):
	"""Display name per designation, the most recognizable SIMBAD has: a proper name in words ("Barnard's star"; the shortest, where spellings vary in length), else a Bayer, Flamsteed or variable-star name ("α Cen A", "61 Cyg A", "EV Lac"), else a discoverer's catalog ("Ross 248", "Wolf 359"), else a Gliese number ("GJ 1061"), else the designation itself; survey acronyms ("HD", "UCAC4", "LP") read no better than the designation. Fetches only designations missing from the cache at `path`."""
	cache = {}
	if os.path.exists(path):
		cache = {r["id"]: r for r in csv.DictReader(open(path))}
	missing = [i for i in ids if i not in cache]
	for k in range(0, len(missing), CHUNK):
		cache |= fetch(missing[k:k + CHUNK])
	os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
	with open(path, "w", newline="") as f:
		out = csv.DictWriter(f, ["id", "main_id", "proper", "gj"])
		out.writeheader()
		out.writerows(cache.values())
	return {i: display(cache[i]) for i in ids}


def fetch(ids):
	"""SIMBAD rows for these designations; a designation SIMBAD does not know is cached as having no name."""
	quoted = ", ".join("'" + i.replace("'", "''") + "'" for i in ids)
	query = ("SELECT i.id, b.main_id, n.id FROM ident AS i JOIN basic AS b ON b.oid = i.oidref "
			 "LEFT JOIN ident AS n ON n.oidref = i.oidref AND (n.id LIKE 'NAME %' OR n.id LIKE 'GJ %') "
			 f"WHERE i.id IN ({quoted})")
	data = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}).encode()
	with urlopen(TAP, data=data, timeout=600) as r:
		body = r.read().decode()
	rows = list(csv.reader(body.splitlines()))
	if rows[0] != ["id", "main_id", "id"]:
		raise RuntimeError("unexpected SIMBAD response")
	found = {i: {"id": i, "main_id": "", "proper": "", "gj": ""} for i in ids}
	for key, main, other in rows[1:]:
		row = found[key]
		row["main_id"] = main
		if other.startswith("GJ "):
			if not row["gj"] or (len(other), other) < (len(row["gj"]), row["gj"]):
				row["gj"] = other
		elif other:
			row["proper"] = "|".join(sorted({*row["proper"].split("|"), other.removeprefix("NAME ")} - {""}))
	return found


def display(row):
	main = re.sub(r"\s+", " ", row["main_id"])
	# A proper name is words; "StM 187A" and "COCONUTS-2A" are catalog numbers filed as names. Two equally short spellings ("Rigil Kentaurus", "Rigel Kentaurus") leave no one name, and the designation decides.
	words = [n for n in row["proper"].split("|") + [main.removeprefix("NAME ")] * main.startswith("NAME ") if re.fullmatch(r"[A-Za-z' ]+", n)]
	shortest = [n for n in set(words) if len(n) == min(map(len, words))] if words else []
	if main.startswith("NAME ") and main.removeprefix("NAME ") in words:
		return main.removeprefix("NAME ")
	if len(shortest) == 1:
		return shortest[0]
	if re.match(r"[A-Z][a-z]+ \d", main):
		return main
	# "* " and "V* " mark Bayer, Flamsteed and variable-star names; "** " a double-star catalog number
	if not main.startswith(("* ", "V* ")):
		return row["gj"] or row["id"]
	main = main.removeprefix("V* ").removeprefix("* ")
	# Bayer letters as Greek: "alf Cen A" -> "α Cen A", "alf02 Lib" -> "α² Lib"
	return re.sub(r"^(\w\w\w?\.?)(\d*)(?= )",
				  lambda m: GREEK.get(m[1], m[1]) + m[2].lstrip("0").translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")) if m[1] in GREEK else m[0],
				  main)
