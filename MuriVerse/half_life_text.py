"""Read the qualifier off PEPlife2's free-text half-life value.

`Half Life` is free text: "97.3 +/- 3.2 (T1/2 Elimination)", "~4", "<5", "1 to 4", "2-3". The
numeric value in hours (`half_life_h`) is parsed upstream; this module only labels how that
number was obtained, so the training table says which values are exact, which are the stated
value of an approximation or one-sided bound, and which are the median (midpoint) of a range.

Units: PEPlife2 writes "Hours-1" for a handful of rows whose value is a T1/2 in hours (the
unit is a typo for hours, not a rate constant); every such row is `T1/2 Absorption` and is
removed by the elimination-half-life gate, but the unit is normalised here so nothing
downstream can mistake it for a per-hour rate.
"""

from __future__ import annotations

import re

import numpy as np

NUM = r"(\d+(?:\.\d+)?)"

# Unit text (lower-cased) -> the one canonical spelling for it. "hours-1" is a typo for
# hours. No multiplier to hours lives here: `half_life_h` is converted upstream, and this
# module only labels the value it is handed.
UNIT_NAMES = [
    (r"^s(ec(ond)?s?)?$", "seconds"),
    (r"^min(ute)?s?$", "minutes"),
    (r"^h(ou)?rs?(-1)?$", "hours"),
    (r"^d(ay)?s?$", "days"),
    (r"^w(eek)?s?$", "weeks"),
]


def normalise_unit(units: object) -> str:
    """'Hours-1' is hours; every other recognised spelling maps to one canonical unit."""
    text = str(units).strip().lower()
    for pattern, name in UNIT_NAMES:
        if re.match(pattern, text):
            return name
    return str(units)


def half_life_qualifier(text: object, units: object) -> str:
    """How the numeric half-life was obtained from `text`.

    exact | approximate | upper bound | lower bound | range (median of the endpoints)
    """
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return "exact"
    t = re.sub(r"\([^)]*\)\s*$", "", str(text)).strip()
    t = t.replace("\u2013", "-").replace("\u2212", "-").replace(",", "")
    t = re.sub(r"\u20ac\u201c", "-", t)                   # mojibake en dash
    if re.search(rf"{NUM}\s*(?:to|-)\s*{NUM}", t):
        return "range (median of endpoints)"
    # An Excel date written over a range ("02-Mar" for 2-3) is still a range.
    if re.fullmatch(r"\d{1,2}-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)", t, re.I):
        return "range (median of endpoints)"
    match = re.search(NUM, t)
    head = t[:match.start()].strip() if match else t
    if head.startswith("~") or "approx" in head.lower():
        return "approximate"
    if head.startswith(("<", "\u2264")):
        return "upper bound"
    if head.startswith((">", "\u2265")):
        return "lower bound"
    return "exact"
