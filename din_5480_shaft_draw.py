"""
DIN 5480 shaft cross-section DXF generator using ezdxf.

Target designation:
    DIN 5480 - W 30 x 1.5 x 18 x f8

Interpretation:
    W   = external spline shaft / Welle
    30  = reference diameter, mm
    1.5 = module, mm
    18  = number of teeth
    f8  = shaft spline tolerance class/deviation

Important:
    This creates a clean 2D cross-section draft/profile. It does not embed
    copyrighted DIN 5480 tolerance tables. The tip/root diameters below are
    editable defaults and should be verified against your DIN 5480 sheet or
    inspection/calculation source before manufacturing.

Install:
    pip install ezdxf

Run:
    python din5480_spline_shaft_ezdxf.py
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, radians, sin, sqrt
from pathlib import Path
from typing import List, Tuple

import ezdxf
from ezdxf import units

Point = Tuple[float, float]


@dataclass(frozen=True)
class DIN5480ExternalSpline:
    designation: str = "DIN 5480 - W 30 x 1.5 x 18 x f8"
    reference_diameter: float = 30.0
    module: float = 1.5
    teeth: int = 18
    pressure_angle_deg: float = 30.0
    tolerance_class: str = "f8"

    # Editable geometry defaults for drafting.
    # Common DIN 5480 external spline references often use tip dia near d - 0.2*m.
    # Root diameter is set as a conservative visual draft approximation only.
    tip_diameter: float = 29.70
    root_diameter: float = 25.80

    points_per_flank: int = 12
    tooth_thickness_fraction: float = 0.42

    def validate(self) -> None:
        if self.teeth < 3:
            raise ValueError("teeth must be >= 3")
        if not (0 < self.root_diameter < self.reference_diameter < self.tip_diameter + self.module):
            raise ValueError("Check root/reference/tip diameters")
        if self.module <= 0:
            raise ValueError("module must be positive")


def involute_xy(base_radius: float, t: float) -> Point:
    return (
        base_radius * (cos(t) + t * sin(t)),
        base_radius * (sin(t) - t * cos(t)),
    )


def rotate(p: Point, angle: float) -> Point:
    x, y = p
    c, s = cos(angle), sin(angle)
    return (x * c - y * s, x * s + y * c)


def polar(radius: float, angle: float) -> Point:
    return (radius * cos(angle), radius * sin(angle))


def angle_of(p: Point) -> float:
    from math import atan2

    return atan2(p[1], p[0])


def involute_t_for_radius(base_radius: float, radius: float) -> float:
    if radius <= base_radius:
        return 0.0
    return sqrt((radius / base_radius) ** 2 - 1.0)


def one_tooth(params: DIN5480ExternalSpline) -> List[Point]:
    z = params.teeth
    alpha = radians(params.pressure_angle_deg)

    rp = params.reference_diameter / 2.0
    rb = rp * cos(alpha)
    ra = params.tip_diameter / 2.0
    rf = params.root_diameter / 2.0

    half_pitch_angle = pi / z
    half_tooth_angle_at_reference = half_pitch_angle * params.tooth_thickness_fraction

    tp = involute_t_for_radius(rb, rp)
    pitch_involute_angle = angle_of(involute_xy(rb, tp))

    left_rotation = -half_tooth_angle_at_reference - pitch_involute_angle
    right_rotation = half_tooth_angle_at_reference + pitch_involute_angle

    ta = involute_t_for_radius(rb, ra)
    ts = [ta * i / max(1, params.points_per_flank - 1) for i in range(params.points_per_flank)]

    left_flank = [rotate(involute_xy(rb, t), left_rotation) for t in ts]
    right_flank = [
        rotate((involute_xy(rb, t)[0], -involute_xy(rb, t)[1]), right_rotation)
        for t in reversed(ts)
    ]

    left_root = polar(rf, angle_of(left_flank[0]))
    right_root = polar(rf, angle_of(right_flank[-1]))
    return [left_root, *left_flank, *right_flank, right_root]


def spline_outline(params: DIN5480ExternalSpline) -> List[Point]:
    tooth = one_tooth(params)
    pts: List[Point] = []
    for i in range(params.teeth):
        a = 2.0 * pi * i / params.teeth
        pts.extend(rotate(p, a) for p in tooth)
    pts.append(pts[0])
    return pts


def setup_doc() -> ezdxf.document.Drawing:
    doc = ezdxf.new("R2010", setup=True)
    doc.units = units.MM
    doc.header["$INSUNITS"] = units.MM

    layers = {
        "PROFILE": 7,
        "CENTER": 1,
        "REFERENCE": 8,
        "TEXT": 2,
    }
    for name, color in layers.items():
        if name not in doc.layers:
            doc.layers.add(name, color=color)

    if "CENTER" not in doc.linetypes:
        doc.linetypes.add("CENTER", pattern=[1.25, 0.75, -0.125, 0.125, -0.125])
    if "DASHED" not in doc.linetypes:
        doc.linetypes.add("DASHED", pattern=[0.6, 0.35, -0.25])

    return doc


def add_text(msp, text: str, xy: Point, height: float = 2.0) -> None:
    msp.add_text(text, dxfattribs={"height": height, "layer": "TEXT"}).set_placement(xy)


def add_cross_section(msp, params: DIN5480ExternalSpline) -> None:
    outline = spline_outline(params)
    msp.add_lwpolyline(outline, close=True, dxfattribs={"layer": "PROFILE", "lineweight": 35})

    # Reference, tip, and root circles for checking/inspection.
    circles = [
        (params.tip_diameter, "tip dia da"),
        (params.reference_diameter, "reference dia dB"),
        (params.root_diameter, "root dia df"),
    ]
    for diameter, label in circles:
        msp.add_circle(
            (0, 0),
            diameter / 2.0,
            dxfattribs={"layer": "REFERENCE", "linetype": "DASHED"},
        )
        add_text(msp, f"{label}: Ø{diameter:.2f}", (20, diameter / 2.0 - 3), 1.8)

    r = params.tip_diameter / 2.0 + 5
    msp.add_line((-r, 0), (r, 0), dxfattribs={"layer": "CENTER", "linetype": "CENTER"})
    msp.add_line((0, -r), (0, r), dxfattribs={"layer": "CENTER", "linetype": "CENTER"})

    add_text(msp, params.designation, (-24, -25), 2.5)
    add_text(msp, f"z={params.teeth}, m={params.module}, alpha={params.pressure_angle_deg:g} deg, class={params.tolerance_class}", (-24, -29), 2.0)
    add_text(msp, "Cross-section only. Verify da/df and tolerance values before manufacturing.", (-24, -33), 1.8)


def build_dxf(
    params: DIN5480ExternalSpline = DIN5480ExternalSpline(),
    output: str | Path = "DIN5480_W30x1_5x18xf8_cross_section.dxf",
) -> Path:
    params.validate()
    doc = setup_doc()
    msp = doc.modelspace()
    add_cross_section(msp, params)
    output = Path(output)
    doc.saveas(output)
    return output


if __name__ == "__main__":
    out = build_dxf()
    print(f"Wrote {out.resolve()}")
