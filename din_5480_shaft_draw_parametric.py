"""
Parametric DIN 5480 external spline shaft cross-section DXF generator.

Purpose:
    Build an external spline shaft cross-section from calculator/DIN geometry
    values, for example values from:
        https://ondrives.com/spline-calculator

Important:
    This script does NOT contain DIN 5480 tolerance tables.
    Enter the calculated values explicitly:
        - number of teeth z
        - module m
        - pressure angle alpha
        - shaft tip diameter da
        - shaft root diameter df
        - shaft root form diameter dFf
        - actual tooth thickness on pitch circle

Example used during validation:
    Shaft designation: DIN 5480 - W 30 x 1.5 x 18 x f8
    z  = 18
    m  = 1.5
    alpha = 30 deg
    d  = m*z = 27.00
    da = 29.70
    df = 26.49
    dFf = 26.73
    tooth thickness actual on pitch circle = 3.099

Output:
    Cross-section only.
    Thin profile line.
    Reference circles and two main centerlines.
    No text, no dimensions, no radial tooth-center lines.

Install:
    pip install ezdxf

Run default example:
    python din_5480_shaft_draw.py

Run with custom parameters:
    python din_5480_shaft_draw.py --nominal 30 --module 1.5 --teeth 18 \
        --tip 29.70 --root 26.49 --root-form 26.73 --tooth-thickness 3.099
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import acos, atan2, cos, pi, radians, sin, tan
from pathlib import Path
from typing import List, Tuple

import ezdxf
from ezdxf import units

Point = Tuple[float, float]


@dataclass(frozen=True)
class SplineShaftGeometry:
    """All dimensions are millimetres; angles are degrees except internal radians."""

    designation: str = "DIN 5480 - W 30 x 1.5 x 18 x f8"
    nominal_diameter: float = 30.0           # dB from DIN designation / calculator nominal diameter
    module: float = 1.5                      # m
    teeth: int = 18                          # z
    pressure_angle_deg: float = 30.0         # alpha

    # Shaft external spline values from calculator.
    tip_diameter: float = 29.70              # da, shaft max/min selected by user
    root_diameter: float = 26.49             # df
    root_form_diameter: float = 26.73        # dFf
    tooth_thickness_pitch_actual: float = 3.099

    # The working flank shape that matched the calculator during tuning.
    # This affects only the lower non-working start area near dFf and fades out at da.
    # Usually keep 0.020 unless you are calibrating to a different calculator view.
    lower_flank_relief_angle: float = 0.020

    # Fallback only if the exact tangent fillet construction fails.
    fallback_root_under_form_angle: float = 0.030

    # Drawing fidelity.
    points_per_flank: int = 56
    points_per_root_fillet: int = 24
    points_per_root_arc: int = 16
    points_per_tip_arc: int = 16

    # Drawing options.
    draw_reference_circles: bool = True
    draw_axes: bool = True

    @property
    def pitch_diameter(self) -> float:
        return self.module * self.teeth

    @property
    def base_diameter(self) -> float:
        return self.pitch_diameter * cos(radians(self.pressure_angle_deg))

    @property
    def circular_pitch(self) -> float:
        return pi * self.pitch_diameter / self.teeth

    @property
    def profile_shift_estimate(self) -> float:
        # Useful for checking DIN 5480 calculator logic; not required for drawing.
        return (self.nominal_diameter - self.pitch_diameter - 1.1 * self.module) / (2.0 * self.module)

    def validate(self) -> None:
        if self.teeth < 3:
            raise ValueError("teeth must be >= 3")
        if self.module <= 0:
            raise ValueError("module must be positive")
        if self.pressure_angle_deg <= 0 or self.pressure_angle_deg >= 60:
            raise ValueError("pressure angle looks invalid")
        if not (0 < self.base_diameter < self.root_form_diameter < self.pitch_diameter < self.tip_diameter):
            raise ValueError(
                "Expected base_diameter < root_form_diameter < pitch_diameter < tip_diameter. "
                "Check module, teeth, da, and dFf."
            )
        if not (0 < self.root_diameter < self.root_form_diameter):
            raise ValueError("Expected root_diameter < root_form_diameter")
        if not (0 < self.tooth_thickness_pitch_actual < self.circular_pitch):
            raise ValueError("tooth thickness must be between 0 and circular pitch")


def polar(radius: float, angle: float) -> Point:
    return radius * cos(angle), radius * sin(angle)


def angle_of(p: Point) -> float:
    return atan2(p[1], p[0])


def dist(a: Point, b: Point) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return (dx * dx + dy * dy) ** 0.5


def normalize(v: Point) -> Point:
    x, y = v
    length = (x * x + y * y) ** 0.5
    if length == 0:
        return 1.0, 0.0
    return x / length, y / length


def inv_alpha_from_radius(base_radius: float, radius: float) -> float:
    """Involute function inv(alpha)=tan(alpha)-alpha."""
    if radius <= base_radius:
        return 0.0
    alpha = acos(base_radius / radius)
    return tan(alpha) - alpha


def sample_arc_ccw(radius: float, a1: float, a2: float, n: int) -> List[Point]:
    while a2 < a1:
        a2 += 2.0 * pi
    if n <= 1:
        return [polar(radius, a2)]
    return [polar(radius, a1 + (a2 - a1) * i / (n - 1)) for i in range(n)]


def sample_arc_about_center(c: Point, radius: float, a1: float, a2: float, ccw: bool, n: int) -> List[Point]:
    if ccw:
        while a2 < a1:
            a2 += 2.0 * pi
    else:
        while a2 > a1:
            a2 -= 2.0 * pi
    return [
        (
            c[0] + radius * cos(a1 + (a2 - a1) * i / max(1, n - 1)),
            c[1] + radius * sin(a1 + (a2 - a1) * i / max(1, n - 1)),
        )
        for i in range(n)
    ]


def interp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def nominal_half_tooth_angle(g: SplineShaftGeometry, radius: float) -> float:
    """
    External involute tooth half-angle at radius.

    half_angle(r) = s/(2*rp) + inv(alpha_p) - inv(alpha_r)
    where s is actual circular tooth thickness at pitch diameter d=m*z.
    """
    rb = g.base_diameter / 2.0
    rp = g.pitch_diameter / 2.0
    half_pitch = g.tooth_thickness_pitch_actual / (2.0 * rp)
    return half_pitch + inv_alpha_from_radius(rb, rp) - inv_alpha_from_radius(rb, radius)


def relieved_half_tooth_angle(g: SplineShaftGeometry, radius: float) -> float:
    """
    Lower-flank relief used to match the calculator view.

    It fades from lower_flank_relief_angle at dFf to zero at da.
    """
    rff = g.root_form_diameter / 2.0
    ra = g.tip_diameter / 2.0
    t = (radius - rff) / (ra - rff)
    fade = 1.0 - smoothstep(t)
    return nominal_half_tooth_angle(g, radius) - g.lower_flank_relief_angle * fade


def tooth_angles(g: SplineShaftGeometry, radius: float, center_angle: float) -> Tuple[float, float]:
    half = relieved_half_tooth_angle(g, radius)
    if half <= 0:
        raise ValueError("lower_flank_relief_angle is too large; tooth half-angle became negative")
    return center_angle - half, center_angle + half


def tooth_points(g: SplineShaftGeometry, center_angle: float) -> Tuple[List[Point], float, float]:
    """Return one tooth profile from left root-form to right root-form."""
    rff = g.root_form_diameter / 2.0
    ra = g.tip_diameter / 2.0

    radii = [interp(rff, ra, i / (g.points_per_flank - 1)) for i in range(g.points_per_flank)]

    left_flank: List[Point] = []
    for r in radii:
        left, _ = tooth_angles(g, r, center_angle)
        left_flank.append(polar(r, left))

    right_flank: List[Point] = []
    for r in reversed(radii):
        _, right = tooth_angles(g, r, center_angle)
        right_flank.append(polar(r, right))

    left_tip, _ = tooth_angles(g, ra, center_angle)
    _, right_tip = tooth_angles(g, ra, center_angle)
    tip_arc = sample_arc_ccw(ra, left_tip, right_tip, g.points_per_tip_arc)[1:-1]

    left_root_form, _ = tooth_angles(g, rff, center_angle)
    _, right_root_form = tooth_angles(g, rff, center_angle)

    return [*left_flank, *tip_arc, *right_flank], left_root_form, right_root_form


def flank_point(g: SplineShaftGeometry, radius: float, center_angle: float, side: str) -> Point:
    left, right = tooth_angles(g, radius, center_angle)
    return polar(radius, left if side == "left" else right)


def flank_tangent_at_form(g: SplineShaftGeometry, center_angle: float, side: str) -> Point:
    """Numerical tangent of the working flank at dFf, toward the tip."""
    rff = g.root_form_diameter / 2.0
    eps = min(0.03, (g.tip_diameter - g.root_form_diameter) / 30.0)
    p0 = flank_point(g, rff, center_angle, side)
    p1 = flank_point(g, rff + eps, center_angle, side)
    return normalize((p1[0] - p0[0], p1[1] - p0[1]))


def fallback_root_transition(g: SplineShaftGeometry, center_angle: float, side: str) -> List[Point]:
    rf = g.root_diameter / 2.0
    rff = g.root_form_diameter / 2.0
    left_form, right_form = tooth_angles(g, rff, center_angle)
    form_angle = left_form if side == "left" else right_form
    half_form = right_form - center_angle
    half_root = half_form + g.fallback_root_under_form_angle
    root_angle = center_angle - half_root if side == "left" else center_angle + half_root
    pts = []
    for i in range(g.points_per_root_fillet):
        s = smoothstep(i / max(1, g.points_per_root_fillet - 1))
        pts.append(polar(interp(rf, rff, s), interp(root_angle, form_angle, s)))
    if side == "right":
        pts.reverse()
    return pts


def circular_root_fillet(g: SplineShaftGeometry, center_angle: float, side: str) -> List[Point]:
    """
    True circular root fillet.

    The fillet circle is tangent to:
      1) the root circle df;
      2) the working flank at root-form diameter dFf.

    The selected solution touches the root circle outside the tooth, toward the
    gap. This avoids cutting into the tooth base.
    """
    rf = g.root_diameter / 2.0
    rff = g.root_form_diameter / 2.0
    p = flank_point(g, rff, center_angle, side)
    tangent = flank_tangent_at_form(g, center_angle, side)

    candidates = []
    for normal in [(-tangent[1], tangent[0]), (tangent[1], -tangent[0])]:
        pdotn = p[0] * normal[0] + p[1] * normal[1]
        denom = 2.0 * (pdotn - rf)
        if abs(denom) < 1e-9:
            continue

        # Solve |p + rho*n| = rf + rho.
        rho = (rf * rf - rff * rff) / denom
        if rho <= 0:
            continue

        c = (p[0] + rho * normal[0], p[1] + rho * normal[1])
        dc = dist(c, (0.0, 0.0))
        if dc <= rf:
            continue

        q = (c[0] * rf / dc, c[1] * rf / dc)
        root_angle = angle_of(q)
        form_angle = angle_of(p)

        # Root tangency must lie outside the tooth, toward the gap.
        delta = ((root_angle - form_angle + pi) % (2.0 * pi)) - pi
        if side == "left" and delta >= 0:
            continue
        if side == "right" and delta <= 0:
            continue

        local_score = abs(delta)
        candidates.append((local_score, rho, c, q, p))

    if not candidates:
        return fallback_root_transition(g, center_angle, side)

    _, rho, c, q, p = min(candidates, key=lambda item: item[0])
    a_q = atan2(q[1] - c[1], q[0] - c[0])
    a_p = atan2(p[1] - c[1], p[0] - c[0])

    diff_ccw = (a_p - a_q) % (2.0 * pi)
    ccw = diff_ccw <= pi
    pts = sample_arc_about_center(c, rho, a_q, a_p, ccw, g.points_per_root_fillet)

    # Outline direction: left side root -> form; right side form -> root.
    if side == "right":
        pts.reverse()
    return pts


def spline_outline(g: SplineShaftGeometry) -> List[Point]:
    g.validate()

    z = g.teeth
    rf = g.root_diameter / 2.0

    data = []
    for i in range(z):
        center = 2.0 * pi * i / z
        pts, left_form, right_form = tooth_points(g, center)
        data.append((center, pts, left_form, right_form))

    outline: List[Point] = []

    for i in range(z):
        center, pts, _, _ = data[i]
        next_center, _, _, _ = data[(i + 1) % z]

        left_fillet = circular_root_fillet(g, center, "left")
        right_fillet = circular_root_fillet(g, center, "right")
        next_left_fillet = circular_root_fillet(g, next_center, "left")

        right_root_angle = angle_of(right_fillet[-1])
        next_left_root_angle = angle_of(next_left_fillet[0])

        outline.extend(left_fillet)
        outline.extend(pts[1:-1])
        outline.extend(right_fillet)
        outline.extend(sample_arc_ccw(rf, right_root_angle, next_left_root_angle, g.points_per_root_arc)[1:])

    outline.append(outline[0])
    return outline


def setup_doc() -> ezdxf.document.Drawing:
    doc = ezdxf.new("R2010", setup=True)
    doc.units = units.MM
    doc.header["$INSUNITS"] = units.MM

    if "PROFILE" not in doc.layers:
        doc.layers.add("PROFILE", color=7)
    if "REFERENCE" not in doc.layers:
        doc.layers.add("REFERENCE", color=8)
    if "AXIS" not in doc.layers:
        doc.layers.add("AXIS", color=8)

    if "CENTER" not in doc.linetypes:
        doc.linetypes.add("CENTER", pattern=[1.25, 0.75, -0.125, 0.125, -0.125])
    if "DASHED" not in doc.linetypes:
        doc.linetypes.add("DASHED", pattern=[0.6, 0.35, -0.25])

    return doc


def add_cross_section(msp, g: SplineShaftGeometry) -> None:
    msp.add_lwpolyline(spline_outline(g), close=True, dxfattribs={"layer": "PROFILE"})

    r_outer = g.tip_diameter / 2.0 + 2.0

    if g.draw_reference_circles:
        for dia in (
            g.tip_diameter,
            g.nominal_diameter,
            g.pitch_diameter,
            g.base_diameter,
            g.root_form_diameter,
            g.root_diameter,
        ):
            msp.add_circle(
                (0, 0),
                dia / 2.0,
                dxfattribs={"layer": "REFERENCE", "linetype": "DASHED"},
            )

    if g.draw_axes:
        msp.add_line((-r_outer, 0), (r_outer, 0), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
        msp.add_line((0, -r_outer), (0, r_outer), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})


def build_dxf(
    geometry: SplineShaftGeometry = SplineShaftGeometry(),
    output: str | Path = "DIN5480_spline_shaft_cross_section.dxf",
) -> Path:
    doc = setup_doc()
    add_cross_section(doc.modelspace(), geometry)
    output = Path(output)
    doc.saveas(output)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build DIN 5480 external spline shaft cross-section DXF.")
    parser.add_argument("--designation", default=SplineShaftGeometry.designation)
    parser.add_argument("--nominal", type=float, default=SplineShaftGeometry.nominal_diameter, help="Nominal/reference diameter dB")
    parser.add_argument("--module", type=float, default=SplineShaftGeometry.module, help="Module m")
    parser.add_argument("--teeth", type=int, default=SplineShaftGeometry.teeth, help="Number of teeth z")
    parser.add_argument("--pressure-angle", type=float, default=SplineShaftGeometry.pressure_angle_deg)
    parser.add_argument("--tip", type=float, default=SplineShaftGeometry.tip_diameter, help="Shaft tip diameter da")
    parser.add_argument("--root", type=float, default=SplineShaftGeometry.root_diameter, help="Shaft root diameter df")
    parser.add_argument("--root-form", type=float, default=SplineShaftGeometry.root_form_diameter, help="Shaft root form diameter dFf")
    parser.add_argument("--tooth-thickness", type=float, default=SplineShaftGeometry.tooth_thickness_pitch_actual, help="Actual shaft tooth thickness at pitch diameter")
    parser.add_argument("--lower-flank-relief", type=float, default=SplineShaftGeometry.lower_flank_relief_angle)
    parser.add_argument("--no-reference-circles", action="store_true")
    parser.add_argument("--no-axes", action="store_true")
    parser.add_argument("--output", default="DIN5480_spline_shaft_cross_section.dxf")
    return parser.parse_args()


def geometry_from_args(args: argparse.Namespace) -> SplineShaftGeometry:
    return SplineShaftGeometry(
        designation=args.designation,
        nominal_diameter=args.nominal,
        module=args.module,
        teeth=args.teeth,
        pressure_angle_deg=args.pressure_angle,
        tip_diameter=args.tip,
        root_diameter=args.root,
        root_form_diameter=args.root_form,
        tooth_thickness_pitch_actual=args.tooth_thickness,
        lower_flank_relief_angle=args.lower_flank_relief,
        draw_reference_circles=not args.no_reference_circles,
        draw_axes=not args.no_axes,
    )


if __name__ == "__main__":
    args = parse_args()
    geom = geometry_from_args(args)
    out = build_dxf(geom, args.output)
    print(f"Wrote {out.resolve()}")
    print(f"pitch diameter d = {geom.pitch_diameter:.6f}")
    print(f"base diameter db = {geom.base_diameter:.6f}")
    print(f"circular pitch = {geom.circular_pitch:.6f}")
    print(f"profile shift estimate x = {geom.profile_shift_estimate:.6f}")
