"""
DIN 5480 external spline shaft cross-section DXF generator.

Target:
    DIN 5480 - W 30 x 1.5 x 18 x f8

Calculator values used, shaft MAX actual:
    z  = 18
    m  = 1.5
    alpha = 30 deg
    d  = m*z = 27.00
    da = 29.70
    df = 26.49
    dFf = 26.73
    tooth thickness actual on pitch circle = 3.099

This version keeps the working tooth flank unchanged. Only the non-working
root transition between df and dFf is changed: it is built as a TRUE circular
fillet tangent to the root circle and tangent to the working flank at dFf.

Output:
    Cross-section only.
    Thin profile line.
    Reference circles and two main centerlines.
    No text, no dimensions, no radial tooth-center lines.

Install:
    pip install ezdxf
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, atan2, cos, pi, radians, sin, tan
from pathlib import Path
from typing import List, Tuple

import ezdxf
from ezdxf import units

Point = Tuple[float, float]


@dataclass(frozen=True)
class DIN5480ExternalSpline:
    reference_diameter: float = 30.0
    module: float = 1.5
    teeth: int = 18
    pressure_angle_deg: float = 30.0

    # Ondrives shaft MAX actual values
    tip_diameter: float = 29.70
    root_diameter: float = 26.49
    root_form_diameter: float = 26.73
    tooth_thickness_pitch_actual: float = 3.099

    # Fallback only. The normal path uses tangent circular fillets.
    root_under_form_angle: float = 0.030

    points_per_flank: int = 48
    points_per_root_transition: int = 20
    points_per_arc: int = 14
    draw_reference_circles: bool = True
    draw_axes: bool = True

    @property
    def pitch_diameter(self) -> float:
        return self.module * self.teeth

    @property
    def base_diameter(self) -> float:
        return self.pitch_diameter * cos(radians(self.pressure_angle_deg))

    def validate(self) -> None:
        if self.teeth < 3:
            raise ValueError("teeth must be >= 3")
        if self.module <= 0:
            raise ValueError("module must be positive")
        if not (0 < self.base_diameter < self.root_form_diameter < self.pitch_diameter < self.tip_diameter):
            raise ValueError("Expected db < dFf < d < da")
        if not (0 < self.root_diameter < self.root_form_diameter):
            raise ValueError("Expected df < dFf")
        circular_pitch = pi * self.pitch_diameter / self.teeth
        if not (0 < self.tooth_thickness_pitch_actual < circular_pitch):
            raise ValueError("tooth thickness must be less than circular pitch")


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


def nominal_half_tooth_angle(params: DIN5480ExternalSpline, radius: float) -> float:
    """
    External involute tooth half-angle at radius.

    half_angle(r) = s/(2*rp) + inv(alpha_p) - inv(alpha_r)
    where s is actual circular tooth thickness at pitch diameter d=m*z.
    """
    rb = params.base_diameter / 2.0
    rp = params.pitch_diameter / 2.0
    half_pitch = params.tooth_thickness_pitch_actual / (2.0 * rp)
    return half_pitch + inv_alpha_from_radius(rb, rp) - inv_alpha_from_radius(rb, radius)


def relieved_half_tooth_angle(params: DIN5480ExternalSpline, radius: float) -> float:
    """
    This is the previously accepted working tooth shape.

    It narrows only the non-working lower flank region near dFf and fades to
    the nominal involute at the tip. Do not change unless the tooth flank itself
    needs to be recalibrated.
    """
    root_relief_angle = 0.020
    rff = params.root_form_diameter / 2.0
    ra = params.tip_diameter / 2.0
    t = (radius - rff) / (ra - rff)
    fade = 1.0 - smoothstep(t)
    return nominal_half_tooth_angle(params, radius) - root_relief_angle * fade


def tooth_angles(params: DIN5480ExternalSpline, radius: float, center_angle: float) -> Tuple[float, float]:
    half = relieved_half_tooth_angle(params, radius)
    if half <= 0:
        raise ValueError("Root relief is too large; tooth half-angle became negative")
    return center_angle - half, center_angle + half


def tooth_points(params: DIN5480ExternalSpline, center_angle: float) -> Tuple[List[Point], float, float]:
    """Return one tooth profile from left root-form to right root-form."""
    rff = params.root_form_diameter / 2.0
    ra = params.tip_diameter / 2.0

    radii = [interp(rff, ra, i / (params.points_per_flank - 1)) for i in range(params.points_per_flank)]

    left_flank: List[Point] = []
    for r in radii:
        left, _ = tooth_angles(params, r, center_angle)
        left_flank.append(polar(r, left))

    right_flank: List[Point] = []
    for r in reversed(radii):
        _, right = tooth_angles(params, r, center_angle)
        right_flank.append(polar(r, right))

    left_tip, _ = tooth_angles(params, ra, center_angle)
    _, right_tip = tooth_angles(params, ra, center_angle)
    tip_arc = sample_arc_ccw(ra, left_tip, right_tip, params.points_per_arc)[1:-1]

    left_root_form, _ = tooth_angles(params, rff, center_angle)
    _, right_root_form = tooth_angles(params, rff, center_angle)

    return [*left_flank, *tip_arc, *right_flank], left_root_form, right_root_form


def flank_point(params: DIN5480ExternalSpline, radius: float, center_angle: float, side: str) -> Point:
    left, right = tooth_angles(params, radius, center_angle)
    return polar(radius, left if side == "left" else right)


def flank_tangent_at_form(params: DIN5480ExternalSpline, center_angle: float, side: str) -> Point:
    """Numerical tangent of the unchanged working flank at dFf, toward the tip."""
    rff = params.root_form_diameter / 2.0
    eps = min(0.03, (params.tip_diameter - params.root_form_diameter) / 30.0)
    p0 = flank_point(params, rff, center_angle, side)
    p1 = flank_point(params, rff + eps, center_angle, side)
    return normalize((p1[0] - p0[0], p1[1] - p0[1]))


def circular_root_fillet(
    params: DIN5480ExternalSpline,
    center_angle: float,
    side: str,
) -> List[Point]:
    """
    True circular root fillet.

    The fillet circle is tangent to:
      1) the root circle df;
      2) the unchanged working flank at root-form diameter dFf.

    This removes the spline-like double inflection while preserving the tooth
    flank above dFf.
    """
    rf = params.root_diameter / 2.0
    rff = params.root_form_diameter / 2.0
    p = flank_point(params, rff, center_angle, side)
    tangent = flank_tangent_at_form(params, center_angle, side)

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

        q = (c[0] * rf / dc, c[1] * rf / dc)  # tangent point on df
        root_angle = angle_of(q)
        form_angle = angle_of(p)

        # Critical topology rule:
        # the root tangency must lie OUTSIDE the tooth, toward the gap.
        # Left flank: root angle must be smaller than form angle.
        # Right flank: root angle must be larger than form angle.
        # Otherwise the fillet cuts into the tooth.
        delta = ((root_angle - form_angle + pi) % (2.0 * pi)) - pi
        if side == "left" and delta >= 0:
            continue
        if side == "right" and delta <= 0:
            continue

        # Prefer the smaller local outside fillet.
        local_score = abs(delta)
        candidates.append((local_score, rho, c, q, p))

    if not candidates:
        # Fallback: smooth polar transition, only if tangent construction fails.
        left_form, right_form = tooth_angles(params, rff, center_angle)
        form_angle = left_form if side == "left" else right_form
        half_form = right_form - center_angle
        half_root = half_form + params.root_under_form_angle
        root_angle = center_angle - half_root if side == "left" else center_angle + half_root
        pts = []
        for i in range(params.points_per_root_transition):
            s = smoothstep(i / max(1, params.points_per_root_transition - 1))
            pts.append(polar(interp(rf, rff, s), interp(root_angle, form_angle, s)))
        if side == "right":
            pts.reverse()
        return pts

    _, rho, c, q, p = min(candidates, key=lambda item: item[0])
    a_q = atan2(q[1] - c[1], q[0] - c[0])
    a_p = atan2(p[1] - c[1], p[0] - c[0])

    diff_ccw = (a_p - a_q) % (2.0 * pi)
    ccw = diff_ccw <= pi
    pts = sample_arc_about_center(c, rho, a_q, a_p, ccw, params.points_per_root_transition)

    # Outline direction: left side root -> form; right side form -> root.
    if side == "right":
        pts.reverse()
    return pts


def spline_outline(params: DIN5480ExternalSpline) -> List[Point]:
    params.validate()

    z = params.teeth
    rf = params.root_diameter / 2.0

    data = []
    for i in range(z):
        center = 2.0 * pi * i / z
        pts, left_form, right_form = tooth_points(params, center)
        data.append((center, pts, left_form, right_form))

    outline: List[Point] = []

    for i in range(z):
        center, pts, left_form, right_form = data[i]
        next_center, _, _, _ = data[(i + 1) % z]

        left_fillet = circular_root_fillet(params, center, "left")
        right_fillet = circular_root_fillet(params, center, "right")
        next_left_fillet = circular_root_fillet(params, next_center, "left")

        right_root_angle = angle_of(right_fillet[-1])
        next_left_root_angle = angle_of(next_left_fillet[0])

        outline.extend(left_fillet)
        outline.extend(pts[1:-1])
        outline.extend(right_fillet)
        outline.extend(sample_arc_ccw(rf, right_root_angle, next_left_root_angle, params.points_per_arc)[1:])

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


def add_cross_section(msp, params: DIN5480ExternalSpline) -> None:
    msp.add_lwpolyline(spline_outline(params), close=True, dxfattribs={"layer": "PROFILE"})

    r_outer = params.tip_diameter / 2.0 + 2.0

    if params.draw_reference_circles:
        for dia in (
            params.tip_diameter,
            params.reference_diameter,
            params.pitch_diameter,
            params.base_diameter,
            params.root_form_diameter,
            params.root_diameter,
        ):
            msp.add_circle(
                (0, 0),
                dia / 2.0,
                dxfattribs={"layer": "REFERENCE", "linetype": "DASHED"},
            )

    if params.draw_axes:
        msp.add_line((-r_outer, 0), (r_outer, 0), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
        msp.add_line((0, -r_outer), (0, r_outer), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})


def build_dxf(
    params: DIN5480ExternalSpline = DIN5480ExternalSpline(),
    output: str | Path = "DIN5480_W30x1_5x18xf8_cross_section.dxf",
) -> Path:
    doc = setup_doc()
    add_cross_section(doc.modelspace(), params)
    output = Path(output)
    doc.saveas(output)
    return output


if __name__ == "__main__":
    out = build_dxf()
    print(f"Wrote {out.resolve()}")
