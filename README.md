# DIN 5480 Parametric Spline Shaft DXF Generator

Parametric external spline shaft cross-section generator for **DIN 5480** splines using **Python + ezdxf**.

The script generates a 2D DXF cross-section of an external involute spline shaft.

The geometry was tuned and verified visually against the Ondrives spline calculator:

https://ondrives.com/spline-calculator

---

## Features

* External spline shaft geometry (`W` profile)
* Parametric DIN 5480 input
* DXF output via `ezdxf`
* Involute tooth profile
* Tangent circular root fillets
* Reference circles
* Center lines
* Thin CAD-friendly geometry
* No dimensions or annotations
* Suitable for visual fit/interference checks

---

## Example

Validated example:

```text
DIN 5480 - W 42 x 2 x 20 x f8
```

Parameters:

```text
z   = 20
m   = 2.0
α   = 30°
d   = 40.00
da  = 41.60
df  = 37.30
dFf = 37.64
```

Run:

```bash
python din_5480_shaft_draw_parametric.py \
    --nominal 42 \
    --module 2 \
    --teeth 20 \
    --tip 41.60 \
    --root 37.30 \
    --root-form 37.64 \
    --tooth-thickness 2.9
```

Output:

```text
DIN5480_spline_shaft_cross_section.dxf
```

---

## Installation

Create environment:

```bash
python -m venv .venv
```

Activate:

### Windows PowerShell

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
pip install ezdxf
```

---

## Command Line Parameters

| Parameter           | Description                              |
| ------------------- | ---------------------------------------- |
| `--nominal`         | Nominal/reference diameter               |
| `--module`          | Module `m`                               |
| `--teeth`           | Number of teeth `z`                      |
| `--pressure-angle`  | Pressure angle                           |
| `--tip`             | Shaft tip diameter `da`                  |
| `--root`            | Shaft root diameter `df`                 |
| `--root-form`       | Shaft root-form diameter `dFf`           |
| `--tooth-thickness` | Actual tooth thickness at pitch diameter |
| `--output`          | Output DXF filename                      |

---

## Notes

This project does **not** contain DIN 5480 standard tables.

All dimensions should be taken from:

* DIN 5480 documentation
* your CAD system
* inspection reports
* spline calculators such as Ondrives

The script is intended for:

* CAD automation
* spline visualization
* geometry experiments
* interference checking
* DXF generation

It is **not** a certified manufacturing calculation tool.

---

## Geometry Notes

The generated spline uses:

* involute working flanks
* tangent circular root fillets
* root-form diameter transition
* actual tooth thickness values

The root fillet geometry was iteratively tuned to avoid:

* tooth overbuild
* impossible sharp corners
* inward-cutting fillets
* non-tangent spline transitions

---

## Dependencies

* Python 3.10+
* ezdxf

https://github.com/mozman/ezdxf

---

## License

MIT License
