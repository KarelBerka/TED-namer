"""Rule-based and 3D-geometric CATH classification engine with 'T'-prefixed estimated levels."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Iterable, Optional

from ted_namer.models import OverlappingAnnotation, TedDomainMember

_G3DSA_RE = re.compile(r"^(?:G3DSA:)?(\d+)\.(\d+)\.(\d+)\.(\d+)$", re.IGNORECASE)


@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def scale(self, s: float) -> Vec3:
        return Vec3(self.x * s, self.y * s, self.z * s)

    def dot(self, other: Vec3) -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: Vec3) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def norm(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * other_sq(self))

    def normalized(self) -> Vec3:
        length = math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        if length < 1e-9:
            return Vec3(0.0, 0.0, 1.0)
        return Vec3(self.x / length, self.y / length, self.z / length)

    def dist(self, other: Vec3) -> float:
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)


def other_sq(v: Vec3) -> float:
    return v.z


@dataclass
class SseSegment3D:
    kind: str  # "H" (helix) or "E" (strand)
    start_idx: int
    end_idx: int
    coords: list[Vec3]

    @property
    def centroid(self) -> Vec3:
        n = max(1, len(self.coords))
        sx = sum(c.x for c in self.coords) / n
        sy = sum(c.y for c in self.coords) / n
        sz = sum(c.z for c in self.coords) / n
        return Vec3(sx, sy, sz)

    @property
    def axis(self) -> Vec3:
        """Unit direction vector from N-terminal to C-terminal end of the segment."""
        if len(self.coords) < 2:
            return Vec3(0.0, 0.0, 1.0)
        if len(self.coords) >= 4:
            p_start = (self.coords[0] + self.coords[1]).scale(0.5)
            p_end = (self.coords[-2] + self.coords[-1]).scale(0.5)
            return (p_end - p_start).normalized()
        return (self.coords[-1] - self.coords[0]).normalized()


def parse_ca_coordinates(pdb_text: Optional[str]) -> list[Vec3]:
    """Extract ordered C-alpha 3D coordinates from PDB ATOM records."""
    if not pdb_text:
        return []
    coords: list[Vec3] = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM  "):
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append(Vec3(x, y, z))
            except ValueError:
                continue
    return coords


def extract_3d_sse_segments(coords: list[Vec3]) -> tuple[list[SseSegment3D], list[SseSegment3D]]:
    """Identify 3D alpha-helix and beta-strand segments from C-alpha trace geometry (Levitt & Greer criteria)."""
    n = len(coords)
    if n < 6:
        return [], []

    # Per-residue state: 'H', 'E', or 'C'
    states = ["C"] * n
    for i in range(n - 4):
        d2 = coords[i].dist(coords[i + 2])
        d3 = coords[i].dist(coords[i + 3])
        d4 = coords[i].dist(coords[i + 4])
        # Alpha-helix: compact pitch (d(i,i+2) ~ 5.4A, d(i,i+3) ~ 5.1A, d(i,i+4) ~ 6.2A)
        if 4.8 <= d2 <= 6.2 and 4.4 <= d3 <= 6.0 and 5.2 <= d4 <= 7.2:
            for k in range(i, i + 4):
                states[k] = "H"
        # Beta-strand: extended chain (d(i,i+2) ~ 6.7A, d(i,i+3) ~ 9.9A, d(i,i+4) >= 11.5A)
        elif d2 >= 6.1 and d3 >= 8.8 and d4 >= 11.2:
            for k in range(i, i + 4):
                if states[k] == "C":
                    states[k] = "E"

    helices: list[SseSegment3D] = []
    strands: list[SseSegment3D] = []
    idx = 0
    while idx < n:
        st = states[idx]
        if st == "C":
            idx += 1
            continue
        end = idx
        while end + 1 < n and states[end + 1] == st:
            end += 1
        seg_len = end - idx + 1
        if st == "H" and seg_len >= 5:
            helices.append(SseSegment3D(kind="H", start_idx=idx, end_idx=end, coords=coords[idx : end + 1]))
        elif st == "E" and seg_len >= 3:
            strands.append(SseSegment3D(kind="E", start_idx=idx, end_idx=end, coords=coords[idx : end + 1]))
        idx = end + 1

    return helices, strands


def strip_t_prefixes(cath_code: str) -> str:
    """Convert a 'T'-prefixed CATH code (e.g. '3.T40.T50.T10') into a numeric CATH code ('3.40.50.10')."""
    if not cath_code or cath_code == "-":
        return "-"
    primary = cath_code.split(",")[0].strip()
    parts = [p.lstrip("Tt").strip() for p in primary.split(".") if p.strip()]
    return ".".join(p for p in parts if p)


def is_estimated_cath_code(cath_code: str) -> bool:
    """Return True if any level of the CATH code carries the 'T' estimation prefix."""
    if not cath_code or cath_code == "-":
        return False
    primary = cath_code.split(",")[0].strip()
    return any(p.strip().upper().startswith("T") for p in primary.split("."))


def format_t_prefixed_cath(
    full_numeric_parts: list[str],
    known_depth: int = 1,
) -> str:
    """Format a 4-part CATH code where levels up to `known_depth` are direct numbers and deeper levels have 'T'."""
    formatted: list[str] = []
    for idx, part in enumerate(full_numeric_parts):
        clean_num = part.lstrip("Tt").strip()
        if idx < known_depth:
            formatted.append(clean_num)
        else:
            formatted.append(f"T{clean_num}")
    return ".".join(formatted)


def classify_cath_with_rules(
    member: TedDomainMember,
    pdb_text: Optional[str] = None,
) -> tuple[list[str], int, str]:
    """Determine CATH [C, A, T, H], `known_depth` (1..3), and rule description from SSEs and 3D coordinates.

    Any level <= `known_depth` is determined by a direct structural rule (no 'T' prefix).
    Any level > `known_depth` is an estimate and receives the 'T' prefix.
    """
    h = max(0, member.num_helix)
    s = max(0, member.num_strand)
    t = max(0, member.num_turn)
    nres = max(1, member.nres_domain)
    total_ss = h + s
    norm_rg = member.norm_rg if member.norm_rg is not None else 0.30
    pack = member.packing_density if member.packing_density is not None else 10.0

    ca_coords = parse_ca_coordinates(pdb_text)
    helices_3d, strands_3d = extract_3d_sse_segments(ca_coords) if ca_coords else ([], [])

    # =========================================================================
    # LEVEL 1: CLASS (C) — Always determined directly by SSE composition rules
    # =========================================================================
    if total_ss <= 1 and (nres < 45 or (total_ss == 0 and t >= 2)):
        # Class 4: Few Secondary Structures -> 4.10 is the sole CATH architecture in Class 4 (direct depth=2)
        return ["4", "10", "10", "10"], 2, "rule_class4_irregular"

    helix_frac = h / max(1, total_ss)
    strand_frac = s / max(1, total_ss)

    # =========================================================================
    # CLASS 1: MAINLY ALPHA
    # =========================================================================
    if s == 0 or (h >= 2 and helix_frac >= 0.80):
        # Rule 1A: Single alpha-helix or extended coiled-coil -> direct 1.20.5 (depth=3)
        if h <= 1 or (h == 2 and norm_rg >= 0.42):
            return ["1", "20", "5", "10"], 3, "rule_alpha_coiled_coil_1.20.5"

        # Rule 1B: 2-helix antiparallel hairpin -> direct 1.10.287 (depth=3)
        if h == 2 and norm_rg < 0.40:
            return ["1", "10", "287", "10"], 3, "rule_alpha_hairpin_1.10.287"

        # Rule 1C: 3D inter-helical axis crossing angle analysis when PDB trace is available
        if len(helices_3d) >= 2:
            angles_deg: list[float] = []
            for i in range(len(helices_3d)):
                for j in range(i + 1, len(helices_3d)):
                    # Check if helices i and j are in spatial contact (centroid or closest CA < 13.5 A)
                    if helices_3d[i].centroid.dist(helices_3d[j].centroid) <= 16.0:
                        cos_val = min(1.0, abs(helices_3d[i].axis.dot(helices_3d[j].axis)))
                        angles_deg.append(math.degrees(math.acos(cos_val)))
            if angles_deg:
                mean_angle = sum(angles_deg) / len(angles_deg)
                if h >= 6 and norm_rg >= 0.33:
                    # Repetitive alpha-solenoid / horseshoe (TPR/Ankyrin) -> direct 1.25.40 (depth=3)
                    return ["1", "25", "40", "10"], 3, "rule_3d_alpha_solenoid_1.25.40"
                if mean_angle < 34.0:
                    # Parallel/antiparallel bundle -> direct Architecture 1.20 (depth=2; depth=3 if 4-helix bundle)
                    depth = 3 if 3 <= h <= 5 else 2
                    return ["1", "20", "120", "10"], depth, "rule_3d_updown_bundle_1.20"
                # Oblique/orthogonal crossing helices -> direct Architecture 1.10 (depth=2)
                return ["1", "10", "10", "10"], 2, "rule_3d_orthogonal_bundle_1.10"

        # Fallback SSE rules for Class 1 when 3D trace is not present
        if h >= 6 and norm_rg >= 0.33:
            return ["1", "25", "40", "10"], 2, "rule_sse_alpha_horseshoe_1.25"
        if 3 <= h <= 5 and pack >= 10.6 and norm_rg < 0.30:
            return ["1", "10", "10", "10"], 2, "rule_sse_orthogonal_bundle_1.10"
        return ["1", "20", "120", "10"], 2, "rule_sse_updown_bundle_1.20"

    # =========================================================================
    # CLASS 2: MAINLY BETA
    # =========================================================================
    if h == 0 or (s >= 3 and strand_frac >= 0.78):
        # Rule 2A: Radial multi-bladed beta-propeller (>= 14 strands, compact Rg) -> direct 2.130.10 (depth=3)
        if s >= 14 and norm_rg < 0.29:
            return ["2", "130", "10", "10"], 3, "rule_beta_propeller_2.130.10"

        # Rule 2B: 3D beta-strand contact graph analysis (closed barrel vs 2-layer sandwich)
        if len(strands_3d) >= 4:
            adj_deg = [0] * len(strands_3d)
            for i in range(len(strands_3d)):
                for j in range(i + 1, len(strands_3d)):
                    # Minimum CA-CA distance between strands i and j
                    min_d = min(
                        c1.dist(c2) for c1 in strands_3d[i].coords for c2 in strands_3d[j].coords
                    )
                    if min_d <= 5.8:
                        adj_deg[i] += 1
                        adj_deg[j] += 1
            closed_fraction = sum(1 for d in adj_deg if d >= 2) / len(strands_3d)
            if closed_fraction >= 0.75:
                # Closed beta-barrel or roll -> direct Architecture 2.30 or 2.40
                if s <= 5 and nres <= 95:
                    return ["2", "30", "30", "10"], 3, "rule_3d_sh3_beta_roll_2.30.30"
                return ["2", "40", "50", "10"], 2, "rule_3d_beta_barrel_2.40"
            if 6 <= s <= 10:
                # 2-sheet beta-sandwich -> direct Architecture 2.60 (depth=3 if canonical Ig-like)
                if 7 <= s <= 9 and 80 <= nres <= 155 and pack >= 10.8:
                    return ["2", "60", "40", "10"], 3, "rule_3d_ig_sandwich_2.60.40"
                return ["2", "60", "40", "10"], 2, "rule_3d_beta_sandwich_2.60"

        # Fallback SSE rules for Class 2
        if s <= 5:
            if nres <= 90:
                return ["2", "30", "30", "10"], 2, "rule_sse_beta_roll_2.30"
            return ["2", "40", "50", "10"], 2, "rule_sse_beta_barrel_2.40"
        if 6 <= s <= 10:
            if pack >= 10.8:
                return ["2", "60", "40", "10"], 2, "rule_sse_beta_sandwich_2.60"
            return ["2", "60", "120", "10"], 2, "rule_sse_beta_sandwich_2.60"
        return ["2", "40", "160", "10"], 2, "rule_sse_beta_barrel_2.40"

    # =========================================================================
    # CLASS 3: ALPHA BETA
    # =========================================================================
    # Rule 3A: (beta/alpha)8 TIM Barrel -> direct 3.20.20 (depth=3)
    if h >= 6 and s >= 6 and 0.70 <= (h / max(1, s)) <= 1.40 and nres >= 175:
        if helices_3d and strands_3d:
            # Verify inner beta-barrel surrounded by outer alpha-helices (R_strand < R_helix)
            all_pts = [seg.centroid for seg in helices_3d + strands_3d]
            center = Vec3(
                sum(p.x for p in all_pts) / len(all_pts),
                sum(p.y for p in all_pts) / len(all_pts),
                sum(p.z for p in all_pts) / len(all_pts),
            )
            r_strand = sum(seg.centroid.dist(center) for seg in strands_3d) / len(strands_3d)
            r_helix = sum(seg.centroid.dist(center) for seg in helices_3d) / len(helices_3d)
            if r_strand < r_helix:
                return ["3", "20", "20", "10"], 3, "rule_3d_tim_barrel_3.20.20"
        else:
            return ["3", "20", "20", "10"], 3, "rule_sse_tim_barrel_3.20.20"

    # Rule 3B: Compact Ubiquitin-like beta-grasp roll -> direct 3.10.20 (depth=3)
    if 1 <= h <= 2 and 3 <= s <= 5 and nres <= 110:
        return ["3", "10", "20", "10"], 3, "rule_ub_roll_3.10.20"

    # Rule 3C: 3D layer analysis for 3-Layer (aba) Sandwich (3.40) vs 2-Layer (ab) Sandwich (3.30)
    if len(helices_3d) >= 2 and len(strands_3d) >= 3:
        s_pts = [seg.centroid for seg in strands_3d]
        s_center = Vec3(
            sum(p.x for p in s_pts) / len(s_pts),
            sum(p.y for p in s_pts) / len(s_pts),
            sum(p.z for p in s_pts) / len(s_pts),
        )
        # Approximate sheet normal from strand axis and sheet spread vector
        v1 = strands_3d[0].axis
        v2 = (strands_3d[-1].centroid - strands_3d[0].centroid).normalized()
        normal = v1.cross(v2).normalized()
        if normal.dot(normal) > 0.1:
            signed_dists = [(h_seg.centroid - s_center).dot(normal) for h_seg in helices_3d]
            pos_side = sum(1 for d in signed_dists if d > 3.0)
            neg_side = sum(1 for d in signed_dists if d < -3.0)
            if pos_side >= 1 and neg_side >= 1:
                # Helices on BOTH faces of the beta-sheet -> direct 3-Layer (aba) Rossmann sandwich (3.40.50, depth=3)
                return ["3", "40", "50", "10"], 3, "rule_3d_3layer_aba_rossmann_3.40.50"
            # Helices on ONE face of the beta-sheet -> direct 2-Layer Sandwich (3.30.70, depth=3)
            return ["3", "30", "70", "10"], 3, "rule_3d_2layer_ab_plait_3.30.70"

    # Fallback SSE rules for Class 3
    if h <= 3 and 3 <= s <= 6 and nres <= 165:
        return ["3", "30", "70", "10"], 2, "rule_sse_2layer_sandwich_3.30"
    if h >= 3 and s >= 4:
        return ["3", "40", "50", "10"], 2, "rule_sse_3layer_rossmann_3.40"

    return ["3", "30", "70", "10"], 2, "rule_sse_2layer_sandwich_3.30"


def classify_cath_from_sse(member: TedDomainMember) -> tuple[str, str, str, str]:
    """Backwards-compatible helper returning the 4 numeric CATH parts (C, A, T, H) from SSE rules."""
    parts, _depth, _rule = classify_cath_with_rules(member, pdb_text=None)
    return (parts[0], parts[1], parts[2], parts[3])


def infer_cath_classification(
    member: TedDomainMember,
    overlapping_annotations: Iterable[OverlappingAnnotation] = (),
    pdb_text: Optional[str] = None,
) -> tuple[str, str, str]:
    """Determine the full 4-level CATH classification for a TED domain.

    Levels determined by official TED assignment, high-IoU G3DSA match, or direct structural/3D rules
    are written as plain numbers; remaining estimated levels are prefixed with 'T'.
    """
    raw_ted = (member.cath_label or "-").split(",")[0].strip()
    ted_parts = [p.strip() for p in raw_ted.split(".") if p.strip() and p.strip() != "-"]

    # Case 1: TED already provides a full 4-level H-code
    if len(ted_parts) >= 4:
        numeric = ".".join(ted_parts[:4])
        return numeric, numeric, "ted_official"

    # Check if there is an overlapping CATH-Gene3D (G3DSA:c.a.t.h) annotation from InterPro
    best_g3dsa_parts: list[str] | None = None
    best_g3dsa_iou = 0.0
    best_g3dsa_ted_cov = 0.0
    best_g3dsa_score = -1.0
    for annot in overlapping_annotations:
        if annot.source_db == "cathgene3d" or annot.accession.upper().startswith("G3DSA:"):
            if annot.iou >= 0.25 or annot.ted_coverage >= 0.45:
                m = _G3DSA_RE.match(annot.accession.strip())
                if m:
                    score = annot.iou + 0.5 * annot.ted_coverage
                    if score > best_g3dsa_score:
                        best_g3dsa_score = score
                        best_g3dsa_iou = annot.iou
                        best_g3dsa_ted_cov = annot.ted_coverage
                        best_g3dsa_parts = [m.group(1), m.group(2), m.group(3), m.group(4)]

    # Case 2: TED provides a 3-level T-code (e.g. '2.60.40' or '1.20.5')
    if len(ted_parts) == 3:
        if best_g3dsa_parts and best_g3dsa_parts[:3] == ted_parts[:3]:
            full_parts = ted_parts[:3] + [best_g3dsa_parts[3]]
            # If G3DSA matches with high boundary IoU (>=0.70), H-level is direct; else T-prefixed
            k_depth = 4 if (best_g3dsa_iou >= 0.70 and best_g3dsa_ted_cov >= 0.80) else 3
            return (
                format_t_prefixed_cath(full_parts, known_depth=k_depth),
                ".".join(full_parts),
                "ted_t_plus_g3dsa",
            )
        full_parts = ted_parts[:3] + ["10"]
        return (
            format_t_prefixed_cath(full_parts, known_depth=3),
            ".".join(full_parts),
            "ted_t_plus_sse",
        )

    # Case 3: TED cath_label is '-' (unassigned) and G3DSA overlaps
    if best_g3dsa_parts:
        # High-IoU G3DSA match directly establishes C.A.T (depth=3), or C.A (depth=2) for partial match
        k_depth = 3 if (best_g3dsa_iou >= 0.65 and best_g3dsa_ted_cov >= 0.75) else 1
        return (
            format_t_prefixed_cath(best_g3dsa_parts, known_depth=k_depth),
            ".".join(best_g3dsa_parts),
            "g3dsa_overlap",
        )

    # Case 4: Direct structural & 3D geometric rules (Class C, Architecture A, Topology T)
    rule_parts, known_depth, rule_name = classify_cath_with_rules(member, pdb_text=pdb_text)
    return (
        format_t_prefixed_cath(rule_parts, known_depth=known_depth),
        ".".join(rule_parts),
        rule_name,
    )
