//! Custom vector-drawn ChiralCraft icons:
//! - 18x2 `Tools` palette icons (`tools-20261118.png`) with blue-gradient active highlight
//! - Unified General + Style top toolbar icons (`p08-opened.png`) with zero emoji dependencies
//! - Chemical formula subscript/superscript formatting helpers

use egui::{Color32, FontId, Painter, Pos2, Rect, Sense, Shape, Stroke, Ui, Vec2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    Lasso,
    Marquee,
    Perspective,
    Fragmentation,
    SolidBond,
    Eraser,
    MultipleBond,
    DoubleBond,
    TripleBond,
    DativeBond,
    QuadrupleBond,
    TextLabel,
    DashedBond,
    PenCurve,
    HashedBond,
    Arrow,
    EquilibriumArrow,
    NoGoArrow,
    RetrosyntheticArrow,
    CurvedArrow,
    HashWedgeBond,
    Orbital,
    OrbitalS,
    OrbitalD,
    BoldBond,
    ShapeRect,
    ShapeOval,
    ShapeLine,
    WedgeBond,
    Bracket,
    RoundBracket,
    CurlyBracket,
    HollowWedgeBond,
    ChemSymbol,
    SymbolMinus,
    SymbolRadical,
    SymbolLonePair,
    WavyBond,
    ReactionMap,
    Table,
    TlcPlate,
    AcyclicChain,
    Templates,
    Cyclopropane,
    Cyclobutane,
    Cyclopentane,
    Cyclohexane,
    Cycloheptane,
    Cyclooctane,
    Chair1,
    Chair2,
    Cyclopentadiene,
    Benzene,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToolbarIcon {
    NewDoc,
    Open,
    Save,
    Print,
    Help,
    Undo,
    Redo,
    Cut,
    Copy,
    Paste,
    ZoomIn,
    ZoomOut,
    CleanUp,
    AlignLeft,
    AlignCenter,
    AlignRight,
    AlignJustify,
}

/// Converts plain ASCII chemical formulas like `C9H8O4`, `CH3`, `NH4+`, `SO4-`
/// into crisp Unicode subscript/superscript display strings (`C₉H₈O₄`, `CH₃`, `NH₄⁺`, `SO₄⁻`).
pub fn format_chemical_formula(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    let chars: Vec<char> = s.chars().collect();
    for (i, &ch) in chars.iter().enumerate() {
        let prev_is_chem = i > 0 && (chars[i - 1].is_ascii_alphabetic() || chars[i - 1] == ')' || chars[i - 1] == ']' || chars[i - 1].is_ascii_digit());
        if ch.is_ascii_digit() && prev_is_chem {
            let sub = match ch {
                '0' => '₀',
                '1' => '₁',
                '2' => '₂',
                '3' => '₃',
                '4' => '₄',
                '5' => '₅',
                '6' => '₆',
                '7' => '₇',
                '8' => '₈',
                '9' => '₉',
                _ => ch,
            };
            out.push(sub);
        } else if ch == '+' && i > 0 && i == chars.len() - 1 {
            out.push('⁺');
        } else if ch == '-' && i > 0 && i == chars.len() - 1 {
            out.push('⁻');
        } else {
            out.push(ch);
        }
    }
    out
}

/// Draws a classic Windows/ChiralCraft steel-blue vertical gradient strip inside `rect`.
pub fn draw_steel_blue_gradient(p: &Painter, rect: Rect) {
    let bands = 10;
    let h = rect.height() / (bands as f32);
    for i in 0..bands {
        let t = (i as f32) / ((bands - 1) as f32);
        let r = (158.0 + t * 44.0) as u8;
        let g = (190.0 + t * 34.0) as u8;
        let b = (224.0 + t * 20.0) as u8;
        let y0 = rect.min.y + (i as f32) * h;
        let band = Rect::from_min_max(Pos2::new(rect.min.x, y0), Pos2::new(rect.max.x, (y0 + h + 0.5).min(rect.max.y)));
        p.rect_filled(band, 0.0, Color32::from_rgb(r, g, b));
    }
}

/// Renders a crisp `22x22` ChiralCraft top toolbar icon button matching `p08-opened.png`.
pub fn toolbar_icon_button(ui: &mut Ui, icon: ToolbarIcon, active: bool, tooltip: &str) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(Vec2::splat(22.0), Sense::click());
    let p = ui.painter();
    if active {
        p.rect_filled(rect, 1.0, Color32::from_rgb(194, 224, 250));
        p.rect_stroke(rect, 1.0, Stroke::new(1.0_f32, Color32::from_rgb(43, 87, 154)), egui::StrokeKind::Inside);
    } else if resp.hovered() {
        p.rect_filled(rect, 1.0, Color32::from_rgb(224, 238, 252));
        p.rect_stroke(rect, 1.0, Stroke::new(1.0_f32, Color32::from_rgb(122, 174, 224)), egui::StrokeKind::Inside);
    }
    let r = rect.shrink(3.5);
    let c = r.center();
    let dark = Color32::from_gray(72);
    let stroke = Stroke::new(1.2_f32, dark);
    let thin = Stroke::new(1.0_f32, dark);

    match icon {
        ToolbarIcon::NewDoc => {
            let doc = Rect::from_min_max(Pos2::new(r.min.x + 2.0, r.min.y + 1.0), Pos2::new(r.max.x - 2.0, r.max.y - 1.0));
            p.rect_filled(doc, 0.0, Color32::from_rgb(106, 92, 150));
            p.rect_stroke(doc, 0.0, Stroke::new(1.0_f32, Color32::from_rgb(68, 58, 102)), egui::StrokeKind::Inside);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(doc.max.x - 4.0, doc.min.y), Pos2::new(doc.max.x, doc.min.y + 4.0), Pos2::new(doc.max.x - 4.0, doc.min.y + 4.0)],
                Color32::from_rgb(220, 215, 240),
                Stroke::NONE,
            ));
        }
        ToolbarIcon::Open => {
            let tab = Rect::from_min_max(Pos2::new(r.min.x + 1.0, r.min.y + 3.0), Pos2::new(r.min.x + 7.0, r.min.y + 6.0));
            let body = Rect::from_min_max(Pos2::new(r.min.x + 1.0, r.min.y + 5.0), Pos2::new(r.max.x - 1.0, r.max.y - 1.0));
            p.rect_filled(tab, 0.0, Color32::from_rgb(206, 178, 118));
            p.rect_filled(body, 0.0, Color32::from_rgb(225, 198, 134));
            p.rect_stroke(body, 0.0, Stroke::new(1.0_f32, Color32::from_rgb(138, 114, 70)), egui::StrokeKind::Inside);
            // Orange arrow
            p.line_segment([Pos2::new(c.x - 3.0, r.min.y + 3.0), Pos2::new(c.x + 3.0, r.min.y + 3.0)], Stroke::new(1.6_f32, Color32::from_rgb(220, 108, 30)));
            p.add(Shape::convex_polygon(
                vec![Pos2::new(c.x + 4.5, r.min.y + 3.0), Pos2::new(c.x + 1.5, r.min.y + 0.5), Pos2::new(c.x + 1.5, r.min.y + 5.5)],
                Color32::from_rgb(220, 108, 30),
                Stroke::NONE,
            ));
        }
        ToolbarIcon::Save => {
            let fr = r.shrink(1.0);
            p.rect_filled(fr, 0.0, Color32::from_rgb(136, 146, 158));
            p.rect_stroke(fr, 0.0, Stroke::new(1.0_f32, Color32::from_gray(80)), egui::StrokeKind::Inside);
            p.rect_filled(Rect::from_min_max(Pos2::new(fr.min.x + 3.0, fr.min.y + 1.0), Pos2::new(fr.max.x - 3.0, fr.min.y + 6.0)), 0.0, Color32::from_gray(235));
            p.rect_filled(Rect::from_min_max(Pos2::new(fr.min.x + 3.0, fr.max.y - 5.0), Pos2::new(fr.max.x - 3.0, fr.max.y - 1.0)), 0.0, Color32::from_gray(220));
        }
        ToolbarIcon::Print => {
            let top = Rect::from_min_max(Pos2::new(r.min.x + 3.0, r.min.y + 1.0), Pos2::new(r.max.x - 3.0, r.min.y + 6.0));
            let mid = Rect::from_min_max(Pos2::new(r.min.x + 1.0, r.min.y + 5.0), Pos2::new(r.max.x - 1.0, r.max.y - 3.0));
            p.rect_filled(top, 0.0, Color32::WHITE);
            p.rect_stroke(top, 0.0, thin, egui::StrokeKind::Inside);
            p.rect_filled(mid, 1.0, Color32::from_gray(165));
            p.rect_stroke(mid, 1.0, thin, egui::StrokeKind::Inside);
            p.rect_filled(Rect::from_min_max(Pos2::new(r.min.x + 3.0, r.max.y - 4.0), Pos2::new(r.max.x - 3.0, r.max.y - 1.0)), 0.0, Color32::WHITE);
        }
        ToolbarIcon::Help => {
            p.circle_filled(c, 6.8, Color32::from_rgb(18, 90, 170));
            p.text(c, egui::Align2::CENTER_CENTER, "?", FontId::proportional(11.0), Color32::WHITE);
        }
        ToolbarIcon::Undo => {
            let blue = Color32::from_rgb(76, 122, 178);
            let bs = Stroke::new(2.0_f32, blue);
            p.line_segment([Pos2::new(c.x - 4.0, c.y - 1.0), Pos2::new(c.x + 3.0, c.y - 1.0)], bs);
            p.line_segment([Pos2::new(c.x + 3.0, c.y - 1.0), Pos2::new(c.x + 4.0, c.y + 4.0)], bs);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(r.min.x, c.y - 1.0), Pos2::new(c.x - 1.5, c.y - 4.8), Pos2::new(c.x - 1.5, c.y + 2.8)],
                blue,
                Stroke::NONE,
            ));
        }
        ToolbarIcon::Redo => {
            let blue = Color32::from_rgb(122, 146, 174);
            let bs = Stroke::new(2.0_f32, blue);
            p.line_segment([Pos2::new(c.x - 3.0, c.y - 1.0), Pos2::new(c.x + 4.0, c.y - 1.0)], bs);
            p.line_segment([Pos2::new(c.x - 3.0, c.y - 1.0), Pos2::new(c.x - 4.0, c.y + 4.0)], bs);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(r.max.x, c.y - 1.0), Pos2::new(c.x + 1.5, c.y - 4.8), Pos2::new(c.x + 1.5, c.y + 2.8)],
                blue,
                Stroke::NONE,
            ));
        }
        ToolbarIcon::Cut => {
            p.line_segment([Pos2::new(r.min.x + 3.0, r.min.y + 1.0), Pos2::new(r.max.x - 3.0, r.max.y - 4.0)], stroke);
            p.line_segment([Pos2::new(r.max.x - 3.0, r.min.y + 1.0), Pos2::new(r.min.x + 3.0, r.max.y - 4.0)], stroke);
            p.circle_stroke(Pos2::new(r.min.x + 3.0, r.max.y - 2.5), 2.3, stroke);
            p.circle_stroke(Pos2::new(r.max.x - 3.0, r.max.y - 2.5), 2.3, stroke);
        }
        ToolbarIcon::Copy => {
            let r1 = Rect::from_min_max(Pos2::new(r.min.x + 1.0, r.min.y + 1.0), Pos2::new(r.max.x - 4.0, r.max.y - 4.0));
            let r2 = Rect::from_min_max(Pos2::new(r.min.x + 5.0, r.min.y + 4.0), Pos2::new(r.max.x - 1.0, r.max.y - 1.0));
            p.rect_filled(r1, 0.0, Color32::WHITE);
            p.rect_stroke(r1, 0.0, thin, egui::StrokeKind::Inside);
            p.rect_filled(r2, 0.0, Color32::from_rgb(230, 240, 250));
            p.rect_stroke(r2, 0.0, Stroke::new(1.0_f32, Color32::from_rgb(65, 110, 165)), egui::StrokeKind::Inside);
        }
        ToolbarIcon::Paste => {
            let brd = Rect::from_min_max(Pos2::new(r.min.x + 1.0, r.min.y + 2.0), Pos2::new(r.max.x - 4.0, r.max.y - 2.0));
            let doc = Rect::from_min_max(Pos2::new(r.min.x + 5.0, r.min.y + 5.0), Pos2::new(r.max.x - 1.0, r.max.y - 1.0));
            p.rect_filled(brd, 0.0, Color32::from_rgb(176, 152, 118));
            p.rect_stroke(brd, 0.0, thin, egui::StrokeKind::Inside);
            p.rect_filled(doc, 0.0, Color32::from_rgb(80, 134, 194));
            p.rect_stroke(doc, 0.0, Stroke::new(1.0_f32, Color32::from_rgb(45, 85, 135)), egui::StrokeKind::Inside);
        }
        ToolbarIcon::ZoomIn | ToolbarIcon::ZoomOut => {
            let gc = Pos2::new(c.x - 1.5, c.y - 1.5);
            p.circle_filled(gc, 4.5, Color32::WHITE);
            p.circle_stroke(gc, 4.5, stroke);
            p.line_segment([Pos2::new(gc.x + 3.2, gc.y + 3.2), Pos2::new(r.max.x - 1.0, r.max.y - 1.0)], Stroke::new(2.0_f32, dark));
            let orange = Stroke::new(1.3_f32, Color32::from_rgb(215, 100, 20));
            p.line_segment([Pos2::new(gc.x - 2.5, gc.y), Pos2::new(gc.x + 2.5, gc.y)], orange);
            if icon == ToolbarIcon::ZoomIn {
                p.line_segment([Pos2::new(gc.x, gc.y - 2.5), Pos2::new(gc.x, gc.y + 2.5)], orange);
            }
        }
        ToolbarIcon::CleanUp => {
            let pts = poly_pts(c, 5.5, 6, -90.0);
            draw_polygon(p, &pts, Stroke::new(1.1_f32, Color32::from_rgb(45, 95, 165)));
            p.circle_filled(pts[0], 1.6, Color32::from_rgb(200, 50, 50));
            p.circle_filled(pts[3], 1.6, Color32::from_rgb(45, 95, 165));
        }
        ToolbarIcon::AlignLeft | ToolbarIcon::AlignCenter | ToolbarIcon::AlignRight | ToolbarIcon::AlignJustify => {
            let widths = [1.0_f32, 0.65, 0.9, 0.7];
            for (i, &w) in widths.iter().enumerate() {
                let y = r.min.y + 2.0 + (i as f32) * 3.5;
                let lw = if icon == ToolbarIcon::AlignJustify { r.width() } else { r.width() * w };
                let x0 = match icon {
                    ToolbarIcon::AlignLeft | ToolbarIcon::AlignJustify => r.min.x,
                    ToolbarIcon::AlignCenter => c.x - lw * 0.5,
                    ToolbarIcon::AlignRight => r.max.x - lw,
                    _ => r.min.x,
                };
                p.line_segment([Pos2::new(x0, y), Pos2::new(x0 + lw, y)], Stroke::new(1.1_f32, Color32::from_gray(90)));
            }
        }
    }
    resp.on_hover_text(tooltip)
}

/// Renders a single `26x26` ChiralCraft iconic tool button with exact blue-gradient
/// active selection highlight (`tools-20261118.png`) and optional flyout indicator.
pub fn tool_icon_button(
    ui: &mut Ui,
    tool: Tool,
    active: bool,
    has_flyout: bool,
    tooltip: &str,
) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(Vec2::splat(26.0), Sense::click());
    let p = ui.painter();

    if active {
        // Exact ChiralCraft blue gradient active selection (`tools-20261118.png`)
        let bands = 8;
        let h = rect.height() / (bands as f32);
        for i in 0..bands {
            let t = (i as f32) / ((bands - 1) as f32);
            let cr = (130.0 + t * 88.0) as u8;
            let cg = (178.0 + t * 58.0) as u8;
            let cb = (232.0 + t * 20.0) as u8;
            let y0 = rect.min.y + (i as f32) * h;
            p.rect_filled(
                Rect::from_min_max(Pos2::new(rect.min.x, y0), Pos2::new(rect.max.x, (y0 + h + 0.5).min(rect.max.y))),
                0.0,
                Color32::from_rgb(cr, cg, cb),
            );
        }
        p.rect_stroke(
            rect,
            0.0,
            Stroke::new(1.2_f32, Color32::from_rgb(20, 58, 138)),
            egui::StrokeKind::Inside,
        );
    } else if resp.hovered() {
        p.rect_filled(rect, 0.0, Color32::from_rgb(222, 238, 252));
        p.rect_stroke(
            rect,
            0.0,
            Stroke::new(1.0_f32, Color32::from_rgb(122, 174, 224)),
            egui::StrokeKind::Inside,
        );
    }

    draw_tool_glyph(p, rect.shrink(4.0), tool);

    if has_flyout {
        let br = rect.max - Vec2::new(2.0, 2.0);
        p.add(Shape::convex_polygon(
            vec![br, br - Vec2::new(4.2, 0.0), br - Vec2::new(0.0, 4.2)],
            Color32::BLACK,
            Stroke::NONE,
        ));
    }

    resp.on_hover_text(tooltip)
}

fn poly_pts(c: Pos2, r: f32, n: usize, rot_deg: f32) -> Vec<Pos2> {
    let rot = rot_deg.to_radians();
    (0..n)
        .map(|i| {
            let a = rot + (i as f32) * std::f32::consts::TAU / (n as f32);
            Pos2::new(c.x + r * a.cos(), c.y + r * a.sin())
        })
        .collect()
}

fn draw_polygon(p: &Painter, pts: &[Pos2], stroke: Stroke) {
    for i in 0..pts.len() {
        p.line_segment([pts[i], pts[(i + 1) % pts.len()]], stroke);
    }
}

fn draw_tool_glyph(p: &Painter, r: Rect, tool: Tool) {
    let c = r.center();
    let s = Stroke::new(1.3_f32, Color32::BLACK);
    let thin = Stroke::new(1.0_f32, Color32::BLACK);
    let bold = Stroke::new(2.6_f32, Color32::BLACK);
    let x0 = r.min.x;
    let y0 = r.min.y;
    let x1 = r.max.x;
    let y1 = r.max.y;

    match tool {
        Tool::Lasso => {
            p.circle_stroke(Pos2::new(c.x + 1.0, c.y - 1.5), 5.8, thin);
            p.line_segment([Pos2::new(x0 + 2.0, y1 - 1.5), Pos2::new(c.x - 3.0, c.y + 2.5)], thin);
            p.circle_filled(Pos2::new(c.x - 3.2, c.y + 2.4), 1.4, Color32::BLACK);
        }
        Tool::Marquee => {
            let mr = Rect::from_min_max(Pos2::new(x0 + 1.0, y0 + 1.0), Pos2::new(x1 - 4.0, y1 - 4.0));
            p.rect_stroke(mr, 0.0, thin, egui::StrokeKind::Inside);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(c.x + 1.0, c.y + 1.0), Pos2::new(x1, c.y + 4.5), Pos2::new(c.x + 4.5, y1)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::Perspective => {
            p.circle_stroke(c, 6.5, thin);
            p.line_segment([Pos2::new(c.x - 5.0, c.y), Pos2::new(c.x + 4.0, c.y)], s);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(c.x + 4.5, c.y), Pos2::new(c.x + 0.5, c.y - 3.2), Pos2::new(c.x + 0.5, c.y + 3.2)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::Fragmentation => {
            p.line_segment([Pos2::new(x0 + 3.0, y0 + 1.5), Pos2::new(c.x, y0 + 4.5)], s);
            p.line_segment([Pos2::new(c.x, y0 + 4.5), Pos2::new(c.x, y1 - 4.5)], s);
            p.line_segment([Pos2::new(c.x, y1 - 4.5), Pos2::new(x1 - 3.0, y1 - 1.5)], s);
            for k in 0..4 {
                let t0 = k as f32 / 4.0;
                let t1 = (k as f32 + 0.55) / 4.0;
                p.line_segment(
                    [Pos2::new(x0 + t0 * r.width(), c.y - 2.0 + t0 * 4.0), Pos2::new(x0 + t1 * r.width(), c.y - 2.0 + t1 * 4.0)],
                    Stroke::new(1.2_f32, Color32::RED),
                );
            }
        }
        Tool::SolidBond => {
            p.line_segment([Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x1 - 2.0, y1 - 2.0)], s);
        }
        Tool::Eraser => {
            let er = Rect::from_center_size(c, Vec2::new(13.0, 8.0));
            p.rect_filled(er, 0.0, Color32::WHITE);
            p.rect_stroke(er, 0.0, thin, egui::StrokeKind::Inside);
            p.line_segment([Pos2::new(c.x, er.min.y), Pos2::new(c.x - 2.0, er.max.y)], thin);
        }
        Tool::MultipleBond | Tool::DoubleBond => {
            p.line_segment([Pos2::new(x0 + 3.2, y0 + 1.2), Pos2::new(x1 - 1.2, y1 - 3.2)], thin);
            p.line_segment([Pos2::new(x0 + 1.2, y0 + 3.2), Pos2::new(x1 - 3.2, y1 - 1.2)], thin);
        }
        Tool::TripleBond => {
            for d in [-2.6_f32, 0.0, 2.6] {
                p.line_segment([Pos2::new(x0 + 2.0 + d, y0 + 2.0 - d), Pos2::new(x1 - 2.0 + d, y1 - 2.0 - d)], thin);
            }
        }
        Tool::QuadrupleBond => {
            for d in [-3.3_f32, -1.1, 1.1, 3.3] {
                p.line_segment([Pos2::new(x0 + 2.5 + d, y0 + 2.5 - d), Pos2::new(x1 - 2.5 + d, y1 - 2.5 - d)], thin);
            }
        }
        Tool::DativeBond => {
            p.line_segment([Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x1 - 2.0, y1 - 2.0)], s);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(x1 - 1.0, y1 - 1.0), Pos2::new(x1 - 6.0, y1 - 2.0), Pos2::new(x1 - 2.0, y1 - 6.0)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::TextLabel => {
            p.text(c, egui::Align2::CENTER_CENTER, "A", FontId::proportional(15.5), Color32::BLACK);
        }
        Tool::DashedBond => {
            for k in 0..5 {
                let t = k as f32 / 4.0;
                p.circle_filled(Pos2::new(x0 + 2.0 + t * (r.width() - 4.0), y0 + 2.0 + t * (r.height() - 4.0)), 1.05, Color32::BLACK);
            }
        }
        Tool::PenCurve => {
            p.line_segment([Pos2::new(x0 + 2.0, y1 - 2.0), Pos2::new(c.x - 2.0, y0 + 4.0)], thin);
            p.line_segment([Pos2::new(c.x - 2.0, y0 + 4.0), Pos2::new(x1 - 2.0, y0 + 2.0)], thin);
            p.circle_filled(Pos2::new(x1 - 3.5, y1 - 3.5), 2.6, Color32::BLACK);
        }
        Tool::HashedBond => {
            for k in 0..5 {
                let t = k as f32 / 4.0;
                let mx = x0 + 3.0 + t * (r.width() - 6.0);
                let my = y0 + 3.0 + t * (r.height() - 6.0);
                p.line_segment([Pos2::new(mx - 2.5, my + 2.5), Pos2::new(mx + 2.5, my - 2.5)], thin);
            }
        }
        Tool::Arrow => {
            p.line_segment([Pos2::new(x0 + 1.0, c.y), Pos2::new(x1 - 2.0, c.y)], s);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(x1, c.y), Pos2::new(x1 - 5.5, c.y - 3.0), Pos2::new(x1 - 5.5, c.y + 3.0)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::EquilibriumArrow => {
            p.line_segment([Pos2::new(x0 + 1.0, c.y - 2.0), Pos2::new(x1 - 1.0, c.y - 2.0)], thin);
            p.line_segment([Pos2::new(x1 - 4.5, c.y - 4.5), Pos2::new(x1 - 1.0, c.y - 2.0)], thin);
            p.line_segment([Pos2::new(x0 + 1.0, c.y + 2.0), Pos2::new(x1 - 1.0, c.y + 2.0)], thin);
            p.line_segment([Pos2::new(x0 + 4.5, c.y + 4.5), Pos2::new(x0 + 1.0, c.y + 2.0)], thin);
        }
        Tool::NoGoArrow => {
            p.line_segment([Pos2::new(x0 + 1.0, c.y), Pos2::new(x1 - 2.0, c.y)], s);
            p.add(Shape::convex_polygon(
                vec![Pos2::new(x1, c.y), Pos2::new(x1 - 5.0, c.y - 2.8), Pos2::new(x1 - 5.0, c.y + 2.8)],
                Color32::BLACK,
                Stroke::NONE,
            ));
            p.line_segment([Pos2::new(c.x - 3.0, c.y - 3.0), Pos2::new(c.x + 2.0, c.y + 3.0)], Stroke::new(1.2_f32, Color32::RED));
            p.line_segment([Pos2::new(c.x - 3.0, c.y + 3.0), Pos2::new(c.x + 2.0, c.y - 3.0)], Stroke::new(1.2_f32, Color32::RED));
        }
        Tool::RetrosyntheticArrow => {
            p.line_segment([Pos2::new(x0 + 1.0, c.y - 2.0), Pos2::new(x1 - 3.0, c.y - 2.0)], thin);
            p.line_segment([Pos2::new(x0 + 1.0, c.y + 2.0), Pos2::new(x1 - 3.0, c.y + 2.0)], thin);
            p.line_segment([Pos2::new(x1 - 4.5, c.y - 4.5), Pos2::new(x1 - 0.5, c.y)], thin);
            p.line_segment([Pos2::new(x1 - 4.5, c.y + 4.5), Pos2::new(x1 - 0.5, c.y)], thin);
        }
        Tool::CurvedArrow => {
            let pts = [Pos2::new(x0 + 2.0, y1 - 3.0), Pos2::new(x0 + 4.0, y0 + 3.0), Pos2::new(x1 - 2.0, c.y)];
            for w in pts.windows(2) { p.line_segment([w[0], w[1]], thin); }
            p.add(Shape::convex_polygon(
                vec![Pos2::new(x1 - 1.0, c.y + 1.0), Pos2::new(x1 - 5.0, c.y - 2.5), Pos2::new(x1 - 2.5, c.y - 4.0)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::HashWedgeBond => {
            for k in 1..=5 {
                let t = k as f32 / 5.0;
                let mx = x0 + 2.0 + t * (r.width() - 4.0);
                let my = y0 + 2.0 + t * (r.height() - 4.0);
                let hw = 3.2 * t;
                p.line_segment([Pos2::new(mx - hw, my + hw), Pos2::new(mx + hw, my - hw)], thin);
            }
        }
        Tool::Orbital => {
            p.circle_filled(Pos2::new(c.x, c.y - 3.5), 3.1, Color32::BLACK);
            p.circle_stroke(Pos2::new(c.x, c.y + 3.5), 3.1, thin);
        }
        Tool::OrbitalS => {
            p.circle_stroke(c, 5.5, thin);
        }
        Tool::OrbitalD => {
            for (dx, dy, filled) in [(-3.2, -3.2, true), (3.2, 3.2, true), (-3.2, 3.2, false), (3.2, -3.2, false)] {
                let pc = Pos2::new(c.x + dx, c.y + dy);
                if filled { p.circle_filled(pc, 2.5, Color32::BLACK); } else { p.circle_stroke(pc, 2.5, thin); }
            }
        }
        Tool::BoldBond => {
            p.line_segment([Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x1 - 2.0, y1 - 2.0)], bold);
        }
        Tool::ShapeRect => {
            let sr = Rect::from_min_max(Pos2::new(x0 + 1.0, y0 + 3.0), Pos2::new(x1 - 2.0, y1 - 3.0));
            p.rect_filled(sr.translate(Vec2::splat(1.5)), 0.0, Color32::from_gray(160));
            p.rect_filled(sr, 0.0, Color32::WHITE);
            p.rect_stroke(sr, 0.0, thin, egui::StrokeKind::Inside);
        }
        Tool::ShapeOval => {
            p.circle_stroke(c, 6.2, thin);
        }
        Tool::ShapeLine => {
            p.line_segment([Pos2::new(x0 + 1.5, y1 - 2.5), Pos2::new(x1 - 1.5, y0 + 2.5)], s);
        }
        Tool::WedgeBond => {
            p.add(Shape::convex_polygon(
                vec![Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x1 - 1.0, y1 - 4.0), Pos2::new(x1 - 4.0, y1 - 1.0)],
                Color32::BLACK,
                Stroke::NONE,
            ));
        }
        Tool::Bracket => {
            p.line_segment([Pos2::new(x0 + 5.0, y0 + 2.0), Pos2::new(x0 + 2.0, y0 + 2.0)], s);
            p.line_segment([Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x0 + 2.0, y1 - 2.0)], s);
            p.line_segment([Pos2::new(x0 + 2.0, y1 - 2.0), Pos2::new(x0 + 5.0, y1 - 2.0)], s);
            p.line_segment([Pos2::new(x1 - 5.0, y0 + 2.0), Pos2::new(x1 - 2.0, y0 + 2.0)], s);
            p.line_segment([Pos2::new(x1 - 2.0, y0 + 2.0), Pos2::new(x1 - 2.0, y1 - 2.0)], s);
            p.line_segment([Pos2::new(x1 - 2.0, y1 - 2.0), Pos2::new(x1 - 5.0, y1 - 2.0)], s);
        }
        Tool::RoundBracket => {
            p.text(c, egui::Align2::CENTER_CENTER, "( )", FontId::proportional(13.0), Color32::BLACK);
        }
        Tool::CurlyBracket => {
            p.text(c, egui::Align2::CENTER_CENTER, "{ }", FontId::proportional(13.0), Color32::BLACK);
        }
        Tool::HollowWedgeBond => {
            let pts = [Pos2::new(x0 + 2.0, y0 + 2.0), Pos2::new(x1 - 1.0, y1 - 4.0), Pos2::new(x1 - 4.0, y1 - 1.0)];
            draw_polygon(p, &pts, thin);
        }
        Tool::ChemSymbol => {
            p.circle_stroke(c, 5.5, thin);
            p.line_segment([Pos2::new(c.x - 3.0, c.y), Pos2::new(c.x + 3.0, c.y)], thin);
            p.line_segment([Pos2::new(c.x, c.y - 3.0), Pos2::new(c.x, c.y + 3.0)], thin);
        }
        Tool::SymbolMinus => {
            p.circle_stroke(c, 5.5, thin);
            p.line_segment([Pos2::new(c.x - 3.0, c.y), Pos2::new(c.x + 3.0, c.y)], thin);
        }
        Tool::SymbolRadical => {
            p.circle_filled(c, 2.5, Color32::BLACK);
        }
        Tool::SymbolLonePair => {
            p.circle_filled(Pos2::new(c.x - 2.8, c.y), 1.8, Color32::BLACK);
            p.circle_filled(Pos2::new(c.x + 2.8, c.y), 1.8, Color32::BLACK);
        }
        Tool::WavyBond => {
            let pts = [
                Pos2::new(x0 + 2.0, y0 + 3.0),
                Pos2::new(x0 + 6.0, y0 + 2.0),
                Pos2::new(c.x - 1.0, c.y + 1.0),
                Pos2::new(c.x + 3.0, c.y),
                Pos2::new(x1 - 2.0, y1 - 2.0),
            ];
            for w in pts.windows(2) {
                p.line_segment([w[0], w[1]], thin);
            }
        }
        Tool::ReactionMap => {
            p.text(c, egui::Align2::CENTER_CENTER, "A→A", FontId::proportional(8.5), Color32::BLACK);
        }
        Tool::Table => {
            let tr = r.shrink(2.0);
            p.rect_stroke(tr, 0.0, thin, egui::StrokeKind::Inside);
            p.line_segment([Pos2::new(c.x, tr.min.y), Pos2::new(c.x, tr.max.y)], thin);
            p.line_segment([Pos2::new(tr.min.x, c.y), Pos2::new(tr.max.x, c.y)], thin);
        }
        Tool::TlcPlate => {
            let tr = Rect::from_center_size(c, Vec2::new(11.0, 15.0));
            p.rect_stroke(tr, 0.0, thin, egui::StrokeKind::Inside);
            p.line_segment([Pos2::new(tr.min.x, tr.max.y - 3.0), Pos2::new(tr.max.x, tr.max.y - 3.0)], thin);
            p.circle_filled(Pos2::new(c.x - 2.0, c.y), 1.2, Color32::BLACK);
            p.circle_filled(Pos2::new(c.x + 2.0, c.y - 3.0), 1.2, Color32::BLACK);
        }
        Tool::AcyclicChain => {
            let pts = [
                Pos2::new(x0 + 1.0, c.y + 1.0),
                Pos2::new(x0 + 5.0, c.y - 3.0),
                Pos2::new(x0 + 9.0, c.y + 1.0),
                Pos2::new(x0 + 13.0, c.y - 3.0),
            ];
            for w in pts.windows(2) {
                p.line_segment([w[0], w[1]], thin);
            }
            p.text(Pos2::new(x1 - 2.0, y1 - 2.0), egui::Align2::RIGHT_BOTTOM, "n", FontId::proportional(8.5), Color32::BLACK);
        }
        Tool::Templates => {
            p.circle_filled(Pos2::new(c.x, y0 + 4.0), 2.5, Color32::BLACK);
            p.line_segment([Pos2::new(c.x, y0 + 6.0), Pos2::new(c.x, c.y + 2.0)], bold);
            p.rect_stroke(
                Rect::from_min_max(Pos2::new(x0 + 2.0, c.y + 2.0), Pos2::new(x1 - 2.0, y1 - 2.0)),
                1.0,
                s,
                egui::StrokeKind::Inside,
            );
        }
        Tool::Cyclopropane => draw_polygon(p, &poly_pts(c, 6.5, 3, -90.0), thin),
        Tool::Cyclobutane => draw_polygon(p, &poly_pts(c, 7.0, 4, 45.0), thin),
        Tool::Cyclopentane => draw_polygon(p, &poly_pts(c, 7.0, 5, -90.0), thin),
        Tool::Cyclohexane => draw_polygon(p, &poly_pts(c, 7.0, 6, -90.0), thin),
        Tool::Cycloheptane => draw_polygon(p, &poly_pts(c, 7.2, 7, -90.0), thin),
        Tool::Cyclooctane => draw_polygon(p, &poly_pts(c, 7.2, 8, 22.5), thin),
        Tool::Chair1 => {
            let pts = [
                Pos2::new(x0 + 1.0, c.y - 1.0),
                Pos2::new(c.x, c.y - 3.0),
                Pos2::new(x1 - 1.0, c.y - 1.0),
                Pos2::new(x1 - 3.0, c.y + 3.0),
                Pos2::new(c.x - 1.0, c.y + 4.0),
                Pos2::new(x0 + 3.0, c.y + 2.0),
            ];
            draw_polygon(p, &pts, thin);
        }
        Tool::Chair2 => {
            let pts = [
                Pos2::new(x0 + 1.0, c.y + 1.0),
                Pos2::new(c.x, c.y - 3.0),
                Pos2::new(x1 - 1.0, c.y + 1.0),
                Pos2::new(x1 - 3.0, c.y + 4.0),
                Pos2::new(c.x, c.y + 2.0),
                Pos2::new(x0 + 3.0, c.y + 4.0),
            ];
            draw_polygon(p, &pts, thin);
        }
        Tool::Cyclopentadiene => {
            let pts = poly_pts(c, 7.0, 5, -90.0);
            draw_polygon(p, &pts, thin);
            let inner = poly_pts(c, 4.8, 5, -90.0);
            p.line_segment([inner[1], inner[2]], thin);
            p.line_segment([inner[3], inner[4]], thin);
        }
        Tool::Benzene => {
            let pts = poly_pts(c, 7.2, 6, -90.0);
            draw_polygon(p, &pts, thin);
            let inner = poly_pts(c, 5.0, 6, -90.0);
            for i in [0, 2, 4] {
                p.line_segment([inner[i], inner[(i + 1) % 6]], thin);
            }
        }
    }
}
