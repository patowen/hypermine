use yakui::{
    align, colored_box, colored_box_container, label, pad,
    paint::{PaintMesh, Vertex},
    widget::{LayoutContext, PaintContext, Widget},
    widgets::Pad,
    Alignment, Color, Constraints, Rect, TextureId, Vec2,
};

use crate::Sim;

pub struct GuiState {
    show_gui: bool,
}

impl GuiState {
    pub fn new() -> Self {
        GuiState { show_gui: true }
    }

    /// Toggles whether the GUI is shown
    pub fn toggle_gui(&mut self) {
        self.show_gui = !self.show_gui;
    }

    /// Prepare the GUI for rendering. This should be called between
    /// Yakui::start and yakui::finish.
    pub fn prepare_gui(&self, sim: &Sim) {
        if !self.show_gui {
            return;
        }

        align(Alignment::CENTER, || {
            colored_box(Color::BLACK.with_alpha(0.9), [5.0, 5.0]);
        });

        align(Alignment::TOP_LEFT, || {
            pad(Pad::all(8.0), || {
                colored_box_container(Color::BLACK.with_alpha(0.7), || {
                    label(format!("Selected material: {:?}", sim.selected_material()));
                });
            });
        });
    }
}

#[derive(Debug)]
struct ItemIconWidget(Option<TextureId>);

impl ItemIconWidget {
    const CUBE_VERTICES: [([f32; 2], [f32; 2], f32); 12] = {
        let r3 = 1.7320508f32; // Square root of 3
        [
            ([0.0, -2.0], [0.0, 0.0], 1.0),
            ([-r3, -1.0], [0.0, 1.0], 1.0),
            ([0.0, 0.0], [1.0, 1.0], 1.0),
            ([r3, -1.0], [1.0, 0.0], 1.0),
            ([-r3, -1.0], [0.0, 0.0], 0.9),
            ([-r3, 1.0], [0.0, 1.0], 0.9),
            ([0.0, 2.0], [1.0, 1.0], 0.9),
            ([0.0, 0.0], [1.0, 0.0], 0.9),
            ([0.0, 0.0], [0.0, 0.0], 0.8),
            ([0.0, 2.0], [0.0, 1.0], 0.8),
            ([r3, 1.0], [1.0, 1.0], 0.8),
            ([r3, -1.0], [1.0, 0.0], 0.8),
        ]
    };

    const CUBE_INDICES: [u16; 18] = [0, 1, 2, 3, 0, 2, 4, 5, 6, 7, 4, 6, 8, 9, 10, 11, 8, 10];
}

impl Widget for ItemIconWidget {
    type Props<'a> = Option<TextureId>;

    type Response = ();

    fn new() -> Self {
        Self(None)
    }

    fn update(&mut self, props: Self::Props<'_>) -> Self::Response {
        self.0 = props;
    }

    fn layout(&self, _ctx: LayoutContext<'_>, input: Constraints) -> Vec2 {
        input.constrain_min(Vec2::splat(200.0))
    }

    fn paint(&self, ctx: PaintContext<'_>) {
        let layout_node = ctx.layout.get(ctx.dom.current()).unwrap();

        if let Some(texture) = self.0 {
            let vertices = Self::CUBE_VERTICES.map(|(p, t, c)| {
                Vertex::new(
                    Vec2::new(p[0] + 2.0, p[1] + 2.0) * layout_node.rect.size() / 4.0
                        + layout_node.rect.pos(),
                    t,
                    [c, c, c, 1.0],
                )
            });
            let mut hexagon = PaintMesh::new(vertices, Self::CUBE_INDICES);
            hexagon.texture = Some((texture, Rect::ONE));
            ctx.paint.add_mesh(hexagon);
        }
    }
}
