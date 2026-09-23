use common::math::MVector;
use yakui::{
    Alignment, Color, Response, align, canvas, colored_box, colored_box_container, label, pad,
    widgets::{CanvasResponse, Pad},
};

use crate::{Sim, graphics::Frustum};

pub struct GuiState {
    show_gui: bool,
    show_home_waypoint: bool,
}

impl GuiState {
    pub fn new() -> Self {
        GuiState {
            show_gui: true,
            show_home_waypoint: false,
        }
    }

    /// Toggles whether the GUI is shown
    pub fn toggle_gui(&mut self) {
        self.show_gui = !self.show_gui;
    }

    /// Toggles whether the home waypoint is shown
    pub fn toggle_show_home_waypoint(&mut self) {
        self.show_home_waypoint = !self.show_home_waypoint;
    }

    /// Prepare the GUI for rendering. This should be called between
    /// Yakui::start and Yakui::finish.
    pub fn run(&self, sim: &Sim, frustum: Frustum) {
        if self.show_home_waypoint {
            // Choose a waypoint that appears to be at the origin of the world. 10 steps is sufficient for
            // the waypoint icon's location to be indistinguishable from the correct location.
            let waypoint_location = sim.view_relative_origin(10);

            pad(Pad::all(8.0), || {
                waypoint_icon(waypoint_location, frustum);
            });
        }

        if !self.show_gui {
            return;
        }

        align(Alignment::CENTER, || {
            colored_box(Color::BLACK.with_alpha(0.9), [5.0, 5.0]);
        });

        align(Alignment::TOP_LEFT, || {
            pad(Pad::all(8.0), || {
                colored_box_container(Color::BLACK.with_alpha(0.7), || {
                    let material_count_string = if sim.cfg.gameplay_enabled {
                        sim.count_inventory_entities_matching_material(sim.selected_material())
                            .to_string()
                    } else {
                        "∞".to_string()
                    };
                    label(format!(
                        "Selected material: {:?} (×{})",
                        sim.selected_material(),
                        material_count_string
                    ));
                });
            });
        });
    }
}

fn waypoint_icon(world_vector: MVector<f32>, frustum: Frustum) -> Response<CanvasResponse> {
    canvas(move |ctx| {
        let viewport = ctx.layout.get(ctx.dom.root()).unwrap().rect;
        let panel_region = ctx.layout.get(ctx.dom.current()).unwrap().rect;

        let c1 = [1.0, 1.0, 1.0, 1.0]; // Fill color
        let c2 = [0.0, 0.0, 0.0, 1.0]; // Outline color
        let r1 = 8.0; // Fill radius
        let r2 = 12.0; // Outline radius

        let bounds = yakui::Rect::from_pos_size(
            panel_region.pos() + yakui::Vec2::new(r2, r2),
            panel_region.size() - yakui::Vec2::new(r2, r2) * 2.0,
        );

        let target = get_layout_coords_of_world_vector(world_vector, frustum, &viewport, &bounds);
        let vertices = [
            ([r1, 0.0], c1),
            ([0.0, r1], c1),
            ([-r1, 0.0], c1),
            ([0.0, -r1], c1),
            ([r2, 0.0], c2),
            ([0.0, r2], c2),
            ([-r2, 0.0], c2),
            ([0.0, -r2], c2),
        ]
        .map(|(v, c)| {
            new_yakui_vertex(
                na::Vector2::from(v) + target.layout_center + target.capped_offset,
                [0.0, 0.0],
                c,
            )
        });

        let indices = [4, 5, 6, 4, 6, 7, 0, 1, 2, 0, 2, 3];
        let mesh = yakui::paint::PaintMesh::new(vertices, indices);
        ctx.paint.add_mesh(mesh);
    })
}

/// Returns the coordinates, along with whether the coordinates were clamped to the fit in the bounds
fn get_layout_coords_of_world_vector(
    point: MVector<f32>,
    frustum: Frustum,
    viewport: &yakui::Rect,
    bounds: &yakui::Rect,
) -> WorldPointInLayout {
    let viewport_center: na::Vector2<f32> =
        (viewport.pos() + viewport.size() * 0.5).to_array().into();
    let viewport_scale: na::Vector2<f32> = (viewport.size() * 0.5).to_array().into();
    let bounds_center: na::Vector2<f32> = (bounds.pos() + bounds.size() * 0.5).to_array().into();
    let bounds_scale: na::Vector2<f32> = (bounds.size() * 0.5).to_array().into();

    let raw_projected_point = frustum.projection().matrix() * na::Vector4::from(point);

    // Projected point location in viewport coordinates relative to `bounds_center`
    let projected_point = na::Vector3::new(
        raw_projected_point.x * viewport_scale.x
            + (viewport_center.x - bounds_center.x) * raw_projected_point.w,
        raw_projected_point.y * viewport_scale.y
            + (viewport_center.y - bounds_center.y) * raw_projected_point.w,
        raw_projected_point.w,
    );

    // Relative to `bounds_center`, how much we need to scale a point on the boundary to reach `projected_point.xy()`
    let projected_point_xy_to_boundary_factor = projected_point
        .xy()
        .component_div(&bounds_scale)
        .abs()
        .max();
    if projected_point.z <= 0.0 || projected_point_xy_to_boundary_factor > projected_point.z {
        let offset = if projected_point.xy().norm_squared() < 1e-16 {
            // Choose an arbitrary vector on the boundary if too close to 0
            na::Vector2::new(0.0, bounds_scale.y)
        } else {
            projected_point.xy() / projected_point_xy_to_boundary_factor
        };
        WorldPointInLayout {
            layout_center: bounds_center,
            capped_offset: offset,
            in_bounds: false,
        }
    } else {
        WorldPointInLayout {
            layout_center: bounds_center,
            capped_offset: projected_point.xy() / projected_point.z,
            in_bounds: true,
        }
    }
}

struct WorldPointInLayout {
    layout_center: na::Vector2<f32>,
    capped_offset: na::Vector2<f32>,

    #[expect(unused)] // We want to keep track of this even if we're not currently using it.
    in_bounds: bool,
}

fn new_yakui_vertex(
    position: impl Into<na::Vector2<f32>>,
    texcoord: impl Into<na::Vector2<f32>>,
    color: impl Into<na::Vector4<f32>>,
) -> yakui::paint::Vertex {
    let position = position.into();
    let texcoord = texcoord.into();
    let color = color.into();
    yakui::paint::Vertex::new(
        yakui::Vec2::new(position.x, position.y),
        yakui::Vec2::new(texcoord.x, texcoord.y),
        yakui::Vec4::new(color.x, color.y, color.z, color.w),
    )
}
