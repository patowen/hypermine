use yakui::{
    align, colored_box, colored_box_container, label, pad, widgets::Pad, Alignment, Color,
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
    /// Yakui::start and Yakui::finish.
    pub fn run(&self, sim: &Sim) {
        if !self.show_gui {
            return;
        }

        align(Alignment::CENTER, || {
            colored_box(Color::BLACK.with_alpha(0.9), [5.0, 5.0]);
        });

        align(Alignment::TOP_LEFT, || {
            pad(Pad::all(8.0), || {
                colored_box_container(Color::BLACK.with_alpha(0.7), || {
                    let selected_material_string = if let Some(material) = sim.selected_material() {
                        format!(
                            "Selected material: {:?} (\u{00D7}{})",
                            material,
                            sim.inventory_contents_matching_material(material).len()
                        )
                    } else {
                        "Selected material: None".to_string()
                    };
                    label(selected_material_string);
                });
            });
        });
    }
}
