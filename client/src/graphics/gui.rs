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
                    let material_count_string = if sim.cfg.gameplay_enabled {
                        sim.inventory_contents_matching_material(sim.selected_material())
                            .len()
                            .to_string()
                    } else {
                        "\u{221E}".to_string() // \u{221E} is the infinity synbol
                    };
                    label(format!(
                        "Selected material: {:?} (\u{00D7}{})", // \u{00D7} is the multiplication synbol
                        sim.selected_material(),
                        material_count_string
                    ));
                });
            });
        });
    }
}
