use common::{math, proto::Position};

pub struct LocalCharacterController {
    /// The last extrapolated inter-frame view position, used for rendering and gravity-specific
    /// orientation computations
    position: Position,

    /// The quaternion adjustment to the character position to represent its actual apparent orientation
    orientation: na::UnitQuaternion<f32>,

    /// The up vector relative to both position and orientation
    up: na::UnitVector3<f32>,
}

impl LocalCharacterController {
    pub fn new() -> Self {
        LocalCharacterController {
            position: Position::origin(),
            orientation: na::UnitQuaternion::identity(),
            up: na::Vector::z_axis(),
        }
    }

    /// Get the current position with orientation applied to it
    pub fn oriented_position(&self) -> Position {
        Position {
            node: self.position.node,
            local: self.position.local * self.orientation.to_homogeneous(),
        }
    }

    pub fn orientation(&self) -> na::UnitQuaternion<f32> {
        self.orientation
    }

    /// Updates the LocalCharacter based on outside information. Note that the `up` parameter is relative
    /// only to `position`, not the character's orientation.
    pub fn update_position(
        &mut self,
        position: Position,
        up: na::UnitVector3<f32>,
        preserve_up_alignment: bool,
    ) {
        if preserve_up_alignment {
            // Rotate the character orientation to stay consistent with changes in gravity
            self.orientation = math::rotation_between_axis(&self.up, &up, 1e-5)
                .unwrap_or(na::UnitQuaternion::identity())
                * self.orientation;
        }

        self.position = position;

        // Change of coordinates
        self.up = self.orientation.conjugate() * up;
    }

    /// Rotates the camera's view by locally adding pitch and yaw.
    pub fn look_free(&mut self, delta_yaw: f32, delta_pitch: f32, delta_roll: f32) {
        self.orientation *= na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw)
            * na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), delta_pitch)
            * na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), delta_roll);
    }

    /// Rotates the camera's view with standard first-person walking simulator mouse controls. This function
    /// is designed to be flexible enough to work with any starting orientation, but it works best when the
    /// camera is level, not rolled to the left or right.
    pub fn look_level(&mut self, delta_yaw: f32, delta_pitch: f32) {
        // Handle yaw. This is as simple as rotating the view about the up vector
        self.orientation *= na::UnitQuaternion::from_axis_angle(&self.up, delta_yaw);

        // Handling pitch is more compicated because the view angle needs to be capped. The rotation axis
        // is the camera's local x-axis (left-right axis). If the camera is level, this axis is perpendicular
        // to the up vector.

        // We need to know the current pitch to properly cap pitch changes, and this is only well-defined
        // if the pitch axis is not too similar to the up vector, so we skip applying pitch changes if this
        // isn't the case.
        if self.up.x.abs() < 0.9 {
            // Compute the current pitch by ignoring the x-component of the up vector and assuming the camera
            // is level.
            let current_pitch = -self.up.z.atan2(self.up.y);
            let mut target_pitch = current_pitch + delta_pitch;
            if delta_pitch > 0.0 {
                target_pitch = target_pitch
                    .min(std::f32::consts::FRAC_PI_2) // Cap the view angle at looking straight up
                    .max(current_pitch); // But if already upside-down, don't make any corrections.
            } else {
                target_pitch = target_pitch
                    .max(-std::f32::consts::FRAC_PI_2) // Cap the view angle at looking straight down
                    .min(current_pitch); // But if already upside-down, don't make any corrections.
            }

            self.orientation *= na::UnitQuaternion::from_axis_angle(
                &na::Vector3::x_axis(),
                target_pitch - current_pitch,
            );
        }
    }

    /// Instantly updates the current orientation quaternion to make the camera level. This function
    /// is designed to be numerically stable for any camera orientation.
    pub fn align_to_gravity(&mut self) {
        if self.up.z.abs() < 0.9 {
            // If facing not too vertically, roll the camera to make it level.
            let delta_roll = -self.up.x.atan2(self.up.y);
            self.orientation *=
                na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), delta_roll);
        } else if self.up.y > 0.0 {
            // Otherwise, if not upside-down, pan the camera to make it level.
            let delta_yaw = (self.up.x / self.up.z).atan();
            self.orientation *=
                na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw);
        } else {
            // Otherwise, rotate the camera to look straight up or down.
            self.orientation *= na::UnitQuaternion::rotation_between(
                &(na::Vector3::z() * self.up.z.signum()),
                &self.up,
            )
            .unwrap();
        }
    }

    /// Returns an orientation quaternion that is as faithful as possible to the current orientation quaternion
    /// while being restricted to ensuring the view is level and does not look up or down. This function's main
    /// purpose is to figure out what direction the character should go when a movement key is pressed.
    pub fn get_horizontal_orientation(&mut self) -> na::UnitQuaternion<f32> {
        let forward = if self.up.x.abs() < 0.9 {
            // Rotate the local forward vector about the locally horizontal axis until it is horizontal
            na::Vector3::new(0.0, -self.up.z, self.up.y)
        } else {
            // Project the local forward vector to the level plane
            na::Vector3::z() - self.up.into_inner() * self.up.z
        };

        self.orientation * na::UnitQuaternion::face_towards(&forward, &self.up)
    }

    pub fn renormalize_orientation(&mut self) {
        self.orientation.renormalize_fast();
    }
}
