/// Rotates the camera's view by locally adding pitch and yaw.
fn look_free(
    orientation: &mut na::UnitQuaternion<f32>,
    delta_yaw: f32,
    delta_pitch: f32,
    delta_roll: f32,
) {
    *orientation *= na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw)
        * na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), delta_pitch)
        * na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), delta_roll);
}

/// Rotates the camera's view with standard first-person walking simulator mouse controls. This function
/// is designed to be flexible enough to work with any starting orientation, but it works best when the
/// camera is level, not rolled to the left or right.
fn look_level(
    orientation: &mut na::UnitQuaternion<f32>,
    up: &na::UnitVector3<f32>,
    delta_yaw: f32,
    delta_pitch: f32,
) {
    let up = orientation.conjugate() * up;

    // Handle yaw. This is as simple as rotating the view about the up vector
    *orientation *= na::UnitQuaternion::from_axis_angle(&up, delta_yaw);

    // Handling pitch is more compicated because the view angle needs to be capped. The rotation axis
    // is the camera's local x-axis (left-right axis). If the camera is level, this axis is perpendicular
    // to the up vector.

    // We need to know the current pitch to properly cap pitch changes, and this is only well-defined
    // if the pitch axis is not too similar to the up vector, so we skip applying pitch changes if this
    // isn't the case.
    if up.x.abs() < 0.9 {
        // Compute the current pitch by ignoring the x-component of the up vector and assuming the camera
        // is level.
        let current_pitch = -up.z.atan2(up.y);
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

        *orientation *= na::UnitQuaternion::from_axis_angle(
            &na::Vector3::x_axis(),
            target_pitch - current_pitch,
        );
    }
}

/// Instantly updates the current orientation quaternion to make the camera level. This function
/// is designed to be numerically stable for any camera orientation.
fn align_to_gravity(orientation: &mut na::UnitQuaternion<f32>, up: &na::UnitVector3<f32>) {
    let up = orientation.conjugate() * up;

    if up.z.abs() < 0.9 {
        // If facing not too vertically, roll the camera to make it level.
        let delta_roll = -up.x.atan2(up.y);
        *orientation *= na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), delta_roll);
    } else if up.y > 0.0 {
        // Otherwise, if not upside-down, pan the camera to make it level.
        let delta_yaw = (up.x / up.z).atan();
        *orientation *= na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw);
    } else {
        // Otherwise, rotate the camera to look straight up or down.
        *orientation *=
            na::UnitQuaternion::rotation_between(&(na::Vector3::z() * up.z.signum()), &up).unwrap();
    }
}

/// Returns an orientation quaternion that is as faithful as possible to the current orientation quaternion
/// while being restricted to ensuring the view is level and does not look up or down.
fn get_horizontal_orientation(
    orientation: &mut na::UnitQuaternion<f32>,
    up: &na::UnitVector3<f32>,
) -> na::UnitQuaternion<f32> {
    let up = orientation.conjugate() * up;

    let forward = if up.x.abs() < 0.9 {
        // Rotate the local forward vector about the locally horizontal axis until it is horizontal
        na::Vector3::new(0.0, -up.z, up.y)
    } else {
        // Project the local forward vector to the level plane
        na::Vector3::z() - up.into_inner() * up.z
    };

    *orientation * na::UnitQuaternion::face_towards(&forward, &up)
}
