use std::{
    fs::{self, File},
    io::BufReader,
    path::PathBuf,
};

use anyhow::{Context, anyhow, bail, ensure};
use common::Anonymize;

use crate::graphics::{
    asset_loader::AssetLoadContext,
    meshes::{MeshMaterial, MeshMaterialDefinition},
};

pub struct PngArray {
    pub path: PathBuf,
    pub size: usize,
}

impl PngArray {
    async fn load_inner(self, context: &skid_steer::Context<'_>) -> anyhow::Result<MeshMaterial> {
        tracing::trace!("Started loading png array");
        let ctx: &AssetLoadContext = context.get().unwrap();
        let full_path = ctx
            .find_asset(&self.path)
            .ok_or_else(|| anyhow!("{} not found", self.path.anonymize().display()))?;
        let mut paths = fs::read_dir(&full_path)
            .with_context(|| format!("reading {}", full_path.anonymize().display()))?
            .map(|x| x.map(|x| x.path()))
            .collect::<Result<Vec<_>, _>>()
            .with_context(|| format!("reading {}", full_path.anonymize().display()))?;
        if paths.len() < self.size {
            bail!(
                "{}: expected {} textures, found {}",
                full_path.anonymize().display(),
                self.size,
                paths.len()
            );
        }
        paths.sort();
        paths.truncate(self.size);
        let mut dims: Option<(u32, u32)> = None;
        let mut image_data: Vec<u8> = Vec::new();
        for (i, path) in paths.iter().enumerate() {
            tracing::trace!(layer=i, path=%path.anonymize().display(), "loading");
            let file = File::open(path)
                .with_context(|| format!("reading {}", path.anonymize().display()))?;
            let decoder = png::Decoder::new(BufReader::new(file));
            let mut reader = decoder
                .read_info()
                .with_context(|| format!("decoding {}", path.anonymize().display()))?;
            let info = reader.info();
            let step_size = info.width as usize * info.height as usize * 4;
            ensure!(info.color_type == png::ColorType::Rgba);
            ensure!(info.bit_depth == png::BitDepth::Eight);
            ensure!(reader.output_buffer_size() == Some(step_size));
            if let Some(dims) = dims {
                if dims != (info.width, info.height) {
                    bail!(
                        "inconsistent dimensions: expected {}x{}, got {}x{}",
                        dims.0,
                        dims.1,
                        info.width,
                        info.height
                    );
                }
            } else {
                dims = Some((info.width, info.height));
                image_data.resize(step_size * self.size, 0);
            }
            reader
                .next_frame(&mut image_data[i * step_size..(i + 1) * step_size])
                .with_context(|| format!("decoding {}", path.anonymize().display()))?;
        }
        let (width, height) = dims.unwrap();
        Ok(MeshMaterial::from_definition(
            ctx,
            MeshMaterialDefinition {
                width,
                height,
                array_layers: self.size as u32,
                srgb_rgba_color_data: image_data,
            },
        )
        .await)
    }
}

impl skid_steer::Source for PngArray {
    type Output = MeshMaterial;

    async fn load(self, context: &skid_steer::Context<'_>) -> Option<MeshMaterial> {
        self.load_inner(context)
            .await
            .inspect_err(|e| tracing::error!("{}", e))
            .ok()
    }

    fn free(mut output: Self::Output, context: &skid_steer::Context) {
        let ctx: &AssetLoadContext = context.get().unwrap();
        unsafe { output.destroy(ctx.device()) };
    }
}
