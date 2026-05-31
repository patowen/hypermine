use std::{
    fs::{self, File},
    io::BufReader,
    path::PathBuf,
    sync::Arc,
};

use anyhow::{Context, anyhow, bail};
use ash::vk;
use common::Anonymize;
use lahar::{DedicatedImage, parallel_queue};
use tracing::trace;

use crate::{
    Config,
    asset_loader::{Allocation, StagingRing},
    graphics::Base,
    loader::{LoadCtx, LoadFuture, Loadable},
};

pub struct PngArray {
    pub path: PathBuf,
    pub size: usize,
}

impl PngArray {
    async fn load_inner(self, context: &skid_steer::Context<'_>) -> anyhow::Result<DedicatedImage> {
        let cfg: &Config = context.get().unwrap();
        let handle: &parallel_queue::Handle = context.get().unwrap();
        let gfx: &Arc<Base> = context.get().unwrap();
        let staging_buffer: &StagingRing = context.get().unwrap();
        let full_path = cfg
            .find_asset(&self.path)
            .ok_or_else(|| anyhow!("{} not found", self.path.anonymize().display()))
            .inspect_err(|e| tracing::error!("{}", e))?;
        let mut paths = fs::read_dir(&full_path)
            .with_context(|| format!("reading {}", full_path.anonymize().display()))?
            .map(|x| x.map(|x| x.path()))
            .collect::<Result<Vec<_>, _>>()
            .with_context(|| format!("reading {}", full_path.anonymize().display()))?;
        if paths.is_empty() {
            bail!("{} is empty", full_path.anonymize().display());
        }
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
        let work = unsafe { handle.begin(&gfx.device) };
        let mut mem2: Option<Allocation> = None;
        for (i, path) in paths.iter().enumerate() {
            tracing::trace!(layer=i, path=%path.anonymize().display(), "loading");
            let file = File::open(path)
                .with_context(|| format!("reading {}", path.anonymize().display()))?;
            let decoder = png::Decoder::new(BufReader::new(file));
            let mut reader = decoder
                .read_info()
                .with_context(|| format!("decoding {}", path.anonymize().display()))?;
            let info = reader.info();
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
                mem2 = Some(unsafe {
                    staging_buffer
                        .alloc(
                            info.width as usize * info.height as usize * 4 * self.size,
                            1, /* TODO: Is an alignment of 1 safe? */
                            work.time().get(),
                        )
                        .expect("TODO")
                });
            }
            let mem2 = mem2.as_mut().unwrap();
            let step_size = info.width as usize * info.height as usize * 4;
            reader
                .next_frame(&mut mem2.bytes[i * step_size..(i + 1) * step_size])
                .with_context(|| format!("decoding {}", path.anonymize().display()))?;
        }
        let (width, height) = dims.unwrap();
        let mem2 = mem2.unwrap();
        unsafe {
            let image = DedicatedImage::new(
                &gfx.device,
                &gfx.memory_properties,
                &vk::ImageCreateInfo::default()
                    .image_type(vk::ImageType::TYPE_2D)
                    .format(vk::Format::R8G8B8A8_SRGB)
                    .extent(vk::Extent3D {
                        width,
                        height,
                        depth: 1,
                    })
                    .mip_levels(1)
                    .array_layers(self.size as u32)
                    .samples(vk::SampleCountFlags::TYPE_1)
                    .usage(vk::ImageUsageFlags::SAMPLED | vk::ImageUsageFlags::TRANSFER_DST),
            );

            let range = vk::ImageSubresourceRange {
                aspect_mask: vk::ImageAspectFlags::COLOR,
                base_mip_level: 0,
                level_count: 1,
                base_array_layer: 0,
                layer_count: self.size as u32,
            };
            let src = staging_buffer.buffer();
            let buffer_offset = mem2.offset;
            let dst = image.handle;

            gfx.device.cmd_pipeline_barrier(
                work.cmd(),
                vk::PipelineStageFlags::TOP_OF_PIPE,
                vk::PipelineStageFlags::TRANSFER,
                vk::DependencyFlags::default(),
                &[],
                &[],
                &[vk::ImageMemoryBarrier::default()
                    .dst_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                    .src_queue_family_index(vk::QUEUE_FAMILY_IGNORED)
                    .dst_queue_family_index(vk::QUEUE_FAMILY_IGNORED)
                    .old_layout(vk::ImageLayout::UNDEFINED)
                    .new_layout(vk::ImageLayout::TRANSFER_DST_OPTIMAL)
                    .image(dst)
                    .subresource_range(range)],
            );
            gfx.device.cmd_copy_buffer_to_image(
                work.cmd(),
                src,
                dst,
                vk::ImageLayout::TRANSFER_DST_OPTIMAL,
                &[vk::BufferImageCopy {
                    buffer_offset,
                    image_subresource: vk::ImageSubresourceLayers {
                        aspect_mask: vk::ImageAspectFlags::COLOR,
                        mip_level: 0,
                        base_array_layer: 0,
                        layer_count: range.layer_count,
                    },
                    image_extent: vk::Extent3D {
                        width,
                        height,
                        depth: 1,
                    },
                    ..Default::default()
                }],
            );
            gfx.device.cmd_pipeline_barrier(
                work.cmd(),
                vk::PipelineStageFlags::TRANSFER,
                vk::PipelineStageFlags::FRAGMENT_SHADER,
                vk::DependencyFlags::default(),
                &[],
                &[],
                &[vk::ImageMemoryBarrier::default()
                    .src_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                    .dst_access_mask(vk::AccessFlags::SHADER_READ)
                    .src_queue_family_index(gfx.queue_family)
                    .dst_queue_family_index(gfx.queue_family)
                    .old_layout(vk::ImageLayout::TRANSFER_DST_OPTIMAL)
                    .new_layout(vk::ImageLayout::SHADER_READ_ONLY_OPTIMAL)
                    .image(dst)
                    .subresource_range(range)],
            );
            work.end();
            //TODO: Await work completion

            trace!(
                width = width,
                height = height,
                path = %full_path.anonymize().display(),
                "loaded array"
            );
            Ok(image)
        }
    }
}

impl Loadable for PngArray {
    type Output = DedicatedImage;

    fn load(self, handle: &LoadCtx) -> LoadFuture<'_, Self::Output> {
        Box::pin(async move {
            let full_path = handle
                .cfg
                .find_asset(&self.path)
                .ok_or_else(|| anyhow!("{} not found", self.path.anonymize().display()))?;
            let mut paths = fs::read_dir(&full_path)
                .with_context(|| format!("reading {}", full_path.anonymize().display()))?
                .map(|x| x.map(|x| x.path()))
                .collect::<Result<Vec<_>, _>>()
                .with_context(|| format!("reading {}", full_path.anonymize().display()))?;
            if paths.is_empty() {
                bail!("{} is empty", full_path.anonymize().display());
            }
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
            let mut mem = None;
            for (i, path) in paths.iter().enumerate() {
                trace!(layer=i, path=%path.anonymize().display(), "loading");
                let file = File::open(path)
                    .with_context(|| format!("reading {}", path.anonymize().display()))?;
                let decoder = png::Decoder::new(BufReader::new(file));
                let mut reader = decoder
                    .read_info()
                    .with_context(|| format!("decoding {}", path.anonymize().display()))?;
                let info = reader.info();
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
                    mem = Some(
                        handle
                            .staging
                            .alloc(info.width as usize * info.height as usize * 4 * self.size)
                            .await
                            .ok_or_else(|| {
                                anyhow!(
                                    "{}: image array too large",
                                    full_path.anonymize().display()
                                )
                            })?,
                    );
                }
                let mem = mem.as_mut().unwrap();
                let step_size = info.width as usize * info.height as usize * 4;
                reader
                    .next_frame(&mut mem[i * step_size..(i + 1) * step_size])
                    .with_context(|| format!("decoding {}", path.anonymize().display()))?;
            }
            let (width, height) = dims.unwrap();
            let mem = mem.unwrap();
            unsafe {
                let image = DedicatedImage::new(
                    &handle.gfx.device,
                    &handle.gfx.memory_properties,
                    &vk::ImageCreateInfo::default()
                        .image_type(vk::ImageType::TYPE_2D)
                        .format(vk::Format::R8G8B8A8_SRGB)
                        .extent(vk::Extent3D {
                            width,
                            height,
                            depth: 1,
                        })
                        .mip_levels(1)
                        .array_layers(self.size as u32)
                        .samples(vk::SampleCountFlags::TYPE_1)
                        .usage(vk::ImageUsageFlags::SAMPLED | vk::ImageUsageFlags::TRANSFER_DST),
                );

                let range = vk::ImageSubresourceRange {
                    aspect_mask: vk::ImageAspectFlags::COLOR,
                    base_mip_level: 0,
                    level_count: 1,
                    base_array_layer: 0,
                    layer_count: self.size as u32,
                };
                let src = handle.staging.buffer();
                let buffer_offset = mem.offset();
                let dst = image.handle;

                handle
                    .transfer
                    .run(move |xf, cmd| {
                        xf.device.cmd_pipeline_barrier(
                            cmd,
                            vk::PipelineStageFlags::TOP_OF_PIPE,
                            vk::PipelineStageFlags::TRANSFER,
                            vk::DependencyFlags::default(),
                            &[],
                            &[],
                            &[vk::ImageMemoryBarrier::default()
                                .dst_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                                .src_queue_family_index(vk::QUEUE_FAMILY_IGNORED)
                                .dst_queue_family_index(vk::QUEUE_FAMILY_IGNORED)
                                .old_layout(vk::ImageLayout::UNDEFINED)
                                .new_layout(vk::ImageLayout::TRANSFER_DST_OPTIMAL)
                                .image(dst)
                                .subresource_range(range)],
                        );
                        xf.device.cmd_copy_buffer_to_image(
                            cmd,
                            src,
                            dst,
                            vk::ImageLayout::TRANSFER_DST_OPTIMAL,
                            &[vk::BufferImageCopy {
                                buffer_offset,
                                image_subresource: vk::ImageSubresourceLayers {
                                    aspect_mask: vk::ImageAspectFlags::COLOR,
                                    mip_level: 0,
                                    base_array_layer: 0,
                                    layer_count: range.layer_count,
                                },
                                image_extent: vk::Extent3D {
                                    width,
                                    height,
                                    depth: 1,
                                },
                                ..Default::default()
                            }],
                        );
                        xf.stages |= vk::PipelineStageFlags::FRAGMENT_SHADER;
                        xf.image_barriers.push(
                            vk::ImageMemoryBarrier::default()
                                .src_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                                .dst_access_mask(vk::AccessFlags::SHADER_READ)
                                .src_queue_family_index(xf.queue_family)
                                .dst_queue_family_index(xf.dst_queue_family)
                                .old_layout(vk::ImageLayout::TRANSFER_DST_OPTIMAL)
                                .new_layout(vk::ImageLayout::SHADER_READ_ONLY_OPTIMAL)
                                .image(dst)
                                .subresource_range(range),
                        );
                    })
                    .await?;

                trace!(
                    width = width,
                    height = height,
                    path = %full_path.anonymize().display(),
                    "loaded array"
                );
                Ok(image)
            }
        })
    }
}

impl skid_steer::Source for PngArray {
    type Output = DedicatedImage; // TODO: We may want a dedicated struct here.

    async fn load(self, context: &skid_steer::Context<'_>) -> Option<DedicatedImage> {
        self.load_inner(context)
            .await
            .inspect_err(|e| tracing::error!("{}", e))
            .ok()
    }
}
