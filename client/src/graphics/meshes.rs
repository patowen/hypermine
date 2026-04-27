use std::mem;

use anyhow::{Ok, Result, anyhow};
use ash::{Device, vk};
use lahar::{BufferRegionAlloc, DedicatedImage};
use memoffset::offset_of;
use vk_shader_macros::include_glsl;

use crate::loader::LoadCtx;

use super::Base;
use common::{defer, math};

const VERT: &[u32] = include_glsl!("shaders/mesh.vert");
const FRAG: &[u32] = include_glsl!("shaders/mesh.frag");

pub struct Meshes {
    pipeline_layout: vk::PipelineLayout,
    pipeline: vk::Pipeline,
}

impl Meshes {
    pub fn new(gfx: &Base, ds_layout: vk::DescriptorSetLayout) -> Self {
        let device = &*gfx.device;
        unsafe {
            // Construct the shader modules
            let vert = device
                .create_shader_module(&vk::ShaderModuleCreateInfo::default().code(VERT), None)
                .unwrap();
            // Note that these only need to live until the pipeline itself is constructed
            let v_guard = defer(|| device.destroy_shader_module(vert, None));

            let frag = device
                .create_shader_module(&vk::ShaderModuleCreateInfo::default().code(FRAG), None)
                .unwrap();
            let f_guard = defer(|| device.destroy_shader_module(frag, None));

            // Define the outward-facing interface of the shaders, incl. uniforms, samplers, etc.
            let pipeline_layout = device
                .create_pipeline_layout(
                    &vk::PipelineLayoutCreateInfo::default()
                        .set_layouts(&[gfx.common_layout, ds_layout])
                        .push_constant_ranges(&[vk::PushConstantRange {
                            stage_flags: vk::ShaderStageFlags::VERTEX,
                            offset: 0,
                            size: 64,
                        }]),
                    None,
                )
                .unwrap();

            let entry_point = cstr!("main").as_ptr();
            let mut pipelines = device
                .create_graphics_pipelines(
                    gfx.pipeline_cache,
                    &[vk::GraphicsPipelineCreateInfo::default()
                        .stages(&[
                            vk::PipelineShaderStageCreateInfo {
                                stage: vk::ShaderStageFlags::VERTEX,
                                module: vert,
                                p_name: entry_point,
                                ..Default::default()
                            },
                            vk::PipelineShaderStageCreateInfo {
                                stage: vk::ShaderStageFlags::FRAGMENT,
                                module: frag,
                                p_name: entry_point,
                                ..Default::default()
                            },
                        ])
                        .vertex_input_state(
                            &vk::PipelineVertexInputStateCreateInfo::default()
                                .vertex_binding_descriptions(&[vk::VertexInputBindingDescription {
                                    binding: 0,
                                    stride: mem::size_of::<Vertex>() as u32,
                                    input_rate: vk::VertexInputRate::VERTEX,
                                }])
                                .vertex_attribute_descriptions(&[
                                    vk::VertexInputAttributeDescription {
                                        location: 0,
                                        binding: 0,
                                        format: vk::Format::R32G32B32_SFLOAT,
                                        offset: offset_of!(Vertex, position) as u32,
                                    },
                                    vk::VertexInputAttributeDescription {
                                        location: 1,
                                        binding: 0,
                                        format: vk::Format::R32G32B32_SFLOAT,
                                        offset: offset_of!(Vertex, texcoords) as u32,
                                    },
                                    vk::VertexInputAttributeDescription {
                                        location: 2,
                                        binding: 0,
                                        format: vk::Format::R32G32B32_SFLOAT,
                                        offset: offset_of!(Vertex, normal) as u32,
                                    },
                                ]),
                        )
                        .input_assembly_state(
                            &vk::PipelineInputAssemblyStateCreateInfo::default()
                                .topology(vk::PrimitiveTopology::TRIANGLE_LIST),
                        )
                        .viewport_state(
                            &vk::PipelineViewportStateCreateInfo::default()
                                .scissor_count(1)
                                .viewport_count(1),
                        )
                        .rasterization_state(
                            &vk::PipelineRasterizationStateCreateInfo::default()
                                .cull_mode(vk::CullModeFlags::BACK)
                                .front_face(vk::FrontFace::COUNTER_CLOCKWISE)
                                .polygon_mode(vk::PolygonMode::FILL)
                                .line_width(1.0),
                        )
                        .multisample_state(
                            &vk::PipelineMultisampleStateCreateInfo::default()
                                .rasterization_samples(vk::SampleCountFlags::TYPE_1),
                        )
                        .depth_stencil_state(
                            &vk::PipelineDepthStencilStateCreateInfo::default()
                                .depth_test_enable(true)
                                .depth_write_enable(true)
                                .depth_compare_op(vk::CompareOp::GREATER),
                        )
                        .color_blend_state(
                            &vk::PipelineColorBlendStateCreateInfo::default().attachments(&[
                                vk::PipelineColorBlendAttachmentState {
                                    blend_enable: vk::TRUE,
                                    src_color_blend_factor: vk::BlendFactor::ONE,
                                    dst_color_blend_factor: vk::BlendFactor::ZERO,
                                    color_blend_op: vk::BlendOp::ADD,
                                    color_write_mask: vk::ColorComponentFlags::R
                                        | vk::ColorComponentFlags::G
                                        | vk::ColorComponentFlags::B,
                                    ..Default::default()
                                },
                            ]),
                        )
                        .dynamic_state(
                            &vk::PipelineDynamicStateCreateInfo::default().dynamic_states(&[
                                vk::DynamicState::VIEWPORT,
                                vk::DynamicState::SCISSOR,
                            ]),
                        )
                        .layout(pipeline_layout)
                        .render_pass(gfx.render_pass)
                        .subpass(0)],
                    None,
                )
                .unwrap()
                .into_iter();

            let pipeline = pipelines.next().unwrap();
            gfx.set_name(pipeline, cstr!("meshes"));

            // Clean up the shaders explicitly, so the defer guards don't hold onto references we're
            // moving into `Self` to be returned
            v_guard.invoke();
            f_guard.invoke();

            Self {
                pipeline_layout,
                pipeline,
            }
        }
    }

    pub unsafe fn draw(
        &mut self,
        device: &Device,
        common_ds: vk::DescriptorSet,
        cmd: vk::CommandBuffer,
        mesh: &Mesh,
        transform: &na::Matrix4<f32>,
    ) {
        unsafe {
            device.cmd_bind_pipeline(cmd, vk::PipelineBindPoint::GRAPHICS, self.pipeline);
            device.cmd_bind_descriptor_sets(
                cmd,
                vk::PipelineBindPoint::GRAPHICS,
                self.pipeline_layout,
                0,
                &[common_ds, mesh.ds],
                &[],
            );
            device.cmd_push_constants(
                cmd,
                self.pipeline_layout,
                vk::ShaderStageFlags::VERTEX,
                0,
                &mem::transmute::<na::Matrix4<f32>, [u8; 64]>(*transform),
            );
            device.cmd_bind_vertex_buffers(
                cmd,
                0,
                &[mesh.geom.vertices.buffer],
                &[mesh.geom.vertices.offset],
            );
            device.cmd_bind_index_buffer(
                cmd,
                mesh.geom.indices.buffer,
                mesh.geom.indices.offset,
                vk::IndexType::UINT32,
            );
            device.cmd_draw_indexed(cmd, mesh.geom.index_count, 1, 0, 0, 0);
        }
    }

    pub unsafe fn destroy(&mut self, device: &Device) {
        unsafe {
            device.destroy_pipeline(self.pipeline, None);
            device.destroy_pipeline_layout(self.pipeline_layout, None);
        }
    }
}

#[repr(C)]
pub struct Vertex {
    pub position: math::MPoint<f32>,
    pub texcoords: na::Vector3<f32>,
    pub normal: math::MDirection<f32>,
}

#[derive(Copy, Clone)]
pub struct Geometry {
    vertices: BufferRegionAlloc,
    indices: BufferRegionAlloc,
    index_count: u32,
}

impl Geometry {
    pub async unsafe fn new(
        ctx: &LoadCtx,
        vertices: impl ExactSizeIterator<Item = Vertex>,
        indices: impl ExactSizeIterator<Item = u32>,
    ) -> Result<Geometry> {
        let vertex_byte_size = vertices.len() * mem::size_of::<Vertex>();
        let mut v_staging = ctx
            .staging
            .alloc(vertex_byte_size)
            .await
            .ok_or_else(|| anyhow!("too large"))?;
        for (vertex, storage) in vertices.zip(v_staging.chunks_exact_mut(mem::size_of::<Vertex>()))
        {
            unsafe {
                std::ptr::write_unaligned(storage.as_ptr() as *mut Vertex, vertex);
            }
        }

        let num_indices = indices.len();
        let index_byte_size = num_indices * mem::size_of::<u32>();
        let mut i_staging = ctx
            .staging
            .alloc(indices.len() * mem::size_of::<u32>())
            .await
            .ok_or_else(|| anyhow!("too large"))?;
        for (idx, storage) in indices.zip(i_staging.chunks_exact_mut(mem::size_of::<u32>())) {
            storage.copy_from_slice(&idx.to_ne_bytes());
        }

        let vert_alloc = ctx.vertex_alloc.lock().unwrap().alloc(
            &ctx.gfx.device,
            vertex_byte_size as vk::DeviceSize,
            4,
        );
        let staging_buffer = ctx.staging.buffer();
        let vert_buffer = vert_alloc.buffer;
        let vert_src_offset = v_staging.offset();
        let vert_dst_offset = vert_alloc.offset;
        let vertex_upload = unsafe {
            ctx.transfer.run(move |xf, cmd| {
                xf.device.cmd_copy_buffer(
                    cmd,
                    staging_buffer,
                    vert_buffer,
                    &[vk::BufferCopy {
                        src_offset: vert_src_offset,
                        dst_offset: vert_dst_offset,
                        size: vertex_byte_size as vk::DeviceSize,
                    }],
                );
                xf.stages |= vk::PipelineStageFlags::VERTEX_INPUT;
                xf.buffer_barriers.push(
                    vk::BufferMemoryBarrier::default()
                        .src_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                        .dst_access_mask(vk::AccessFlags::VERTEX_ATTRIBUTE_READ)
                        .src_queue_family_index(xf.queue_family)
                        .dst_queue_family_index(xf.dst_queue_family)
                        .buffer(vert_buffer)
                        .offset(vert_dst_offset)
                        .size(vertex_byte_size as vk::DeviceSize),
                );
            })
        };

        let idx_alloc = ctx.index_alloc.lock().unwrap().alloc(
            &ctx.gfx.device,
            index_byte_size as vk::DeviceSize,
            4,
        );
        let idx_buffer = idx_alloc.buffer;
        let idx_src_offset = i_staging.offset();
        let idx_dst_offset = idx_alloc.offset;
        let index_upload = unsafe {
            ctx.transfer.run(move |xf, cmd| {
                xf.device.cmd_copy_buffer(
                    cmd,
                    staging_buffer,
                    idx_buffer,
                    &[vk::BufferCopy {
                        src_offset: idx_src_offset,
                        dst_offset: idx_dst_offset,
                        size: index_byte_size as vk::DeviceSize,
                    }],
                );
                xf.stages |= vk::PipelineStageFlags::VERTEX_INPUT;
                xf.buffer_barriers.push(
                    vk::BufferMemoryBarrier::default()
                        .src_access_mask(vk::AccessFlags::TRANSFER_WRITE)
                        .dst_access_mask(vk::AccessFlags::INDEX_READ)
                        .src_queue_family_index(xf.queue_family)
                        .dst_queue_family_index(xf.dst_queue_family)
                        .buffer(idx_buffer)
                        .offset(idx_dst_offset)
                        .size(index_byte_size as vk::DeviceSize),
                );
            })
        };

        // Upload concurrently
        let (r1, r2) = tokio::join!(vertex_upload, index_upload);
        r1?;
        r2?;

        Ok(Geometry {
            vertices: vert_alloc,
            indices: idx_alloc,
            index_count: num_indices as u32,
        })
    }
}

#[derive(Copy, Clone)]
pub struct Mesh {
    pub geom: Geometry,
    pub pool: vk::DescriptorPool,
    pub ds: vk::DescriptorSet,
    // TODO: Make shareable
    pub color: DedicatedImage,
    pub color_view: vk::ImageView,
}

impl crate::loader::Cleanup for Mesh {
    unsafe fn cleanup(mut self, gfx: &Base) {
        unsafe {
            let device = &*gfx.device;
            device.destroy_descriptor_pool(self.pool, None);
            device.destroy_image_view(self.color_view, None);
            self.color.destroy(device);
        }
    }
}
