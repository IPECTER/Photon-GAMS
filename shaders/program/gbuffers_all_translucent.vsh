/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/gbuffers_all_translucent:
  Handle translucent terrain, translucent entities (Iris), translucent handheld
  items and gbuffers_textured

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

out vec2 uv;
out vec2 light_levels;
out vec3 position_view;
out vec3 position_scene;
out vec4 tint;

flat out vec3 light_color;
flat out vec3 ambient_color;
flat out uint material_mask;
flat out mat3 tbn;

#if defined PROGRAM_GBUFFERS_WATER
out vec2 atlas_tile_coord;
out vec3 position_tangent;
flat out vec2 atlas_tile_offset;
flat out vec2 atlas_tile_scale;
#endif

#if defined WORLD_OVERWORLD
#include "/include/fog/overworld/parameters.glsl"
flat out OverworldFogParameters fog_params;
#endif

// --------------
//   Attributes
// --------------

attribute vec4 at_tangent;
attribute vec3 mc_Entity;
attribute vec2 mc_midTexCoord;

// ------------
//   Uniforms
// ------------

uniform sampler2D noisetex;

uniform sampler2D colortex4; // Sky map, lighting colors
uniform sampler2D colortex9; // Sky SH

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferProjectionInverse;

uniform vec3 cameraPosition;

uniform float near;
uniform float far;

uniform ivec2 atlasSize;
uniform int renderStage;

uniform int worldTime;
uniform int worldDay;
uniform int frameCounter;
uniform float frameTimeCounter;

uniform float sunAngle;
uniform float rainStrength;
uniform float wetness;

uniform vec2 view_res;
uniform vec2 view_pixel_size;
uniform vec2 taa_offset;

uniform vec3 light_dir;
uniform vec3 sun_dir;

uniform float eye_skylight;

uniform float biome_temperate;
uniform float biome_arid;
uniform float biome_snowy;
uniform float biome_taiga;
uniform float biome_jungle;
uniform float biome_swamp;
uniform float biome_may_rain;
uniform float biome_may_snow;
uniform float biome_temperature;
uniform float biome_humidity;

uniform float world_age;
uniform float time_sunrise;
uniform float time_noon;
uniform float time_sunset;
uniform float time_midnight;

uniform float desert_sandstorm;

#if defined PROGRAM_GBUFFERS_ENTITIES_TRANSLUCENT
uniform int entityId;
#endif

#if defined PROGRAM_GBUFFERS_BLOCK_TRANSLUCENT
uniform int blockEntityId;
#endif

#if (defined PROGRAM_GBUFFERS_ENTITIES_TRANSLUCENT || defined PROGRAM_GBUFFERS_HAND_WATER) && defined IS_IRIS
uniform int currentRenderedItemId;
#endif

#if defined (PHYSICS_MOD_OCEAN) && defined (PHYSICS_OCEAN)
#include "/include/misc/oceans.glsl"
#endif

#include "/include/utility/space_conversion.glsl"
#include "/include/vertex/displacement.glsl"
#include "/include/vertex/utility.glsl"

#if defined WORLD_OVERWORLD
#include "/include/misc/weather.glsl"
#endif

void main() {
	uv            = mat2(gl_TextureMatrix[0]) * gl_MultiTexCoord0.xy + gl_TextureMatrix[0][3].xy;
	light_levels  = clamp01(gl_MultiTexCoord1.xy * rcp(240.0));
	tint          = gl_Color;
	material_mask = get_material_mask();
	tbn           = get_tbn_matrix();

	int lighting_color_x = SKY_MAP_LIGHT_X;
	light_color   = texelFetch(colortex4, ivec2(lighting_color_x, 0), 0).rgb;
#if defined WORLD_OVERWORLD && defined SH_SKYLIGHT
	ambient_color = texelFetch(colortex9, ivec2(9, 0), 0).rgb;
#else	
	ambient_color = texelFetch(colortex4, ivec2(lighting_color_x, 1), 0).rgb;
#endif

	bool is_top_vertex = uv.y < mc_midTexCoord.y;

	vec4 vert = gl_Vertex;

#if defined (PHYSICS_MOD_OCEAN) && defined (PHYSICS_OCEAN)
	if(physics_iterationsNormal >= 1.0 && material_mask == 1) {
		// basic texture to determine how shallow/far away from the shore the water is
		physics_localWaviness = texelFetch(physics_waviness, ivec2(gl_Vertex.xz) - physics_textureOffset, 0).r;
		// transform gl_Vertex (since it is the raw mesh, i.e. not transformed yet)
		vec4 finalPosition = vec4(gl_Vertex.x, gl_Vertex.y + physics_waveHeight(gl_Vertex.xz, PHYSICS_ITERATIONS_OFFSET, physics_localWaviness, physics_gameTime), gl_Vertex.z, gl_Vertex.w);
		// pass this to the fragment shader to fetch the texture there for per fragment normals
		physics_localPosition = finalPosition.xyz;

		// now use finalPosition instead of gl_Vertex
		vert.xyz = finalPosition.xyz;
	}
#endif

	position_scene = transform(gl_ModelViewMatrix, vert.xyz);                                 // To view space
	position_scene = view_to_scene_space(position_scene);                                          // To scene space
	position_scene = position_scene + cameraPosition;                                              // To world space
	position_scene = animate_vertex(position_scene, is_top_vertex, light_levels.y, material_mask); // Apply vertex animations
	position_scene = position_scene - cameraPosition;                                              // Back to scene space
	position_scene = world_curvature(position_scene);                                              // Apply world curvature

#if defined PROGRAM_GBUFFERS_WATER
	tint.a = 1.0;

	if (material_mask == 62) {
		// Nether portal
		position_tangent = (position_scene - gbufferModelViewInverse[3].xyz) * tbn;

		// (from fayer3)
		vec2 uv_minus_mid = uv - mc_midTexCoord;
		atlas_tile_offset = min(uv, mc_midTexCoord - uv_minus_mid);
		atlas_tile_scale = abs(uv_minus_mid) * 2.0;
		atlas_tile_coord = sign(uv_minus_mid) * 0.5 + 0.5;
	}
#endif

#if defined PROGRAM_GBUFFERS_TEXTURED
	// Make world border emissive
	if (renderStage == MC_RENDER_STAGE_WORLD_BORDER) material_mask = 4;
#endif

#if defined PROGRAM_GBUFFERS_TEXTURED && !defined IS_IRIS
	// Make enderman/nether portal particles glow
	if (gl_Color.r > gl_Color.g && gl_Color.g < 0.6 && gl_Color.b > 0.4) material_mask = 47;
#endif

#if defined PROGRAM_GBUFFERS_ENTITIES_TRANSLUCENT
 	if (material_mask == 102) material_mask = 80;
 #endif

#if defined PROGRAM_GBUFFERS_WATER
	// Fix issue where the normal of the bottom of the water surface is flipped
	if (dot(position_scene, tbn[2]) > 0.0) tbn[2] = -tbn[2];
#endif

#if defined WORLD_OVERWORLD
	fog_params = get_fog_parameters(get_weather());
#endif

	position_view = scene_to_view_space(position_scene);
	vec4 position_clip = project(gl_ProjectionMatrix, position_view);

#if   defined TAA && defined TAAU
	position_clip.xy  = position_clip.xy * taau_render_scale + position_clip.w * (taau_render_scale - 1.0);
	position_clip.xy += taa_offset * position_clip.w;
#elif defined TAA
	position_clip.xy += taa_offset * position_clip.w * 0.66;
#endif

	gl_Position = position_clip;
}