/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/c2_dof
  Calculate depth of field

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

layout (location = 0) out vec3 scene_color;

/* RENDERTARGETS: 0 */

in vec2 uv;

// ------------
//   Uniforms
// ------------

uniform sampler2D noisetex;

uniform sampler2D colortex0;

uniform sampler2D depthtex0;

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferProjectionInverse;

uniform float near, far;

uniform float aspectRatio;
uniform float centerDepthSmooth;

uniform int frameCounter;

uniform vec2 view_pixel_size;
uniform vec2 taa_offset;

#include "/include/misc/distant_horizons.glsl"
#include "/include/utility/random.glsl"
#include "/include/utility/sampling.glsl"
#include "/include/utility/fast_math.glsl"
#include "/include/utility/space_conversion.glsl"

void main() {
	ivec2 texel = ivec2(gl_FragCoord.xy);

	float depth = texelFetch(depthtex0, texel, 0).x;

#ifdef DISTANT_HORIZONS
	float depth_dh = texelFetch(dhDepthTex, texel, 0).x;

	if (is_distant_horizons_terrain(depth, depth_dh)) {
		depth = view_to_screen_space_depth(
			gbufferProjection,
			screen_to_view_space_depth(dhProjectionInverse, depth_dh)
		);
	}
#endif

	if (depth < hand_depth) {
		scene_color = texelFetch(colortex0, texel, 0).rgb;
		return;
	};

	// Calculate vogel disk rotation
	float theta  = texelFetch(noisetex, texel & 511, 0).b;
	      theta  = r1(frameCounter, theta);
	      theta *= tau;

	// Calculate circle of confusion
	float focus = DOF_FOCUS < 0.0 ? centerDepthSmooth : view_to_screen_space_depth(gbufferProjection, DOF_FOCUS);
	vec2 CoC = min(abs(depth - focus), 0.1) * (DOF_INTENSITY * 0.2 / 1.37) * vec2(DOF_SIZE_X, aspectRatio * DOF_SIZE_Y) * gbufferProjection[1][1];

	scene_color = vec3(0.0);

#ifdef DOF_CA                                                                                     // TODO
    float caDist = abs(depth - view_to_screen_space_depth(gbufferProjection, 1.8)) * 2.75 - 0.1; // Magic numbers ;-;
    caDist = max0(1.0 - caDist - (1.0 - rcp(16.0))) * 16.0;                                     //
    caDist = clamp01(smoothstep(0.135, 0.147, caDist) * 5.0);                                  // More magic numbers ._.
    caDist *= clamp01(smoothstep(0.0, 0.1, (CoC.x + CoC.y) * 10.0));                          // Focused
    caDist *= (1.0 - step(1.0 - eps, depth)) * 0.75 + 0.25;                                  // Sky
    caDist *= step(focus, depth) * 0.5 + 0.5;                                               // Closer than focus point
    caDist *= DOF_CA_INTENSITY * 0.01;                                                     // Adjust strength

	const mat3x2 ca = mat3x2(
		polar_to_cartesian2(DOF_CA_R_OFFSET, DOF_CA_R_ANGLE * degree),
		polar_to_cartesian2(DOF_CA_G_OFFSET, DOF_CA_G_ANGLE * degree),
		polar_to_cartesian2(DOF_CA_B_OFFSET, DOF_CA_B_ANGLE * degree)
	);
#endif

	vec2 m = vec2(1.0 - 2.0 * view_pixel_size * rcp(taau_render_scale));
	for (int i = 0; i < DOF_SAMPLES; ++i) {
		vec2 offset = vogel_disk_sample(i, DOF_SAMPLES, theta) * CoC * ((1.0 - step(1.0 - eps, depth)) * 0.6 + 0.4);
		offset = vec2(uv + offset);
#ifdef DOF_CA
		scene_color.r += textureLod(colortex0, clamp(offset + caDist * ca[0], vec2(0.0), m) * taau_render_scale, 0).r;
		scene_color.g += textureLod(colortex0, clamp(offset + caDist * ca[1], vec2(0.0), m) * taau_render_scale, 0).g;
		scene_color.b += textureLod(colortex0, clamp(offset + caDist * ca[2], vec2(0.0), m) * taau_render_scale, 0).b;
#else
		scene_color += textureLod(colortex0, clamp(offset, vec2(0.0), m) * taau_render_scale, 0).rgb;
#endif
	}
	//scene_color = vec3(caDist);
	//scene_color = vec3(smoothstep(0.0, 0.2, (CoC.x + CoC.y) * 10.0));
	scene_color *= rcp(DOF_SAMPLES);
}

