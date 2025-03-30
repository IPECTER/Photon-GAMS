#if !defined INCLUDE_SKY_CLOUDS_TOWERING_CUMULUS
#define INCLUDE_SKY_CLOUDS_TOWERING_CUMULUS

// Volumetric towering cumulus clouds - more vertically developed than regular cumulus

#include "common.glsl"

// altitude_fraction := 0 at the bottom of the cloud layer and 1 at the top
float clouds_towering_cumulus_altitude_shaping(float density, float altitude_fraction) {
	// Carve egg shape to make the cloud more vertical
	density -= smoothstep(0.0, 8.0, altitude_fraction) * 0.6;

	// Reduce density at the top and bottom of the cloud
	density *= smoothstep(0.0, 0.2, altitude_fraction);

	return density;
}

float clouds_towering_cumulus_density(vec3 pos) {
	const float wind_angle = CLOUDS_TOWERING_CUMULUS_WIND_ANGLE * degree;
	const vec2 wind_velocity = CLOUDS_TOWERING_CUMULUS_WIND_SPEED * vec2(cos(wind_angle), sin(wind_angle));

	float r = length(pos);
	if (r < clouds_towering_cumulus_radius || r > clouds_towering_cumulus_top_radius) return 0.0;

	float altitude_fraction = (r - clouds_towering_cumulus_radius) * clouds_params.towering_cumulus_altitude_scale;

	pos.xz += cameraPosition.xz * CLOUDS_SCALE + wind_velocity * world_age;

	// 2D noise for base shape and coverage
	vec3 noise = vec3(
		texture(noisetex, (0.0000001 / CLOUDS_TOWERING_CUMULUS_SIZE) * pos.xz).x, // cloud coverage
		texture(noisetex, (0.000002 / CLOUDS_TOWERING_CUMULUS_SIZE) * pos.xz).x, // cloud coverage
		texture(noisetex, (0.000008 / CLOUDS_TOWERING_CUMULUS_SIZE) * pos.xz).w  // cloud shape
	);

	float density = mix(clouds_towering_cumulus_coverage.x, clouds_towering_cumulus_coverage.y, pow(noise.x, 1.0));
		  density = linear_step(1.0 - density, 1.0, noise.y);
	      density = clouds_towering_cumulus_altitude_shaping(density, altitude_fraction);

	if (density < eps) return 0.0;

#ifndef PROGRAM_PREPARE
	// Curl noise used to warp the 3D noise into swirling shapes
	vec3 curl = (0.01 * CLOUDS_TOWERING_CUMULUS_CURL_STRENGTH) * texture(colortex7, 0.001 * pos).xyz * smoothstep(0.4, 1.0, 1.0 - altitude_fraction);
	vec3 wind = vec3(wind_velocity * world_age, 0.0).xzy;

	// 3D worley noise for detail
	float worley_0 = texture(colortex6, (pos + 0.2 * wind) * 0.00015 + curl * 1.0).x;
	float worley_1 = texture(colortex6, (pos + 0.4 * wind) * 0.00095 + curl * 3.0).x;
#else
	const float worley_0 = 0.5;
	const float worley_1 = 0.5;
#endif // !PROGRAM_PREPARE

	float detail_fade = 0.02 * smoothstep(0.85, 1.0, 1.0 - altitude_fraction)
	                  - 0.35 * smoothstep(0.05, 0.5, altitude_fraction) + 0.6;

	density -= clouds_params.towering_cumulus_detail_weights.x * sqr(worley_0) * dampen(clamp01(1.0 - density));
	density -= clouds_params.towering_cumulus_detail_weights.y * sqr(worley_1) * dampen(clamp01(1.0 - density)) * detail_fade;

	// Adjust density so that the clouds are wispy at the bottom and hard at the top
	density  = max0(density);
	density  = 1.0 - pow(
		1.0 - density, mix(
		clouds_params.towering_cumulus_edge_sharpening.x,
		clouds_params.towering_cumulus_edge_sharpening.y, 
		altitude_fraction
		)
	);
	density *= CLOUDS_ROUGHNESS + 0.9 * smoothstep(0.2, 0.7, altitude_fraction);

	return density;
}

float clouds_towering_cumulus_optical_depth(
	vec3 ray_origin,
	vec3 ray_dir,
	float dither,
	const uint step_count,
	bool include_towering
) {
	const float step_growth = 4.0;

	float step_length = 0.1 * clouds_towering_cumulus_thickness / float(step_count);

	vec3 ray_pos = ray_origin;
	vec4 ray_step = vec4(ray_dir, 1.0) * step_length;

	float optical_depth = 0.0;

	for (uint i = 0u; i < step_count; ++i, ray_pos += ray_step.xyz) {
		ray_step *= step_growth;
		optical_depth += clouds_towering_cumulus_density(ray_pos + ray_step.xyz * dither) * ray_step.w;
	}

	return optical_depth;
}

vec2 clouds_towering_cumulus_scattering(
	float density,
	float light_optical_depth,
	float sky_optical_depth,
	float ground_optical_depth,
	float step_transmittance,
	float cos_theta,
	float bounced_light
) {
	vec2 scattering = vec2(0.0);

	float scatter_amount = clouds_params.l0_scattering_coeff;
	float extinct_amount = clouds_params.l0_extinction_coeff;

	float scattering_integral_times_density = (1.0 - step_transmittance) / clouds_params.l0_extinction_coeff;

	float powder_effect = clouds_powder_effect(density, cos_theta);

	float phase = clouds_phase_single(cos_theta);
	vec3 phase_g = pow(vec3(0.6, 0.9, 0.3), vec3(1.0 + light_optical_depth));

	for (uint i = 0u; i < 8u; ++i) {
		scattering.x += scatter_amount * exp(-extinct_amount *  light_optical_depth) * phase * (1.0 - 0.8 * clouds_params.l0_shadow);
		scattering.x += scatter_amount * exp(-extinct_amount * ground_optical_depth) * isotropic_phase * bounced_light;
		scattering.x += scatter_amount * exp(-extinct_amount *    sky_optical_depth) * isotropic_phase * clouds_params.l0_shadow * 0.2;
		scattering.y += scatter_amount * exp(-extinct_amount *    sky_optical_depth) * isotropic_phase;

		scatter_amount *= 0.55 * mix(lift(clamp01(clouds_params.l0_scattering_coeff / 0.1), 0.33), 1.0, cos_theta * 0.5 + 0.5) * powder_effect;
		extinct_amount *= 0.4;
		phase_g *= 0.8;

		powder_effect = mix(powder_effect, sqrt(powder_effect), 0.5);

		phase = clouds_phase_multi(cos_theta, phase_g);
	}

	return scattering * scattering_integral_times_density;
}

CloudsResult draw_towering_cumulus_clouds(
	vec3 air_viewer_pos,
	vec3 ray_dir,
	vec3 clear_sky,
	float distance_to_terrain,
	float dither
) {
	// ---------------------
	//   Raymarching Setup
	// ---------------------

#if defined PROGRAM_DEFERRED0
	const uint  primary_steps_horizon = CLOUDS_TOWERING_CUMULUS_PRIMARY_STEPS_H / 2;
	const uint  primary_steps_zenith  = CLOUDS_TOWERING_CUMULUS_PRIMARY_STEPS_Z / 2;
#else
	const uint  primary_steps_horizon = CLOUDS_TOWERING_CUMULUS_PRIMARY_STEPS_H;
	const uint  primary_steps_zenith  = CLOUDS_TOWERING_CUMULUS_PRIMARY_STEPS_Z;
#endif // PROGRAM_DEFERRED0
	const uint  lighting_steps        = CLOUDS_TOWERING_CUMULUS_LIGHTING_STEPS;
	const uint  ambient_steps         = CLOUDS_TOWERING_CUMULUS_AMBIENT_STEPS;
	const float max_ray_length        = 100e4;
	const float min_transmittance     = 0.075;
	const float planet_albedo         = 0.4;
	const vec3  sky_dir               = vec3(0.0, 1.0, 0.0);

	uint primary_steps = uint(mix(primary_steps_horizon, primary_steps_zenith, abs(ray_dir.y)));

	float r = length(air_viewer_pos);

	vec2 dists = intersect_spherical_shell(air_viewer_pos, ray_dir, clouds_towering_cumulus_radius, clouds_towering_cumulus_top_radius);
	bool planet_intersected = intersect_sphere(air_viewer_pos, ray_dir, min(r - 10.0, planet_radius)).y >= 0.0;
	bool terrain_intersected = distance_to_terrain >= 0.0 && r < clouds_towering_cumulus_radius && distance_to_terrain * CLOUDS_SCALE < dists.y;

	if (dists.y < 0.0                                   // volume not intersected
	 || planet_intersected && r < clouds_towering_cumulus_radius // planet blocking clouds
	 || terrain_intersected                             // terrain blocking clouds
	) { return clouds_not_hit; }

	float ray_length = (distance_to_terrain >= 0.0) ? distance_to_terrain : dists.y;
	      ray_length = clamp(ray_length - dists.x, 0.0, max_ray_length);
	float step_length = ray_length * rcp(float(primary_steps));

	vec3 ray_step = ray_dir * step_length;
	vec3 ray_origin = air_viewer_pos + ray_dir * (dists.x + step_length * dither);

	vec2 scattering = vec2(0.0); // x: direct light, y: skylight
	float transmittance = 1.0;

	float distance_sum = 0.0;
	float distance_weight_sum = 0.0;

	// ------------------
	//   Lighting Setup
	// ------------------

	float altocumulus_shadow = linear_step(0.5, 0.6, clouds_params.l1_coverage.x) * dampen(day_factor);

	bool  moonlit            = sun_dir.y < -0.04;
	vec3  light_dir          = moonlit ? moon_dir : sun_dir;
	float cos_theta          = dot(ray_dir, light_dir);
	float bounced_light      = planet_albedo * light_dir.y * rcp_pi;

	float extinction_coeff   = mix(0.07, 0.1, smoothstep(0.0, 0.3, abs(sun_dir.y))) * (1.0 - 0.33 * rainStrength) * (1.0 - 0.6 * altocumulus_shadow) * CLOUDS_TOWERING_CUMULUS_DENSITY;
	float scattering_coeff   = extinction_coeff * mix(1.00, 1.6, rainStrength);

	float dynamic_thickness  = mix(0.5, 1.0, smoothstep(0.4, 0.6, clouds_towering_cumulus_coverage.y));
	vec2  detail_weights     = vec2(0.33, 0.40) * CLOUDS_TOWERING_CUMULUS_DETAIL_STRENGTH;
	vec2  edge_sharpening    = vec2(3.0, 8.0);

	// --------------------
	//   Raymarching Loop
	// --------------------

	for (uint i = 0u; i < primary_steps; ++i) {
		if (transmittance < min_transmittance) break;

		vec3 ray_pos = ray_origin + ray_step * i;
		float r_sample = length(ray_pos);

		float density = clouds_towering_cumulus_density(ray_pos);

		if (density < eps) continue;

		float distance_to_sample = distance(ray_origin, ray_pos);
		float distance_fade = smoothstep(0.95, 1.0, distance_to_sample * rcp(max_ray_length));

		density *= 1.0 - distance_fade;

		float step_optical_depth = density * extinction_coeff * step_length;
		float step_transmittance = exp(-step_optical_depth);

#if defined PROGRAM_DEFERRED0
		vec2 hash = vec2(0.0);
#else
		vec2 hash = hash2(fract(ray_pos));
#endif // PROGRAM_DEFERRED0

		float light_optical_depth  = clouds_towering_cumulus_optical_depth(ray_pos, light_dir, hash.x, lighting_steps, false);
		float sky_optical_depth    = clouds_towering_cumulus_optical_depth(ray_pos, sky_dir, hash.y, ambient_steps, false);

		float ground_optical_depth = mix(
			density, 
			1.0, 
			clamp01(r_sample/clouds_towering_cumulus_top_radius * 2.0 - 1.0)
		) * (r_sample - clouds_towering_cumulus_radius) / clouds_towering_cumulus_thickness;

		scattering += clouds_towering_cumulus_scattering(
			density,
			light_optical_depth,
			sky_optical_depth,
			ground_optical_depth,
			step_transmittance,
			cos_theta,
			bounced_light
		) * transmittance;

		transmittance *= step_transmittance;

		distance_sum += distance_to_sample * density;
		distance_weight_sum += density;
	}

	vec3 light_color  = sunlight_color * atmosphere_transmittance(ray_origin, light_dir);
		 light_color  = atmosphere_post_processing(light_color);
	     light_color *= moonlit ? moon_color : sun_color;

	float clouds_transmittance = linear_step(min_transmittance, 1.0, transmittance);

	vec3 clouds_scattering = scattering.x * light_color + scattering.y * sky_color;
	if (distance_to_terrain < 0.0) clouds_scattering = clouds_aerial_perspective(clouds_scattering, clouds_transmittance, air_viewer_pos, ray_origin, ray_dir, clear_sky);
	
	float apparent_distance = (distance_weight_sum == 0.0)
		? 1e6
		: (distance_sum / distance_weight_sum) + distance(air_viewer_pos, ray_origin);

	return CloudsResult(
		vec4(clouds_scattering, scattering.y),
		clouds_transmittance,
		apparent_distance
	);
}
#endif // INCLUDE_SKY_CLOUDS_TOWERING_CUMULUS