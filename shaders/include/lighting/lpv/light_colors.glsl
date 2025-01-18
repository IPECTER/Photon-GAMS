#ifndef INCLUDE_LIGHTING_LPV_LIGHT_COLORS
#define INCLUDE_LIGHTING_LPV_LIGHT_COLORS

#if (!defined(WORLD_END) && defined(END_COLORED_LIGHTING) || !defined(END_COLORED_LIGHTING))
const vec3[64] light_color = vec3[64](
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) * 12.0, // Strong white light
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) *  6.0, // Medium white light
	vec3(1.00, 1.00, 1.00) * vec3(1.0, 1.0, 1.0) * 1.0, // Weak white light
	vec3(1.00, 0.65, 0.25) * vec3(1.00, 0.65, 0.25) * 12.0, // Strong golden light
	vec3(MEDIUM_GOLDEN_LIGHT_R, MEDIUM_GOLDEN_LIGHT_G, MEDIUM_GOLDEN_LIGHT_B) * vec3(MEDIUM_GOLDEN_LIGHT_R, MEDIUM_GOLDEN_LIGHT_G, MEDIUM_GOLDEN_LIGHT_B) * MEDIUM_GOLDEN_LIGHT_I, // Medium golden light
	vec3(1.00, 0.57, 0.30) * vec3(1.00, 0.57, 0.30) *  2.2, // Weak golden light
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) *  4.0, // Redstone components
	vec3(LAVA_LIGHT_R, LAVA_LIGHT_G, LAVA_LIGHT_B) * vec3(LAVA_LIGHT_R, LAVA_LIGHT_G, LAVA_LIGHT_B) * LAVA_LIGHT_I, // Lava
	vec3(MEDIUM_ORANGE_LIGHT_R, MEDIUM_ORANGE_LIGHT_G, MEDIUM_ORANGE_LIGHT_B) * vec3(MEDIUM_ORANGE_LIGHT_R, MEDIUM_ORANGE_LIGHT_G, MEDIUM_ORANGE_LIGHT_B) * MEDIUM_ORANGE_LIGHT_I, // Medium orange light
	vec3(1.00, 0.63, 0.15) * vec3(1.00, 0.63, 0.15) * 4.0, // Brewing stand
	vec3(1.00, 0.57, 0.30) * vec3(1.00, 0.57, 0.30) * 12.0, // Medium golden light (Jack o' Lantern)
	vec3(0.45, 0.73, 1.00) * vec3(0.45, 0.73, 1.00) *  6.0, // Soul lights
	vec3(0.45, 0.73, 1.00) * vec3(0.45, 0.73, 1.00) * 14.0, // Beacon
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * 1.5, // End portal frame
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * SCULK_I, // Sculk
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * 5.0, // Pink glow
	vec3(0.75, 1.00, 0.50) * vec3(0.75, 1.00, 0.50) * 1.0, // Sea pickle
	vec3(1.00, 0.50, 0.25) * vec3(1.00, 0.50, 0.25) * 1.8, // Nether plants
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) * 1.0, // Restone wire
	vec3(0.97, 0.87, 0.57) * vec3(0.97, 0.87, 0.57) * 8.0, // Ochre froglight
	vec3(0.58, 0.77, 0.52) * vec3(0.58, 0.77, 0.52) * 8.0, // Verdant froglight
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 8.0, // Pearlescent froglight
	vec3(0.10, 0.86, 1.0) * vec3(0.10, 0.86, 1.0) * 2.0, // Enchanting table
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 4.0, // Amethyst cluster
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 4.0, // Calibrated sculk sensor
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * 6.0, // Active sculk sensor
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) * 3.5, // Redstone block
	vec3(1.00, 0.10, 1.00) * vec3(1.00, 0.10, 1.00) * 3.0, // Purple weak light
	vec3(0.10, 0.10, 1.00) *  3.3, // Lapis block
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) * 64.0, // Lightning rod
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * 12.0, // Nether portal
	vec3(0.0),  // End portal
	vec3(1.0, 0.1, 0.1) * vec3(1.0, 0.1, 0.1) *  1.5, // Red
	vec3(1.0, 0.5, 0.1) * vec3(1.0, 0.5, 0.1) *  1.5, // Orange
	vec3(1.0, 0.8, 0.2) * vec3(1.0, 0.8, 0.2) * 1.0, // Yellow
	vec3(0.7, 0.7, 0.0) * vec3(0.7, 0.7, 0.0) *  1.0, // Brown
	vec3(0.1, 1.0, 0.1) * vec3(0.1, 1.0, 0.1) * 1.0, // Green
	vec3(0.5, 1.0, 0.5) * vec3(0.5, 1.0, 0.5) * 1.0, // Lime
	vec3(0.1, 0.1, 1.0) * vec3(0.1, 0.1, 1.0) * 1.5, // Blue
	vec3(0.5, 0.5, 1.0) * vec3(0.5, 0.5, 1.0) * 1.0, // Light blue
	vec3(0.1, 1.0, 1.0) * vec3(0.1, 1.0, 1.0) * 1.0, // Cyan
	vec3(0.7, 0.1, 1.0) * vec3(0.7, 0.1, 1.0) * 1.0, // Purple
	vec3(1.0, 0.1, 1.0) * vec3(1.0, 0.1, 1.0) * 1.0, // Magenta
	vec3(1.0, 0.5, 1.0) * vec3(1.0, 0.5, 1.0) * 1.0, // Pink
	vec3(0.1, 0.1, 0.1) * vec3(0.1, 0.1, 0.1) * 1.0, // Black
	vec3(0.9, 0.9, 0.9) * vec3(0.9, 0.9, 0.9) * 1.0, // White
	vec3(0.3, 0.3, 0.3) * vec3(0.3, 0.3, 0.3) * 1.0, // Gray
	vec3(0.7, 0.7, 0.7) * vec3(0.7, 0.7, 0.7) * 1.0,  // Light gray
	vec3(1.00, 1.00, 1.00) * 256.0, // Lightning rod
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0)   // Unused
);
#endif

#if (defined(WORLD_END) && defined(END_COLORED_LIGHTING))
const vec3[64] light_color = vec3[64](
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) * 12.0, // Strong white light
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) *  6.0, // Medium white light
	vec3(1.00, 1.00, 1.00) * vec3(1.0, 1.0, 1.0) * 1.0, // Weak white light
	vec3(1.00, 0.65, 0.25) * vec3(1.00, 0.65, 0.25) * 12.0, // Strong golden light
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) *  MEDIUM_GOLDEN_LIGHT_I, // Medium golden light
	vec3(1.00, 0.57, 0.30) * vec3(1.00, 0.57, 0.30) *  2.2, // Weak golden light
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) *  4.0, // Redstone components
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * LAVA_LIGHT_I, // Lava
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * MEDIUM_ORANGE_LIGHT_I, // Medium orange light
	vec3(1.00, 0.63, 0.15) * vec3(1.00, 0.63, 0.15) * 4.0, // Brewing stand
	vec3(1.00, 0.57, 0.30) * vec3(1.00, 0.57, 0.30) * 12.0, // Medium golden light (Jack o' Lantern)
	vec3(0.45, 0.73, 1.00) * vec3(0.45, 0.73, 1.00) *  6.0, // Soul lights
	vec3(0.45, 0.73, 1.00) * vec3(0.45, 0.73, 1.00) * 14.0, // Beacon
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * 1.5, // End portal frame
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * SCULK_I, // Sculk
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * 5.0, // Pink glow
	vec3(0.75, 1.00, 0.50) * vec3(0.75, 1.00, 0.50) * 1.0, // Sea pickle
	vec3(1.00, 0.50, 0.25) * vec3(1.00, 0.50, 0.25) * 1.8, // Nether plants
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) * 1.0, // Restone wire
	vec3(0.97, 0.87, 0.57) * vec3(0.97, 0.87, 0.57) * 8.0, // Ochre froglight
	vec3(0.58, 0.77, 0.52) * vec3(0.58, 0.77, 0.52) * 8.0, // Verdant froglight
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 8.0, // Pearlescent froglight
	vec3(0.10, 0.86, 1.0) * vec3(0.10, 0.86, 1.0) * 2.0, // Enchanting table
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 4.0, // Amethyst cluster
	vec3(0.75, 0.44, 1.00) * vec3(0.75, 0.44, 1.00) * 4.0, // Calibrated sculk sensor
	vec3(0.75, 1.00, 0.83) * vec3(0.75, 1.00, 0.83) * 6.0, // Active sculk sensor
	vec3(1.00, 0.18, 0.10) * vec3(1.00, 0.18, 0.10) * 3.5, // Redstone block
	vec3(1.00, 0.10, 1.00) * vec3(1.00, 0.10, 1.00) * 3.0, // Purple weak light
	vec3(0.10, 0.10, 1.00) *  3.3, // Lapis block
	vec3(1.00, 1.00, 1.00) * vec3(1.00, 1.00, 1.00) * 64.0, // Lightning rod
	vec3(0.60, 0.10, 1.00) * vec3(0.60, 0.10, 1.00) * 12.0, // Nether portal
	vec3(0.0),  // End portal
	vec3(1.0, 0.1, 0.1) * vec3(1.0, 0.1, 0.1) *  1.5, // Red
	vec3(1.0, 0.5, 0.1) * vec3(1.0, 0.5, 0.1) *  1.5, // Orange
	vec3(1.0, 0.8, 0.2) * vec3(1.0, 0.8, 0.2) * 1.0, // Yellow
	vec3(0.7, 0.7, 0.0) * vec3(0.7, 0.7, 0.0) *  1.0, // Brown
	vec3(0.1, 1.0, 0.1) * vec3(0.1, 1.0, 0.1) * 1.0, // Green
	vec3(0.5, 1.0, 0.5) * vec3(0.5, 1.0, 0.5) * 1.0, // Lime
	vec3(0.1, 0.1, 1.0) * vec3(0.1, 0.1, 1.0) * 1.5, // Blue
	vec3(0.5, 0.5, 1.0) * vec3(0.5, 0.5, 1.0) * 1.0, // Light blue
	vec3(0.1, 1.0, 1.0) * vec3(0.1, 1.0, 1.0) * 1.0, // Cyan
	vec3(0.7, 0.1, 1.0) * vec3(0.7, 0.1, 1.0) * 1.0, // Purple
	vec3(1.0, 0.1, 1.0) * vec3(1.0, 0.1, 1.0) * 1.0, // Magenta
	vec3(1.0, 0.5, 1.0) * vec3(1.0, 0.5, 1.0) * 1.0, // Pink
	vec3(0.1, 0.1, 0.1) * vec3(0.1, 0.1, 0.1) * 1.0, // Black
	vec3(0.9, 0.9, 0.9) * vec3(0.9, 0.9, 0.9) * 1.0, // White
	vec3(0.3, 0.3, 0.3) * vec3(0.3, 0.3, 0.3) * 1.0, // Gray
	vec3(0.7, 0.7, 0.7) * vec3(0.7, 0.7, 0.7) * 1.0,  // Light gray
	vec3(1.00, 1.00, 1.00) * 256.0, // Lightning rod
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0),  // Unused
	vec3(0.0)   // Unused
);
#endif

const vec3[16] tint_color = vec3[16](
	vec3(1.0, 0.1, 0.1), // Red
	vec3(1.0, 0.5, 0.1), // Orange
	vec3(1.0, 1.0, 0.1), // Yellow
	vec3(0.7, 0.7, 0.0), // Brown
	vec3(0.1, 1.0, 0.1), // Green
	vec3(0.5, 1.0, 0.5), // Lime
	vec3(0.1, 0.1, 1.0), // Blue
	vec3(0.5, 0.5, 1.0), // Light blue
	vec3(0.1, 1.0, 1.0), // Cyan
	vec3(0.7, 0.1, 1.0), // Purple
	vec3(1.0, 0.1, 1.0), // Magenta
	vec3(1.0, 0.5, 1.0), // Pink
	vec3(0.1, 0.1, 0.1), // Black
	vec3(0.9, 0.9, 0.9), // White
	vec3(0.3, 0.3, 0.3), // Gray
	vec3(0.7, 0.7, 0.7)  // Light gray
);

#endif // INCLUDE_LIGHTING_LPV_LIGHT_COLORS
