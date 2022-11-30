shader_type spatial;
render_mode depth_draw_alpha_prepass;

uniform sampler2D albedoTexture : hint_albedo;
uniform sampler2D normalTexture;

uniform sampler2D dissolveTexture : hint_albedo;
uniform float dissolveValue : hint_range(0, 1);

uniform float outlineWidth;
uniform vec4 growColor : hint_color;

void fragment() {
    vec4 mainTex = texture(albedoTexture, UV);
    vec4 normalMap = texture(normalTexture, UV);
    vec4 dissolveTex = texture(dissolveTexture, UV);

    float cutoff = dot(dissolveTex.rgb, vec3(0.3, 0.3, 0.3)) -
        float(-0.4 + clamp(dissolveValue, 0, 1));

    vec3 dissolveOutline = vec3(round(1.0 - float(cutoff - outlineWidth))) *
        growColor.rgb;

    ALBEDO = mainTex.rgb;
    NORMALMAP = normalMap.xyz;
    METALLIC = 0.5-(ALBEDO.r); // This should use a map
    ROUGHNESS = 0.6+(ALBEDO.r); // This should use a map
    ALPHA = round(cutoff) * mainTex.a;
    EMISSION = dissolveOutline;
}
