#version 150 core
uniform vec4 Lpos = vec4(0.0, 0.0, 5.0, 1.0);
uniform vec3 Lamb = vec3(0.2);
uniform vec3 Ldiff = vec3(1.0);
uniform vec3 Lspec = vec3(1.0);
layout (std140) uniform Material {
	vec3 Kamb;
	vec3 Kdiff;
	vec3 Kspec;
	float Kshi;
};
in vec4 P;
in vec3 N;
out vec4 fragment;
void main() {
   vec3 V = -normalize(P.xyz);
   vec3 L = normalize((Lpos * P.w - P * Lpos.w).xyz);
   vec3 Iamb = Kamb * Lamb;
   vec3 Idiff = max(dot(N, L), 0.0) * Kdiff * Ldiff + Iamb;
   vec3 H = normalize(L + V);
   vec3 Ispec = pow(max(dot(normalize(N), H), 0.0), Kshi) * Kspec * Lspec;
   fragment = vec4(Idiff + Ispec, 1.0);
}
