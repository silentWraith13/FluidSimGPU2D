cbuffer CameraConstants : register(b2)
{
	float4x4 ProjectionMatrix;
	float4x4 ViewMatrix;
}

cbuffer ParticleSizeConstants : register(b5)
{
	float particleSize;
	int  debugType;
	float2 padding;
}

cbuffer ModelConstants : register(b3)
{
	float4x4 ModelMatrix;
	float4 ModelColor;
}

struct Particle 
{
	float2 position;
	float2 velocity;
	float type;
	float3 padding;
};

struct ParticleDensity 
{
	float density;
	float3 padding2;
};


StructuredBuffer<ParticleDensity> ParticleDensityRO : register(t1);
SamplerState diffuseSampler : register(s0);

struct VSParticleOut
{
	float2 position : POSITION;
	float density : DENSITY;
};

struct PSParticle
{
	float4 position : SV_POSITION;
	float4 color : COLOR;
	float2 uv : TEXCOORD;
};

//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------

VSParticleOut VertexMain(uint ID : SV_VertexID)
{
	VSParticleOut Out;
	Out.position = float2(0,0);
	Out.density = ParticleDensityRO[ID].density;

	return Out;
}

//--------------------------------------------------------------------------------------
// Pixel Shader
//--------------------------------------------------------------------------------------

float4 PixelMain(VSParticleOut input) : SV_Target
{	
	 float normalizedDensity = saturate(input.density / 1000.f);

	float4 color = float4(1,1,1,1);

	return color;
}