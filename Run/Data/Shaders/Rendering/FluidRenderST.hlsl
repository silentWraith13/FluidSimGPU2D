cbuffer CameraConstants : register(b2)
{
	float4x4 ProjectionMatrix;
	float4x4 ViewMatrix;
}

cbuffer ParticleSizeConstants : register(b5)
{
	float particleSize;
	float3 padding;
}

cbuffer ModelConstants : register(b3)
{
	float4x4 ModelMatrix;
	float4 ModelColor;
}

cbuffer Solver : register(b6)
{
	float2 position;
	float density;
	float padding1;
}

struct VSParticleOut
{
	float2 position : POSITION;
	float4 color : COLOR;
	float2 uv : TEXCOORD;
};

struct GSParticleOut
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
	VSParticleOut Out = (VSParticleOut)0;
	Out.position = position;
	Out.color = float4(1, 0, 0, 1);
	return Out;
}

//--------------------------------------------------------------------------------------
// Particle Geometry Shader
//--------------------------------------------------------------------------------------
static const float2 g_positions[4] = { float2(-1, 1), float2(1, 1), float2(-1, -1), float2(1, -1) };
static const float2 g_texcoords[4] = { float2(0, 1), float2(1, 1), float2(0, 0), float2(1, 0) };

[maxvertexcount(4)]
void GeometryMain(point VSParticleOut In[1], inout TriangleStream<GSParticleOut> SpriteStream)
{
	[unroll]
	for (int i = 0; i < 4; i++)
	{
		GSParticleOut Out = (GSParticleOut)0;
		float4 position = float4(In[0].position, 0, 1) + particleSize * float4(g_positions[i], 0, 0);
		float4 viewPosition = mul(ViewMatrix, position);
		float4 clipPosition = mul(ProjectionMatrix, viewPosition);
		Out.position = clipPosition;
		Out.uv = g_texcoords[i];
		Out.color = In[0].color;
		SpriteStream.Append(Out);
	}
	SpriteStream.RestartStrip();
}


float4 PixelMain(GSParticleOut input) : SV_Target
{	
	return input.color;
}