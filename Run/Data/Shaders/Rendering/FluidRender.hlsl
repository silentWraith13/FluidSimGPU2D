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

cbuffer ScreenSizeConstants : register(b4)
{
    float screenWidth;
    float screenHeight;
	float time;
    float padding3;
};

StructuredBuffer<Particle> ParticlesRO : register(t0);
StructuredBuffer<ParticleDensity> ParticleDensityRO : register(t1);
SamplerState diffuseSampler : register(s0);
Texture2D fluidTexture : register(t0);
Texture2D noiseTexture : register(t1);

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
// Visualization Helper Functions
//--------------------------------------------------------------------------------------
static const float4 Rainbow[5] = 
{
	float4(1, 0, 0, 1), // red
	float4(1, 1, 0, 1), // orange
	float4(0, 1, 0, 1), // green
	float4(0, 1, 1, 1), // teal
	float4(0, 0, 1, 1), // blue
};

float4 VisualizeNumber(float n)
{
	return lerp(Rainbow[floor(n * 4.0f)], Rainbow[ceil(n * 4.0f)], frac(n * 4.0f));
}

float4 VisualizeNumber(float n, float lower, float upper)
{
	return VisualizeNumber(saturate((n - lower) / (upper - lower)));
}

//static const float4 WaterPalette[5] =
//{
//		float4(0.0f, 0.0f, 1.0f, 1.0f), // deep blue
//	float4(0.2f, 0.5f, 1.0f, 1.0f), // medium blue
//	float4(0.5f, 0.8f, 1.0f, 1.0f), // light blue
//	float4(0.7f, 0.9f, 1.0f, 1.0f), // very light blue
//	float4(0.9f, 1.0f, 1.0f, 1.0f), // almost white-blue
//};
static const float4 WaterPalette1[5] =
{
	float4(0.9f, 0.95f, 1.0f, 1.0f),          // White (surface)
	float4(0.8f, 0.9f, 0.95f, 1.0f),         // Very light blue
	float4(0.6f, 0.8f, 0.9f, 1.0f),          // Lighter blue
	float4(44.0f / 255.0f, 86.0f / 255.0f, 114.0f / 255.0f, 1.0f), // Light blue
	float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f)   // Start of medium blue (end of WaterPalette1)
};

// Light blue to medium blue transition
static const float4 WaterPalette2[5] =
{
	float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f),  // Start of medium blue (end of WaterPalette1)
	float4(7.0f / 255.0f, 65.0f / 255.0f, 95.0f / 255.0f, 1.0f),   // Medium blue
	float4(5.0f / 255.0f, 55.0f / 255.0f, 90.0f / 255.0f, 1.0f),   // Deeper blue
	float4(2.0f / 255.0f, 45.0f / 255.0f, 85.0f / 255.0f, 1.0f),   // Darker blue
	float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f)    // Start of navy blue (end of WaterPalette2)
};

// Medium blue to navy blue transition
static const float4 WaterPalette3[5] =
{
	float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f),    // Start of navy blue (end of WaterPalette2)
	float4(0.0f / 255.0f, 30.0f / 255.0f, 50.0f / 255.0f, 1.0f),    // Navy blue
	float4(0.0f / 255.0f, 25.0f / 255.0f, 45.0f / 255.0f, 1.0f),    // Darker navy blue
	float4(0.0f / 255.0f, 20.0f / 255.0f, 40.0f / 255.0f, 1.0f),    // Even darker navy blue
	float4(0.0f / 255.0f, 15.0f / 255.0f, 35.0f / 255.0f, 1.0f)     // Dark navy blue (high density)
};



float4 VisualizeNumberW1(float n)
{
	return lerp(WaterPalette1[floor(n * 4.0f)], WaterPalette1[ceil(n * 4.0f)], frac(n * 4.0f));
}

float4 VisualizeNumberW1(float n, float lower, float upper)
{
	return VisualizeNumberW1(saturate((n - lower) / (upper - lower)));
}

float4 VisualizeNumberW2(float n)
{
	return lerp(WaterPalette2[floor(n * 4.0f)], WaterPalette2[ceil(n * 4.0f)], frac(n * 4.0f));
}

float4 VisualizeNumberW2(float n, float lower, float upper)
{
	return VisualizeNumberW2(saturate((n - lower) / (upper - lower)));
}

float4 VisualizeNumberW3(float n)
{
	return lerp(WaterPalette3[floor(n * 4.0f)], WaterPalette3[ceil(n * 4.0f)], frac(n * 4.0f));
}

float4 VisualizeNumberW3(float n, float lower, float upper)
{
	return VisualizeNumberW3(saturate((n - lower) / (upper - lower)));
}
//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------

VSParticleOut VertexMain(uint ID : SV_VertexID)
{
	VSParticleOut Out = (VSParticleOut)0;
	Out.position = ParticlesRO[ID].position;
	
	if (debugType == 0)
	{
		Out.color = VisualizeNumber(ParticleDensityRO[ID].density, 1000.0f, 2000.0f);
	}
	
	if (debugType == 1)
	{
		if (ParticleDensityRO[ID].density >= 100.f && ParticleDensityRO[ID].density <= 1000.f)
		{
			Out.color = VisualizeNumberW1(ParticleDensityRO[ID].density,500.0f,1000.0f);
		}

		if (ParticleDensityRO[ID].density >= 1000.f && ParticleDensityRO[ID].density <= 1500.f)
		{
			Out.color = VisualizeNumberW2(ParticleDensityRO[ID].density, 1000.0f, 1500.0f);
		}

		if (ParticleDensityRO[ID].density >= 1500.f && ParticleDensityRO[ID].density <=4000.f)
		{
			Out.color = VisualizeNumberW3(ParticleDensityRO[ID].density, 1500.0f, 2000.0f);
		}
	}

	if (debugType == 2)
	{
		if (ParticleDensityRO[ID].density >= 100.f && ParticleDensityRO[ID].density <= 1000.f)
		{
			Out.color = VisualizeNumberW1(ParticleDensityRO[ID].density, 500.0f, 1000.0f);
		}

		if (ParticleDensityRO[ID].density >= 1000.f && ParticleDensityRO[ID].density <= 1500.f)
		{
			Out.color = VisualizeNumberW2(ParticleDensityRO[ID].density, 1000.0f, 1500.0f);
		}

		if (ParticleDensityRO[ID].density >= 1500.f && ParticleDensityRO[ID].density <= 4000.f)
		{
			Out.color = VisualizeNumberW3(ParticleDensityRO[ID].density, 1500.0f, 2000.0f);
		}
	}

	if (debugType == 3)
	{
		if (ParticleDensityRO[ID].density >= 100.f && ParticleDensityRO[ID].density <= 1000.f)
		{
			Out.color = VisualizeNumberW1(ParticleDensityRO[ID].density, 500.0f, 1000.0f);
		}

		if (ParticleDensityRO[ID].density >= 1000.f && ParticleDensityRO[ID].density <= 1500.f)
		{
			Out.color = VisualizeNumberW2(ParticleDensityRO[ID].density, 1000.0f, 1500.0f);
		}

		if (ParticleDensityRO[ID].density >= 1500.f && ParticleDensityRO[ID].density <= 4000.f)
		{
			Out.color = VisualizeNumberW3(ParticleDensityRO[ID].density, 1500.0f, 2000.0f);
		}
	}

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

//--------------------------------------------------------------------------------------
// Pixel Shader
//--------------------------------------------------------------------------------------

float4 PixelMain(GSParticleOut input) : SV_Target
{	
	float4 color = float4(0, 0, 0, 0);
	float4 noiseTextureColor = noiseTexture.Sample(diffuseSampler, input.uv);

	float2 animatedUV = input.uv + noiseTextureColor.xy * 0.1 + float2(time * 0.1, time * 0.1);

	if (debugType == 0)
	{
		color = input.color; // Density Debug
	}
	else if (debugType == 1)
	{
		color = input.color;  // Color palette
	}
	else if (debugType == 2)
	{
		color = fluidTexture.Sample(diffuseSampler, input.uv) * input.color * ModelColor; //Additive and gaussian blurred texture of a blue disc with varying opacity and a black background
	}
	if (debugType == 3)
	{
		color = fluidTexture.Sample(diffuseSampler, input.uv) * input.color * ModelColor;
	}

	return color;
}