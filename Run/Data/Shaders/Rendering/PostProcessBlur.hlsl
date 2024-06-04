//---------------------------------------------------------------------------------------
cbuffer CameraConstants : register(b2)
{
	float4x4 ProjectionMatrix;
	float4x4 ViewMatrix;
}
//---------------------------------------------------------------------------------------
cbuffer ModelConstants : register(b3)
{
	float4x4 ModelMatrix;
	float4 ModelColor;
}
//---------------------------------------------------------------------------------------
cbuffer ScreenSizeConstants : register(b4)
{
	float screenWidth;
	float screenHeight;
	float2 padding3;
};
//---------------------------------------------------------------------------------------
struct ParticleDensity
{
    float density;
    float3 padding2;
};
//---------------------------------------------------------------------------------------
struct vs_input_t
{
	float3 localPosition : POSITION;
	float4 color : COLOR;
	float2 uv : TEXCOORD;
};
//---------------------------------------------------------------------------------------
struct v2p_t
{
	float4 position : SV_POSITION;
	float4 color : COLOR;
	float2 uv : TEXCOORD;
};
//---------------------------------------------------------------------------------------
SamplerState diffuseSampler : register(s0);
//---------------------------------------------------------------------------------------
Texture2D diffuseTexture : register(t0);
//---------------------------------------------------------------------------------------
StructuredBuffer<ParticleDensity> ParticleDensityRO : register(t1);
//---------------------------------------------------------------------------------------
v2p_t VertexMain(vs_input_t input)
{
	float4 localPosition = float4(input.localPosition, 1);
	float4 worldPosition = mul(ModelMatrix, localPosition);
	float4 viewPosition = mul(ViewMatrix, worldPosition);
	float4 clipPosition = mul(ProjectionMatrix, viewPosition);
	v2p_t v2p;
	v2p.position = clipPosition;
	v2p.color = input.color;
	v2p.uv = input.uv;
	return v2p;
}
//---------------------------------------------------------------------------------------
float RangeMap(float inValue, float inStart, float inEnd, float outStart, float outEnd)
{
	if (inStart == inEnd)
	{
		return outStart + (outEnd - outStart) * 0.5;
	}
	float normalizedValue = (inValue - inStart) / (inEnd - inStart);

	return outStart + (outEnd - outStart) * normalizedValue;
}
//---------------------------------------------------------------------------------------
static const float4 WaterPalette1[5] =
{
	float4(0.9f, 0.95f, 1.0f, 1.0f),          // White (surface)
	float4(0.8f, 0.9f, 0.95f, 1.0f),         // Very light blue
	float4(0.6f, 0.8f, 0.9f, 1.0f),          // Lighter blue
	float4(44.0f / 255.0f, 86.0f / 255.0f, 114.0f / 255.0f, 1.0f), // Light blue
	float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f)   // Start of medium blue (end of WaterPalette1)
};
//---------------------------------------------------------------------------------------
// Light blue to medium blue transition
static const float4 WaterPalette2[5] =
{
	float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f),  // Start of medium blue (end of WaterPalette1)
	float4(7.0f / 255.0f, 65.0f / 255.0f, 95.0f / 255.0f, 1.0f),   // Medium blue
	float4(5.0f / 255.0f, 55.0f / 255.0f, 90.0f / 255.0f, 1.0f),   // Deeper blue
	float4(2.0f / 255.0f, 45.0f / 255.0f, 85.0f / 255.0f, 1.0f),   // Darker blue
	float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f)    // Start of navy blue (end of WaterPalette2)
};
//---------------------------------------------------------------------------------------
// Medium blue to navy blue transition
static const float4 WaterPalette3[5] =
{
	float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f),    // Start of navy blue (end of WaterPalette2)
	float4(0.0f / 255.0f, 30.0f / 255.0f, 50.0f / 255.0f, 1.0f),    // Navy blue
	float4(0.0f / 255.0f, 25.0f / 255.0f, 45.0f / 255.0f, 1.0f),    // Darker navy blue
	float4(0.0f / 255.0f, 20.0f / 255.0f, 40.0f / 255.0f, 1.0f),    // Even darker navy blue
	float4(0.0f / 255.0f, 15.0f / 255.0f, 35.0f / 255.0f, 1.0f)     // Dark navy blue (high density)
};
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW1(float n)
{
	return lerp(WaterPalette1[floor(n * 4.0f)], WaterPalette1[ceil(n * 4.0f)], frac(n * 4.0f));
}
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW1(float n, float lower, float upper)
{
	return VisualizeNumberW1(saturate((n - lower) / (upper - lower)));
}
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW2(float n)
{
	return lerp(WaterPalette2[floor(n * 4.0f)], WaterPalette2[ceil(n * 4.0f)], frac(n * 4.0f));
}
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW2(float n, float lower, float upper)
{
	return VisualizeNumberW2(saturate((n - lower) / (upper - lower)));
}
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW3(float n)
{
	return lerp(WaterPalette3[floor(n * 4.0f)], WaterPalette3[ceil(n * 4.0f)], frac(n * 4.0f));
}
//---------------------------------------------------------------------------------------
float4 VisualizeNumberW3(float n, float lower, float upper)
{
	return VisualizeNumberW3(saturate((n - lower) / (upper - lower)));
}
//---------------------------------------------------------------------------------------
float4 InterpolateColor(float value)
{
	const int TOTAL_COLORS = 15;
    const int PALETTE_SIZE = 5;
    const float4 colorPalette[TOTAL_COLORS] =
	{
        // WaterPalette1
        float4(0.9f, 0.95f, 1.0f, 1.0f),
        float4(0.8f, 0.9f, 0.95f, 1.0f),
        float4(0.6f, 0.8f, 0.9f, 1.0f),
        float4(44.0f / 255.0f, 86.0f / 255.0f, 114.0f / 255.0f, 1.0f),
        float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f),

        // WaterPalette2
        float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f),
        float4(7.0f / 255.0f, 65.0f / 255.0f, 95.0f / 255.0f, 1.0f),
        float4(5.0f / 255.0f, 55.0f / 255.0f, 90.0f / 255.0f, 1.0f),
        float4(2.0f / 255.0f, 45.0f / 255.0f, 85.0f / 255.0f, 1.0f),
        float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f),

        // WaterPalette3
		float4(10.0f / 255.0f, 75.0f / 255.0f, 99.0f / 255.0f, 1.0f),
		float4(7.0f / 255.0f, 65.0f / 255.0f, 95.0f / 255.0f, 1.0f),
		float4(5.0f / 255.0f, 55.0f / 255.0f, 90.0f / 255.0f, 1.0f),
		float4(2.0f / 255.0f, 45.0f / 255.0f, 85.0f / 255.0f, 1.0f),
		float4(0.0f / 255.0f, 34.0f / 255.0f, 57.0f / 255.0f, 1.0f),
    };

	int numColors = 15;

	// Map 'value' to the range of indices in the color palette array.
	value = saturate(value) * (numColors - 1);
	int idx = (int)value;
	float lerpFactor = frac(value);

	// Get the two nearest colors from the palette to interpolate between.
	float4 color1 = colorPalette[max(0, idx)];
	float4 color2 = colorPalette[min(idx + 1, numColors - 1)];

	// Return the interpolated color.
	return lerp(color1, color2, lerpFactor);
}
//---------------------------------------------------------------------------------------
float4 PixelMain(v2p_t input) : SV_Target0
{
	 float4 sampledColor = diffuseTexture.Sample(diffuseSampler, input.uv);

	 // If the color is black, return black
	 if (sampledColor.r == 0 && sampledColor.g == 0 && sampledColor.b == 0)
	 {
		 return float4(0, 0, 0, 1);
	 }

	 // Use the blue channel to determine the interpolated color
	 float4 color = InterpolateColor(sampledColor.b);

	 return color;
}
//---------------------------------------------------------------------------------------