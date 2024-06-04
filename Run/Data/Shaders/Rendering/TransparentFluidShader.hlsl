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
	float time;
	float padding3;
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
Texture2D noiseTexture : register(t1);
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
float4 PixelMain(v2p_t input) : SV_Target0
{
	float4 noiseTextureColor = noiseTexture.Sample(diffuseSampler, input.uv);
	float2 animatedUV = input.uv + noiseTextureColor.xy * 0.1 + float2(time * 0.1, time * 0.1);
    
	float4 color = diffuseTexture.Sample(diffuseSampler, input.uv);
    float4 white = float4(1, 1, 1, 1);
    float4 vibrantBlue = float4(0.1, 0.6, 0.9, 1); 
    float4 navyBlue = float4(0.0, 0.2, 0.4, 1);   
    float minWaterThreshold = 0.2f;
    float midWaterThreshold = 0.4f;
    float maxWaterThreshold = 0.9f;
	
	if (color.r == 0 && color.g == 0 && color.b == 0)
	{
		return float4(0, 0, 0, 0);
	}
    
	if (color.r == 0 && color.g == 0)
    {
        if (color.b <= minWaterThreshold && color.b >= 0.1)
        {
            color = white;
        }
        else if (color.b >= maxWaterThreshold)
        {
            color = navyBlue;
        }
        else if (color.b <= midWaterThreshold)
        {
            float edgeness = RangeMap(color.b, minWaterThreshold, midWaterThreshold, 0, 1);
            edgeness = saturate(edgeness);
            color = lerp(white, vibrantBlue, smoothstep(0, 1, edgeness));
        }
        else
        {
            float edgeness = RangeMap(color.b, midWaterThreshold, maxWaterThreshold, 0, 1);
            edgeness = saturate(edgeness);
            color = lerp(vibrantBlue, navyBlue, smoothstep(0, 1, edgeness));
        }
    }
	color.a = 0.7;
    return color;
}
//---------------------------------------------------------------------------------------