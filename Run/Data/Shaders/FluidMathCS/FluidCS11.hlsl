//------------------------------------------------------------------------------------------------------
struct Particle
{
    float2 position;
    float2 velocity;
    float type;
    float3 padding;
};
//------------------------------------------------------------------------------------------------------
struct ParticleForces
{
    float2 acceleration;
};
//------------------------------------------------------------------------------------------------------
struct ParticleDensity
{
    float density;
    float3 padding;
};
//------------------------------------------------------------------------------------------------------
cbuffer cbSimulationConstants : register( b0 )
{
    uint          g_iNumParticles;
    float         g_fTimeStep;
    float         g_fSmoothlen;
    float         g_fPressureStiffness;
                  
    float         g_fRestDensity;
    float         g_fDensityCoef;
    float         g_fGradPressureCoef;
    float         g_fLapViscosityCoef;
                  
    float         g_surfaceTensionCoef;
    float         g_fWallStiffness;
                  
    float4        g_vGravity;
                  
    float4        g_vGridDim;
                  
    float3        g_vPlanes1;
                  
    float3        g_vPlanes2;
                  
    float3        g_vPlanes3;
                  
    float3        g_vPlanes4;
                  
    float2        mousePos;
    int           isHoldingLeftClick;
                  
    int           isHoldingRightClick;

    float2        tumblingConfigOBB1HalfDimensions;
    float2        tumblingConfigOBB1Center;
    float4x4      tumblingConfigOBB1LocalToWorldMatrix;
    float4x4      tumblingConfigOBB1WorldToLocalMatrix;

    float2        tumblingConfigOBB2HalfDimensions;
    float2        tumblingConfigOBB2Center;
    float4x4      tumblingConfigOBB2LocalToWorldMatrix;
    float4x4      tumblingConfigOBB2WorldToLocalMatrix;
    
    float2        tumblingConfigOBB3HalfDimensions;
    float2        tumblingConfigOBB3Center;
    float4x4      tumblingConfigOBB3LocalToWorldMatrix;
    float4x4      tumblingConfigOBB3WorldToLocalMatrix;
    
    float2        tumblingConfigOBB4HalfDimensions;
    float2        tumblingConfigOBB4Center;
    float4x4      tumblingConfigOBB4LocalToWorldMatrix;
    float4x4      tumblingConfigOBB4WorldToLocalMatrix;
   
    float2        tumblingConfigOBB5HalfDimensions;
    float2        tumblingConfigOBB5Center;
    float4x4      tumblingConfigOBB5LocalToWorldMatrix;
    float4x4      tumblingConfigOBB5WorldToLocalMatrix;

    float2        tumblingConfigCapsule1Start;
    float2        tumblingConfigCapsule1End;

    float2        tumblingConfigCapsule2Start;
    float2        tumblingConfigCapsule2End;
    
    float         tumblingConfigCapsule1Radius;
    float         tumblingConfigCapsule2Radius;
    float2        padding;

    float2		  tumblingConfigHorizontalMovingAABB1Size;
    float2		  tumblingConfigHorizontalMovingAABB1Center;
                  
    float2		  tumblingDiscCenter;
    float		  tumblingDiscRadius;
                  
    float2		  windmillDisc1Center;
    float		  windmillDisc1Radius;

    float2		  windmillConfigOBB1HalfDimensions;
    float2		  windmillConfigOBB1Center;

    float4x4	  windmillConfigOBB1LocalToWorldMatrix;
    float4x4	  windmillConfigOBB1WorldToLocalMatrix;

    float2		  windmillConfigOBB2HalfDimensions;
    float2		  windmillConfigOBB2Center;

    float4x4	  windmillConfigOBB2LocalToWorldMatrix;
    float4x4	  windmillConfigOBB2WorldToLocalMatrix;

    float2		windmillConfigBottomBorderCapsuleStart;
    float2		windmillConfigBottomBorderCapsuleEnd;

    float		windmillConfigBottomBorderCapsuleRadius;
    float3       windmillBorderCapsulePAdding;

    float2		windmillConfigLeftBorderCapsuleStart;
    float2		windmillConfigLeftBorderCapsuleEnd;

    float		windmillConfigLeftBorderCapsuleRadius;
    float3		windmillConfigLeftBorderCapsulePadding;
   
    float2		windmillConfigRightBorderCapsuleStart;
    float2		windmillConfigRightBorderCapsuleEnd;

    float		windmillConfigRightBorderCapsuleRadius;
    float3		windmillConfigRightBorderCapsulePadding;

    float2		windmillConfigTopBorderCapsule1Start;
    float2		windmillConfigTopBorderCapsule1End;
    
    float		windmillConfigTopBorderCapsule1Radius;
    float3		windmillConfigTopBorderCapsule1Padding;

    float2		windmillConfigTopBorderCapsule2Start;
    float2		windmillConfigTopBorderCapsule2End;
   
    float		windmillConfigTopBorderCapsule2Radius;
    float3		windmillConfigTopBorderCapsule2Padding;

    float2		windmillConfigTopBorderCapsule3Start;
    float2		windmillConfigTopBorderCapsule3End;

    float		windmillConfigTopBorderCapsule3Radius;
    float3		windmillConfigTopBorderCapsule3Padding;

    float2		hourglassConfigCapsule1Start;
    float2		hourglassConfigCapsule1End;

    float		hourglassConfigCapsule1Radius;
    float3		hourglassConfigCapsule1Padding;

    float2		hourglassConfigCapsule2Start;
    float2		hourglassConfigCapsule2End;

    float		hourglassConfigCapsule2Radius;
    float3		hourglassConfigCapsule2Padding;

    float2		hourglassConfigCapsule3Start;
    float2		hourglassConfigCapsule3End;

    float		hourglassConfigCapsule3Radius;
    float3		hourglassConfigCapsule3Padding;

    float2		hourglassConfigCapsule4Start;
    float2		hourglassConfigCapsule4End;

    float		hourglassConfigCapsule4Radius;
    float3		hourglassConfigCapsule4Padding;

    float2		hourglassConfigCapsule5Start;
    float2		hourglassConfigCapsule5End;

    float		hourglassConfigCapsule5Radius;
    float3		hourglassConfigCapsule5Padding;

    float2		hourglassConfigCapsule6Start;
    float2		hourglassConfigCapsule6End;

    float		hourglassConfigCapsule6Radius;
    float3		hourglassConfigCapsule6Padding;

    float2		hourglassConfigMiddleCapsule1Start;
    float2		hourglassConfigMiddleCapsule1End;

    float		hourglassConfigMiddleCapsule1Radius;
    float3		hourglassConfigMiddleCapsule1Padding;

    float2		hourglassConfigMiddleCapsule2Start;
    float2		hourglassConfigMiddleCapsule2End;

    float		hourglassConfigMiddleCapsule2Radius;
    float3		hourglassConfigMiddleCapsule2Padding;

    float2		waterElevatorConfigVerticalMovingAABBSize;
    float2		waterElevatorConfigVerticalMovingAABBCenter;

    float2		waterElevatorDisc1Center;
    float		waterElevatorDisc1Radius;
    float		waterElevator1padding;

    float2		waterElevatorDisc2Center;
    float		waterElevatorDisc2Radius;
    float		waterElevator2padding;

    float2		waterElevatorDisc3Center;
    float		waterElevatorDisc3Radius;
    float		waterElevator3padding;
};
//------------------------------------------------------------------------------------------------------
#define SIMULATION_BLOCK_SIZE 256
//--------------------------------------------------------------------------------------
// Structured Buffers
//--------------------------------------------------------------------------------------
RWStructuredBuffer<Particle> ParticlesRW : register(u0);
StructuredBuffer<Particle> ParticlesRO : register(t0);

RWStructuredBuffer<ParticleDensity> ParticlesDensityRW : register(u0);
StructuredBuffer<ParticleDensity> ParticlesDensityRO : register(t1);

RWStructuredBuffer<ParticleForces> ParticlesForcesRW : register(u0);
StructuredBuffer<ParticleForces> ParticlesForcesRO : register(t2);

RWStructuredBuffer<unsigned int> GridRW : register(u0);
StructuredBuffer<unsigned int> GridRO : register(t3);

RWStructuredBuffer<uint2> GridIndicesRW : register(u0);
StructuredBuffer<uint2> GridIndicesRO : register(t4);


//--------------------------------------------------------------------------------------
// Grid Construction
//--------------------------------------------------------------------------------------

float2 GridCalculateCell(float2 position)
{
    return clamp(position * g_vGridDim.xy + g_vGridDim.zw, float2(0, 0), float2(255, 255));
}

unsigned int GridConstructKey(uint2 xy)
{
    return dot(xy.yx, uint2(256, 1));
}

unsigned int GridConstructKeyValuePair(uint2 xy, uint value)
{
    return dot(uint3(xy.yx, value), uint3(256 * 256 * 256, 256 * 256, 1));
}

unsigned int GridGetKey(unsigned int keyvaluepair)
{
    return (keyvaluepair >> 16);
}

unsigned int GridGetValue(unsigned int keyvaluepair)
{
    return (keyvaluepair & 0xFFFF);
}


//--------------------------------------------------------------------------------------
// Build Grid
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void BuildGridCS(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int P_ID = DTid.x;

    float2 position = ParticlesRO[P_ID].position;
    float2 grid_xy = GridCalculateCell(position);

    GridRW[P_ID] = GridConstructKeyValuePair((uint2)grid_xy, P_ID);
}


//--------------------------------------------------------------------------------------
// Build Grid Indices
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void ClearGridIndicesCS(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    GridIndicesRW[DTid.x] = uint2(0, 0);
}

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void BuildGridIndicesCS(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int G_ID = DTid.x; // Grid ID to operate on
    unsigned int G_ID_PREV = (G_ID == 0) ? g_iNumParticles : G_ID; G_ID_PREV--;
    unsigned int G_ID_NEXT = G_ID + 1; if (G_ID_NEXT == g_iNumParticles) { G_ID_NEXT = 0; }

    unsigned int cell = GridGetKey(GridRO[G_ID]);
    unsigned int cell_prev = GridGetKey(GridRO[G_ID_PREV]);
    unsigned int cell_next = GridGetKey(GridRO[G_ID_NEXT]);
    if (cell != cell_prev)
    {
        // I'm the start of a cell
        GridIndicesRW[cell].x = G_ID;
    }
    if (cell != cell_next)
    {
        // I'm the end of a cell
        GridIndicesRW[cell].y = G_ID + 1;
    }
}


//--------------------------------------------------------------------------------------
// Rearrange Particles
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void RearrangeParticlesCS(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int ID = DTid.x;
    const unsigned int G_ID = GridGetValue(GridRO[ID]);
    ParticlesRW[ID] = ParticlesRO[G_ID];
}


//--------------------------------------------------------------------------------------
// Density Calculation
//--------------------------------------------------------------------------------------

float CalculateDensity(float r_sq)
{
    const float h_sq = g_fSmoothlen * g_fSmoothlen;
    // Implements this equation:
    // W_poly6(r, h) = 315 / (64 * pi * h^9) * (h^2 - r^2)^3
    // g_fDensityCoef = fParticleMass * 315.0f / (64.0f * PI * fSmoothlen^9)
    return g_fDensityCoef * (h_sq - r_sq) * (h_sq - r_sq) * (h_sq - r_sq);
}


//--------------------------------------------------------------------------------------
// Simple N^2 Algorithm
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void DensityCS_Simple(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{

    const unsigned int P_ID = DTid.x;
    const float h_sq = g_fSmoothlen * g_fSmoothlen;
    float2 P_position = ParticlesRO[P_ID].position;

    float density = 0;

    // Calculate the density based on all neighbors
    for (uint N_ID = 0; N_ID < g_iNumParticles; N_ID++)
    {
        float2 N_position = ParticlesRO[N_ID].position;

        float2 diff = N_position - P_position;
        float r_sq = dot(diff, diff);
        if (r_sq < h_sq)
        {
            density += CalculateDensity(r_sq);
        }
    }

    ParticlesDensityRW[P_ID].density = density;
}


//--------------------------------------------------------------------------------------
// Optimized Grid + Sort Algorithm
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void DensityCS_Grid(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int P_ID = DTid.x;
    const float h_sq = g_fSmoothlen * g_fSmoothlen;
    float2 P_position = ParticlesRO[P_ID].position;
    float type = ParticlesRO[P_ID].type;
    float density = 0;

    // Calculate the density based on neighbors from the 8 adjacent cells + current cell
    int2 G_XY = (int2)GridCalculateCell(P_position);
    for (int Y = max(G_XY.y - 1, 0); Y <= min(G_XY.y + 1, 255); Y++)
    {
        for (int X = max(G_XY.x - 1, 0); X <= min(G_XY.x + 1, 255); X++)
        {
            unsigned int G_CELL = GridConstructKey(uint2(X, Y));
            uint2 G_START_END = GridIndicesRO[G_CELL];
            for (unsigned int N_ID = G_START_END.x; N_ID < G_START_END.y; N_ID++)
            {
                float2 N_position = ParticlesRO[N_ID].position;

                float2 diff = N_position - P_position;
                float r_sq = dot(diff, diff);
                if (r_sq < h_sq)
                {

                    density += CalculateDensity(r_sq);


                }
            }
        }
    }

    ParticlesDensityRW[P_ID].density = density;
}


//--------------------------------------------------------------------------------------
// Force Calculation
//--------------------------------------------------------------------------------------

float CalculatePressure(float density)
{
    // Implements this equation:
    // Pressure = B * ((rho / rho_0)^y  - 1)
    return g_fPressureStiffness * max(pow(density / g_fRestDensity, 3) - 1, 0);
}

float2 CalculateGradPressure(float r, float P_pressure, float N_pressure, float N_density, float2 diff)
{
    const float h = g_fSmoothlen;
    float avg_pressure = 0.5f * (N_pressure + P_pressure);

    if (diff.x == 0 && diff.y==0)
    {
        diff = float2(0.001, 0);
    }
    // Implements this equation:
    // W_spkiey(r, h) = 15 / (pi * h^6) * (h - r)^3
    // GRAD( W_spikey(r, h) ) = -45 / (pi * h^6) * (h - r)^2
    // g_fGradPressureCoef = fParticleMass * -45.0f / (PI * fSmoothlen^6)
    return g_fGradPressureCoef * avg_pressure / N_density * (h - r) * (h - r) / r * (diff);
}

float2 CalculateLapVelocity(float r, float2 P_velocity, float2 N_velocity, float N_density)
{
    const float h = g_fSmoothlen;
    float2 vel_diff = (N_velocity - P_velocity);
    // Implements this equation:
    // W_viscosity(r, h) = 15 / (2 * pi * h^3) * (-r^3 / (2 * h^3) + r^2 / h^2 + h / (2 * r) - 1)
    // LAPLACIAN( W_viscosity(r, h) ) = 45 / (pi * h^6) * (h - r)
    // g_fLapViscosityCoef = fParticleMass * fViscosity * 45.0f / (PI * fSmoothlen^6)
    return g_fLapViscosityCoef / N_density * (h - r) * vel_diff;
}

float2 CalculateSurfaceTension(float r, float2 diff, float P_density, float N_density)
{
    const float h = g_fSmoothlen;
    float r_norm = r / h;
    float2 surfaceTensionForce = g_surfaceTensionCoef * (h - r) * normalize(diff) / (P_density + N_density);

    return surfaceTensionForce;
}

//--------------------------------------------------------------------------------------
// Simple N^2 Algorithm
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void ForceCS_Simple(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int P_ID = DTid.x; // Particle ID to operate on

    float2 P_position = ParticlesRO[P_ID].position;
    float2 P_velocity = ParticlesRO[P_ID].velocity;
    float P_density = ParticlesDensityRO[P_ID].density;
    float P_pressure = CalculatePressure(P_density);

    const float h_sq = g_fSmoothlen * g_fSmoothlen;

    float2 acceleration = float2(0, 0);

    // Calculate the acceleration based on all neighbors
    for (uint N_ID = 0; N_ID < g_iNumParticles; N_ID++)
    {
        float2 N_position = ParticlesRO[N_ID].position;

        float2 diff = N_position - P_position;
        float r_sq = dot(diff, diff);
        if (r_sq < h_sq && P_ID != N_ID)
        {
            float2 N_velocity = ParticlesRO[N_ID].velocity;
            float N_density = ParticlesDensityRO[N_ID].density;
            float N_pressure = CalculatePressure(N_density);
            float r = sqrt(r_sq);

            // Pressure Term
            acceleration += CalculateGradPressure(r, P_pressure, N_pressure, N_density, diff);

            // Viscosity Term
            acceleration += CalculateLapVelocity(r, P_velocity, N_velocity, N_density);
        }
    }

    ParticlesForcesRW[P_ID].acceleration = acceleration / P_density;
}
//------------------------------------------------------------------------------------------------------
float2 ApplyMouseForce(float2 mousePos, float radius, float strength, float2 position, float2 velocity)
{
    float2 acc = float2(0, 0);
    float2 offset = mousePos - position;
    float sqrDist = dot(offset, offset);

    if (sqrDist < radius * radius)
    {
        float dist = sqrt(sqrDist);
        float2 dirToInputPoint = dist <= 0.0001 ? float2(0, 0) : offset / dist;
        float centreT = 1 - dist / radius;

        float gravityWeight = 1 - (centreT * saturate(strength / 10));
        acc = g_vGravity.xy * gravityWeight + dirToInputPoint * centreT * strength;
        acc -= velocity * centreT;
    }

    return acc;
}
//------------------------------------------------------------------------------------------------------
float2 ApplyExternalForces(float2 position, float2 velocity, float2 acceleration)
{
    float2 gravAcceleration = g_vGravity.xy;
    if (isHoldingLeftClick == 1)
    {
        acceleration += ApplyMouseForce(mousePos, 0.3, 30, position, velocity);

    }
    if (isHoldingRightClick == 1)
    {
        acceleration += ApplyMouseForce(mousePos, 0.3, -5, position, velocity);

    }
    acceleration += gravAcceleration;
    return acceleration;

}
//-----------------------------------------------------------------------------------
float2 ComputeWallCollisionForce(float2 position)
{
    float2 wallAccumulatedForce = float2(0, 0);

    float dist1 = dot(float3(position, 1), g_vPlanes1);
    wallAccumulatedForce += min(dist1, 0) * -g_fWallStiffness * g_vPlanes1.xy;

    float dist2 = dot(float3(position, 1), g_vPlanes2);
    wallAccumulatedForce += min(dist2, 0) * -g_fWallStiffness * g_vPlanes2.xy;

    float dist3 = dot(float3(position, 1), g_vPlanes3);
    wallAccumulatedForce += min(dist3, 0) * -g_fWallStiffness * g_vPlanes3.xy;

    float dist4 = dot(float3(position, 1), g_vPlanes4);
    wallAccumulatedForce += min(dist4, 0) * -g_fWallStiffness * g_vPlanes4.xy;

    return wallAccumulatedForce;
}
//-----------------------------------------------------------------------------------
void ProcessDiscCollision(inout float2 position, inout float2 velocity, float2 discCenter, float discRadius)
{
    float2 toDisc = position - discCenter;
    float distanceToDisc = length(toDisc);

    if (distanceToDisc < discRadius)
    {
        float overlap = discRadius - distanceToDisc;
        toDisc = normalize(toDisc) * overlap;
        position += toDisc;
        float distanceToDiscSq = dot(toDisc, toDisc);
        velocity -= 2.0f * dot(velocity, toDisc) * toDisc / distanceToDiscSq;
    }
}
//-----------------------------------------------------------------------------------
void ProcessAABBCollision(inout float2 position, inout float2 velocity, float2 AABBcenter, float2 AABBsize)
{
    const float2 halfSize = AABBsize * 0.5;
    float2 edgeDistance = halfSize - abs(position - AABBcenter);

    if (edgeDistance.x >= 0 && edgeDistance.y >= 0)
    {
        if (edgeDistance.x < edgeDistance.y)
        {
            position.x = halfSize.x * sign(position.x - AABBcenter.x) + AABBcenter.x;
            velocity.x *= -1;
        }
        else
        {
            position.y = halfSize.y * sign(position.y - AABBcenter.y) + AABBcenter.y;
            velocity.y *= -1;
        }
    }
}
//-----------------------------------------------------------------------------------
void ProcessOBBCollision(inout float2 position, inout float2 velocity, float4x4 worldToLocal, float4x4 localToWorld, float2 OBBcenter, float2 OBBsize)
{
    float4 posLocal4 = mul(worldToLocal, float4(position, 0, 1));
    float2 posLocal = posLocal4.xy;

    // Collision detection in local space
    const float2 halfSize = OBBsize * 0.5;
    float2 edgeDst = halfSize - abs(posLocal);

    if (edgeDst.x >= 0 && edgeDst.y >= 0)
    {
        // Transform velocity to local space
        float2 velLocal = mul(worldToLocal, float4(velocity, 0, 0)).xy;

        // Collision response in local space
        if (edgeDst.x < edgeDst.y)
        {
            posLocal.x = halfSize.x * sign(posLocal.x);
            velLocal.x *= -1;
        }
        else
        {
            posLocal.y = halfSize.y * sign(posLocal.y);
            velLocal.y *= -1;
        }

        // Transform the updated position and velocity back to world space
        position = mul(localToWorld, float4(posLocal, 0, 1)).xy;
        velocity = mul(localToWorld, float4(velLocal, 0, 0)).xy;
    }
}
//----------------------------------------------------------------------------
float2 GetNearestPointOnLineSegment2D(float2 referencePos, float2 lineSegStart, float2 lineSegEnd)
{
    float2 result = lineSegStart;

    float2 SE = lineSegEnd - lineSegStart;
    float2 SP = referencePos - lineSegStart;
    float2 EP = referencePos - lineSegEnd;

    if (dot(SP, SE) <= 0.0f)
    {
        result = lineSegStart;
    }
   
    else if (dot(EP, SE) >= 0.0f)
    {
        result = lineSegEnd;
    }
    else
    {
        float2 projection = dot(SP, SE) / dot(SE, SE) * SE;
        result = lineSegStart + projection;
    }

    return result;
}
//----------------------------------------------------------------------------
float2 GetNearestPointInCapsule2D(float2 referencePos, float2 capsuleStart, float2 capsuleEnd, float radius)
{
    float2 nearestPointOnLine = GetNearestPointOnLineSegment2D(referencePos, capsuleStart, capsuleEnd);
    float2 toReference = referencePos - nearestPointOnLine;
    float distanceToReference = length(toReference);

    if (distanceToReference <= radius)
    {
        return nearestPointOnLine;
    }
    else
    {
        float2 clampedDisplacement = normalize(toReference) * radius;
        return nearestPointOnLine;
    }
}
//----------------------------------------------------------------------------
void ProcessCapsuleCollision(inout float2 position, inout float2 velocity, float2 capsuleStart, float2 capsuleEnd, float capsuleRadius)
{
    float2 nearestPoint = GetNearestPointInCapsule2D(position, capsuleStart, capsuleEnd, capsuleRadius);
    float2 toNearest = position - nearestPoint;
    float distSquared = dot(toNearest, toNearest);

    if (distSquared <= (capsuleRadius * capsuleRadius))
    {
        float distance = sqrt(distSquared);
        float overlap = capsuleRadius - distance;

        if (distance > 0.0f)
        {
            float2 collisionNormal = toNearest / distance;
            position += overlap * collisionNormal;
            velocity -= 1.0f * dot(velocity, collisionNormal) * collisionNormal;
        }
    }
}


//--------------------------------------------------------------------------------------
// Optimized Grid + Sort Algorithm
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void ForceCS_Grid(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int P_ID = DTid.x; // Particle ID to operate on

    float2 P_position = ParticlesRO[P_ID].position;
    float2 P_velocity = ParticlesRO[P_ID].velocity;
    float P_density = ParticlesDensityRO[P_ID].density;
    float P_pressure = CalculatePressure(P_density);
    float type = ParticlesRO[P_ID].type;
    const float h_sq = g_fSmoothlen * g_fSmoothlen;

    float2 acceleration = float2(0, 0);

    // Calculate the acceleration based on neighbors from the 8 adjacent cells + current cell
    int2 G_XY = (int2)GridCalculateCell(P_position);
    for (int Y = max(G_XY.y - 1, 0); Y <= min(G_XY.y + 1, 255); Y++)
    {
        for (int X = max(G_XY.x - 1, 0); X <= min(G_XY.x + 1, 255); X++)
        {
            unsigned int G_CELL = GridConstructKey(uint2(X, Y));
            uint2 G_START_END = GridIndicesRO[G_CELL];
            for (unsigned int N_ID = G_START_END.x; N_ID < G_START_END.y; N_ID++)
            {
                float2 N_position = ParticlesRO[N_ID].position;

                float2 diff = N_position - P_position;
                float r_sq = dot(diff, diff);
                if (r_sq < h_sq && P_ID != N_ID)
                {
                    float2 N_velocity = ParticlesRO[N_ID].velocity;
                    float N_density = ParticlesDensityRO[N_ID].density;
                    float N_pressure = CalculatePressure(N_density);
                    float r = sqrt(r_sq);

                    // Pressure Term
                    acceleration += CalculateGradPressure(r, P_pressure, N_pressure, N_density, diff);

                    // Viscosity Term
                    acceleration += CalculateLapVelocity(r, P_velocity, N_velocity, N_density);                    
                }
            }
        }
    }
    ParticlesForcesRW[P_ID].acceleration = acceleration / P_density;
}

//--------------------------------------------------------------------------------------
// Integration
//--------------------------------------------------------------------------------------

[numthreads(SIMULATION_BLOCK_SIZE, 1, 1)]
void IntegrateCS(uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex)
{
    const unsigned int P_ID = DTid.x;

    float2 position = ParticlesRO[P_ID].position;
    float2 velocity = ParticlesRO[P_ID].velocity;
    float2 acceleration = ParticlesForcesRO[P_ID].acceleration;
    float type = ParticlesRO[P_ID].type;

    // Apply the forces from the map walls
    acceleration += ComputeWallCollisionForce(position);

    //Tumbling config collisions
    ProcessOBBCollision(position, velocity, tumblingConfigOBB1WorldToLocalMatrix, tumblingConfigOBB1LocalToWorldMatrix, tumblingConfigOBB1Center, tumblingConfigOBB1HalfDimensions * 2);
    ProcessOBBCollision(position, velocity, tumblingConfigOBB2WorldToLocalMatrix, tumblingConfigOBB2LocalToWorldMatrix, tumblingConfigOBB2Center, tumblingConfigOBB2HalfDimensions * 2);
    ProcessOBBCollision(position, velocity, tumblingConfigOBB3WorldToLocalMatrix, tumblingConfigOBB3LocalToWorldMatrix, tumblingConfigOBB3Center, tumblingConfigOBB3HalfDimensions * 2);
    ProcessOBBCollision(position, velocity, tumblingConfigOBB4WorldToLocalMatrix, tumblingConfigOBB4LocalToWorldMatrix, tumblingConfigOBB4Center, tumblingConfigOBB4HalfDimensions * 2);
    ProcessOBBCollision(position, velocity, tumblingConfigOBB5WorldToLocalMatrix, tumblingConfigOBB5LocalToWorldMatrix, tumblingConfigOBB5Center, tumblingConfigOBB5HalfDimensions * 2);
    ProcessCapsuleCollision(position, velocity, tumblingConfigCapsule1Start, tumblingConfigCapsule1End, tumblingConfigCapsule1Radius);
    ProcessCapsuleCollision(position, velocity, tumblingConfigCapsule2Start, tumblingConfigCapsule2End, tumblingConfigCapsule2Radius);
    ProcessAABBCollision( position,  velocity, tumblingConfigHorizontalMovingAABB1Center, tumblingConfigHorizontalMovingAABB1Size);
    ProcessDiscCollision(position, velocity, tumblingDiscCenter, tumblingDiscRadius);

    //Windmill config  collisions
    ProcessDiscCollision(position, velocity, windmillDisc1Center, windmillDisc1Radius);
    ProcessOBBCollision(position, velocity, windmillConfigOBB1WorldToLocalMatrix, windmillConfigOBB1LocalToWorldMatrix, windmillConfigOBB1Center, windmillConfigOBB1HalfDimensions * 2);
    ProcessOBBCollision(position, velocity, windmillConfigOBB2WorldToLocalMatrix, windmillConfigOBB2LocalToWorldMatrix, windmillConfigOBB2Center, windmillConfigOBB2HalfDimensions * 2);
    ProcessCapsuleCollision(position, velocity, windmillConfigBottomBorderCapsuleStart, windmillConfigBottomBorderCapsuleEnd, windmillConfigBottomBorderCapsuleRadius);
    ProcessCapsuleCollision(position, velocity, windmillConfigLeftBorderCapsuleStart, windmillConfigLeftBorderCapsuleEnd, windmillConfigLeftBorderCapsuleRadius);
    ProcessCapsuleCollision(position, velocity, windmillConfigRightBorderCapsuleStart, windmillConfigRightBorderCapsuleEnd, windmillConfigRightBorderCapsuleRadius);
    ProcessCapsuleCollision(position, velocity, windmillConfigTopBorderCapsule1Start, windmillConfigTopBorderCapsule1End, windmillConfigTopBorderCapsule1Radius);
    ProcessCapsuleCollision(position, velocity, windmillConfigTopBorderCapsule2Start, windmillConfigTopBorderCapsule2End, windmillConfigTopBorderCapsule2Radius);
    ProcessCapsuleCollision(position, velocity, windmillConfigTopBorderCapsule3Start, windmillConfigTopBorderCapsule3End, windmillConfigTopBorderCapsule3Radius);

    //Hourglass Config collisions
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule1Start, hourglassConfigCapsule1End, hourglassConfigCapsule1Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule2Start, hourglassConfigCapsule2End, hourglassConfigCapsule2Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule3Start, hourglassConfigCapsule3End, hourglassConfigCapsule3Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule4Start, hourglassConfigCapsule4End, hourglassConfigCapsule4Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule5Start, hourglassConfigCapsule5End, hourglassConfigCapsule5Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigCapsule6Start, hourglassConfigCapsule6End, hourglassConfigCapsule6Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigMiddleCapsule1Start, hourglassConfigMiddleCapsule1End, hourglassConfigMiddleCapsule1Radius);
    ProcessCapsuleCollision(position, velocity, hourglassConfigMiddleCapsule2Start, hourglassConfigMiddleCapsule2End, hourglassConfigMiddleCapsule2Radius);

    //Water Elevator Config Collisions
    ProcessAABBCollision(position, velocity, waterElevatorConfigVerticalMovingAABBCenter, waterElevatorConfigVerticalMovingAABBSize);
    ProcessDiscCollision(position, velocity, waterElevatorDisc1Center, waterElevatorDisc1Radius);
    ProcessDiscCollision(position, velocity, waterElevatorDisc2Center, waterElevatorDisc2Radius);
    ProcessDiscCollision(position, velocity, waterElevatorDisc3Center, waterElevatorDisc3Radius);
 
    //Apply gravity and mouse forces(LMB and RMB)
    acceleration += ApplyExternalForces(position, velocity, acceleration);

    // Integrate
    velocity += g_fTimeStep * acceleration;
    position += g_fTimeStep * velocity;

    // Update
    ParticlesRW[P_ID].position = position;
    ParticlesRW[P_ID].velocity = velocity;
}
//-----------------------------------------------------------------------------------