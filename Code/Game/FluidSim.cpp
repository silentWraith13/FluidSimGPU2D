#include "Game/FluidSim.hpp"
#include "Game/App.hpp"
#include "Engine/Renderer/Camera.hpp"
#include "Engine/Renderer/UnorderedAccessBuffer.hpp"
#include "Engine/Core/Clock.hpp"
#include "Engine/Renderer/Shader.hpp"
#include "Engine/Renderer/ConstantBuffer.hpp"
#include "Engine/Renderer/VertexBuffer.hpp"
#include "Engine/Renderer/BitmapFont.hpp"
#include "Engine/Math/LineSegment2.hpp"
#include "Engine/Core/Time.hpp"
#include "ThirdParty/ImGUI/imgui.h"
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

//--------------------------------------------------------------------------------------------------------------------------------------------------------
FluidSim::FluidSim()
{
	m_rng = RandomNumberGenerator((unsigned int)std::chrono::system_clock::now().time_since_epoch().count());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
FluidSim::~FluidSim()
{
	delete	m_forceSimpleCS;
	m_forceSimpleCS = nullptr;

	delete m_densitySimpleCS;
	m_densitySimpleCS = nullptr;

	delete m_integrateCS;
	m_integrateCS = nullptr;

	delete m_particleUAV;
	m_particleUAV = nullptr;

	delete m_particleDensityUAV;
	m_particleDensityUAV = nullptr;

	delete m_particleForcesUAV;
	m_particleForcesUAV = nullptr;

	delete m_sortedParticlesUAV;
	m_sortedParticlesUAV = nullptr;

	delete m_simulationConstantsCBO;
	m_simulationConstantsCBO = nullptr;

	delete m_particleSizeCBO;
	m_particleSizeCBO = nullptr;

	delete m_gridUAV;
	m_gridUAV = nullptr;

	delete m_gridPingPongUAV;
	m_gridPingPongUAV = nullptr;
	delete m_gridIndicesUAV;
	m_gridIndicesUAV = nullptr;

	delete m_densityGridCS;
	m_densityGridCS = nullptr;

	delete m_forcesGridCS;
	m_forcesGridCS = nullptr;

	delete m_buildGridCS;
	m_buildGridCS = nullptr;

	delete m_ClearGridIndicesGridCS;
	m_ClearGridIndicesGridCS = nullptr;

	delete m_buildGridIndicesCS;
	m_buildGridIndicesCS = nullptr;

	delete m_bitonicSortCS;
	m_bitonicSortCS = nullptr;

	delete m_matrixTransposeCS;
	m_matrixTransposeCS = nullptr;

	delete m_rearrangeParticlesCS;
	m_rearrangeParticlesCS = nullptr;

	delete m_sortCBO;
	m_sortCBO = nullptr;

	delete m_postProcessRenderTexture;
	m_postProcessRenderTexture = nullptr;

	delete m_screenSizeCBO;
	m_screenSizeCBO = nullptr;

	for (int i = 0; i < m_tumblingConfigShapes.size(); i++)
	{
		if (m_tumblingConfigShapes[i])
		{
			delete m_tumblingConfigShapes[i];
			m_tumblingConfigShapes[i] = nullptr;
		}
	}
	m_tumblingConfigShapes.clear();
	m_tumblingConfigShapesVerts.clear();
	m_windmillConfigShapesVerts.clear();

	delete m_particleDensityRenderTexture;
	m_particleDensityRenderTexture = nullptr;

	delete m_transparentFluidRenderTexture;
	m_transparentFluidRenderTexture = nullptr;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::Startup(unsigned int numParticles, InitialParticleConfigurationGPU initialParticleConfig)
{
	InitializeWorldCamera();

	if (initialParticleConfig == DAM)
	{
		InitializeParticlesInADam( numParticles);
	}

	if (initialParticleConfig == RANDOM_SHAPES)
	{
		InitializeParticlesInARandomShapeConfig(numParticles);
	}
	
	if (initialParticleConfig == WINDMILL)
	{
		InitializeParticlesInAWindmillConfig(numParticles);
	}

	if (initialParticleConfig == HOURGLASS)
	{
		InitializeParticlesInAnHourglassConfig(numParticles);
	}

	if (initialParticleConfig == TRANSPARENT_FLUID)
	{
		InitializeBackgroundTextureForTransparentFluidConfig();
		InitializeParticlesForTransparentFluidConfig(numParticles);
	}

	if (initialParticleConfig == WATER_ELEVATOR)
	{
		InitializeParticlesForWaterElevatorConfig(numParticles);
	}

	InitializeParticleTexturesAndShaders();
	CreateComputeShaders();
	CreateParticleSizeConstantBuffer();
	CreateSimulationConstantBuffer();
	CreateSortingConstantBuffer();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::Update(float deltaSeconds)
{
	SimulateFluid(deltaSeconds);

	HandleKeyboardInputs();
	
	if (g_theGame->m_initialParticleConfigGPU == WINDMILL)
	{
		UpdateRotationOfWindmillMesh(deltaSeconds);
	}

	if (g_theGame->m_initialParticleConfigGPU == HOURGLASS)
	{
		UpdateHourglassMeshConfig(deltaSeconds);
	}

	if (g_theGame->m_initialParticleConfigGPU == WATER_ELEVATOR)
	{
		UpdateShapesForWaterElevatorConfig(deltaSeconds);
	}

	if (g_theGame->m_initialParticleConfigGPU == TRANSPARENT_FLUID)
	{
		UpdateShapesForWaterElevatorConfig(deltaSeconds);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::Render()
{
	g_theRenderer->ClearScreen(m_clearColor);
	g_theRenderer->BeginCamera(m_worldCamera);


	if (m_particleRenderType == CARTOON)
	{
		RenderCartoonFluid();
	}

	else if (g_theGame->m_initialParticleConfigGPU == TRANSPARENT_FLUID)
	{
		RenderTransparentFluidForTheTransparentConfig();
		RenderShapesForWaterElevatorConfig();
	}

	else
	{
		RenderFluidParticles();
		RenderShapesBasedOnConfig();
	}

	ImGuiRenderFunctions();
	g_theRenderer->EndCamera(m_worldCamera);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::Shutdown()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::EndFrame()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeWorldCamera()
{
	m_worldCamera.m_mode = Camera::eMode_Orthographic;
	m_worldCamera.SetOrthographicView(Vec2(0.f, 0.f), Vec2(MAP_WIDTH, MAP_HEIGHT));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::HandleKeyboardInputs()
{
	if (g_theInput->WasKeyJustReleased(KEYCODE_LEFT_MOUSE))
	{
		m_isHoldingLeftClick = 0;
	}

	if (g_theInput->IsKeyDown(KEYCODE_LEFT_MOUSE))
	{
		m_isHoldingLeftClick = 1;
	}

	if (g_theInput->IsKeyDown(KEYCODE_RIGHT_MOUSE))
	{
		m_isHoldingRightClick = 1;
	}
	if (g_theInput->WasKeyJustReleased(KEYCODE_RIGHT_MOUSE))
	{
		m_isHoldingRightClick = 0;
	}

	if (g_theInput->WasKeyJustPressed('1'))
	{
		m_showTumblingConfigMovingAABB = !m_showTumblingConfigMovingAABB;
	}

	if (g_theInput->WasKeyJustPressed('W'))
	{
		m_shouldRotateHourglassMesh = !m_shouldRotateHourglassMesh;
	}

	if (g_theInput->WasKeyJustPressed('R'))
	{
		m_renderTheTopCapsuleForHourglassMesh = !m_renderTheTopCapsuleForHourglassMesh;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticleTexturesAndShaders()
{
	AddVertsForAABB2D(m_screenVerts, m_screenBounds, Rgba8::WHITE, AABB2(0.f, 1.f, 1.f, 0.f));
	m_particleTexture = g_theRenderer->CreateOrGetTextureFromFile("Data/Images/Water_Additive.png");
	m_postProcessRenderTexture = g_theRenderer->CreateRenderTexture(g_theWindow->GetClientDimensions().x, g_theWindow->GetClientDimensions().y);
	m_particleDensityRenderTexture = g_theRenderer->CreateRenderTexture(g_theWindow->GetClientDimensions().x, g_theWindow->GetClientDimensions().y);
	m_postProcessShader = g_theRenderer->CreateShaderOrGetFromFile("Data/Shaders/Rendering/PostProcessBlur");
	m_transparentFluidShader = g_theRenderer->CreateShaderOrGetFromFile("Data/Shaders/Rendering/TransparentFluidShader");
	m_noiseTexture = g_theRenderer->CreateOrGetTextureFromFile("Data/Images/Noise2D.png");
	InitializeParticleRenderingShader();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesInADam(unsigned int numParticles)
{
	m_numParticles = numParticles;
	int startingWidth = (int)sqrt(numParticles);
	std::vector<Particle> particles(numParticles);

	for (int i = 0; i < (int)numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		particles[i].m_position = Vec2(0.0f + m_initialParticleSpacing * (float)x , m_initialParticleSpacing * (float)y + 0.f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_renderType = 1;
	m_particleRenderType = BLUE_PALLETE;
	m_particleSize = 0.003f;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesInARandomShapeConfig(unsigned int numParticles)
{
	m_numParticles = numParticles;
	int startingWidth = (int)sqrt(numParticles);
	std::vector<Particle> particles(numParticles);

	for (int i = 0; i < (int)numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		particles[i].m_position = Vec2(0.8f + m_initialParticleSpacing * (float)x, m_initialParticleSpacing * (float)y + 2.f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_renderType = 1;
	m_particleRenderType = BLUE_PALLETE;
	m_particleSize = 0.003f;
	InitializeShapesForRandomShapesConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesInAWindmillConfig(unsigned int numParticles)
{
	m_numParticles = numParticles;
	int startingWidth = (int)sqrt(numParticles);
	std::vector<Particle> particles(numParticles);

	for (int i = 0; i < (int)numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		particles[i].m_position = Vec2(0.8f + m_initialParticleSpacing * (float)x, m_initialParticleSpacing * (float)y + 2.f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_renderType = 1;
	m_particleRenderType = BLUE_PALLETE;
	m_particleSize = 0.003f;
	InitializeShapesForWindmillConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesInAnHourglassConfig(unsigned int numParticles)
{
	m_numParticles = numParticles;
	int startingWidth = (int)sqrt(numParticles);
	std::vector<Particle> particles(numParticles);

	for (int i = 0; i < (int)numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		particles[i].m_position = Vec2(0.8f + m_initialParticleSpacing * (float)x, m_initialParticleSpacing * (float)y + 2.f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_renderType = 1;
	m_particleRenderType = BLUE_PALLETE;
	m_particleSize = 0.003f;
	InitializeShapesForHourglassConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderFluidParticles()
{
	g_theRenderer->SetRasterizerModes(RasterizerMode::SOLID_CULL_BACK_FRONT_CW);
	
	BindCBForParticleSize();
	
	g_theRenderer->BindTexture(m_particleTexture, 0);
	g_theRenderer->BindTexture(m_noiseTexture, 1);
	g_theRenderer->BindShader(m_renderingShader);
	g_theRenderer->VSSetShaderResource(0, 1, m_particleUAV);
	g_theRenderer->VSSetShaderResource(1, 1, m_particleDensityUAV);
	g_theRenderer->IASetVertexBuffers(nullptr);
	g_theRenderer->IASetPrimitiveTopology();
	
	g_theRenderer->Draw(m_numParticles);
	
	g_theRenderer->VSSetShaderResource(0, 1, m_nullUAV);
	g_theRenderer->VSSetShaderResource(1, 1, m_nullUAV);
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->BindTexture(nullptr, 0);
	g_theRenderer->BindTexture(nullptr, 1);

	g_theRenderer->SetRasterizerModes(RasterizerMode::SOLID_CULL_BACK);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticleRenderingShader()
{
	std::string shaderPath = "Data/Shaders/Rendering/FluidRender";
	m_renderingShader = g_theRenderer->CreateShaderOrGetFromFile(shaderPath.c_str());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::CreateParticleSizeConstantBuffer()
{
	m_particleSizeCBO = g_theRenderer->CreateConstantBuffer(sizeof(ParticleSize));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::BindCBForParticleSize()
{
	ParticleSize constant = {};
	constant.particleSize = m_particleSize;
	constant.debugType = m_renderType;
	g_theRenderer->CopyCPUtoGPU(&constant, sizeof(ParticleSize), m_particleSizeCBO);
	g_theRenderer->BindConstantBufferWithGS(k_particleSizeConstantsSlot, m_particleSizeCBO);

	ScreenSizeConstants constants = {};
	constants.screenHeight = (float)g_theWindow->GetClientDimensions().y;
	constants.screenWidth = (float)g_theWindow->GetClientDimensions().x;
	constants.time = (float)GetCurrentTimeSeconds();
	g_theRenderer->CopyCPUtoGPU(&constants, sizeof(ScreenSizeConstants), m_screenSizeCBO);
	g_theRenderer->BindConstantBufferWithGS(4, m_screenSizeCBO);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::CreateSimulationConstantBuffer()
{
	m_screenSizeCBO = g_theRenderer->CreateConstantBuffer(sizeof(ScreenSizeConstants));
	m_simulationConstantsCBO = g_theRenderer->CreateConstantBuffer(sizeof(SimulationConstants));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::CreateSortingConstantBuffer()
{
	m_sortCBO = g_theRenderer->CreateConstantBuffer(sizeof(Sort));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::BindConstantBufferForSimulationConstants(float deltaSeconds)
{
	SimulationConstants constants = {};

	constants.numFluidParticles = m_numParticles;
	constants.timeStep =  std::min(m_maxAllowableTimeStep, deltaSeconds);
	constants.smoothingLength = m_smoothingLength;
	constants.pressureStiffness = m_pressureStiffness;
	constants.restDensity = m_restDensity;
	constants.densityCoefficient = m_particleMass * static_cast < float>(315.0f / (64.0f * M_PI * powf(m_smoothingLength, 9)));
	constants.gradPressureCoefficient = m_particleMass * static_cast<float>(-45.0f / (M_PI * powf(m_smoothingLength, 6)));
	constants.laplcaianViscosityCoefficient = m_particleMass * m_viscosity * static_cast < float>(45.0f / (M_PI * powf(m_smoothingLength, 6)));
	constants.surfaceTensionCoefficient = m_particleMass * static_cast<float>(1.0f / (M_PI * powf(m_smoothingLength, 6)));
	constants.wallStiffness = m_wallStiffness;

	constants.vGravity = m_gravityDirection;
	constants.vGridDim.x = 1.0f / m_smoothingLength;
	constants.vGridDim.y = 1.0f / m_smoothingLength;
	constants.vGridDim.z = 0;
	constants.vGridDim.w = 0;
	
	constants.leftMapBounds = Vec3(1, 0, 0);
	constants.bottomMapBounds = Vec3(0, 1, 0);
	constants.rightMapBounds = Vec3(-1, 0, MAP_WIDTH);

	if (g_theGame->m_initialParticleConfigGPU == RANDOM_SHAPES)
	{
		constants.topMapBounds = Vec3(0, -1, MAP_HEIGHT + 15.f);
		UpdateSimulationConstantBufferFromRandomShapesConfig(constants);
	}

	else if(g_theGame->m_initialParticleConfigGPU == WINDMILL)
	{
		constants.topMapBounds = Vec3(0, -1, MAP_HEIGHT + 5.f);
		constants.bottomMapBounds = Vec3(0, 1, MAP_HEIGHT + 15.f);
		UpdateSimulationConstantBufferForWindmillConfig(constants);
	}
	
	else if (g_theGame->m_initialParticleConfigGPU == HOURGLASS)
	{
		constants.topMapBounds = Vec3(0, -1, MAP_HEIGHT + 15.f);
		UpdateSimulationConstantBufferForHourglassConfig(constants);
	}

	else if (g_theGame->m_initialParticleConfigGPU == TRANSPARENT_FLUID)
	{
		constants.bottomMapBounds = Vec3(0, 1, MAP_HEIGHT + 12.f);
		UpdateSimulationConstantBufferForWaterElevatorConfig(constants);
	}

	else if (g_theGame->m_initialParticleConfigGPU == WATER_ELEVATOR)
	{
		constants.bottomMapBounds = Vec3(0, 1, MAP_HEIGHT + 12.f);
		UpdateSimulationConstantBufferForWaterElevatorConfig(constants);
	}

	else
	{
		constants.topMapBounds = Vec3(0, -1, MAP_HEIGHT);
	}

	constants.mousePos = GetMouseCursorPos();
	constants.isHoldingLeftClick = m_isHoldingLeftClick;
	constants.isHoldingRightClick = m_isHoldingRightClick;

	g_theRenderer->CopyCPUtoGPU(&constants, sizeof(SimulationConstants), m_simulationConstantsCBO);
	g_theRenderer->CSSetConstantBuffers(0, 1, m_simulationConstantsCBO);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::SimulateFluid(float deltaSeconds)
{
	BindConstantBufferForSimulationConstants(deltaSeconds);
	
	switch (m_simMode)
	{
	case SIMPLE:
	{
		SimulateFluidSimple();
	}
	break;

	case GRID:
	{
		SimulateFluidGrid();
	}
	break;
	}

	UnsetAllResources();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::SimulateFluidSimple()
{
	// Setup
	g_theRenderer->CSSetShaderResources(0, 1, m_particleUAV);
	
	// Density
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleDensityUAV);
	g_theRenderer->CSSetShader(m_densitySimpleCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	// Force
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleForcesUAV);
	g_theRenderer->CSSetShaderResources(1, 1, m_particleDensityUAV);
	g_theRenderer->CSSetShader(m_forceSimpleCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	// Integrate
	g_theRenderer->CopyResource(m_sortedParticlesUAV, m_particleUAV);
	g_theRenderer->CSSetShaderResources(0, 1, m_sortedParticlesUAV);
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleUAV);
	g_theRenderer->CSSetShaderResources(2, 1, m_particleForcesUAV);
	g_theRenderer->CSSetShader(m_integrateCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::SimulateFluidGrid()
{
	//Setup
	g_theRenderer->CSSetConstantBuffers(0, 1, m_simulationConstantsCBO);
	g_theRenderer->CSSetUnorderedAccessViews(0, m_gridUAV);
	g_theRenderer->CSSetShaderResources(0, 1, m_particleUAV);

	//Build Grid
	g_theRenderer->CSSetShader(m_buildGridCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	//Sort grid
	GPUSort(m_gridUAV, m_gridPingPongUAV);

	//Setup
	g_theRenderer->CSSetConstantBuffers(0, 1, m_simulationConstantsCBO);
	g_theRenderer->CSSetUnorderedAccessViews(0, m_gridIndicesUAV);
	g_theRenderer->CSSetShaderResources(3, 1, m_gridUAV);

	//Build Grid Indices
	g_theRenderer->CSSetShader(m_ClearGridIndicesGridCS);
	g_theRenderer->DispatchCS(NUM_GRID_INDICES / SIMULATION_BLOCK_SIZE, 1, 1);
	g_theRenderer->CSSetShader(m_buildGridIndicesCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	//Setup
	g_theRenderer->CSSetUnorderedAccessViews(0, m_sortedParticlesUAV);
	g_theRenderer->CSSetShaderResources(0, 1, m_particleUAV);
	g_theRenderer->CSSetShaderResources(3, 1, m_gridUAV);

	//Rearrange
	g_theRenderer->CSSetShader(m_rearrangeParticlesCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	//Setup
	g_theRenderer->CSSetUnorderedAccessViews(0, m_nullUAV);
	g_theRenderer->CSSetShaderResources(0, 1, m_nullUAV);
	g_theRenderer->CSSetShaderResources(0, 1, m_sortedParticlesUAV);
	g_theRenderer->CSSetShaderResources(3, 1, m_gridUAV);
	g_theRenderer->CSSetShaderResources(4, 1, m_gridIndicesUAV);
	
	//Density
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleDensityUAV);
	g_theRenderer->CSSetShader(m_densityGridCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	//Force
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleForcesUAV);
	g_theRenderer->CSSetShaderResources(1, 1, m_particleDensityUAV);
	g_theRenderer->CSSetShader(m_forcesGridCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);

	//Integrate
	g_theRenderer->CSSetUnorderedAccessViews(0, m_particleUAV);
	g_theRenderer->CSSetShaderResources(2, 1, m_particleForcesUAV);
	g_theRenderer->CSSetShader(m_integrateCS);
	g_theRenderer->DispatchCS(m_numParticles / SIMULATION_BLOCK_SIZE, 1, 1);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::CreateComputeShaders()
{
	std::string shaderPath = "Data/Shaders/FluidMathCS/FluidCS11";
	m_integrateCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "IntegrateCS");
	m_densitySimpleCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "DensityCS_Simple");
	m_forceSimpleCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "ForceCS_Simple");
	m_densityGridCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "DensityCS_Grid");
	m_forcesGridCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "ForceCS_Grid");
	m_ClearGridIndicesGridCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "ClearGridIndicesCS");
	m_buildGridCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "BuildGridCS");
	m_buildGridIndicesCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "BuildGridIndicesCS");
	m_rearrangeParticlesCS = g_theRenderer->CreateComputeShaderOnly(shaderPath.c_str(), "RearrangeParticlesCS");
	
	std::string shaderPath2 = "Data/Shaders/FluidMathCS/ComputeShaderSort11";
	m_bitonicSortCS = g_theRenderer->CreateComputeShaderOnly(shaderPath2.c_str(), "BitonicSort");
	m_matrixTransposeCS = g_theRenderer->CreateComputeShaderOnly(shaderPath2.c_str(), "MatrixTranspose");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiGravityButtons()
{
	if (ImGui::Button("Open Gravity Controls"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Gravity Control");
	}

	if (ImGui::BeginPopup("Gravity Control"))
	{
		int gravityDirection = 0;
		if (m_gravityDirection == GRAVITY_UP) gravityDirection = 0;
		else if (m_gravityDirection == GRAVITY_DOWN) gravityDirection = 1;
		else if (m_gravityDirection == GRAVITY_LEFT) gravityDirection = 2;
		else if (m_gravityDirection == GRAVITY_RIGHT) gravityDirection = 3;
		else if (m_gravityDirection == GRAVITY_DOWN_FAST) gravityDirection = 4;
		const char* items[] = { "Up", "Down", "Left", "Right", "Down Fast"};

		ImGui::Combo("Direction", &gravityDirection, items, IM_ARRAYSIZE(items));
		switch (gravityDirection)
		{
		case 0: m_gravityDirection = GRAVITY_UP; break;
		case 1: m_gravityDirection = GRAVITY_DOWN; break;
		case 2: m_gravityDirection = GRAVITY_LEFT; break;
		case 3: m_gravityDirection = GRAVITY_RIGHT; break;
		case 4: m_gravityDirection = GRAVITY_DOWN_FAST; break;
		}
		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiFpsText()
{
 	float fps = ImGui::GetIO().Framerate;
 
 	int windowWidth = g_theWindow->GetClientDimensions().x;
 	ImVec2 windowPos = ImVec2((float)windowWidth - 80.f, 8.f);
 	ImVec2 windowPosPivot = ImVec2(0.5f, 0.5f); 
 
 	ImGui::SetNextWindowPos(windowPos, ImGuiCond_Always, windowPosPivot);
 	ImGui::SetNextWindowBgAlpha(0.0f);
 
 	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoMove |
 		ImGuiWindowFlags_NoResize |
 		ImGuiWindowFlags_NoSavedSettings |
 		ImGuiWindowFlags_NoTitleBar |
 		ImGuiWindowFlags_NoCollapse |
 		ImGuiWindowFlags_NoScrollbar |
 		ImGuiWindowFlags_NoBackground;
 
 	if (ImGui::Begin("Performance", NULL, windowFlags))
 	{
 		ImGui::SetWindowFontScale(1.5f);
 		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
 		ImGui::Text("FPS: %.1f", fps);
 		ImGui::PopStyleColor();
 	}
 	ImGui::End();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiSimulationMode()
{
	if (ImGui::Button("Open Simulation Mode"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Simulation Mode");
	}

	if (ImGui::BeginPopup("Simulation Mode"))
	{
		static int selectedSimMode = m_simMode; 
		const char* items[] = { "Simple", "Grid" };

		if (ImGui::Combo("Simulation Mode", &selectedSimMode, items, IM_ARRAYSIZE(items)))
		{
			m_simMode = static_cast<SimulationMode>(selectedSimMode);
		}

		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiRenderTypeForParticles()
{
	if (ImGui::Button("Open Particle Render Mode"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Particle Render Mode");
	}

	if (ImGui::BeginPopup("Particle Render Mode"))
	{
		int selectedmode = m_particleRenderType;

		if (g_theGame->m_initialParticleConfigGPU == TRANSPARENT_FLUID)
		{
			const char* transparentFluidItems[] = { "Density Debug", "Transparent Fluid" };
			int transparentItemsCount = IM_ARRAYSIZE(transparentFluidItems);

			// Adjust index for transparent fluid items
			int transparentModeIndex = m_particleRenderType == DENSITY_DEBUG ? 0 : 1;

			if (ImGui::Combo("Particle Render Mode", &transparentModeIndex, transparentFluidItems, transparentItemsCount))
			{
				m_particleRenderType = (ParticleRenderType)((transparentModeIndex == 0) ? DENSITY_DEBUG : TRANSPARENT);
			}
		}
		else
		{
			const char* allItems[] = { "Density Debug", "Water Pallet", "Cartoon"};
			int allItemsCount = IM_ARRAYSIZE(allItems);

			if (ImGui::Combo("Particle Render Mode", &selectedmode, allItems, allItemsCount))
			{
				m_particleRenderType = static_cast<ParticleRenderType>(selectedmode);
			}
		}

		switch (m_particleRenderType)
		{
		case DENSITY_DEBUG:
			m_renderType = 0;
			m_particleSize = 0.003f;
			break;
		case BLUE_PALLETE:
			m_renderType = 1;
			m_particleSize = 0.003f;
			break;
		case CARTOON:
			m_renderType = 2;
			m_particleSize = 0.0271f;
			break;
		case TRANSPARENT:
			m_renderType = 3;
			m_particleSize = 0.03f;
			break;
		}

		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiViscosity()
{
	ImGui::SliderFloat("Viscosity", &m_viscosity, 0.01f, 1.0f, "%.3f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiRestDensity()
{
	ImGui::SliderFloat("RestDensity", &m_restDensity, 1000.f, 3000.f, "%.1f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiPressureStiffness()
{
	ImGui::SliderFloat("Pressure Stiffness", &m_pressureStiffness, 100.0f, 300.0f, "%.1f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiWallStiffness()
{
	ImGui::SliderFloat("Wall Stiffness", &m_wallStiffness, 100.0f, 5000.0f, "%.1f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiSmoothingLength()
{
	ImGui::SliderFloat("Smoothing Length", &m_smoothingLength, 0.01f, 0.05f, "%.4f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiTimestep()
{
	ImGui::SliderFloat("Time Step", &m_maxAllowableTimeStep, 0.001f, 0.1f, "%.4f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiParticleMass()
{
	ImGui::SliderFloat("Particle Mass", &m_particleMass, 0.0001f, 0.0009f, "%.4f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiParticleRenderSize()
{
	ImGui::SliderFloat("Particle Size", &m_particleSize, 0.003f, 0.05f, "%.4f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiAddGeometry()
{
	ImGui::Text("Press 2 to spawn falling disc");
	ImGui::Text("Press 3 to spawn falling AABB");
	ImGui::Text("Press 4 to spawn horizontally moving AABB");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiChangeBackgroundColor()
{
	if (ImGui::Button("Background Color"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Background Color");
	}

	if (ImGui::BeginPopup("Background Color"))
	{
		static int selectedmode = 0;
		const char* items[] = { "WHITE", "SEA_BLUE", "BLACK" };
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiMillisecondText()
{
	ImGui::Text("m/s - %.2f", m_lastSimulationTimeMs);
	ImGui::Text("m/s - %.2f", m_lastRenderTimeMs);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UnsetAllResources()
{
	g_theRenderer->CSSetUnorderedAccessViews(0, nullptr);
	g_theRenderer->CSSetShaderResources(0, 1, nullptr);
	g_theRenderer->CSSetShaderResources(1, 1, nullptr);
	g_theRenderer->CSSetShaderResources(2, 1, nullptr);
 	g_theRenderer->CSSetShaderResources(3, 1, nullptr);
 	g_theRenderer->CSSetShaderResources(4, 1, nullptr);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::GPUSort(UnorderedAccessBuffer* inUAV, UnorderedAccessBuffer* tempUAV)
{
	g_theRenderer->CSSetConstantBuffers(0, 1, m_sortCBO);

	const unsigned int NUM_ELEMENTS = m_numParticles;
	const unsigned int MATRIX_WIDTH = BITONIC_BLOCK_SIZE;
	const unsigned int MATRIX_HEIGHT = NUM_ELEMENTS / BITONIC_BLOCK_SIZE;

	//Sort the data
	//First sort the rows for the levels <= to the block size
	for (unsigned int level = 2; level <= BITONIC_BLOCK_SIZE; level <<= 1)
	{
		Sort constants = { level, level, MATRIX_HEIGHT, MATRIX_WIDTH };
		g_theRenderer->CopyCPUtoGPU(&constants, sizeof(Sort), m_sortCBO);

		//Sort the row data
		g_theRenderer->CSSetUnorderedAccessViews(0, inUAV);
		g_theRenderer->CSSetShader(m_bitonicSortCS);
		g_theRenderer->DispatchCS(NUM_ELEMENTS / BITONIC_BLOCK_SIZE, 1, 1);
	}

	// Then sort the rows and columns for the levels > than the block size
    // Transpose. Sort the Columns. Transpose. Sort the Rows.
	for (unsigned int level = (BITONIC_BLOCK_SIZE << 1); level <= NUM_ELEMENTS; level <<= 1)
	{
		Sort constants1 = { (level / BITONIC_BLOCK_SIZE), (level & ~NUM_ELEMENTS) / BITONIC_BLOCK_SIZE, MATRIX_WIDTH, MATRIX_HEIGHT };
		g_theRenderer->CopyCPUtoGPU(&constants1, sizeof(Sort), m_sortCBO);

		//Transpose the data from buffer 1 into buffer2
		UnorderedAccessBuffer* nullUAV = nullptr;
		g_theRenderer->CSSetShaderResources(0, 1, nullUAV);
		g_theRenderer->CSSetUnorderedAccessViews(0, tempUAV);
		g_theRenderer->CSSetShaderResources(0, 1, inUAV);
		g_theRenderer->CSSetShader(m_matrixTransposeCS);
		g_theRenderer->DispatchCS(MATRIX_WIDTH / TRANSPOSE_BLOCK_SIZE, MATRIX_HEIGHT / TRANSPOSE_BLOCK_SIZE, 1);


		//Sort the transposed column data
		g_theRenderer->CSSetShader(m_bitonicSortCS);
		g_theRenderer->DispatchCS(NUM_ELEMENTS / BITONIC_BLOCK_SIZE, 1, 1);

		Sort constants2 = { BITONIC_BLOCK_SIZE, level, MATRIX_HEIGHT, MATRIX_WIDTH };
		g_theRenderer->CopyCPUtoGPU(&constants2, sizeof(Sort), m_sortCBO);

		//Transpose the data from buffer 2 back to buffer1
		g_theRenderer->CSSetShaderResources(0, 1, nullUAV);
		g_theRenderer->CSSetUnorderedAccessViews(0, inUAV);
		g_theRenderer->CSSetShaderResources(0, 1, tempUAV);
		g_theRenderer->CSSetShader(m_matrixTransposeCS);
		g_theRenderer->DispatchCS(MATRIX_HEIGHT / TRANSPOSE_BLOCK_SIZE, MATRIX_WIDTH / TRANSPOSE_BLOCK_SIZE, 1);

		//Sort the row data
		g_theRenderer->CSSetShader(m_bitonicSortCS);
		g_theRenderer->DispatchCS(NUM_ELEMENTS / BITONIC_BLOCK_SIZE, 1, 1);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSim::GetMouseCursorPos()
{
	Vec2 normalizedMousePos = g_theInput->GetCursorNormalizedPosition();
	Vec2 mousePos = Vec2(RangeMap(normalizedMousePos.x, 0.f, 1.f, 0.f, MAP_WIDTH), RangeMap(normalizedMousePos.y, 0.f, 1.f, 0.f, MAP_HEIGHT));
	return mousePos;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSim::GetMouseCursorPosition()
{
	Vec2 normalizedMousePos = g_theInput->GetCursorNormalizedPosition();
	Vec2 cursorPos = Vec2(RangeMap(normalizedMousePos.x, 0.f, 1.f, 0.f, MAP_WIDTH), RangeMap(normalizedMousePos.y, 0.f, 1.f, 0.f, MAP_HEIGHT));
	return cursorPos;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Shape2D* FluidSim::GetShapeUnderMouse()
{
	AABB2 worldBoundingBox(Vec2::ZERO, Vec2(MAP_WIDTH, MAP_HEIGHT));

	Vec2 cursorPosition = GetMouseCursorPosition();


	for (int i = 0; i < m_tumblingConfigShapes.size(); i++)
	{
		OBB2D* shape = static_cast<OBB2D*>(m_tumblingConfigShapes[i]);
		if (!shape) continue;

		bool isMouseInside = IsPointInsideOBB2D(cursorPosition, *shape);
		if (isMouseInside)
		{
			return shape;
		}

	}

	return nullptr;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::ImGuiRenderFunctions()
{
	ImGuiGravityButtons();
	ImGuiSimulationMode();
	ImGuiRenderTypeForParticles();
	ImGuiViscosity();
	ImGuiTimestep();
	ImGuiRestDensity();
	ImGuiParticleMass();
	ImGuiPressureStiffness();
	ImGuiWallStiffness();
	ImGuiParticleRenderSize();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderTextureToScreenForCartoonFluid()
{
	ScreenSizeConstants constants = {};
	constants.screenHeight = (float)g_theWindow->GetClientDimensions().y;
	constants.screenWidth = (float)g_theWindow->GetClientDimensions().x;

	g_theRenderer->CopyCPUtoGPU(&constants, sizeof(ScreenSizeConstants), m_screenSizeCBO);
	g_theRenderer->BindConstantBufferWithGS(4, m_screenSizeCBO);

	g_theRenderer->BindRenderTexture(m_postProcessRenderTexture);

	g_theRenderer->VSSetShaderResource(1, 1, m_particleDensityUAV);
	g_theRenderer->IASetVertexBuffers(nullptr);

	g_theRenderer->BindShader(m_postProcessShader);
	g_theRenderer->DrawVertexArray((int)m_screenVerts.size(), m_screenVerts.data());
	g_theRenderer->BindTexture(nullptr);
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->VSSetShaderResource(1, 1, m_nullUAV);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderTextureToScreenForTransparentFluid()
{
	ScreenSizeConstants constants = {};
	constants.screenHeight = (float)g_theWindow->GetClientDimensions().y;
	constants.screenWidth = (float)g_theWindow->GetClientDimensions().x;
	constants.time = (float)GetCurrentTimeSeconds();
	g_theRenderer->CopyCPUtoGPU(&constants, sizeof(ScreenSizeConstants), m_screenSizeCBO);
	g_theRenderer->BindConstantBufferWithGS(4, m_screenSizeCBO);

	g_theRenderer->BindRenderTexture(m_transparentFluidRenderTexture);
	g_theRenderer->BindTexture(m_noiseTexture, 1);
	g_theRenderer->BindShader(m_transparentFluidShader);
	g_theRenderer->DrawVertexArray((int)m_screenVerts.size(), m_screenVerts.data());
	g_theRenderer->BindRenderTexture(nullptr);
	g_theRenderer->BindTexture(nullptr, 1);
	g_theRenderer->BindShader(nullptr);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeShapesForRandomShapesConfig()
{
	m_tumblingConfigShapesVerts.clear();

	for (int i = 0; i < 5; i++)
	{
		float randomAngle = m_rng.RollRandomFloatInRange(0.f, 360.f);
		Vec2 iBasis = Vec2::MakeFromPolarDegrees(randomAngle);
		Vec2 randPos = GetClampedRandomPositionInWorld();
		OBB2D* obb = new OBB2D();
		if (obb)
		{
			obb->m_center = randPos;
			obb->m_iBasisNormal = iBasis;
			obb->m_iBasisNormal.GetNormalized();
			obb->m_halfDimensions = Vec2(0.05f, 0.1f);
			obb->m_localToWorldMatrix = CreateLocalToWorldMatrix(*obb, obb->m_center);
			obb->m_worldToLocalMatrix = CreateWorldToLocalMatrix(obb->m_localToWorldMatrix);
			AddVertsForOBB2D(m_tumblingConfigShapesVerts, *obb, Rgba8::WHITE);
			m_tumblingConfigShapes.push_back(obb);
		}
	}

	Capsule2* capsule = new Capsule2();
	capsule->radius = m_rng.RollRandomFloatInRange(0.05f, 0.07f);
	capsule->m_bone.m_start = GetClampedRandomPositionInWorld();
	capsule->m_bone.m_end = GetClampedRandomPositionInWorld();
	if (capsule)
	{
		AddVertsForCapsule2D(m_tumblingConfigShapesVerts, capsule->m_bone.m_start, capsule->m_bone.m_end, capsule->radius, Rgba8::WHITE);
		m_tumblingConfigShapes.push_back(capsule);
	}
	
	Capsule2* capsule1 = new Capsule2();
	capsule1->radius = m_rng.RollRandomFloatInRange(0.02f, 0.07f);  
	capsule1->m_bone.m_start = GetClampedRandomPositionInWorld();
	capsule1->m_bone.m_end = GetClampedRandomPositionInWorld();
	if (capsule1)
	{
		AddVertsForCapsule2D(m_tumblingConfigShapesVerts, capsule1->m_bone.m_start, capsule1->m_bone.m_end, capsule1->radius, Rgba8::WHITE);
		m_tumblingConfigShapes.push_back(capsule1);
	}

	m_tumblingConfigMovingAABBVelocity = Vec2(1.0f, 0.0f) * 0.2f;

	m_tumblingDisc2D = new Disc2D();
	m_tumblingDisc2D->m_center = Vec2(0.5, 1.0f);
	m_tumblingDisc2D->m_radius = 0.2f;

	AddVertsForDisc2D(m_tumblingConfigShapesVerts, m_tumblingDisc2D->m_center, m_tumblingDisc2D->m_radius, Rgba8::WHITE);
	
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateMovingAABB2ForRandomShapesConfig(float deltaSeconds)
{
	if (g_theInput->WasKeyJustPressed('1'))
	{
		m_tumblingConfigMovingAABB2 = new AABB2(Vec2(2.2f, 0.05f), Vec2(2.3f, 0.6f));
	}

	if (m_tumblingConfigMovingAABB2)
	{
		Vec2 translation = m_tumblingConfigMovingAABBVelocity * deltaSeconds;
		m_tumblingConfigMovingAABB2->Translate(translation);

		float leftBoundary = 0.1f;
		float rightBoundary = MAP_WIDTH - 0.1f;

		// Check for collision with the left and right walls
		if (m_tumblingConfigMovingAABB2->m_mins.x <= leftBoundary || m_tumblingConfigMovingAABB2->m_maxs.x >= rightBoundary)
		{
			// Reverse the horizontal velocity
			m_tumblingConfigMovingAABBVelocity.x = -m_tumblingConfigMovingAABBVelocity.x;

			// Correct the position if AABB goes beyond bounds
			if (m_tumblingConfigMovingAABB2->m_mins.x < leftBoundary)
			{
				float overlap = m_tumblingConfigMovingAABB2->m_mins.x - leftBoundary;
				m_tumblingConfigMovingAABB2->SetCenter(Vec2(m_tumblingConfigMovingAABB2->GetCenter().x - overlap, m_tumblingConfigMovingAABB2->GetCenter().y));
			}
			else if (m_tumblingConfigMovingAABB2->m_maxs.x > rightBoundary)
			{
				float overlap = m_tumblingConfigMovingAABB2->m_maxs.x - rightBoundary;
				m_tumblingConfigMovingAABB2->SetCenter(Vec2(m_tumblingConfigMovingAABB2->GetCenter().x - overlap, m_tumblingConfigMovingAABB2->GetCenter().y));
			}
		}
		m_tumblingConfigMovingAABB2Verts.clear();
		AddVertsForAABB2D(m_tumblingConfigMovingAABB2Verts, *m_tumblingConfigMovingAABB2, Rgba8::WHITE);
	}
	
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderMovingAABB2ForRandomShapesConfig()
{
	g_theRenderer->SetRasterizerModes(RasterizerMode::SOLID_CULL_BACK);
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->BindTexture(nullptr);
	g_theRenderer->DrawVertexArray((int)m_tumblingConfigMovingAABB2Verts.size(), m_tumblingConfigMovingAABB2Verts.data());
	g_theRenderer->SetRasterizerModes(RasterizerMode::SOLID_CULL_BACK);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderShapesForRandomShapesConfig()
{
	g_theRenderer->SetBlendMode(BlendMode::ADDITIVE);
	g_theRenderer->SetRasterizerModes(RasterizerMode::SOLID_CULL_BACK);
	g_theRenderer->BindTexture(nullptr);
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->DrawVertexArray((int)m_tumblingConfigShapesVerts.size(), m_tumblingConfigShapesVerts.data());
	g_theRenderer->SetBlendMode(BlendMode::ALPHA);

	if (m_tumblingConfigMovingAABB2)
	{
		RenderMovingAABB2ForRandomShapesConfig();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Mat44 FluidSim::CreateWorldToLocalMatrix(Mat44& localToWorldMatrix)
{
	Mat44 worldToLocal = localToWorldMatrix.GetOrthonormalInverse();
	return worldToLocal;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Mat44 FluidSim::CreateLocalToWorldMatrix(OBB2D& box, Vec2& obbCenter)
{
	Mat44 localToWorld;
	Vec2 jBasisNormal = box.m_iBasisNormal.GetRotated90Degrees();

	localToWorld.SetIJT2D(box.m_iBasisNormal, jBasisNormal, obbCenter);
	return localToWorld;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateSimulationConstantBufferFromRandomShapesConfig(SimulationConstants& constants)
{
	if (m_tumblingConfigShapes[0])
	{
		OBB2D* obb1 = static_cast<OBB2D*>(m_tumblingConfigShapes[0]);
		constants.tumblingConfigOBB1HalfDimensions = obb1->m_halfDimensions;
		constants.tumblingConfigOBB1Center = obb1->m_center;
		constants.tumblingConfigOBB1LocalToWorldMatrix = obb1->m_localToWorldMatrix;
		constants.tumblingConfigOBB1WorldToLocalMatrix = obb1->m_worldToLocalMatrix;
	}

	if (m_tumblingConfigShapes[1])
	{
		OBB2D* obb2 = static_cast<OBB2D*>(m_tumblingConfigShapes[1]);
		constants.tumblingConfigOBB2HalfDimensions = obb2->m_halfDimensions;
		constants.tumblingConfigOBB2Center = obb2->m_center;
		constants.tumblingConfigOBB2LocalToWorldMatrix = obb2->m_localToWorldMatrix;
		constants.tumblingConfigOBB2WorldToLocalMatrix = obb2->m_worldToLocalMatrix;
	}

	if (m_tumblingConfigShapes[2])
	{
		OBB2D* obb3 = static_cast<OBB2D*>(m_tumblingConfigShapes[2]);
		constants.tumblingConfigOBB3HalfDimensions = obb3->m_halfDimensions;
		constants.tumblingConfigOBB3Center = obb3->m_center;
		constants.tumblingConfigOBB3LocalToWorldMatrix = obb3->m_localToWorldMatrix;
		constants.tumblingConfigOBB3WorldToLocalMatrix = obb3->m_worldToLocalMatrix;
	}

	if (m_tumblingConfigShapes[3])
	{
		OBB2D* obb4 = static_cast<OBB2D*>(m_tumblingConfigShapes[3]);
		constants.tumblingConfigOBB4HalfDimensions = obb4->m_halfDimensions;
		constants.tumblingConfigOBB4Center = obb4->m_center;
		constants.tumblingConfigOBB4LocalToWorldMatrix = obb4->m_localToWorldMatrix;
		constants.tumblingConfigOBB4WorldToLocalMatrix = obb4->m_worldToLocalMatrix;
	}

	if (m_tumblingConfigShapes[4])
	{
		OBB2D* obb5 = static_cast<OBB2D*>(m_tumblingConfigShapes[4]);
		constants.tumblingConfigOBB5HalfDimensions = obb5->m_halfDimensions;
		constants.tumblingConfigOBB5Center = obb5->m_center;
		constants.tumblingConfigOBB5LocalToWorldMatrix = obb5->m_localToWorldMatrix;
		constants.tumblingConfigOBB5WorldToLocalMatrix = obb5->m_worldToLocalMatrix;
	}

	if (m_tumblingConfigShapes[5])
	{
		Capsule2* capsule = static_cast<Capsule2*>(m_tumblingConfigShapes[5]);
		{
			constants.tumblingConfigCapsule1Start = capsule->m_bone.m_start;
			constants.tumblingConfigCapsule1End = capsule->m_bone.m_end;
			constants.tumblingConfigCapsule1Radius = capsule->radius;
		}
	}

	if (m_tumblingConfigShapes[6])
	{
		Capsule2* capsule1 = static_cast<Capsule2*>(m_tumblingConfigShapes[6]);
		{
			constants.tumblingConfigCapsule2Start = capsule1->m_bone.m_start;
			constants.tumblingConfigCapsule2End = capsule1->m_bone.m_end;
			constants.tumblingConfigCapsule2Radius = capsule1->radius;
		}
	}

	if (m_tumblingConfigMovingAABB2)
	{
		constants.tumblingConfigHorizontalMovingAABB1Size = Vec2(m_tumblingConfigMovingAABB2->m_maxs.x - m_tumblingConfigMovingAABB2->m_mins.x, m_tumblingConfigMovingAABB2->m_maxs.y - m_tumblingConfigMovingAABB2->m_mins.y);
		constants.tumblingConfigHorizontalMovingAABB1Center = Vec2((m_tumblingConfigMovingAABB2->m_maxs.x + m_tumblingConfigMovingAABB2->m_mins.x) / 2, (m_tumblingConfigMovingAABB2->m_maxs.y + m_tumblingConfigMovingAABB2->m_mins.y) / 2);
	}
	else
	{
		constants.tumblingConfigHorizontalMovingAABB1Size = Vec2(0.f, 0.f);
		constants.tumblingConfigHorizontalMovingAABB1Center = Vec2(0.f, 0.f);
	}

	if (m_tumblingDisc2D)
	{
		constants.tumblingDiscCenter = m_tumblingDisc2D->m_center;
		constants.tumblingDiscRadius = m_tumblingDisc2D->m_radius;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeShapesForWindmillConfig()
{
	m_windmillDisc1 = new Disc2D();
	m_windmillDisc1->m_center = Vec2(MAP_WIDTH / 2, MAP_HEIGHT / 2);
	m_windmillDisc1->m_radius = 0.1f;
	AddVertsForDisc2D(m_windmillConfigShapesVerts, m_windmillDisc1->m_center, m_windmillDisc1->m_radius, Rgba8::YELLOW);

	m_windmillOBB1 = new OBB2D();
	m_windmillOBB1->m_center = Vec2(MAP_WIDTH / 2, MAP_HEIGHT / 2);
	m_windmillOBB1->m_halfDimensions = Vec2(0.3f, 0.05f);
	float windmillOBB1Orientation = 45.f;
	m_windmillOBB1->m_iBasisNormal = Vec2::MakeFromPolarDegrees(windmillOBB1Orientation);
	m_windmillOBB1->m_localToWorldMatrix = CreateLocalToWorldMatrix(*m_windmillOBB1, m_windmillOBB1->m_center);
	m_windmillOBB1->m_worldToLocalMatrix = CreateWorldToLocalMatrix(m_windmillOBB1->m_localToWorldMatrix);
	AddVertsForOBB2D(m_windmillConfigShapesVerts, *m_windmillOBB1, Rgba8::YELLOW);

	m_windmillOBB2 = new OBB2D();
	m_windmillOBB2->m_center = Vec2(MAP_WIDTH / 2, MAP_HEIGHT / 2);
	m_windmillOBB2->m_halfDimensions = Vec2(0.3f, 0.05f);
	float windmillOBB2Orientation = 135.f;
	m_windmillOBB2->m_iBasisNormal = Vec2::MakeFromPolarDegrees(windmillOBB2Orientation);
	m_windmillOBB2->m_localToWorldMatrix = CreateLocalToWorldMatrix(*m_windmillOBB2, m_windmillOBB2->m_center);
	m_windmillOBB2->m_worldToLocalMatrix = CreateWorldToLocalMatrix(m_windmillOBB2->m_localToWorldMatrix);
	AddVertsForOBB2D(m_windmillConfigShapesVerts, *m_windmillOBB2, Rgba8::YELLOW);

	m_windmillLeftBorderCapsule = new Capsule2();
	m_windmillLeftBorderCapsule->radius = 0.01f;
	m_windmillLeftBorderCapsule->m_bone.m_start = Vec2(0.f, 0.f);
	m_windmillLeftBorderCapsule->m_bone.m_end = Vec2(0.f, MAP_HEIGHT);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillLeftBorderCapsule->m_bone.m_start, m_windmillLeftBorderCapsule->m_bone.m_end, m_windmillLeftBorderCapsule->radius, Rgba8::YELLOW);

	m_windmillRightBorderCapsule = new Capsule2();
	m_windmillRightBorderCapsule->radius = 0.01f;
	m_windmillRightBorderCapsule->m_bone.m_start = Vec2(MAP_WIDTH, 0.f);
	m_windmillRightBorderCapsule->m_bone.m_end = Vec2(MAP_WIDTH, MAP_HEIGHT);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillRightBorderCapsule->m_bone.m_start, m_windmillRightBorderCapsule->m_bone.m_end, m_windmillRightBorderCapsule->radius, Rgba8::YELLOW);

	m_windmillBottomBorderCapsule = new Capsule2();
	m_windmillBottomBorderCapsule->radius = 0.01f;
	m_windmillBottomBorderCapsule->m_bone.m_start = Vec2(0.15f, 0.f);
	m_windmillBottomBorderCapsule->m_bone.m_end = Vec2(MAP_WIDTH, 0.f);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillBottomBorderCapsule->m_bone.m_start, m_windmillBottomBorderCapsule->m_bone.m_end, m_windmillBottomBorderCapsule->radius, Rgba8::YELLOW);

	m_windmillTopBorderCapsule2 = new Capsule2();
	m_windmillTopBorderCapsule2->radius = 0.01f;
	m_windmillTopBorderCapsule2->m_bone.m_start = Vec2(1.3f, MAP_HEIGHT);
	m_windmillTopBorderCapsule2->m_bone.m_end = Vec2(1.7f, MAP_HEIGHT);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule2->m_bone.m_start, m_windmillTopBorderCapsule2->m_bone.m_end, m_windmillTopBorderCapsule2->radius, Rgba8::YELLOW);

	m_windmillTopBorderCapsule1 = new Capsule2();
	m_windmillTopBorderCapsule1->radius = 0.01f;
	m_windmillTopBorderCapsule1->m_bone.m_start = Vec2(0.f , MAP_HEIGHT);
	m_windmillTopBorderCapsule1->m_bone.m_end = Vec2(1.1f, MAP_HEIGHT);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule1->m_bone.m_start, m_windmillTopBorderCapsule1->m_bone.m_end, m_windmillTopBorderCapsule1->radius, Rgba8::YELLOW);

	m_windmillTopBorderCapsule3 = new Capsule2();
	m_windmillTopBorderCapsule3->radius = 0.01f;
	m_windmillTopBorderCapsule3->m_bone.m_start = Vec2(1.9f, MAP_HEIGHT);
	m_windmillTopBorderCapsule3->m_bone.m_end = Vec2(MAP_WIDTH, MAP_HEIGHT);
	AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule3->m_bone.m_start, m_windmillTopBorderCapsule3->m_bone.m_end, m_windmillTopBorderCapsule3->radius, Rgba8::YELLOW);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderShapesForWindwillConfig()
{
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->DrawVertexArray((int)m_windmillConfigShapesVerts.size(), m_windmillConfigShapesVerts.data());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateSimulationConstantBufferForHourglassConfig(SimulationConstants& constants)
{
	if (m_hourglassConfigCapsule1 && m_renderTheTopCapsuleForHourglassMesh)
	{
		constants.hourglassConfigCapsule1Start = m_hourglassConfigCapsule1->m_bone.m_start;
		constants.hourglassConfigCapsule1End = m_hourglassConfigCapsule1->m_bone.m_end;
		constants.hourglassConfigCapsule1Radius = m_hourglassConfigCapsule1->radius;
	}
	else
	{
		constants.hourglassConfigCapsule1Start = Vec2(-10.f, -10.f);
		constants.hourglassConfigCapsule1End = Vec2(-10.f, -10.f);
		constants.hourglassConfigCapsule1Radius = 0.f;
	}
	if (m_hourglassConfigCapsule2)
	{
		constants.hourglassConfigCapsule2Start = m_hourglassConfigCapsule2->m_bone.m_start;
		constants.hourglassConfigCapsule2End = m_hourglassConfigCapsule2->m_bone.m_end;
		constants.hourglassConfigCapsule2Radius = m_hourglassConfigCapsule2->radius;
	}

	if (m_hourglassConfigCapsule3)
	{
		constants.hourglassConfigCapsule3Start = m_hourglassConfigCapsule3->m_bone.m_start;
		constants.hourglassConfigCapsule3End = m_hourglassConfigCapsule3->m_bone.m_end;
		constants.hourglassConfigCapsule3Radius = m_hourglassConfigCapsule3->radius;
	}

	if (m_hourglassConfigCapsule4)
	{
		constants.hourglassConfigCapsule4Start = m_hourglassConfigCapsule4->m_bone.m_start;
		constants.hourglassConfigCapsule4End = m_hourglassConfigCapsule4->m_bone.m_end;
		constants.hourglassConfigCapsule4Radius = m_hourglassConfigCapsule4->radius;
	}

	if (m_hourglassConfigCapsule5)
	{
		constants.hourglassConfigCapsule5Start = m_hourglassConfigCapsule5->m_bone.m_start;
		constants.hourglassConfigCapsule5End = m_hourglassConfigCapsule5->m_bone.m_end;
		constants.hourglassConfigCapsule5Radius = m_hourglassConfigCapsule5->radius;
	}

	if (m_hourglassConfigCapsule6)
	{
		constants.hourglassConfigCapsule6Start = m_hourglassConfigCapsule6->m_bone.m_start;
		constants.hourglassConfigCapsule6End = m_hourglassConfigCapsule6->m_bone.m_end;
		constants.hourglassConfigCapsule6Radius = m_hourglassConfigCapsule6->radius;
	}

	if (m_hourglassConfigMiddleCapsule1)
	{
		constants.hourglassConfigMiddleCapsuleStart1 = m_hourglassConfigMiddleCapsule1->m_bone.m_start;
		constants.hourglassConfigMiddleCapsuleEnd1 = m_hourglassConfigMiddleCapsule1->m_bone.m_end;
		constants.hourglassConfigMiddleCapsule1Radius = m_hourglassConfigMiddleCapsule1->radius;
	}

	if (m_hourglassConfigMiddleCapsule2)
	{
		constants.hourglassConfigMiddleCapsuleStart2 = m_hourglassConfigMiddleCapsule2->m_bone.m_start;
		constants.hourglassConfigMiddleCapsuleEnd2 = m_hourglassConfigMiddleCapsule2->m_bone.m_end;
		constants.hourglassConfigMiddleCapsule2Radius = m_hourglassConfigMiddleCapsule2->radius;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Mat44 FluidSim::GetModelMatrixForHourglassConfig(float deltaSeconds)
{
	Mat44 matrix;
	float rotationSpeed = 0.001f;
	m_totalRotationDegrees += rotationSpeed * deltaSeconds;
	m_totalRotationDegrees = fmod(m_totalRotationDegrees, 360.0f);

	if (m_totalRotationDegrees > 0.1f)
	{
		m_totalRotationDegrees = 0.f;
	}

	if (m_shouldRotateHourglassMesh)
	{
		Vec2 meshCenter = GetHourGlassMeshCenter();
		Mat44 translateToOriginMatrix = Mat44::CreateTranslation2D(-meshCenter);
		Mat44 rotationMatrix = Mat44::CreateZRotationDegrees(m_totalRotationDegrees);
		Mat44 translateBackMatrix = Mat44::CreateTranslation2D(meshCenter);
		matrix.Append(translateBackMatrix);
		matrix.Append(rotationMatrix);
		matrix.Append(translateToOriginMatrix);
	}
	
	return matrix;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSim::GetHourGlassMeshCenter()
{
	Vec2 meshCenter = Vec2(MAP_WIDTH / 2, MAP_HEIGHT / 2);
	return meshCenter;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateHourglassMeshConfig(float deltaSeconds)
{
	Mat44 transformationMatrix = GetModelMatrixForHourglassConfig(deltaSeconds);
	m_hourglassConfigShapesVerts.clear();
	m_hourglassTopCapsuleVerts.clear();
	
	if (m_shouldRotateHourglassMesh)
	{
		if (m_hourglassConfigCapsule1)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule1, transformationMatrix, true);
		}

		if (m_hourglassConfigCapsule2)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule2, transformationMatrix);
		}

		if (m_hourglassConfigCapsule3)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule3, transformationMatrix);
		}

		if (m_hourglassConfigCapsule4)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule4, transformationMatrix);
		}

		if (m_hourglassConfigCapsule5)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule5, transformationMatrix);
		}

		if (m_hourglassConfigCapsule6)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigCapsule6, transformationMatrix);
		}

		if (m_hourglassConfigCapsule6)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigMiddleCapsule1, transformationMatrix);
		}

		if (m_hourglassConfigCapsule6)
		{
			TansformAndAddVertsForHourglassConfigShape(m_hourglassConfigMiddleCapsule2, transformationMatrix);
		}
	}

	else
	{
		if (m_hourglassConfigCapsule1)
		{
			AddVertsForCapsule2D(m_hourglassTopCapsuleVerts, m_hourglassConfigCapsule1->m_bone.m_start, m_hourglassConfigCapsule1->m_bone.m_end, m_hourglassConfigCapsule1->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigCapsule2)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule2->m_bone.m_start, m_hourglassConfigCapsule2->m_bone.m_end, m_hourglassConfigCapsule2->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigCapsule3)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule3->m_bone.m_start, m_hourglassConfigCapsule3->m_bone.m_end, m_hourglassConfigCapsule3->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigCapsule4)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule4->m_bone.m_start, m_hourglassConfigCapsule4->m_bone.m_end, m_hourglassConfigCapsule4->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigCapsule5)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule5->m_bone.m_start, m_hourglassConfigCapsule5->m_bone.m_end, m_hourglassConfigCapsule5->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigCapsule6)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule6->m_bone.m_start, m_hourglassConfigCapsule6->m_bone.m_end, m_hourglassConfigCapsule6->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigMiddleCapsule1)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigMiddleCapsule1->m_bone.m_start, m_hourglassConfigMiddleCapsule1->m_bone.m_end, m_hourglassConfigMiddleCapsule1->radius, Rgba8::WHITE);
		}

		if (m_hourglassConfigMiddleCapsule2)
		{
			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigMiddleCapsule2->m_bone.m_start, m_hourglassConfigMiddleCapsule2->m_bone.m_end, m_hourglassConfigMiddleCapsule2->radius, Rgba8::WHITE);
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::TansformAndAddVertsForHourglassConfigShape(Capsule2* capsule, const Mat44& matrix, bool isTopCapsule /*= false*/)
{
	if (!isTopCapsule)
	{
		if (capsule)
		{
			Vec2 transformedStart = matrix.TransformPosition2D(capsule->m_bone.m_start);
			Vec2 transformedEnd = matrix.TransformPosition2D(capsule->m_bone.m_end);

			capsule->m_bone.m_start = transformedStart;
			capsule->m_bone.m_end = transformedEnd;

			AddVertsForCapsule2D(m_hourglassConfigShapesVerts, transformedStart, transformedEnd, capsule->radius, Rgba8::WHITE);
		}
	}
	else
	{
		if (capsule)
		{
			Vec2 transformedStart = matrix.TransformPosition2D(capsule->m_bone.m_start);
			Vec2 transformedEnd = matrix.TransformPosition2D(capsule->m_bone.m_end);

			capsule->m_bone.m_start = transformedStart;
			capsule->m_bone.m_end = transformedEnd;

			AddVertsForCapsule2D(m_hourglassTopCapsuleVerts, transformedStart, transformedEnd, capsule->radius, Rgba8::WHITE);
		}
	}
	
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesForTransparentFluidConfig(unsigned int numParticles)
{
	m_numParticles = numParticles;
	float areaWidth = MAP_WIDTH;
	int particlesPerRow = static_cast<int>(areaWidth / m_initialParticleSpacing);
	if (particlesPerRow == 0)
	{
		particlesPerRow = 1;
	}

	std::vector<Particle> particles(numParticles);

	for (unsigned int i = 0; i < numParticles; i++)
	{
		int row = i / particlesPerRow;
		int col = i % particlesPerRow;
		particles[i].m_position = Vec2(m_initialParticleSpacing * (float)col, m_initialParticleSpacing * (float)row - 0.3f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_particleRenderType = TRANSPARENT;
	m_renderType = 3;
	m_particleSize = 0.03f;

	InitializeShapesForWaterElevatorConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeBackgroundTextureForTransparentFluidConfig()
{
	m_backgroundTexture = g_theRenderer->CreateOrGetTextureFromFile("Data/Images/BGT.png");
	AddVertsForAABB2D(m_backgroundQuadForTransparentFluidConfigVerts, m_screenBounds, Rgba8::WHITE);
	m_transparentFluidRenderTexture = g_theRenderer->CreateRenderTexture(g_theWindow->GetClientDimensions().x, g_theWindow->GetClientDimensions().y);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderBackgroundTextureForTransparentFluidConfig()
{
	g_theRenderer->SetBlendMode(BlendMode::ALPHA);
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->BindTexture(m_backgroundTexture);
	g_theRenderer->DrawVertexArray((int)m_backgroundQuadForTransparentFluidConfigVerts.size(), m_backgroundQuadForTransparentFluidConfigVerts.data());
	g_theRenderer->BindTexture(nullptr);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateSimulationConstantBufferForWindmillConfig(SimulationConstants& constants)
{
	if (m_windmillDisc1)
	{
		constants.windmillDisc1Center = m_windmillDisc1->m_center;
		constants.windmillDisc1Radius = m_windmillDisc1->m_radius;
	}
	else
	{
		constants.windmillDisc1Center = Vec2(0.f, 0.f);
		constants.windmillDisc1Radius = 0.f;
	}

	if (m_windmillOBB1)
	{
		constants.windmillConfigOBB1HalfDimensions = m_windmillOBB1->m_halfDimensions;
		constants.windmillConfigOBB1Center = m_windmillOBB1->m_center;
		constants.windmillConfigOBB1LocalToWorldMatrix = m_windmillOBB1->m_localToWorldMatrix;
		constants.windmillConfigOBB1WorldToLocalMatrix = m_windmillOBB1->m_worldToLocalMatrix;
	}

	if (m_windmillOBB2)
	{
		constants.windmillConfigOBB2HalfDimensions = m_windmillOBB2->m_halfDimensions;
		constants.windmillConfigOBB2Center = m_windmillOBB2->m_center;
		constants.windmillConfigOBB2LocalToWorldMatrix = m_windmillOBB2->m_localToWorldMatrix;
		constants.windmillConfigOBB2WorldToLocalMatrix = m_windmillOBB2->m_worldToLocalMatrix;
	}

	if (m_windmillBottomBorderCapsule)
	{
		constants.windmillConfigBottomBorderCapsuleStart = m_windmillBottomBorderCapsule->m_bone.m_start;
		constants.windmillConfigBottomBorderCapsuleEnd = m_windmillBottomBorderCapsule->m_bone.m_end;
		constants.windmillConfigBottomBorderCapsuleRadius = m_windmillBottomBorderCapsule->radius;
	}

	if (m_windmillLeftBorderCapsule)
	{
		constants.windmillConfigLeftBorderCapsuleStart = m_windmillLeftBorderCapsule->m_bone.m_start;
		constants.windmillConfigLeftBorderCapsuleEnd = m_windmillLeftBorderCapsule->m_bone.m_end;
		constants.windmillConfigLeftBorderCapsuleRadius = m_windmillLeftBorderCapsule->radius;
	}

	if (m_windmillRightBorderCapsule)
	{
		constants.windmillConfigRightBorderCapsuleStart = m_windmillRightBorderCapsule->m_bone.m_start;
		constants.windmillConfigRightBorderCapsuleEnd = m_windmillRightBorderCapsule->m_bone.m_end;
		constants.windmillConfigRightBorderCapsuleRadius = m_windmillRightBorderCapsule->radius;
	}

	if (m_windmillTopBorderCapsule1)
	{
		constants.windmillConfigTopBorderCapsule1Start = m_windmillTopBorderCapsule1->m_bone.m_start;
		constants.windmillConfigTopBorderCapsule1End = m_windmillTopBorderCapsule1->m_bone.m_end;
		constants.windmillConfigTopBorderCapsule1Radius = m_windmillTopBorderCapsule1->radius;
	}

	if (m_windmillTopBorderCapsule2)
	{
		constants.windmillConfigTopBorderCapsule2Start = m_windmillTopBorderCapsule2->m_bone.m_start;
		constants.windmillConfigTopBorderCapsule2End = m_windmillTopBorderCapsule2->m_bone.m_end;
		constants.windmillConfigTopBorderCapsule2Radius = m_windmillTopBorderCapsule2->radius;
	}

	if (m_windmillTopBorderCapsule3)
	{
		constants.windmillConfigTopBorderCapsule3Start = m_windmillTopBorderCapsule3->m_bone.m_start;
		constants.windmillConfigTopBorderCapsule3End = m_windmillTopBorderCapsule3->m_bone.m_end;
		constants.windmillConfigTopBorderCapsule3Radius = m_windmillTopBorderCapsule3->radius;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateRotationOfWindmillMesh(float deltaSeconds)
{
	m_windmillConfigShapesVerts.clear();

	float rotationSpeed = 45.0f; 

	float rotationDeltaDegrees = rotationSpeed * deltaSeconds;

	if (m_windmillOBB1)
	{
		m_windmillOBB1->RotateAboutCenter(rotationDeltaDegrees);
		m_windmillOBB1->m_localToWorldMatrix = CreateLocalToWorldMatrix(*m_windmillOBB1, m_windmillOBB1->m_center);
		m_windmillOBB1->m_worldToLocalMatrix = CreateWorldToLocalMatrix(m_windmillOBB1->m_localToWorldMatrix);
		AddVertsForOBB2D(m_windmillConfigShapesVerts, *m_windmillOBB1, Rgba8::YELLOW);
	}

	if (m_windmillOBB2)
	{
		m_windmillOBB2->RotateAboutCenter(rotationDeltaDegrees);
		m_windmillOBB2->m_localToWorldMatrix = CreateLocalToWorldMatrix(*m_windmillOBB2, m_windmillOBB2->m_center);
		m_windmillOBB2->m_worldToLocalMatrix = CreateWorldToLocalMatrix(m_windmillOBB2->m_localToWorldMatrix);
		AddVertsForOBB2D(m_windmillConfigShapesVerts, *m_windmillOBB2, Rgba8::YELLOW);
	}

	if (m_windmillDisc1)
	{
		AddVertsForDisc2D(m_windmillConfigShapesVerts, m_windmillDisc1->m_center, m_windmillDisc1->m_radius, Rgba8::YELLOW);
	}

	if (m_windmillLeftBorderCapsule)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillLeftBorderCapsule->m_bone.m_start, m_windmillLeftBorderCapsule->m_bone.m_end, m_windmillLeftBorderCapsule->radius, Rgba8::YELLOW);
	}

	if (m_windmillRightBorderCapsule)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillRightBorderCapsule->m_bone.m_start, m_windmillRightBorderCapsule->m_bone.m_end, m_windmillRightBorderCapsule->radius, Rgba8::YELLOW);
	}

	if (m_windmillBottomBorderCapsule)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillBottomBorderCapsule->m_bone.m_start, m_windmillBottomBorderCapsule->m_bone.m_end, m_windmillBottomBorderCapsule->radius, Rgba8::YELLOW);
	}

	if (m_windmillTopBorderCapsule2)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule2->m_bone.m_start, m_windmillTopBorderCapsule2->m_bone.m_end, m_windmillTopBorderCapsule2->radius, Rgba8::YELLOW);
	}

	if (m_windmillTopBorderCapsule1)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule1->m_bone.m_start, m_windmillTopBorderCapsule1->m_bone.m_end, m_windmillTopBorderCapsule1->radius, Rgba8::YELLOW);
	}

	if (m_windmillTopBorderCapsule3)
	{
		AddVertsForCapsule2D(m_windmillConfigShapesVerts, m_windmillTopBorderCapsule3->m_bone.m_start, m_windmillTopBorderCapsule3->m_bone.m_end, m_windmillTopBorderCapsule3->radius, Rgba8::YELLOW);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeShapesForHourglassConfig()
{
	m_hourglassConfigShapesVerts.clear();
	m_hourglassTopCapsuleVerts.clear();
	
	m_hourglassConfigCapsule1 = new Capsule2();
	m_hourglassConfigCapsule1->m_bone.m_start = Vec2(0.750f, 1.3125f);
	m_hourglassConfigCapsule1->m_bone.m_end = Vec2(2.250f, 1.3125f);
	m_hourglassConfigCapsule1->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassTopCapsuleVerts, m_hourglassConfigCapsule1->m_bone.m_start, m_hourglassConfigCapsule1->m_bone.m_end, m_hourglassConfigCapsule1->radius, Rgba8::WHITE);

	m_hourglassConfigCapsule2 = new Capsule2();
	m_hourglassConfigCapsule2->m_bone.m_start = Vec2(2.250f, 1.3125f);
	m_hourglassConfigCapsule2->m_bone.m_end = Vec2(1.6875f, 0.750f);
	m_hourglassConfigCapsule2->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule2->m_bone.m_start, m_hourglassConfigCapsule2->m_bone.m_end, m_hourglassConfigCapsule2->radius, Rgba8::WHITE);

	m_hourglassConfigCapsule3 = new Capsule2();
	m_hourglassConfigCapsule3->m_bone.m_start = Vec2(1.6875f, 0.750f);
	m_hourglassConfigCapsule3->m_bone.m_end = Vec2(2.250f, 0.1875f);
	m_hourglassConfigCapsule3->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule3->m_bone.m_start, m_hourglassConfigCapsule3->m_bone.m_end, m_hourglassConfigCapsule3->radius, Rgba8::WHITE);

	m_hourglassConfigCapsule4 = new Capsule2();
	m_hourglassConfigCapsule4->m_bone.m_start = Vec2(2.250f, 0.1875f);
	m_hourglassConfigCapsule4->m_bone.m_end = Vec2(0.750f, 0.1875f);
	m_hourglassConfigCapsule4->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule4->m_bone.m_start, m_hourglassConfigCapsule4->m_bone.m_end, m_hourglassConfigCapsule4->radius, Rgba8::WHITE);

	m_hourglassConfigCapsule5 = new Capsule2();
	m_hourglassConfigCapsule5->m_bone.m_start = Vec2(0.750f, 0.1875f);
	m_hourglassConfigCapsule5->m_bone.m_end = Vec2(1.3125f, 0.750f);
	m_hourglassConfigCapsule5->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule5->m_bone.m_start, m_hourglassConfigCapsule5->m_bone.m_end, m_hourglassConfigCapsule5->radius, Rgba8::WHITE);

	m_hourglassConfigCapsule6 = new Capsule2();
	m_hourglassConfigCapsule6->m_bone.m_start = Vec2(1.3125f, 0.750f);
	m_hourglassConfigCapsule6->m_bone.m_end = Vec2(0.750f, 1.3125f);
	m_hourglassConfigCapsule6->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigCapsule6->m_bone.m_start, m_hourglassConfigCapsule6->m_bone.m_end, m_hourglassConfigCapsule6->radius, Rgba8::WHITE);

	m_hourglassConfigMiddleCapsule1 = new Capsule2();
	m_hourglassConfigMiddleCapsule1->m_bone.m_start = Vec2(1.3125f, 0.750f);
	m_hourglassConfigMiddleCapsule1->m_bone.m_end = Vec2(1.40625f, 0.750f);
	m_hourglassConfigMiddleCapsule1->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigMiddleCapsule1->m_bone.m_start, m_hourglassConfigMiddleCapsule1->m_bone.m_end, m_hourglassConfigMiddleCapsule1->radius, Rgba8::WHITE);

	m_hourglassConfigMiddleCapsule2 = new Capsule2();
	m_hourglassConfigMiddleCapsule2->m_bone.m_start = Vec2(1.59375f, 0.750f);
	m_hourglassConfigMiddleCapsule2->m_bone.m_end = Vec2(1.6875f, 0.750f);
	m_hourglassConfigMiddleCapsule2->radius = 0.01f;
	AddVertsForCapsule2D(m_hourglassConfigShapesVerts, m_hourglassConfigMiddleCapsule2->m_bone.m_start, m_hourglassConfigMiddleCapsule2->m_bone.m_end, m_hourglassConfigMiddleCapsule2->radius, Rgba8::WHITE);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderShapesForHourglassConfig()
{
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->DrawVertexArray((int)m_hourglassConfigShapesVerts.size(), m_hourglassConfigShapesVerts.data());

	if (m_renderTheTopCapsuleForHourglassMesh)
	{
		g_theRenderer->BindShader(nullptr);
		g_theRenderer->DrawVertexArray((int)m_hourglassTopCapsuleVerts.size(), m_hourglassTopCapsuleVerts.data());
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSim::GetClampedRandomPositionInWorld()
{
	float randX = m_rng.RollRandomFloatInRange(MAP_WIDTH * 0.1f, MAP_WIDTH * 0.9f);
	float randY = m_rng.RollRandomFloatInRange(MAP_HEIGHT * 0.1f, MAP_HEIGHT * 0.9f);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSim::GetClampedRandomPositionInWorld(float shapeRadius)
{
	float randX = m_rng.RollRandomFloatInRange(shapeRadius, MAP_WIDTH - shapeRadius);
	float randY = m_rng.RollRandomFloatInRange(shapeRadius, MAP_HEIGHT - shapeRadius);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSim::GetRandomPositionInWorld()
{
	float randX = m_rng.RollRandomFloatInRange(0.0f, MAP_WIDTH);
	float randY = m_rng.RollRandomFloatInRange(0.0f, MAP_HEIGHT);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderShapesBasedOnConfig()
{
	if (g_theGame->m_initialParticleConfigGPU == RANDOM_SHAPES)
	{
		RenderShapesForRandomShapesConfig();
	}

	if (g_theGame->m_initialParticleConfigGPU == WINDMILL)
	{
		RenderShapesForWindwillConfig();
	}

	if (g_theGame->m_initialParticleConfigGPU == HOURGLASS)
	{
		g_theRenderer->SetModelConstants(GetModelMatrixForHourglassConfig(g_theApp->m_clock.GetDeltaSeconds()));
		RenderShapesForHourglassConfig();
		g_theRenderer->SetModelConstants();
	}

	if (g_theGame->m_initialParticleConfigGPU == WATER_ELEVATOR)
	{
		RenderShapesForWaterElevatorConfig();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderCartoonFluid()
{
	g_theRenderer->SetBlendMode(BlendMode::ADDITIVE);
	g_theRenderer->SetRenderTarget(m_postProcessRenderTexture);
	g_theRenderer->ClearRenderTarget(m_postProcessRenderTexture, Rgba8::BLACK);

	RenderFluidParticles();

	g_theRenderer->SetBlendMode(BlendMode::ALPHA);
	g_theRenderer->SetBackBufferRenderTarget();

	RenderTextureToScreenForCartoonFluid();
	RenderShapesBasedOnConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderTransparentFluidForTheTransparentConfig()
{
	if (m_particleRenderType == DENSITY_DEBUG)
	{
		RenderBackgroundTextureForTransparentFluidConfig();

		RenderFluidParticles();
	}

	else
	{
		RenderBackgroundTextureForTransparentFluidConfig();

		g_theRenderer->SetBlendMode(BlendMode::ADDITIVE);
		g_theRenderer->SetRenderTarget(m_transparentFluidRenderTexture);
		g_theRenderer->ClearRenderTarget(m_transparentFluidRenderTexture, Rgba8::BLACK);

		RenderFluidParticles();
		g_theRenderer->SetBlendMode(BlendMode::ALPHA);
		g_theRenderer->SetBackBufferRenderTarget();

		RenderTextureToScreenForTransparentFluid();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeParticlesForWaterElevatorConfig(unsigned int numParticles)
{
	m_numParticles = numParticles;
	float areaWidth = MAP_WIDTH; 
	int particlesPerRow = static_cast<int>(areaWidth / m_initialParticleSpacing);
	if (particlesPerRow == 0)
	{
		particlesPerRow = 1;
	}

	std::vector<Particle> particles(numParticles);

	for (unsigned int i = 0; i < numParticles; i++)
	{
		int row = i / particlesPerRow;  
		int col = i % particlesPerRow; 
		particles[i].m_position = Vec2(m_initialParticleSpacing * (float)col, m_initialParticleSpacing * (float)row - 0.3f);
		particles[i].type = 1;
	}

	m_particleUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_sortedParticlesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(Particle), particles.data());
	m_particleDensityUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleDensity), nullptr);
	m_particleForcesUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(ParticleForces), nullptr);
	m_gridUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridPingPongUAV = g_theRenderer->CreateUnorderedAccessBuffer(particles.size(), sizeof(unsigned int), nullptr);
	m_gridIndicesUAV = g_theRenderer->CreateUnorderedAccessBuffer(NUM_GRID_INDICES, sizeof(UINT2), nullptr);
	m_renderType = 1;
	m_particleRenderType = BLUE_PALLETE;
	m_particleSize = 0.003f;

	InitializeShapesForWaterElevatorConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::InitializeShapesForWaterElevatorConfig()
{
	m_waterElevatorConfigDisc1 = new Disc2D();
	m_waterElevatorConfigDisc1->m_center = Vec2(0.75f, 1.125f);
	m_waterElevatorConfigDisc1->m_radius = 0.1f;
	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc1->m_center, m_waterElevatorConfigDisc1->m_radius, Rgba8::WHITE);

	m_waterElevatorConfigDisc2 = new Disc2D();
	m_waterElevatorConfigDisc2->m_center = Vec2(1.5f, 1.125f);
	m_waterElevatorConfigDisc2->m_radius = 0.1f;
	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc2->m_center, m_waterElevatorConfigDisc2->m_radius, Rgba8::WHITE);

	m_waterElevatorConfigDisc3 = new Disc2D();
	m_waterElevatorConfigDisc3->m_center = Vec2(2.25f, 1.125f);
	m_waterElevatorConfigDisc3->m_radius = 0.1f;
	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc3->m_center, m_waterElevatorConfigDisc3->m_radius, Rgba8::WHITE);

	m_waterElevatorAABB2 = new AABB2(Vec2(0.01f, -0.5f), Vec2(2.99f, -0.4f));
	m_tumblingConfigMovingAABBVelocity = Vec2(0.0f, 0.5f) * 0.2f;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::RenderShapesForWaterElevatorConfig()
{
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->DrawVertexArray((int)m_waterElevatorConfigShapesVerts.size(), m_waterElevatorConfigShapesVerts.data());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateShapesForWaterElevatorConfig(float deltaSeconds)
{
	if (m_waterElevatorAABB2)
	{
		Vec2 translation = Vec2(0.f, m_tumblingConfigMovingAABBVelocity.y * deltaSeconds);
		m_waterElevatorAABB2->Translate(translation);

		float bottomBoundary = -0.4f;
		float topBoundary = MAP_HEIGHT - 0.6f;
		float upwardVelocity = 0.2f;  
		float downwardVelocity = 0.4f; 

		if (m_waterElevatorAABB2->m_mins.y < bottomBoundary || m_waterElevatorAABB2->m_maxs.y > topBoundary)
		{
			if (m_tumblingConfigMovingAABBVelocity.y > 0)  
			{
				m_tumblingConfigMovingAABBVelocity.y = -downwardVelocity; 
			}
			else 
			{
				m_tumblingConfigMovingAABBVelocity.y = upwardVelocity; 
			}
		}

		if (m_waterElevatorAABB2->m_mins.y < bottomBoundary)
		{
			float overlap = bottomBoundary - m_waterElevatorAABB2->m_mins.y;
			m_waterElevatorAABB2->Translate(Vec2(0.f, overlap));
		}
		else if (m_waterElevatorAABB2->m_maxs.y > topBoundary)
		{
			float overlap = m_waterElevatorAABB2->m_maxs.y - topBoundary;
			m_waterElevatorAABB2->Translate(Vec2(0.f, -overlap));
		}

		m_waterElevatorConfigShapesVerts.clear();
		AddVertsForAABB2D(m_waterElevatorConfigShapesVerts, *m_waterElevatorAABB2, Rgba8::WHITE);
	}

	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc1->m_center, m_waterElevatorConfigDisc1->m_radius, Rgba8::WHITE);

	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc2->m_center, m_waterElevatorConfigDisc2->m_radius, Rgba8::WHITE);

	AddVertsForDisc2D(m_waterElevatorConfigShapesVerts, m_waterElevatorConfigDisc3->m_center, m_waterElevatorConfigDisc3->m_radius, Rgba8::WHITE);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSim::UpdateSimulationConstantBufferForWaterElevatorConfig(SimulationConstants& constants)
{
	if (m_waterElevatorAABB2)
	{
		constants.waterElevatorConfigVerticalMovingAABB1Size = Vec2(m_waterElevatorAABB2->m_maxs.x - m_waterElevatorAABB2->m_mins.x, m_waterElevatorAABB2->m_maxs.y - m_waterElevatorAABB2->m_mins.y);
		constants.waterElevatorConfigVerticalMovingAABB1Center = Vec2((m_waterElevatorAABB2->m_maxs.x + m_waterElevatorAABB2->m_mins.x) / 2, (m_waterElevatorAABB2->m_maxs.y + m_waterElevatorAABB2->m_mins.y) / 2);
	}

	else
	{
		constants.waterElevatorConfigVerticalMovingAABB1Size = Vec2(0.f, 0.f);
		constants.waterElevatorConfigVerticalMovingAABB1Center = Vec2(0.f, 0.f);
	}

	if (m_waterElevatorConfigDisc1)
	{
		constants.waterElevatorDisc1Center = m_waterElevatorConfigDisc1->m_center;
		constants.waterElevatorDisc1Radius = m_waterElevatorConfigDisc1->m_radius;
	}
	else
	{
		constants.waterElevatorDisc1Center = Vec2(0.f, 0.f);
		constants.waterElevatorDisc1Radius = 0.f;
	}

	if (m_waterElevatorConfigDisc2)
	{
		constants.waterElevatorDisc2Center = m_waterElevatorConfigDisc2->m_center;
		constants.waterElevatorDisc2Radius = m_waterElevatorConfigDisc2->m_radius;
	}
	else
	{
		constants.waterElevatorDisc2Center = Vec2(0.f, 0.f);
		constants.waterElevatorDisc2Radius = 0.f;
	}

	if (m_waterElevatorConfigDisc3)
	{
		constants.waterElevatorDisc3Center = m_waterElevatorConfigDisc3->m_center;
		constants.waterElevatorDisc3Radius = m_waterElevatorConfigDisc3->m_radius;
	}
	else
	{
		constants.waterElevatorDisc3Center = Vec2(0.f, 0.f);
		constants.waterElevatorDisc3Radius = 0.f;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
