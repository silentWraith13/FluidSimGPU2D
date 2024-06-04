#pragma once
#include "Game/GameCommon.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------
class Camera;
class FluidSim;
class FluidSimMT;
//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum SimulationThread
{
	SINGLE_AND_MULTI_THREAD,
	GPU,
	NUM_SIM_MODES
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum InitialParticleConfigurationGPU
{
	DAM,
	RANDOM_SHAPES,
	WINDMILL,
	HOURGLASS,
	TRANSPARENT_FLUID,
	WATER_ELEVATOR,
	NUM_MODES_CONFIG_GPU
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum InitialParticleConfigurationCPU
{
	CPU_DAM,
	RANDOM_SHAPES_CPU,
	NUM_MODES_CONFIG_CPU
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
const static unsigned int   NUM_PARTICLES_4K = 4 * 1024;
const static unsigned int   NUM_PARTICLES_8K = 8 * 1024;
const static unsigned int   NUM_PARTICLES_16K = 16 * 1024;
const static unsigned int   NUM_PARTICLES_32K = 32 * 1024;
const static unsigned int   NUM_PARTICLES_64K = 64 * 1024;
const static unsigned int   NUM_PARTICLES_128K = 128 * 1024;
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class Game
{
public:
	Game();
	~Game();
	void                Startup();
	void                ShutDown();
	void                Update(float deltaSeconds);
	void                Render();
	void				EndFrame();
	void				RenderImGuiUI();
	void				ImGuiResetSim();
	void				ImGuiChangeNumberOfParticlesGPU();
	void				ImGuiInitialParticleConfigGPU();
	void				ImGuiFpsText();
	void				ImGuiInitialParticleConfigCPU();

	void				SwitchToNextMode();
	void				ToggleModes();
	void				SwitchToPreviousMode();
	void				DeleteCurrentGameMode();
	void				CreateCurrentGameMode();
	void				RandomizeGameMode();
	
	void				InitializeGPUFluidSimMode();
	void				InitializeCPUFluidSimMode();

	//Member variables
	AABB2							m_camBounds;
	Camera							m_screenCamera;
	FluidSim*						m_fluidSim = nullptr;
	FluidSimMT*						m_fluidSimCPU = nullptr;
	InitialParticleConfigurationGPU m_initialParticleConfigGPU = NUM_MODES_CONFIG_GPU;
	InitialParticleConfigurationCPU m_initialParticleConfigCPU = NUM_MODES_CONFIG_CPU;
	unsigned int					m_numParticlesGPU = NUM_PARTICLES_64K;
	unsigned int					m_numParticlesCPU = 8000;
	SimulationThread				m_gameMode = SINGLE_AND_MULTI_THREAD;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------