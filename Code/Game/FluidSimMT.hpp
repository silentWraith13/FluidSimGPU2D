#pragma once
#include "Game/Game.hpp"
#include "Game/FluidSim.hpp"
#include "Engine/Core/JobSystem.hpp"
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <cmath>
#include <unordered_map>

//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum CPUMode
{
	SINGLE_THREADED,
	MULTI_THREADED
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct Solver
{
	Vec2 m_particlePosition;
	Vec2 m_particleVelocity;
	Vec2 m_particleAcceleration;
	float m_particleDensity;
	float padding;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
const static float DENSITY_COEFFICIENT = 0.0002f * static_cast<float>(315.0f / (64.0f * M_PI * powf(SMOOTHING_LENGTH, 9)));
const static float GRADUAL_PRESSURE_COEFFICIENT = 0.0002f * static_cast<float>(-45.0f / (M_PI * powf(SMOOTHING_LENGTH, 6)));
const static float LAPLACIAN_VISCOSITY_COEFFICIENT = 0.0002f * 0.1f * static_cast <float>(45.0f / (M_PI * powf(SMOOTHING_LENGTH, 6)));
const static Vec3  LEFT_MAP_BOUNDS = Vec3(1, 0, 0);
const static Vec3  BOTTOM_MAP_BOUNDS = Vec3(0, 1, 0);
const static Vec3  RIGHT_MAP_BOUNDS = Vec3(-1, 0, MAP_WIDTH / 2);
const static Vec3  TOP_MAP_BOUNDS = Vec3(0, -1, MAP_HEIGHT / 2);
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class FluidSimMT
{
public:
	FluidSimMT(); //Constructor
	~FluidSimMT(); //Destructor

	//Basic Frame work
	void				 Startup(unsigned int numParticles, InitialParticleConfigurationCPU initialParticleConfig);
	void				 Render();
	void				 Update(float deltaSeconds);
	void				 InitializeWorldCamera();
	void				 GetMouseCursorPos();
	
	//Particle functions
	void				 RenderParticles();
	void				 DensitySimpleMultiThreaded();
	void				 DensityGridMultiThreaded();
	void				 DensitySimpleSingleThreaded();
	void				 DensityGridSingleThreaded();
	float				 CalculateDensity(float r_sq);
	float				 CalculatePressure(float density);
	Vec2				 CalculateGradPressure(float r, float particlePressure, float N_pressure, float neighborDensity, Vec2 posDiff);
	Vec2				 CalculateLapVelocity(float r, Vec2 particleVelocity, Vec2 neighborVelocity, float neighborDensity);
	void				 ForceSimpleMultiThreaded();
	void				 ForcesGridMultiThreaded();
	void				 ForceSimpleSingleThreaded();
	void				 ForcesGridSingleThreaded();
	void				 IntegrateMultiThreadedSimulation(float deltaSeconds);
	void				 IntegrateSingleThreadedSimulation(float deltaSeconds);
	void				 UpdatePhysics(float deltaSeconds);
	Vec2				 ApplyMouseForces(Vec2 mousePos, float forceRadius, float forceStrength, Vec2 const& position, Vec2 const& velocity);
	void				 UpdateFluidMesh();
	void				 SpawnParticlesAtMousePos();

	//Hash mapping functions
	int					 CalculateHashKey(const Vec2& position) const;
	void				 AddParticleToGrid(Solver* particle);
	std::vector<Solver*> GetParticlesFromGrid(const Vec2& position);
	void				 UpdateGrid();
	std::vector<Solver*> FindNeighbors(Solver* particle) const;
	
	//ImGui
	void				 ImGuiParticleSize();
	void				 ImGuiModeText();
	void				 ImGuiSimulationMode();
	void				 ImGUiNumberOfParticles();
	void				 ImGuiGravityButtons();
	void				 ImGuiMillisecondText();
	void				 ImGuiPhysicsTimeStep();

	//DAM config
	void				 InitializeParticlesInADamConfig(unsigned int numParticles);
	
	//Random Config
	void				 InitializeParticlesInARandomShapeConfig(unsigned int numParticles);
	void				 InitializeShapesForRandomShapeConfig();
	void				 RenderShapesForRandomConfig();
	void				 ProcessCollisionsForRandomShapeConfig();


	//Finding random positions in world
	Vec2 const			 GetClampedRandomPositionInWorld();
	Vec2 const			 GetClampedRandomPositionInWorld(float shapeRadius);
	Vec2 const			 GetRandomPositionInWorld();

	//Member variables
	std::vector<Solver*>				m_solverparticles;
	Camera								m_worldCamera;
	int									m_numParticles = -1;
	std::vector<Vertex_PCU>				m_solverVerts;
	VertexBuffer*						m_solverVertexBuffer = nullptr;
	float								m_frameTimeMs = 0.f;
	float								m_viscosity = 0.2f;
	float								m_restDensity = 1000.f;
	float								m_pressureStiffness = 200.0f;
	float								m_particleMass = 0.0002f;
	float								m_maxAllowableTimeStep = 0.005f;
	float								m_smoothingLength = 0.012f;
	float								m_wallStiffness = 3000.0f;
	float								m_lastSimulationTimeMs = 0.f;
	float								m_lastRenderTimeMs = 0.f;
	float								m_particleSize = 0.005f;
	float								m_initialParticleSpacing = 0.0065f;
	float								m_cellSize = -1.f;
	bool								m_isBufferDirty = false;
	Texture*							m_texture = nullptr;
	std::map<int, std::vector<Solver*>> m_grid;
	SimulationMode						m_simMode = GRID;
	Vec2								m_cursorPos;
	AABB2								m_mapBounds;
	Vec2								m_gravityDirection = GRAVITY_DOWN;
	CPUMode								m_mode =  MULTI_THREADED;
	std::vector<Vertex_PCU>				m_randomShapeConfigVerts;
	std::vector<OBB2D*>					m_randomShapesConfigOBBVector;
	std::vector<Capsule2*>				m_randomShapesConfigCapsuleVector;
	RandomNumberGenerator				m_rng;
	std::vector<Vertex_PCU>				m_windmillConfigVerts;
	std::vector<OBB2D*>					m_windmillConfigOBBVector;
	Disc2D*								m_windmillConfigDisc;
	std::vector<Vertex_PCU>				m_waterElevatorConfigShapesVerts;
	std::vector<Vertex_PCU>				m_waterElevatorConfigDiscShapesVerts;
	std::vector<Disc2D*>				m_waterElevatorConfigDiscVector;
	AABB2*								m_waterElevatorAABB2 = nullptr;
	Vec2								m_tumblingConfigMovingAABBVelocity = Vec2(-1.f, -1.f);
};

//--------------------------------------------------------------------------------------------------------------------------------------------------------
// CLASS DERIVED BY THE JOB CLASS TO MULTI_THREAD
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class CalculateDensityJobSimple : public Job
{
public:
	CalculateDensityJobSimple(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq );

	virtual void Execute() override;
	virtual void OnFinished() override;

	FluidSimMT*			 m_sim;
	float				 m_hSq;
	std::vector<Solver*> m_particleSubset;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class CalculateDensityJobGrid : public Job
{
public:
	CalculateDensityJobGrid(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq);

	virtual void Execute() override;
	virtual void OnFinished() override;

	FluidSimMT*			 m_sim;
	float			     m_hSq;
	std::vector<Solver*> m_particleSubset;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class CalculateForcesJobSimple : public Job
{
public:
	CalculateForcesJobSimple(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq);

	virtual void Execute() override;
	virtual void OnFinished() override;

	FluidSimMT* m_sim;
	float				 m_hSq;
	std::vector<Solver*> m_particleSubset;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class CalculateForcesJobGrid : public Job
{
public:
	CalculateForcesJobGrid(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq);

	virtual void Execute() override;
	virtual void OnFinished() override;

	FluidSimMT*			 m_sim;
	float				 m_hSq;
	std::vector<Solver*> m_particleSubset;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class IntegrateParticlesJob : public Job 
{
public:
	IntegrateParticlesJob(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float deltaSeconds, const Vec3& leftmapBounds, const Vec3& topmapBounds, const Vec3& rightmapBounds, const Vec3& bottommapBounds, float wallstiffness);

	virtual void Execute() override;
	virtual void OnFinished() override;

private:
	std::vector<Solver*> m_particleSubset;
	float				 m_deltaSeconds;
	Vec3				 m_leftMapBounds;
	Vec3				 m_topMapBounds;
	Vec3				 m_rightMapBounds;
	Vec3				 m_bottomMapBounds;
	float				 m_wallStiffness;
	FluidSimMT* m_fluidSim;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------