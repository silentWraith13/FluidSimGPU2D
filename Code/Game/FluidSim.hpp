#pragma once
#include "Game/GameCommon.hpp"
#include "Engine/Renderer/RenderTexture.hpp"
#include "Engine/Math/Capsule2.hpp"
#include "Engine/Math/Disc2D.hpp"
#include "Game/Game.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------
class Camera;
class Shader;
class UnorderedAccessBuffer;
class ConstantBuffer;
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct Particle
{
	Vec2 m_position;
	Vec2 m_velocity;

	float type;
	float padding[3];
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct ParticleDensity
{
	float m_particleDensity;
	float padding[3];
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct ParticleForces
{
	Vec2 m_particleAcceleration;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct ParticleSize
{
	float particleSize;
	int   debugType;
	float padding[2];
};
static const int k_particleSizeConstantsSlot = 5;
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct ScreenSizeConstants 
{
	float screenWidth;
	float screenHeight;
	float time;
	float padding;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
struct Sort
{
	unsigned int Level;
	unsigned int LevelMask;
	unsigned int iWidth;
	unsigned int iHeight;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
 struct SimulationConstants
{
	unsigned int numFluidParticles; 
	float		 timeStep; 
	float		 smoothingLength; 
	float		 pressureStiffness; 
	
	float		 restDensity; 
	float		 densityCoefficient; 
	float		 gradPressureCoefficient; 
	float		 laplcaianViscosityCoefficient; 
	
	float		 surfaceTensionCoefficient;
	float		 wallStiffness; 
	float		 padding[2]; 
	
	Vec2		 vGravity; 
	float		 gravityPadding[2];

	Vec4		 vGridDim; 

	Vec3		 leftMapBounds;
	float		 paddingAfterPlanes1;  

	Vec3		 bottomMapBounds;
	float		 paddingAfterPlanes2;
	
	Vec3		 rightMapBounds;
	float		 paddingAfterPlanes3;  

	Vec3		topMapBounds;
	float		paddingAfterPlanes4;   

	Vec2		mousePos;
	int			isHoldingLeftClick;
	int			isHoldingRightClick;

	Vec2        tumblingConfigOBB1HalfDimensions;
	Vec2        tumblingConfigOBB1Center;
	
	Mat44		tumblingConfigOBB1LocalToWorldMatrix;
	Mat44		tumblingConfigOBB1WorldToLocalMatrix;

	Vec2        tumblingConfigOBB2HalfDimensions;
	Vec2        tumblingConfigOBB2Center;
	
	Mat44		tumblingConfigOBB2LocalToWorldMatrix;
	Mat44		tumblingConfigOBB2WorldToLocalMatrix;

	Vec2        tumblingConfigOBB3HalfDimensions;
	Vec2        tumblingConfigOBB3Center;
	
	Mat44		tumblingConfigOBB3LocalToWorldMatrix;
	Mat44		tumblingConfigOBB3WorldToLocalMatrix;
	
	Vec2        tumblingConfigOBB4HalfDimensions;
	Vec2        tumblingConfigOBB4Center;
	
	Mat44		tumblingConfigOBB4LocalToWorldMatrix;
	Mat44		tumblingConfigOBB4WorldToLocalMatrix;
	
	Vec2        tumblingConfigOBB5HalfDimensions;
	Vec2        tumblingConfigOBB5Center;
	
	Mat44		tumblingConfigOBB5LocalToWorldMatrix;
	Mat44		tumblingConfigOBB5WorldToLocalMatrix;

	Vec2        tumblingConfigCapsule1Start;
	Vec2        tumblingConfigCapsule1End;
	
 	Vec2        tumblingConfigCapsule2Start;
 	Vec2        tumblingConfigCapsule2End;
	
	float       tumblingConfigCapsule1Radius;
 	float       tumblingConfigCapsule2Radius;
 	float		tumblingCapsule2padding[2];

	Vec2		tumblingConfigHorizontalMovingAABB1Size;
	Vec2		tumblingConfigHorizontalMovingAABB1Center;

	Vec2		tumblingDiscCenter;
	float		tumblingDiscRadius;
	float		tumblingpadding;

	Vec2		windmillDisc1Center;
	float		windmillDisc1Radius;
	float		windmillDiscPadding;

	Vec2		windmillConfigOBB1HalfDimensions;
	Vec2		windmillConfigOBB1Center;
	
	Mat44		windmillConfigOBB1LocalToWorldMatrix;
	Mat44		windmillConfigOBB1WorldToLocalMatrix;

	Vec2		windmillConfigOBB2HalfDimensions;
	Vec2		windmillConfigOBB2Center;

	Mat44		windmillConfigOBB2LocalToWorldMatrix;
	Mat44		windmillConfigOBB2WorldToLocalMatrix;

	Vec2		windmillConfigBottomBorderCapsuleStart;
	Vec2		windmillConfigBottomBorderCapsuleEnd;
	
	float		windmillConfigBottomBorderCapsuleRadius;
	float		windmillConfigBottomBorderCapsulePadding[3];

	Vec2		windmillConfigLeftBorderCapsuleStart;
	Vec2		windmillConfigLeftBorderCapsuleEnd;

	float		windmillConfigLeftBorderCapsuleRadius;
	float		windmillConfigLeftBorderCapsulePadding[3];

	Vec2		windmillConfigRightBorderCapsuleStart;
	Vec2		windmillConfigRightBorderCapsuleEnd;

	float		windmillConfigRightBorderCapsuleRadius;
	float		windmillConfigRightBorderCapsulePadding[3];

	Vec2		windmillConfigTopBorderCapsule1Start;
	Vec2		windmillConfigTopBorderCapsule1End;

	float		windmillConfigTopBorderCapsule1Radius;
	float		windmillConfigTopBorderCapsule1Padding[3];

	Vec2		windmillConfigTopBorderCapsule2Start;
	Vec2		windmillConfigTopBorderCapsule2End;
	
	float		windmillConfigTopBorderCapsule2Radius;
	float		windmillConfigTopBorderCapsule2Padding[3];

	Vec2		windmillConfigTopBorderCapsule3Start;
	Vec2		windmillConfigTopBorderCapsule3End;

	float		windmillConfigTopBorderCapsule3Radius;
	float		windmillConfigTopBorderCapsule3Padding[3];

	Vec2		hourglassConfigCapsule1Start;
	Vec2		hourglassConfigCapsule1End;

	float		hourglassConfigCapsule1Radius;
	float		hourglassConfigCapsule1Padding[3];

	Vec2		hourglassConfigCapsule2Start;
	Vec2		hourglassConfigCapsule2End;

	float		hourglassConfigCapsule2Radius;
	float		hourglassConfigCapsule2Padding[3];

	Vec2		hourglassConfigCapsule3Start;
	Vec2		hourglassConfigCapsule3End;

	float		hourglassConfigCapsule3Radius;
	float		hourglassConfigCapsule3Padding[3];

	Vec2		hourglassConfigCapsule4Start;
	Vec2		hourglassConfigCapsule4End;

	float		hourglassConfigCapsule4Radius;
	float		hourglassConfigCapsule4Padding[3];

	Vec2		hourglassConfigCapsule5Start;
	Vec2		hourglassConfigCapsule5End;

	float		hourglassConfigCapsule5Radius;
	float		hourglassConfigCapsule5Padding[3];

	Vec2		hourglassConfigCapsule6Start;
	Vec2		hourglassConfigCapsule6End;

	float		hourglassConfigCapsule6Radius;
	float		hourglassConfigCapsule6Padding[3];

	Vec2		hourglassConfigMiddleCapsuleStart1;
	Vec2		hourglassConfigMiddleCapsuleEnd1;

	float		hourglassConfigMiddleCapsule1Radius;
	float		hourglassConfigMiddleCapsule1Padding[3];

	Vec2		hourglassConfigMiddleCapsuleStart2;
	Vec2		hourglassConfigMiddleCapsuleEnd2;

	float		hourglassConfigMiddleCapsule2Radius;
	float		hourglassConfigMiddleCapsule2Padding[3];

	Vec2		waterElevatorConfigVerticalMovingAABB1Size;
	Vec2		waterElevatorConfigVerticalMovingAABB1Center;

	Vec2		waterElevatorDisc1Center;
	float		waterElevatorDisc1Radius;
	float		waterElevator1padding;

	Vec2		waterElevatorDisc2Center;
	float		waterElevatorDisc2Radius;
	float		waterElevator2padding;

	Vec2		waterElevatorDisc3Center;
	float		waterElevatorDisc3Radius;
	float		waterElevator3padding;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
 struct UINT2
 {
	 unsigned int x;
	 unsigned int y;
 };
 //--------------------------------------------------------------------------------------------------------------------------------------------------------
const static float		  MAP_HEIGHT = 1.5f;
const static float		  MAP_WIDTH = 2 * MAP_HEIGHT;
const static float		  SMOOTHING_LENGTH = 0.012f;
const static unsigned int NUM_GRID_INDICES = 65536;
const static unsigned int SIMULATION_BLOCK_SIZE = 256;
const static unsigned int BITONIC_BLOCK_SIZE = 512;
const static unsigned int TRANSPOSE_BLOCK_SIZE = 16;

const static Vec2		  GRAVITY_DOWN(0, -0.5f);
const static Vec2		  GRAVITY_DOWN_FAST(0, -5.5f);
const static Vec2		  GRAVITY_UP(0, 0.5f);
const static Vec2		  GRAVITY_LEFT(-0.5f, 0);
const static Vec2		  GRAVITY_RIGHT(0.5f, 0.f);
const static Vec2		  GRAVITY = GRAVITY_DOWN;

const static Vec3		  MAP_PLANES[4] =
{
	Vec3(1, 0, 0),
	Vec3(0, 1, 0),
	Vec3(-1, 0, MAP_WIDTH),
	Vec3(0, -1, MAP_HEIGHT)
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum SimulationMode
{
	SIMPLE,
	GRID,
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
enum ParticleRenderType
{
	DENSITY_DEBUG,
	BLUE_PALLETE,
	CARTOON,
	TRANSPARENT,
	COUNT
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------
class FluidSim
{
public:
	FluidSim();
	~FluidSim();

	//Basic
	void Startup(unsigned int numParticles, InitialParticleConfigurationGPU initialParticleConfig);
	void Update(float deltaSeconds);
	void Render();
	void Shutdown();
	void EndFrame();

	//Camera
	void InitializeWorldCamera();

	void HandleKeyboardInputs();

	//Fluid Render
	void InitializeParticleTexturesAndShaders();
	void InitializeParticleRenderingShader();
	void RenderFluidParticles();
	
	//Creation of all Constant buffers
	void CreateParticleSizeConstantBuffer();
	void CreateSimulationConstantBuffer();
	void CreateSortingConstantBuffer();

	//Binding of all Constant buffers
	void BindCBForParticleSize();
	void BindConstantBufferForSimulationConstants(float deltaSeconds);

	//Simulation of fluid(simple and grid)
	void SimulateFluid(float deltaSeconds);
	void SimulateFluidSimple();
	void SimulateFluidGrid();

	//Creation of compute shaders
	void CreateComputeShaders();

	//UI
	void ImGuiGravityButtons();
	void ImGuiFpsText();
	void ImGuiSimulationMode();
	void ImGuiRenderTypeForParticles();
	void ImGuiViscosity();
	void ImGuiRestDensity();
	void ImGuiPressureStiffness();
	void ImGuiWallStiffness();
	void ImGuiSmoothingLength();
	void ImGuiTimestep();
	void ImGuiParticleMass();
	void ImGuiParticleRenderSize();
	void ImGuiAddGeometry();
	void ImGuiChangeBackgroundColor();
	void ImGuiMillisecondText();
	void ImGuiRenderFunctions();

	//Extra GPU functions
	void UnsetAllResources();
	void GPUSort(UnorderedAccessBuffer* inUAV, UnorderedAccessBuffer* tempUAV);

	//Some mouse related functions
	Vec2 GetMouseCursorPos();
	Vec2 GetMouseCursorPosition();
	Shape2D* GetShapeUnderMouse();

	//For post-process effects
	void RenderTextureToScreenForCartoonFluid();
	void RenderTextureToScreenForTransparentFluid();

	//Dam config
	void InitializeParticlesInADam(unsigned int numParticles);

	//Random Config
	void InitializeParticlesInARandomShapeConfig(unsigned int numParticles);
	void InitializeShapesForRandomShapesConfig();
	void UpdateMovingAABB2ForRandomShapesConfig(float deltaSeconds);
	void RenderMovingAABB2ForRandomShapesConfig();
	void RenderShapesForRandomShapesConfig();
	Mat44 CreateWorldToLocalMatrix(Mat44& localToWorldMatrix);
	Mat44 CreateLocalToWorldMatrix(OBB2D& box, Vec2& obbCenter);
	void  UpdateSimulationConstantBufferFromRandomShapesConfig(SimulationConstants& constants);

	//Windmill config
	void InitializeShapesForWindmillConfig();
	void RenderShapesForWindwillConfig();
	void UpdateSimulationConstantBufferForWindmillConfig(SimulationConstants& constants);
	void UpdateRotationOfWindmillMesh(float deltaSeconds);
	void InitializeParticlesInAWindmillConfig(unsigned int numParticles);
	
	//Hourglass Config
	void InitializeParticlesInAnHourglassConfig(unsigned int numParticles);
	void InitializeShapesForHourglassConfig();
	void RenderShapesForHourglassConfig();
	void UpdateSimulationConstantBufferForHourglassConfig(SimulationConstants& constants);
	Mat44 GetModelMatrixForHourglassConfig(float deltaSeconds);
	Vec2 GetHourGlassMeshCenter();
	void UpdateHourglassMeshConfig(float deltaSeconds);
	void TansformAndAddVertsForHourglassConfigShape(Capsule2* capsule, const Mat44& matrix, bool isTopCapsule = false);

	//Transparent fluid config
	void InitializeParticlesForTransparentFluidConfig(unsigned int numParticles);
	void InitializeBackgroundTextureForTransparentFluidConfig();
	void RenderBackgroundTextureForTransparentFluidConfig();
	void RenderTransparentFluidForTheTransparentConfig();
	
	//Finding random positions in world
	Vec2 const GetClampedRandomPositionInWorld();
	Vec2 const GetClampedRandomPositionInWorld(float shapeRadius);
	Vec2 const GetRandomPositionInWorld();

	//Render shapes based on Config
	void RenderShapesBasedOnConfig();
	void RenderCartoonFluid();

	//Water Elevator Config
	void InitializeParticlesForWaterElevatorConfig(unsigned int numParticles);
	void InitializeShapesForWaterElevatorConfig();		
	void RenderShapesForWaterElevatorConfig();
	void UpdateShapesForWaterElevatorConfig(float deltaSeconds);
	void UpdateSimulationConstantBufferForWaterElevatorConfig(SimulationConstants& constants);

	//Member variables
	Camera					m_worldCamera;
	
	UnorderedAccessBuffer*	m_particleUAV = nullptr;
	UnorderedAccessBuffer*	m_nullUAV = nullptr;
	UnorderedAccessBuffer*	m_particleDensityUAV = nullptr;
	UnorderedAccessBuffer*	m_particleForcesUAV = nullptr;
	UnorderedAccessBuffer*	m_sortedParticlesUAV = nullptr;
	UnorderedAccessBuffer*	m_gridUAV = nullptr;
	UnorderedAccessBuffer*	m_gridPingPongUAV = nullptr;
	UnorderedAccessBuffer*	m_gridIndicesUAV = nullptr;
	
	ConstantBuffer*			m_particleSizeCBO = nullptr;
	ConstantBuffer*			m_simulationConstantsCBO = nullptr;
	ConstantBuffer*			m_sortCBO = nullptr;
	ConstantBuffer*			m_screenSizeCBO = nullptr;
	
	RenderTexture*			m_postProcessRenderTexture = nullptr;
	RenderTexture*			m_particleDensityRenderTexture = nullptr;
	
	Texture*				m_particleTexture = nullptr;
	Shape2D*				m_selectedShape = nullptr;
	Shader*					m_postProcessShader = nullptr;
	Shader*					m_renderingShader = nullptr;
	Shader*					m_integrateCS = nullptr;
	Shader*					m_densitySimpleCS = nullptr;
	Shader*					m_forceSimpleCS = nullptr;
	Shader*					m_densityGridCS = nullptr;
	Shader*					m_forcesGridCS = nullptr;
	Shader*					m_buildGridCS = nullptr;
	Shader*					m_ClearGridIndicesGridCS = nullptr;
	Shader*					m_buildGridIndicesCS = nullptr;
	Shader*					m_bitonicSortCS = nullptr;
	Shader*					m_matrixTransposeCS = nullptr;
	Shader*					m_rearrangeParticlesCS = nullptr;
	Shader*					m_densityShader = nullptr;
	Shader*					m_transparentFluidShader = nullptr;
	SimulationMode			m_simMode = GRID;
	
	Vec2					m_cursorPos;
	Vec2					m_gravityDirection = GRAVITY;
	Vec2					m_prevMousePosition = Vec2(0.f, 0.f);
	unsigned int			m_nullUINT = 0;
	unsigned int			m_numParticles = 0;
	
	int						m_isHoldingLeftClick = 0;
	int						m_isHoldingRightClick = 0;
	int						m_renderType = -1;
	
	float					m_viscosity = 0.1f;
	float					m_restDensity = 1000.f;
	float					m_pressureStiffness = 200.0f;
	float					m_particleMass = 0.0002f;
	float					m_maxAllowableTimeStep = 0.005f;
	float					m_smoothingLength = SMOOTHING_LENGTH;
	float					m_wallStiffness = 3000.0f;
	float					m_particleSize = 0.04f;
	float					m_lastSimulationTimeMs = 0.f;
	float					m_lastRenderTimeMs = 0.f;
	float					m_initialParticleSpacing = 0.0045f;
	bool					m_isMovingShape = false;
	
	std::vector<Vertex_PCU> m_screenVerts;
	AABB2					m_screenBounds = AABB2( Vec2(0.f, 0.f), Vec2(MAP_WIDTH, MAP_HEIGHT));
	ParticleRenderType		m_particleRenderType = COUNT;
	Rgba8					m_clearColor = Rgba8::BLACK;

	//Tumbling Config
	std::vector<Shape2D*>   m_tumblingConfigShapes;
	std::vector<Vertex_PCU> m_tumblingConfigShapesVerts;
	std::vector<Vertex_PCU> m_tumblingConfigMovingAABB2Verts;
	RandomNumberGenerator   m_rng;
	Vec2					m_tumblingConfigMovingAABBVelocity;
	AABB2*					m_tumblingConfigMovingAABB2 = nullptr;
	bool					m_showTumblingConfigMovingAABB = false;
	Disc2D*					m_tumblingDisc2D = nullptr;

	//Windmill config
	OBB2D*					m_windmillOBB1 = nullptr;
	OBB2D*					m_windmillOBB2 = nullptr;
	Disc2D*					m_windmillDisc1 = nullptr;
	std::vector<Vertex_PCU> m_windmillConfigShapesVerts;
	Capsule2*				m_windmillLeftBorderCapsule = nullptr;
	Capsule2*				m_windmillRightBorderCapsule = nullptr;
	Capsule2*				m_windmillBottomBorderCapsule = nullptr;
	Capsule2*				m_windmillTopBorderCapsule1 = nullptr;
	Capsule2*				m_windmillTopBorderCapsule2 = nullptr;
	Capsule2*				m_windmillTopBorderCapsule3 = nullptr;

	//Hourglass config
	std::vector<Vertex_PCU> m_hourglassConfigShapesVerts;
	std::vector<Vertex_PCU> m_hourglassTopCapsuleVerts;
	Capsule2*				m_hourglassConfigCapsule1 = nullptr;
	Capsule2*				m_hourglassConfigCapsule2 = nullptr;
	Capsule2*				m_hourglassConfigCapsule3 = nullptr;
	Capsule2*				m_hourglassConfigCapsule4 = nullptr;
	Capsule2*				m_hourglassConfigCapsule5 = nullptr;
	Capsule2*				m_hourglassConfigCapsule6 = nullptr;
	Capsule2*				m_hourglassConfigCapsule7 = nullptr;
	Capsule2*				m_hourglassConfigMiddleCapsule1 = nullptr;
	Capsule2*				m_hourglassConfigMiddleCapsule2 = nullptr;
	float					m_totalRotationDegrees = 0.0f;
	bool					m_shouldRotateHourglassMesh = false;
	bool					m_renderTheTopCapsuleForHourglassMesh = false;

	//Transparent fluid Config
	Texture*				m_backgroundTexture = nullptr;
	std::vector<Vertex_PCU> m_backgroundQuadForTransparentFluidConfigVerts;
	RenderTexture*			m_transparentFluidRenderTexture = nullptr;
	Texture*				m_noiseTexture = nullptr;

	//Water Elevator Config
	Disc2D*					m_waterElevatorConfigDisc1 = nullptr;
	Disc2D*					m_waterElevatorConfigDisc2 = nullptr;
	Disc2D*					m_waterElevatorConfigDisc3 = nullptr;
	AABB2*					m_waterElevatorAABB2 = nullptr;
	std::vector<Vertex_PCU> m_waterElevatorConfigShapesVerts;
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------