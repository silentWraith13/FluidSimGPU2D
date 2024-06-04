#include "Game/FluidSimMT.hpp"
#include "Engine/Renderer/BitmapFont.hpp"
#include "ThirdParty/ImGUI/imgui.h"
#include "Engine/Renderer/VertexBuffer.hpp"
#include "Engine/Math/MathUtils.hpp"
#include <algorithm>
#include <chrono>
//--------------------------------------------------------------------------------------------------------------------------------------------------------
FluidSimMT::FluidSimMT()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
FluidSimMT::~FluidSimMT()
{
	delete m_solverVertexBuffer;
	m_solverVertexBuffer = nullptr;

	m_solverVerts.clear();
	m_solverparticles.clear();

	g_theJobSystem->ClearQueuedJobs();
	g_theJobSystem->WaitUntilCurrentJobsCompletion();
	g_theJobSystem->ClearCompletedJobs();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::Startup(unsigned int numParticles, InitialParticleConfigurationCPU initialParticleConfig)
{
	InitializeWorldCamera();

	m_texture = g_theRenderer->CreateOrGetTextureFromFile("Data/Images/Disc2.png");
	m_cellSize = m_smoothingLength;
	m_mapBounds = AABB2(Vec2(0.f, 0.f), Vec2(RIGHT_MAP_BOUNDS.z, TOP_MAP_BOUNDS.z));
	m_rng = RandomNumberGenerator((unsigned int)std::chrono::system_clock::now().time_since_epoch().count());
	m_numParticles = numParticles;
	
	if (initialParticleConfig == RANDOM_SHAPES_CPU)
	{
		InitializeParticlesInARandomShapeConfig(numParticles);
	}

	if (initialParticleConfig == CPU_DAM)
	{
		InitializeParticlesInADamConfig(numParticles);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::Render()
{
	g_theRenderer->ClearScreen(Rgba8::BLACK);

	g_theRenderer->BeginCamera(m_worldCamera);

 	if (g_theGame->m_initialParticleConfigCPU == RANDOM_SHAPES_CPU)
 	{
 		RenderShapesForRandomConfig();
 	}

	RenderParticles();

	//IMGUI
	ImGuiSimulationMode();
	ImGuiParticleSize();
	ImGuiModeText();
	ImGuiPhysicsTimeStep();
	ImGUiNumberOfParticles();
	ImGuiGravityButtons();

	g_theRenderer->EndCamera(m_worldCamera);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::Update(float deltaSeconds)
{
	GetMouseCursorPos();
	
	if (g_theInput->IsKeyDown('W') && m_mapBounds.IsPointInside(m_cursorPos))
	{
		SpawnParticlesAtMousePos();
	}
	
	UpdatePhysics(deltaSeconds);

	UpdateFluidMesh();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::InitializeWorldCamera()
{
	m_worldCamera.m_mode = Camera::eMode_Orthographic;
	m_worldCamera.SetOrthographicView(Vec2(0.f, 0.f), Vec2(MAP_WIDTH / 2, MAP_HEIGHT / 2));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiModeText()
{
	if (m_mode == MULTI_THREADED)
	{
		ImGui::Text("Multi threaded");
	}

	if (m_mode == SINGLE_THREADED)
	{
		ImGui::Text("Single threaded");
	}

	if (ImGui::Button("Choose CPU Mode"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Choose CPU Mode");
	}

	if (ImGui::BeginPopup("Choose CPU Mode"))
	{
		static int selectedSimMode = m_mode;
		const char* items[] = { "Single Threaded", "Multi-threaded" };

		if (ImGui::Combo("Choose CPU Mode", &selectedSimMode, items, IM_ARRAYSIZE(items)))
		{
			m_mode = static_cast<CPUMode>(selectedSimMode);
		}

		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::RenderParticles()
{
	g_theRenderer->BindTexture(m_texture);
	g_theRenderer->BindShader(nullptr);

	g_theRenderer->DrawVertexBuffer(m_solverVertexBuffer, (int)m_solverVerts.size());
	g_theRenderer->BindTexture(nullptr);
	g_theRenderer->BindShader(nullptr);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::DensitySimpleMultiThreaded()
{
	const float h_sq = m_smoothingLength * m_smoothingLength;
	int numParticles = static_cast<int>(m_solverparticles.size());
	int numJobs = g_theJobSystem->GetNumThreads(); 
	int chunkSize = numParticles / numJobs;

	for (int i = 0; i < numJobs; i++)
	{
		int start = i * chunkSize;
		int end = (i == numJobs - 1) ? numParticles : start + chunkSize;
		std::vector<Solver*> particleSubset(m_solverparticles.begin() + start, m_solverparticles.begin() + end);

		Job* densityJob = new CalculateDensityJobSimple(this, particleSubset, h_sq);
		g_theJobSystem->QueueJob(densityJob);
	}

	while (g_theJobSystem->GetNumCompletedJobs() != numJobs)
	{
		std::this_thread::yield();
	}

	Job* completedJob = g_theJobSystem->RetrieveCompletedJobs();
	while (completedJob != nullptr)
	{
		delete completedJob;
		completedJob = g_theJobSystem->RetrieveCompletedJobs();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::DensityGridMultiThreaded()
{
	const float h_sq = m_smoothingLength * m_smoothingLength;
	int numParticles = static_cast<int>(m_solverparticles.size());
	int numJobs = g_theJobSystem->GetNumThreads();
	int chunkSize = numParticles / numJobs;

	for (int i = 0; i < numJobs; i++)
	{
		int start = i * chunkSize;
		int end = (i == numJobs - 1) ? numParticles : start + chunkSize;
		std::vector<Solver*> particleSubset(m_solverparticles.begin() + start, m_solverparticles.begin() + end);

		Job* densityJobGrid = new CalculateDensityJobGrid(this, particleSubset, h_sq);
		g_theJobSystem->QueueJob(densityJobGrid);
	}

	while (g_theJobSystem->GetNumCompletedJobs() != numJobs)
	{
		std::this_thread::yield();
	}

	Job* completedJob = g_theJobSystem->RetrieveCompletedJobs();
	while (completedJob != nullptr)
	{
		delete completedJob;
		completedJob = g_theJobSystem->RetrieveCompletedJobs();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::DensitySimpleSingleThreaded()
{
	const float h_sq = m_smoothingLength * m_smoothingLength;


	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		float density = 0;
		Solver* particle = m_solverparticles[i];
		if (particle)
		{
			Vec2  particlePos = m_solverparticles[i]->m_particlePosition;

			for (int j = 0; j < m_solverparticles.size(); j++)
			{
				Solver* neighborParticle = m_solverparticles[j];
				if (neighborParticle)
				{
					Vec2 neighborPos = neighborParticle->m_particlePosition;
					Vec2 diff = neighborPos - particlePos;
					float r_sq = DotProduct2D(diff, diff);
					if (r_sq < h_sq)
					{
						density += CalculateDensity(r_sq);
					}
				}
			}
			m_solverparticles[i]->m_particleDensity = density;
		}

	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::DensityGridSingleThreaded()
{
	const float h_sq = m_smoothingLength * m_smoothingLength;

	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		float density = 0;
		Solver* particle = m_solverparticles[i];
		if (particle)
		{
			Vec2  particlePos = m_solverparticles[i]->m_particlePosition;

			std::vector<Solver*> neighbors = FindNeighbors(particle);

			for (int j = 0; j < neighbors.size(); j++)
			{
				Solver* neighborParticle = neighbors[j];
				if (neighborParticle)
				{
					Vec2 neighborPos = neighborParticle->m_particlePosition;
					Vec2 diff = neighborPos - particlePos;
					float r_sq = DotProduct2D(diff, diff);
					if (r_sq < h_sq)
					{
						density += CalculateDensity(r_sq);
					}
				}
			}
			m_solverparticles[i]->m_particleDensity = density;
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
float FluidSimMT::CalculateDensity(float r_sq)
{
	const float h_sq = m_smoothingLength * m_smoothingLength;
	float val = DENSITY_COEFFICIENT * (h_sq - r_sq) * (h_sq - r_sq) * (h_sq - r_sq);
	return val;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
float FluidSimMT::CalculatePressure(float density)
{
	float value = m_pressureStiffness * std::max((float)std::pow(density / m_restDensity, 3) - 1, 0.0f);
	return value;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSimMT::CalculateGradPressure(float r, float particlePressure, float N_pressure, float neighborDensity, Vec2 posDiff)
{
	const float h = m_smoothingLength;
	float avg_pressure = 0.5f * (N_pressure + particlePressure);
	Vec2 gradPressure = GRADUAL_PRESSURE_COEFFICIENT * avg_pressure / neighborDensity * (h - r) * (h - r) / r * (posDiff);
	return gradPressure;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSimMT::CalculateLapVelocity(float r, Vec2 particleVelocity, Vec2 neighborVelocity, float neighborDensity)
{
	const float h = m_smoothingLength;
	Vec2 vel_diff = (neighborVelocity - particleVelocity);
	Vec2 lapVelocity = LAPLACIAN_VISCOSITY_COEFFICIENT / neighborDensity * (h - r) * vel_diff;
	return lapVelocity;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ForceSimpleMultiThreaded()
{
	float h_sq = m_smoothingLength * m_smoothingLength;
	int numParticles = static_cast<int>(m_solverparticles.size());
	int numJobs = g_theJobSystem->GetNumThreads();
	int chunkSize = numParticles / numJobs;

	for (int i = 0; i < numJobs; i++)
	{
		int start = i * chunkSize;
		int end = (i == numJobs - 1) ? numParticles : start + chunkSize;
		std::vector<Solver*> particleSubset(m_solverparticles.begin() + start, m_solverparticles.begin() + end);

		Job* forceJob = new CalculateForcesJobSimple(this, particleSubset, h_sq);
		g_theJobSystem->QueueJob(forceJob);
	}

	while (g_theJobSystem->GetNumCompletedJobs() != numJobs)
	{
		std::this_thread::yield();
	}

	Job* completedJob = g_theJobSystem->RetrieveCompletedJobs();
	while (completedJob != nullptr)
	{
		delete completedJob;
		completedJob = g_theJobSystem->RetrieveCompletedJobs();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ForcesGridMultiThreaded()
{
	float h_sq = m_smoothingLength * m_smoothingLength;
	int numParticles = static_cast<int>(m_solverparticles.size());
	int numJobs = g_theJobSystem->GetNumThreads();
	int chunkSize = numParticles / numJobs;

	for (int i = 0; i < numJobs; i++)
	{
		int start = i * chunkSize;
		int end = (i == numJobs - 1) ? numParticles : start + chunkSize;
		std::vector<Solver*> particleSubset(m_solverparticles.begin() + start, m_solverparticles.begin() + end);

		Job* forceJobGrid = new CalculateForcesJobGrid(this, particleSubset, h_sq);
		g_theJobSystem->QueueJob(forceJobGrid);
	}

	while (g_theJobSystem->GetNumCompletedJobs() != numJobs)
	{
		std::this_thread::yield();
	}

	Job* completedJob = g_theJobSystem->RetrieveCompletedJobs();
	while (completedJob != nullptr)
	{
		delete completedJob;
		completedJob = g_theJobSystem->RetrieveCompletedJobs();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ForceSimpleSingleThreaded()
{
	float h_sq = m_smoothingLength * m_smoothingLength;
	Vec2 gravity = m_gravityDirection;
	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		Solver* particle = m_solverparticles[i];
		if (particle)
		{
			Vec2  particlePos = particle->m_particlePosition;
			Vec2  particleVelocity = particle->m_particleVelocity;
			float particleDensity = particle->m_particleDensity;
			float particlePressure = CalculatePressure(particleDensity);
			Vec2 acceleration = Vec2(0.f, 0.f);

			for (int j = 0; j < m_solverparticles.size(); j++)
			{
				Solver* neighbor = m_solverparticles[j];
				if (neighbor)
				{
					Vec2 neighborPos = neighbor->m_particlePosition;
					Vec2 diff = neighborPos - particlePos;
					float r_sq = diff.GetLengthSquared();

					if (r_sq < h_sq && (particle != neighbor))
					{
						Vec2 neighborVelocity = neighbor->m_particleVelocity;
						float neighborDensity = neighbor->m_particleDensity;
						float neighborPressure = CalculatePressure(neighborDensity);
						float r = sqrt(r_sq);

						acceleration += CalculateGradPressure(r, particlePressure, neighborPressure, neighborDensity, diff);
						acceleration += CalculateLapVelocity(r, particleVelocity, neighborVelocity, neighborDensity);
					}
				}
			}

			particle->m_particleAcceleration = acceleration / particleDensity + gravity;
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ForcesGridSingleThreaded()
{
	float h_sq = m_smoothingLength * m_smoothingLength;
	Vec2 gravity = m_gravityDirection;
	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		Solver* particle = m_solverparticles[i];
		if (particle)
		{
			Vec2  particlePos = particle->m_particlePosition;
			Vec2  particleVelocity = particle->m_particleVelocity;
			float particleDensity = particle->m_particleDensity;
			float particlePressure = CalculatePressure(particleDensity);
			Vec2 acceleration = Vec2(0.f, 0.f);

			std::vector<Solver*> neighbors = FindNeighbors(particle);

			for (int j = 0; j < neighbors.size(); j++)
			{
				Solver* neighbor = neighbors[j];
				if (neighbor)
				{
					Vec2 neighborPos = neighbor->m_particlePosition;
					Vec2 diff = neighborPos - particlePos;
					float r_sq = DotProduct2D(diff, diff);

					if (r_sq < h_sq && (particle != neighbor))
					{
						Vec2 neighborVelocity = neighbor->m_particleVelocity;
						float neighborDensity = neighbor->m_particleDensity;
						float neighborPressure = CalculatePressure(neighborDensity);
						float r = sqrt(r_sq);

						acceleration += CalculateGradPressure(r, particlePressure, neighborPressure, neighborDensity, diff);
						acceleration += CalculateLapVelocity(r, particleVelocity, neighborVelocity, neighborDensity);
					}
				}
			}

			particle->m_particleAcceleration = acceleration / particleDensity + gravity;
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::IntegrateMultiThreadedSimulation(float deltaSeconds)
{
	int numParticles = static_cast<int>(m_solverparticles.size());
	int numJobs = g_theJobSystem->GetNumThreads();
	int chunkSize = numParticles / numJobs;

	for (int i = 0; i < numJobs; i++)
	{
		int start = i * chunkSize;
		int end = (i == numJobs - 1) ? numParticles : start + chunkSize;
		std::vector<Solver*> particleSubset(m_solverparticles.begin() + start, m_solverparticles.begin() + end);

		Job* integrateJob = new IntegrateParticlesJob(this, particleSubset, deltaSeconds, LEFT_MAP_BOUNDS, TOP_MAP_BOUNDS, RIGHT_MAP_BOUNDS, BOTTOM_MAP_BOUNDS, m_wallStiffness);
		g_theJobSystem->QueueJob(integrateJob);
	}

	// Wait for all jobs to complete
	while (g_theJobSystem->GetNumCompletedJobs() != numJobs)
	{
		std::this_thread::yield();
	}

	// Cleanup completed jobs
	Job* completedJob = g_theJobSystem->RetrieveCompletedJobs();
	while (completedJob != nullptr)
	{
		delete completedJob;
		completedJob = g_theJobSystem->RetrieveCompletedJobs();
	}

	m_isBufferDirty = true;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::IntegrateSingleThreadedSimulation(float deltaSeconds)
{
	for (Solver* particle : m_solverparticles)
	{
		if (particle)
		{
			//Wall collisions
			float distanceFromLeftWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), LEFT_MAP_BOUNDS);
			particle->m_particleAcceleration += std::min(distanceFromLeftWall, 0.f) * -m_wallStiffness * Vec2(LEFT_MAP_BOUNDS.x, LEFT_MAP_BOUNDS.y);

			float distanceFromTopWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), TOP_MAP_BOUNDS);
			particle->m_particleAcceleration += std::min(distanceFromTopWall, 0.f) * -m_wallStiffness * Vec2(TOP_MAP_BOUNDS.x, TOP_MAP_BOUNDS.y);

			float distanceFromRightWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), RIGHT_MAP_BOUNDS);
			particle->m_particleAcceleration += std::min(distanceFromRightWall, 0.f) * -m_wallStiffness * Vec2(RIGHT_MAP_BOUNDS.x, RIGHT_MAP_BOUNDS.y);

			float distanceFromBottomWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), BOTTOM_MAP_BOUNDS);
			particle->m_particleAcceleration += std::min(distanceFromBottomWall, 0.f) * -m_wallStiffness * Vec2(BOTTOM_MAP_BOUNDS.x, BOTTOM_MAP_BOUNDS.y);

			//Attract particles to mouse pos
			if (g_theInput->IsKeyDown(KEYCODE_LEFT_MOUSE))
			{
				Vec2 mouseForce = ApplyMouseForces(m_cursorPos, 0.3f, 10.f, Vec2(particle->m_particlePosition), Vec2(particle->m_particleVelocity));
				particle->m_particleAcceleration += mouseForce;
			}

			//Repel particles from mouse pos
			if (g_theInput->IsKeyDown(KEYCODE_RIGHT_MOUSE))
			{
				Vec2 mouseForce = ApplyMouseForces(m_cursorPos, 0.3f, -5.f, Vec2(particle->m_particlePosition), Vec2(particle->m_particleVelocity));
				particle->m_particleAcceleration += mouseForce;
			}

			//Collisions against Capsules and OBB's
 			if (g_theGame->m_initialParticleConfigCPU == RANDOM_SHAPES_CPU)
 			{
 				if (particle)
 				{
 					for (int j = 0; j < m_randomShapesConfigOBBVector.size(); j++)
 					{
 						OBB2D* obb = m_randomShapesConfigOBBVector[j];
 						if (obb)
 						{
 							BounceOffFixedOBB2D(particle->m_particlePosition, m_particleSize + m_particleSize / 10, particle->m_particleVelocity, *obb, 1.f);
 						}
 					}
 
 					for (int k = 0; k < m_randomShapesConfigCapsuleVector.size(); k++)
 					{
 						Capsule2* capsule = m_randomShapesConfigCapsuleVector[k];
 						if (capsule)
 						{
 							BounceOffFixedCapsule2D(particle->m_particlePosition, m_particleSize + m_particleSize / 10, particle->m_particleVelocity, *capsule, 0.5f);
 						}
 					}

					BounceDiscOffFixedDisc2D(particle->m_particlePosition, m_particleSize + m_particleSize / 10, particle->m_particleVelocity, m_windmillConfigDisc->m_center, m_windmillConfigDisc->m_radius, 0.7f);
 				}
 			}

			//Calculate particle pos and velocity
			particle->m_particleVelocity += particle->m_particleAcceleration * deltaSeconds;
			particle->m_particlePosition += particle->m_particleVelocity * deltaSeconds;

			//Reset particles acceleration for next frame
			particle->m_particleAcceleration = Vec2(0.f, 0.f);
		}
	}

	m_isBufferDirty = true;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::UpdatePhysics(float deltaSeconds)
{
	if (m_mode == SINGLE_THREADED)
	{
		if (m_simMode == SIMPLE)
		{
			DensitySimpleSingleThreaded();
			ForceSimpleSingleThreaded();
		}
		if (m_simMode == GRID)
		{
			UpdateGrid();
			DensityGridSingleThreaded();
			ForcesGridSingleThreaded();
		}

		IntegrateSingleThreadedSimulation(std::min(m_maxAllowableTimeStep, deltaSeconds));
	}

	if (m_mode == MULTI_THREADED)
	{
		if (m_simMode == SIMPLE)
		{
			DensitySimpleMultiThreaded();
			ForceSimpleMultiThreaded();
		}
		if (m_simMode == GRID)
		{
			UpdateGrid();
			DensityGridMultiThreaded();
			ForcesGridMultiThreaded();
		}
		IntegrateMultiThreadedSimulation(std::min(m_maxAllowableTimeStep, deltaSeconds));
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::UpdateFluidMesh()
{
	if (m_isBufferDirty)
	{
		m_solverVerts.clear();
		if (m_solverVertexBuffer)
		{
			delete m_solverVertexBuffer;
			m_solverVertexBuffer = nullptr;
		}

		for (int i = 0; i < m_solverparticles.size(); i++)
		{
			Vec2 center = m_solverparticles[i]->m_particlePosition;

			Vec2 bl = center + Vec2(-m_particleSize, -m_particleSize);
			Vec2 br = center + Vec2(m_particleSize, -m_particleSize);
			Vec2 tr = center + Vec2(m_particleSize, m_particleSize);
			Vec2 tl = center + Vec2(-m_particleSize, m_particleSize);

			AABB2 bounds = AABB2(bl, tr);
			AddVertsForAABB2D(m_solverVerts, bounds, Rgba8::WHITE);
		}

		m_solverVertexBuffer = g_theRenderer->CreateVertexBuffer(m_solverVerts.size() * sizeof(m_solverVerts[0]), sizeof(Vertex_PCU));
		g_theRenderer->CopyCPUToGPU(m_solverVerts.data(), m_solverVerts.size() * sizeof(m_solverVerts[0]), m_solverVertexBuffer);
		m_isBufferDirty = false;
	}

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiSimulationMode()
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
int FluidSimMT::CalculateHashKey(const Vec2& position) const
{
	int gridX = static_cast<int>(floor(position.x / m_cellSize));
	int gridY = static_cast<int>(floor(position.y / m_cellSize));
	return gridX * 73856093 ^ gridY * 19349663;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::AddParticleToGrid(Solver* particle)
{
	int hashKey = CalculateHashKey(particle->m_particlePosition);
	m_grid[hashKey].push_back(particle);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<Solver*> FluidSimMT::GetParticlesFromGrid(const Vec2& position)
{
	int hashKey = CalculateHashKey(position);
	return m_grid.count(hashKey) ? m_grid[hashKey] : std::vector<Solver*>();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::UpdateGrid()
{
	m_grid.clear();
	for (Solver* particle : m_solverparticles)
	{
		AddParticleToGrid(particle);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<Solver*> FluidSimMT::FindNeighbors(Solver* particle) const
{
	std::vector<Solver*> neighbors;

	for (int dx = -1; dx <= 1; dx++)
	{
		for (int dy = -1; dy <= 1; dy++)
		{
			int neighborHashKey = CalculateHashKey(Vec2(
				particle->m_particlePosition.x + dx * m_cellSize,
				particle->m_particlePosition.y + dy * m_cellSize));

			if (m_grid.count(neighborHashKey))
			{
				auto pair = m_grid.find(neighborHashKey);
				std::vector<Solver*> cellParticles;
				if (pair != m_grid.end())
				{
					cellParticles = pair->second;
				}
				neighbors.insert(neighbors.end(), cellParticles.begin(), cellParticles.end());
			}
		}
	}
	return neighbors;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::SpawnParticlesAtMousePos()
{
	Solver* solverParticle = new Solver();
	solverParticle->m_particlePosition = m_cursorPos;
	m_solverparticles.push_back(solverParticle);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::GetMouseCursorPos()
{
	Vec2 normalizedMousePos = g_theInput->GetCursorNormalizedPosition();
	m_cursorPos = Vec2(RangeMap(normalizedMousePos.x, 0.f, 1.f, 0.f, MAP_WIDTH / 2), RangeMap(normalizedMousePos.y, 0.f, 1.f, 0.f, MAP_HEIGHT / 2));
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiParticleSize()
{
	ImGui::SliderFloat("Particle Size", &m_particleSize, 0.001f, 0.009f, "%.5f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGUiNumberOfParticles()
{
	ImGui::Text("Num Particles - %i", (int)m_solverparticles.size());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiGravityButtons()
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
		const char* items[] = { "Up", "Down", "Left", "Right" };

		ImGui::Combo("Direction", &gravityDirection, items, IM_ARRAYSIZE(items));
		switch (gravityDirection)
		{
		case 0: m_gravityDirection = GRAVITY_UP; break;
		case 1: m_gravityDirection = GRAVITY_DOWN; break;
		case 2: m_gravityDirection = GRAVITY_LEFT; break;
		case 3: m_gravityDirection = GRAVITY_RIGHT; break;
		}
		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiMillisecondText()
{
	ImGui::Text("m/s - %.2f", m_lastSimulationTimeMs);
	ImGui::Text("m/s - %.2f", m_lastRenderTimeMs);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ImGuiPhysicsTimeStep()
{
	ImGui::SliderFloat("Physics Time Step", &m_maxAllowableTimeStep, 0.001f, 0.009f, "%.5f");
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::InitializeParticlesInADamConfig(unsigned int numParticles)
{
	(void)numParticles;
	int startingWidth = (int)sqrt(m_numParticles);

	for (int i = 0; i < m_numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		Solver* particle = new Solver();
		particle->m_particlePosition = Vec2(0.6f + m_initialParticleSpacing * (float)x, m_initialParticleSpacing * (float)y + 0.1f);
		m_solverparticles.push_back(particle);
		AddParticleToGrid(particle);
	}

	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		Vec2 center = m_solverparticles[i]->m_particlePosition;

		Vec2 bl = center + Vec2(-m_particleSize, -m_particleSize);
		Vec2 br = center + Vec2(m_particleSize, -m_particleSize);
		Vec2 tr = center + Vec2(m_particleSize, m_particleSize);
		Vec2 tl = center + Vec2(-m_particleSize, m_particleSize);

		AABB2 bounds = AABB2(bl, tr);
		AddVertsForAABB2D(m_solverVerts, bounds, Rgba8::WHITE);
	}

	m_solverVertexBuffer = g_theRenderer->CreateVertexBuffer(m_solverVerts.size() * sizeof(m_solverVerts[0]), sizeof(Vertex_PCU));
	g_theRenderer->CopyCPUToGPU(m_solverVerts.data(), m_solverVerts.size() * sizeof(m_solverVerts[0]), m_solverVertexBuffer);
	m_isBufferDirty = false;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 FluidSimMT::ApplyMouseForces(Vec2 mousePos, float forceRadius, float forceStrength, Vec2 const& position, Vec2 const& velocity)
{
	Vec2 totalForce = Vec2(0.f, 0.f);
	Vec2 offset = mousePos - position;
	float sqrDist = DotProduct2D(offset, offset); 

	if (sqrDist < forceRadius * forceRadius)
	{
		float dist = sqrt(sqrDist);
		Vec2 dirToInputPoint = dist <= 0.0001 ? Vec2(0, 0) : offset / dist;
		float centreT = 1 - dist / forceRadius;
		totalForce = m_gravityDirection + dirToInputPoint * centreT * forceStrength;
		totalForce -= velocity * centreT;
	}

	return totalForce;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::InitializeParticlesInARandomShapeConfig(unsigned int numParticles)
{
	(void)numParticles;
	int startingWidth = (int)sqrt(m_numParticles);

	for (int i = 0; i < m_numParticles; i++)
	{
		int x = i % startingWidth;
		int y = i / startingWidth;
		Solver* particle = new Solver();
		particle->m_particlePosition = Vec2(0.0f + m_initialParticleSpacing * (float)x, m_initialParticleSpacing * (float)y + 0.1f);
		m_solverparticles.push_back(particle);
		AddParticleToGrid(particle);
	}

	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		Vec2 center = m_solverparticles[i]->m_particlePosition;

		Vec2 bl = center + Vec2(-m_particleSize, -m_particleSize);
		Vec2 br = center + Vec2(m_particleSize, -m_particleSize);
		Vec2 tr = center + Vec2(m_particleSize, m_particleSize);
		Vec2 tl = center + Vec2(-m_particleSize, m_particleSize);

		AABB2 bounds = AABB2(bl, tr);
		AddVertsForAABB2D(m_solverVerts, bounds, Rgba8::WHITE);
	}

	m_solverVertexBuffer = g_theRenderer->CreateVertexBuffer(m_solverVerts.size() * sizeof(m_solverVerts[0]), sizeof(Vertex_PCU));
	g_theRenderer->CopyCPUToGPU(m_solverVerts.data(), m_solverVerts.size() * sizeof(m_solverVerts[0]), m_solverVertexBuffer);
	m_isBufferDirty = false;

	InitializeShapesForRandomShapeConfig();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::InitializeShapesForRandomShapeConfig()
{
	m_randomShapeConfigVerts.clear();

	m_windmillConfigDisc = new Disc2D();
	m_windmillConfigDisc->m_center = Vec2(MAP_WIDTH / 5, MAP_HEIGHT / 9);
	m_windmillConfigDisc->m_radius = 0.05f;
	AddVertsForDisc2D(m_randomShapeConfigVerts, m_windmillConfigDisc->m_center, m_windmillConfigDisc->m_radius, Rgba8::YELLOW);

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
			obb->m_halfDimensions = Vec2(0.01f, 0.07f);
			AddVertsForOBB2D(m_randomShapeConfigVerts, *obb, Rgba8::WHITE);
			m_randomShapesConfigOBBVector.push_back(obb);
		}
	}

	Capsule2* capsule = new Capsule2();
	capsule->radius = m_rng.RollRandomFloatInRange(0.01f, 0.02f);
	capsule->m_bone.m_start = GetClampedRandomPositionInWorld();
	capsule->m_bone.m_end = GetClampedRandomPositionInWorld();
	if (capsule)
	{
		AddVertsForCapsule2D(m_randomShapeConfigVerts, capsule->m_bone.m_start, capsule->m_bone.m_end, capsule->radius, Rgba8::WHITE);
		m_randomShapesConfigCapsuleVector.push_back(capsule);
	}

	Capsule2* capsule1 = new Capsule2();
	capsule1->radius = m_rng.RollRandomFloatInRange(0.01f, 0.03f);
	capsule1->m_bone.m_start = GetClampedRandomPositionInWorld();
	capsule1->m_bone.m_end = GetClampedRandomPositionInWorld();
	if (capsule1)
	{
		AddVertsForCapsule2D(m_randomShapeConfigVerts, capsule1->m_bone.m_start, capsule1->m_bone.m_end, capsule1->radius, Rgba8::WHITE);
		m_randomShapesConfigCapsuleVector.push_back(capsule1);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::RenderShapesForRandomConfig()
{
	g_theRenderer->BindShader(nullptr);
	g_theRenderer->DrawVertexArray((int)m_randomShapeConfigVerts.size(),m_randomShapeConfigVerts.data());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void FluidSimMT::ProcessCollisionsForRandomShapeConfig()
{
	for (int i = 0; i < m_solverparticles.size(); i++)
	{
		Solver* solver = m_solverparticles[i];

		if (solver)
		{
			for (int j = 0; j < m_randomShapesConfigOBBVector.size(); j++)
			{
				OBB2D* obb = m_randomShapesConfigOBBVector[j];
				if (obb)
				{
					BounceOffFixedOBB2D(solver->m_particlePosition, m_particleSize + m_particleSize / 10, solver->m_particleVelocity, *obb, 1.f);
				}
			}

			for (int k = 0; k < m_randomShapesConfigCapsuleVector.size(); k++)
			{
				Capsule2* capsule = m_randomShapesConfigCapsuleVector[k];
				if (capsule)
				{
					BounceOffFixedCapsule2D(solver->m_particlePosition, m_particleSize + m_particleSize / 10, solver->m_particleVelocity, *capsule, 0.5f);
				}
			}
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSimMT::GetClampedRandomPositionInWorld()
{
	float randX = m_rng.RollRandomFloatInRange(MAP_WIDTH / 2 * 0.4f, MAP_WIDTH / 2 * 0.9f);
	float randY = m_rng.RollRandomFloatInRange(MAP_HEIGHT / 2 * 0.1f, MAP_HEIGHT / 2 * 0.9f);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSimMT::GetClampedRandomPositionInWorld(float shapeRadius)
{
	float randX = m_rng.RollRandomFloatInRange(shapeRadius, MAP_WIDTH / 2 - shapeRadius);
	float randY = m_rng.RollRandomFloatInRange(shapeRadius, MAP_HEIGHT / 2 - shapeRadius);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Vec2 const FluidSimMT::GetRandomPositionInWorld()
{
	float randX = m_rng.RollRandomFloatInRange(0.0f, MAP_WIDTH / 2);
	float randY = m_rng.RollRandomFloatInRange(0.0f, MAP_HEIGHT / 2);

	return Vec2(randX, randY);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
CalculateDensityJobSimple::CalculateDensityJobSimple(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq)
	: Job(1), m_particleSubset(particleSubset), m_hSq(hSq), m_sim(fluidSimMT)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateDensityJobSimple::Execute()
{
	for (Solver* particle : m_particleSubset)
	{
		float density = 0.0f;
		Vec2 particlePos = particle->m_particlePosition;

		for (Solver* neighborParticle : m_sim->m_solverparticles)
		{
			Vec2 neighborPos = neighborParticle->m_particlePosition;
			Vec2 diff = neighborPos - particlePos;
			float r_sq = DotProduct2D(diff, diff);
			if (r_sq < m_hSq)
			{
				density += m_sim->CalculateDensity(r_sq);
			}
		}
		particle->m_particleDensity = density;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateDensityJobSimple::OnFinished()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
CalculateDensityJobGrid::CalculateDensityJobGrid(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq)
	: Job(2), m_particleSubset(particleSubset), m_hSq(hSq), m_sim(fluidSimMT)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateDensityJobGrid::Execute()
{
	for (Solver* particle : m_particleSubset)
	{
		float density = 0.0f;
		Vec2 particlePos = particle->m_particlePosition;

		std::vector<Solver*> neighbors = m_sim->FindNeighbors(particle);

		for (Solver* neighborParticle : neighbors)
		{
			Vec2 neighborPos = neighborParticle->m_particlePosition;
			Vec2 diff = neighborPos - particlePos;
			float r_sq = DotProduct2D(diff, diff);
			if (r_sq < m_hSq)
			{
				density += m_sim->CalculateDensity(r_sq);
			}
		}
		particle->m_particleDensity = density;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateDensityJobGrid::OnFinished()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
CalculateForcesJobSimple::CalculateForcesJobSimple(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq)
	: Job(3), m_particleSubset(particleSubset), m_hSq(hSq), m_sim(fluidSimMT)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateForcesJobSimple::Execute()
{
	Vec2 gravity = m_sim->m_gravityDirection; 

	for (Solver* particle : m_particleSubset)
	{
		Vec2 particlePos = particle->m_particlePosition;
		Vec2 particleVelocity = particle->m_particleVelocity;
		float particleDensity = particle->m_particleDensity;
		float particlePressure = m_sim->CalculatePressure(particleDensity);
		Vec2 acceleration = Vec2(0.f, 0.f);

		for (Solver* neighbor : m_sim->m_solverparticles)
		{
			if (neighbor && neighbor != particle)
			{
				Vec2 neighborPos = neighbor->m_particlePosition;
				Vec2 diff = neighborPos - particlePos;
				float r_sq = DotProduct2D(diff, diff);

				if (r_sq < m_hSq)
				{
					Vec2 neighborVelocity = neighbor->m_particleVelocity;
					float neighborDensity = neighbor->m_particleDensity;
					float neighborPressure = m_sim->CalculatePressure(neighborDensity);
					float r = sqrt(r_sq);

					acceleration += m_sim->CalculateGradPressure(r, particlePressure, neighborPressure, neighborDensity, diff);
					acceleration += m_sim->CalculateLapVelocity(r, particleVelocity, neighborVelocity, neighborDensity);
				}
			}
		}

		particle->m_particleAcceleration = acceleration / particleDensity + gravity;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateForcesJobSimple::OnFinished()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
CalculateForcesJobGrid::CalculateForcesJobGrid(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float hSq)
	: Job(4), m_particleSubset(particleSubset), m_hSq(hSq), m_sim(fluidSimMT)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateForcesJobGrid::Execute()
{
	Vec2 gravity = m_sim->m_gravityDirection;

	for (Solver* particle : m_particleSubset)
	{
		Vec2 particlePos = particle->m_particlePosition;
		Vec2 particleVelocity = particle->m_particleVelocity;
		float particleDensity = particle->m_particleDensity;
		float particlePressure = m_sim->CalculatePressure(particleDensity);
		Vec2 acceleration = Vec2(0.f, 0.f);

		std::vector<Solver*> neighbors = m_sim->FindNeighbors(particle);

		for (Solver* neighbor : neighbors)
		{
			if (neighbor && neighbor != particle)
			{
				Vec2 neighborPos = neighbor->m_particlePosition;
				Vec2 diff = neighborPos - particlePos;
				float r_sq = DotProduct2D(diff, diff);

				if (r_sq < m_hSq)
				{
					Vec2 neighborVelocity = neighbor->m_particleVelocity;
					float neighborDensity = neighbor->m_particleDensity;
					float neighborPressure = m_sim->CalculatePressure(neighborDensity);
					float r = sqrt(r_sq);

					acceleration += m_sim->CalculateGradPressure(r, particlePressure, neighborPressure, neighborDensity, diff);
					acceleration += m_sim->CalculateLapVelocity(r, particleVelocity, neighborVelocity, neighborDensity);
				}
			}
		}

		particle->m_particleAcceleration = acceleration / particleDensity + gravity;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculateForcesJobGrid::OnFinished()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
IntegrateParticlesJob::IntegrateParticlesJob(FluidSimMT* fluidSimMT, std::vector<Solver*> particleSubset, float deltaSeconds, const Vec3& leftmapBounds, const Vec3& topmapBounds, const Vec3& rightmapBounds, const Vec3& bottommapBounds, float wallstiffness)
	: Job(4), m_fluidSim(fluidSimMT), m_particleSubset(particleSubset), m_deltaSeconds(deltaSeconds), m_leftMapBounds(leftmapBounds), m_topMapBounds(topmapBounds), m_rightMapBounds(rightmapBounds), m_bottomMapBounds(bottommapBounds), m_wallStiffness(wallstiffness)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void IntegrateParticlesJob::Execute()
{
	for (Solver* particle : m_particleSubset)
	{
		if (particle)
		{
			//Wall collisions
			float distanceFromLeftWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), m_leftMapBounds);
			particle->m_particleAcceleration += std::min(distanceFromLeftWall, 0.f) * -m_wallStiffness * Vec2(m_leftMapBounds.x, m_leftMapBounds.y);

			float distanceFromTopWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), m_topMapBounds);
			particle->m_particleAcceleration += std::min(distanceFromTopWall, 0.f) * -m_wallStiffness * Vec2(m_topMapBounds.x, m_topMapBounds.y);

			float distanceFromRightWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), m_rightMapBounds);
			particle->m_particleAcceleration += std::min(distanceFromRightWall, 0.f) * -m_wallStiffness * Vec2(m_rightMapBounds.x, m_rightMapBounds.y);

			float distanceFromBottomWall = DotProduct3D(Vec3(particle->m_particlePosition, 1), m_bottomMapBounds);
			particle->m_particleAcceleration += std::min(distanceFromBottomWall, 0.f) * -m_wallStiffness * Vec2(m_bottomMapBounds.x, m_bottomMapBounds.y);

			//Boundary collisions and mouse forces
			if (m_fluidSim)
			{
				//Attract particles to mouse pos
				if (g_theInput->IsKeyDown(KEYCODE_LEFT_MOUSE))
				{
					Vec2 mouseForce = m_fluidSim->ApplyMouseForces(m_fluidSim->m_cursorPos, 0.3f, 10.f, Vec2(particle->m_particlePosition), Vec2(particle->m_particleVelocity));
					particle->m_particleAcceleration += mouseForce;
				}

				//Repel particles from mouse pos
				if (g_theInput->IsKeyDown(KEYCODE_RIGHT_MOUSE))
				{
					Vec2 mouseForce = m_fluidSim->ApplyMouseForces(m_fluidSim->m_cursorPos, 0.3f, -5.f, Vec2(particle->m_particlePosition), Vec2(particle->m_particleVelocity));
					particle->m_particleAcceleration += mouseForce;
				}

				//Collisions against OBB's and Capsules which are randomized
				if (g_theGame->m_initialParticleConfigCPU == RANDOM_SHAPES_CPU)
				{
					for (int j = 0; j < m_fluidSim->m_randomShapesConfigOBBVector.size(); j++)
					{
						OBB2D* obb = m_fluidSim->m_randomShapesConfigOBBVector[j];
						if (obb)
						{
							BounceOffFixedOBB2D(particle->m_particlePosition, m_fluidSim->m_particleSize + m_fluidSim->m_particleSize / 10, particle->m_particleVelocity, *obb, 1.f);
						}
					}

					for (int k = 0; k < m_fluidSim->m_randomShapesConfigCapsuleVector.size(); k++)
					{
						Capsule2* capsule = m_fluidSim->m_randomShapesConfigCapsuleVector[k];
						if (capsule)
						{
							BounceOffFixedCapsule2D(particle->m_particlePosition, m_fluidSim->m_particleSize + m_fluidSim->m_particleSize / 10, particle->m_particleVelocity, *capsule, 0.5f);
						}
					}

					BounceDiscOffFixedDisc2D(particle->m_particlePosition, m_fluidSim->m_particleSize + m_fluidSim->m_particleSize / 10, particle->m_particleVelocity,  m_fluidSim->m_windmillConfigDisc->m_center, m_fluidSim->m_windmillConfigDisc->m_radius, 0.7f);
				}
			}

			//Calculate particle velocity and position
			particle->m_particleVelocity += particle->m_particleAcceleration * m_deltaSeconds;
			particle->m_particlePosition += particle->m_particleVelocity * m_deltaSeconds;

			//Reset acceleration
			particle->m_particleAcceleration = Vec2(0.f, 0.f);
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void IntegrateParticlesJob::OnFinished()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
