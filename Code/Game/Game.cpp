#include "Game/Game.hpp"
#include "Game/App.hpp"
#include "Game/FluidSim.hpp"
#include "Engine/Renderer/Camera.hpp"
#include "Engine/Core/Time.hpp"
#include "Engine/Renderer/BitmapFont.hpp"
#include "ThirdParty/ImGUI/imgui.h"
#include "Game/FluidSimMT.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------
Game::Game()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
Game::~Game()
{
	delete m_fluidSim;
	m_fluidSim = nullptr;

	delete m_fluidSimCPU;
	m_fluidSimCPU = nullptr;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::Startup()
{
	InitializeGPUFluidSimMode();
	InitializeCPUFluidSimMode();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ShutDown()
{
	if (m_fluidSim)
	{
		m_fluidSim->Shutdown();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::Update(float deltaSeconds)
{
	ToggleModes();

	if (g_theInput->WasKeyJustPressed('P'))
	{
		g_theApp->m_clock.TogglePause();
	}
	if (g_theInput->IsKeyDown('T'))
	{
		g_theApp->m_clock.SetTimeScale(0.1f);
	}
	else if(g_theInput->WasKeyJustReleased('T'))
	{
		g_theApp->m_clock.SetTimeScale(1);
	}

	if (m_fluidSim)
	{
		ImGuiChangeNumberOfParticlesGPU();
		ImGuiResetSim();
		ImGuiInitialParticleConfigGPU();
		m_fluidSim->Update(deltaSeconds);
	}
	
	if (m_fluidSimCPU)
	{
		ImGuiResetSim();
		ImGuiInitialParticleConfigCPU();
		m_fluidSimCPU->Update(deltaSeconds);
	}

	if (g_theInput->WasKeyJustPressed(KEYCODE_F8))
	{
		RandomizeGameMode();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::Render()
{	
	if (m_fluidSim)
	{
		m_fluidSim->Render();
	}

	if (m_fluidSimCPU)
	{
		m_fluidSimCPU->Render();
	}

	ImGuiFpsText();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::EndFrame()
{
	if (m_fluidSim)
	{
		m_fluidSim->EndFrame();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::RenderImGuiUI()
{

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ImGuiResetSim()
{
	if (ImGui::Button("Reset Simulation"))
	{
		if (m_gameMode == GPU)
		{
			if (m_fluidSim)
			{
				delete m_fluidSim;
				m_fluidSim = nullptr;
			}

			m_fluidSim = new FluidSim();
			m_fluidSim->Startup(m_numParticlesGPU, m_initialParticleConfigGPU);
		}

		if (m_gameMode == SINGLE_AND_MULTI_THREAD)
		{
			if (m_fluidSimCPU)
			{
				delete m_fluidSimCPU;
				m_fluidSimCPU = nullptr;
			}

			m_fluidSimCPU = new FluidSimMT();
			m_fluidSimCPU->Startup(m_numParticlesCPU, m_initialParticleConfigCPU);
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ImGuiChangeNumberOfParticlesGPU()
{
	if (ImGui::Button("Number of particles"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Number of particles");
	}

	if (ImGui::BeginPopup("Number of particles"))
	{
		int particleOption = 0;
		if (m_numParticlesGPU == NUM_PARTICLES_8K) particleOption = 0;
		else if (m_numParticlesGPU == NUM_PARTICLES_16K) particleOption = 1;
		else if (m_numParticlesGPU == NUM_PARTICLES_32K) particleOption = 2;
		else if (m_numParticlesGPU == NUM_PARTICLES_64K) particleOption = 3;

		const char* items[] = { "8k", "16k", "32k", "64k" };

		if (ImGui::Combo("Num Particles", &particleOption, items, IM_ARRAYSIZE(items)))
		{
			switch (particleOption)
			{
			case 0: m_numParticlesGPU = NUM_PARTICLES_8K; break;
			case 1: m_numParticlesGPU = NUM_PARTICLES_16K; break;
			case 2: m_numParticlesGPU = NUM_PARTICLES_32K; break;
			case 3: m_numParticlesGPU = NUM_PARTICLES_64K; break;
			}
			
			if (m_fluidSim)
			{
				delete m_fluidSim;
				m_fluidSim = nullptr;
			}

			m_fluidSim = new FluidSim();
			m_fluidSim->Startup(m_numParticlesGPU, m_initialParticleConfigGPU);
		}
		ImGui::EndPopup();
	}
	
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ImGuiInitialParticleConfigGPU()
{
	if (ImGui::Button("Initial Particle Config"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Initial Particle Config");
	}

	if (ImGui::BeginPopup("Initial Particle Config"))
	{
		int particleOption = 0;
		if (m_initialParticleConfigGPU == DAM) particleOption = 0;
		else if (m_initialParticleConfigGPU == RANDOM_SHAPES) particleOption = 1;
		else if (m_initialParticleConfigGPU == WINDMILL) particleOption = 2;
		else if (m_initialParticleConfigGPU == HOURGLASS) particleOption = 3;
		else if (m_initialParticleConfigGPU == TRANSPARENT_FLUID) particleOption = 4;
		else if (m_initialParticleConfigGPU == WATER_ELEVATOR) particleOption = 5;
		const char* items[] = { "DAM", "RANDOM SHAPES", "WINDMILL", "HOURGLASS", "TRANSPARENT_FLUID", "WATER_ELEVATOR"};

		if (ImGui::Combo("Num Particles", &particleOption, items, IM_ARRAYSIZE(items)))
		{
			switch (particleOption)
			{
			case 0: m_initialParticleConfigGPU = DAM;
				break;
			case 1: m_initialParticleConfigGPU = RANDOM_SHAPES; 
				break;
			case 2: m_initialParticleConfigGPU = WINDMILL; 
				break;
			case 3: m_initialParticleConfigGPU = HOURGLASS;
				break;
			case 4: m_initialParticleConfigGPU = TRANSPARENT_FLUID;
				break;
			case 5: m_initialParticleConfigGPU = WATER_ELEVATOR;
				break;
			}

			if (m_fluidSim)
			{
				delete m_fluidSim;
				m_fluidSim = nullptr;
			}

			m_fluidSim = new FluidSim();
			m_fluidSim->Startup(m_numParticlesGPU, m_initialParticleConfigGPU);
		}
		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::SwitchToPreviousMode()
{
	m_gameMode = static_cast<SimulationThread>((static_cast<int>(m_gameMode) + NUM_SIM_MODES - 1) % NUM_SIM_MODES);
	CreateCurrentGameMode();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::DeleteCurrentGameMode()
{
	if (m_fluidSim) 
	{
		delete m_fluidSim;
		m_fluidSim = nullptr;
	}

	if (m_fluidSimCPU)
	{
		delete m_fluidSimCPU;
		m_fluidSimCPU = nullptr;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::CreateCurrentGameMode()
{
	DeleteCurrentGameMode();
	
	switch (m_gameMode) 
	{
	case SINGLE_AND_MULTI_THREAD:
		m_fluidSimCPU = new FluidSimMT();
		if (m_fluidSimCPU)
		{
			m_fluidSimCPU->Startup(m_numParticlesCPU, m_initialParticleConfigCPU);
		}
		break;
	case GPU:
		m_fluidSim = new FluidSim();
		if (m_fluidSim)
		{
			m_fluidSim->Startup(m_numParticlesGPU, m_initialParticleConfigGPU);
		}
		break;
	default:
		break;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::SwitchToNextMode()
{
	m_gameMode = static_cast<SimulationThread>((static_cast<int>(m_gameMode) + 1) % NUM_SIM_MODES);
	CreateCurrentGameMode();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ToggleModes()
{
	if (g_theInput->WasKeyJustPressed(KEYCODE_F6))
	{
		SwitchToPreviousMode();
	}

	if (g_theInput->WasKeyJustPressed(KEYCODE_F7))
	{
		SwitchToNextMode();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::ImGuiFpsText()
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
void Game::ImGuiInitialParticleConfigCPU()
{
	if (ImGui::Button("Initial Particle Config"))
	{
		ImVec2 buttonSize = ImGui::GetItemRectSize();
		ImVec2 buttonPos = ImGui::GetItemRectMin();
		ImVec2 popupPos = ImVec2(buttonPos.x, buttonPos.y + buttonSize.y);
		ImGui::SetNextWindowPos(popupPos);
		ImGui::OpenPopup("Initial Particle Config");
	}

	if (ImGui::BeginPopup("Initial Particle Config"))
	{
		int particleOption = 0;
		if (m_initialParticleConfigCPU == CPU_DAM) particleOption = 0;
		else if (m_initialParticleConfigCPU == RANDOM_SHAPES_CPU) particleOption = 1;
		
		const char* items[] = { "DAM",  "RANDOM SHAPES" };

		if (ImGui::Combo("Num Particles", &particleOption, items, IM_ARRAYSIZE(items)))
		{
			switch (particleOption)
			{
			case 0: m_initialParticleConfigCPU = CPU_DAM;
				break;
			case 1: m_initialParticleConfigCPU = RANDOM_SHAPES_CPU;
				break;
			}

			if (m_fluidSimCPU)
			{
				delete m_fluidSimCPU;
				m_fluidSimCPU = nullptr;
			}

			m_fluidSimCPU = new FluidSimMT();
			m_fluidSimCPU->Startup(m_numParticlesCPU, m_initialParticleConfigCPU);
		}
		ImGui::EndPopup();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::InitializeGPUFluidSimMode()
{
	m_initialParticleConfigGPU = DAM;

	if (m_gameMode == GPU)
	{
		if (m_initialParticleConfigGPU == DAM)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, DAM);
			}
		}

		if (m_initialParticleConfigGPU == RANDOM_SHAPES)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, RANDOM_SHAPES);
			}
		}

		if (m_initialParticleConfigGPU == WINDMILL)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, WINDMILL);
			}
		}

		if (m_initialParticleConfigGPU == HOURGLASS)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, HOURGLASS);
			}
		}

		if (m_initialParticleConfigGPU == TRANSPARENT_FLUID)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, TRANSPARENT_FLUID);
			}
		}

		if (m_initialParticleConfigGPU == WATER_ELEVATOR)
		{
			m_fluidSim = new FluidSim();
			if (m_fluidSim)
			{
				m_fluidSim->Startup(m_numParticlesGPU, WATER_ELEVATOR);
			}
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::InitializeCPUFluidSimMode()
{
	m_initialParticleConfigCPU = CPU_DAM;

	if (m_gameMode == SINGLE_AND_MULTI_THREAD)
	{
		m_fluidSimCPU = new FluidSimMT();
		if (m_fluidSimCPU)
		{
			m_fluidSimCPU->Startup(m_numParticlesCPU, m_initialParticleConfigCPU);
		}
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void Game::RandomizeGameMode()
{
	if (m_gameMode == GPU)
	{
		if (m_fluidSim)
		{
			delete m_fluidSim;
			m_fluidSim = nullptr;
		}

		m_fluidSim = new FluidSim();
		m_fluidSim->Startup(m_numParticlesGPU, m_initialParticleConfigGPU);
	}

	if (m_gameMode == SINGLE_AND_MULTI_THREAD)
	{
		if (m_fluidSimCPU)
		{
			delete m_fluidSimCPU;
			m_fluidSimCPU = nullptr;
		}

		m_fluidSimCPU = new FluidSimMT();
		m_fluidSimCPU->Startup(m_numParticlesCPU, m_initialParticleConfigCPU);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
