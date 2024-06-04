#include <windows.h> 
#include "Game/App.hpp"
#include "Engine/Renderer/Renderer.hpp"
#include "Engine/Core/Time.hpp"
#include "Engine/Core/EngineCommon.hpp"
#include "Engine/Audio/AudioSystem.hpp"
#include "Engine/Window/Window.hpp"
#include "Game/GameCommon.hpp"
#include "Game/Game.hpp"
#include "Engine/Core/DevConsole.hpp"
#include "Engine/Input/InputSystem.hpp"
#include "ThirdParty/ImGUI/imgui.h"
#include "ThirdParty/ImGUI/imgui_impl_dx11.h"
#include "ThirdParty/ImGUI/imgui_impl_win32.h"
#include "Engine/Core/JobSystem.hpp"
#include "Engine/Core/NamedProperties.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------
App*          g_theApp = nullptr;				
Renderer*     g_theRenderer = nullptr;
InputSystem*  g_theInput = nullptr;
AudioSystem*  g_theAudio = nullptr;
Window*		  g_theWindow = nullptr;
Game*		  g_theGame = nullptr;
JobSystem*	  g_theJobSystem = nullptr;
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::Startup() 
{
	m_devConsoleCamera.m_mode = Camera::eMode_Orthographic;
	m_devConsoleCamera.SetOrthographicView(Vec2(0.f, 0.f), Vec2(1600.f, 800.f));

	//constructs the engine components
	EventSystemConfig eventSystemConfig;
	g_theEventSystem = new EventSystem(eventSystemConfig);
	g_theEventSystem->SubscribeEventCallbackFunction("Quit", App::Event_Quit);

	InputSystemConfig inputSystemConfig;
	g_theInput = new InputSystem(inputSystemConfig);

	WindowConfig  windowConfig;
	windowConfig.m_windowTitle = "FluidSimGPU2D";
	
	windowConfig.m_windowMode = FULLSCREEN;
	windowConfig.m_inputSystem = g_theInput;
	g_theWindow = new Window(windowConfig);

	RendererConfig rendererConfig;
	//rendererConfig.m_window = g_theWindow;
	g_theRenderer = new Renderer(rendererConfig);

	AudioSystemConfig audioSystemConfig;
	g_theAudio = new AudioSystem(audioSystemConfig);

	DevConsoleConfig devConsoleConfig;
	devConsoleConfig.m_renderer = g_theRenderer;
	devConsoleConfig.m_camera = &m_devConsoleCamera;
	g_theDevConsole = new DevConsole(devConsoleConfig);
	 	
	JobSystemConfig jobSystemConfig;
 	jobSystemConfig.m_numWorkerThreads = std::thread::hardware_concurrency();
 	g_theJobSystem = new JobSystem(jobSystemConfig);

	//Calling startup of all engine components
	g_theEventSystem->Startup();
	g_theWindow->Startup();
	g_theRenderer->Startup();
	g_theDevConsole->Startup();
	g_theInput->Startup();
	g_theAudio->Startup();
	g_theJobSystem->Startup();
 	ImGui::CreateContext();
 	//Gui::StyleColorsDark();
 	ImGui_ImplWin32_Init(g_theWindow->GetHwnd());
 	ImGui_ImplDX11_Init(g_theRenderer->m_device, g_theRenderer->m_deviceContext);

	g_theGame = new Game();
	g_theGame->Startup();

	PrintControlsToDevConsole();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::Shutdown()
{
	//Delete ImGui components
	ImGui_ImplDX11_Shutdown();
	ImGui_ImplWin32_Shutdown();
	ImGui::DestroyContext();

	//calling Shut-down of all engine components
	g_theGame->ShutDown();
	g_theAudio->Shutdown();
	g_theInput->Shutdown();
	g_theDevConsole->Shutdown();
	g_theRenderer->Shutdown();
	g_theWindow->Shutdown();
	g_theEventSystem->Shutdown();
	g_theJobSystem->ShutDown();

	//deleting engine components
	delete g_theGame;
	g_theGame = nullptr;

	delete g_theAudio;
	g_theAudio = nullptr;
	
	delete g_theRenderer;
	g_theRenderer = nullptr;

	delete g_theWindow;
	g_theWindow = nullptr;

	delete g_theDevConsole;
	g_theDevConsole = nullptr;

	delete g_theInput;
	g_theInput = nullptr;
	
	delete g_theEventSystem;
	g_theEventSystem = nullptr;

	 delete g_theJobSystem;
 	g_theJobSystem = nullptr;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
App::App()
{
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
App::~App()
{
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::RunFrame()
{
	float deltaSeconds = m_clock.GetDeltaSeconds();
	BeginFrame();
	Update(deltaSeconds);
	Render();
	EndFrame();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::Run()
{
	while (!IsQuitting() )
	{
		g_theApp->RunFrame();
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
bool App::HandleQuitRequested()
{
	m_isQuitting = true;
	return true;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
bool App::Event_Quit(EventArgs& args)
{
	(void)args;
	
	if (!g_theApp)
	{
		return false;
	}

	g_theApp->HandleQuitRequested();
	return true;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::PrintControlsToDevConsole()
{
	std::string switchSimulationModeText = "F6/F7 - Switch Simulation Modes";
	std::string lmbMouseForceText = "LMB - Attract particles to mouse position";
	std::string rmbMouseForceText = "LMB - Repel particles from mouse center";
	std::string escapeText = "Esc - Quit";
	std::string slowTimeText = "T - Slow time";
	std::string stepTimeText = "O - Step single frame";
	std::string pauseText = "P - Pause";

	g_theDevConsole->AddLine(Rgba8::WHITE, switchSimulationModeText);
	g_theDevConsole->AddLine(Rgba8::WHITE, lmbMouseForceText);
	g_theDevConsole->AddLine(Rgba8::WHITE, rmbMouseForceText);
	g_theDevConsole->AddLine(Rgba8::WHITE, escapeText);
	g_theDevConsole->AddLine(Rgba8::WHITE, slowTimeText);
	g_theDevConsole->AddLine(Rgba8::WHITE, stepTimeText);
	g_theDevConsole->AddLine(Rgba8::WHITE, pauseText);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::BeginFrame()
{
	Clock::TickSytemClock();

	g_theDevConsole->BeginFrame();
	g_theInput->BeginFrame();
	g_theAudio->BeginFrame();
	g_theWindow->BeginFrame();
	g_theRenderer->BeginFrame();
	g_theJobSystem->BeginFrame();
	
	ImGui_ImplDX11_NewFrame();
	ImGui_ImplWin32_NewFrame();
	ImGui::NewFrame();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::Update(float deltaSeconds)
{
	if (g_theInput->WasKeyJustPressed('O'))
	{
		m_clock.StepSingleFrame();
	}

	if (g_theGame != nullptr)
	{
		g_theGame->Update(deltaSeconds);
	}

	if (g_theInput->WasKeyJustPressed(KEYCODE_ESC))
	{
		EventArgs args;
		FireEvent("Quit", args);
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::Render() const
{
	Rgba8 clearColor{255,0,255,255};

	if(g_theGame)
	{
		g_theGame->Render();
	}

	AABB2 bounds(Vec2(0.f,0.f), Vec2(1000.f, 500.f));

	if (g_theDevConsole)
	{
		g_theDevConsole->Render(bounds);
	}

	ImGui::Render();
	ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
void App::EndFrame()
{
	g_theRenderer->Endframe();
	g_theWindow->EndFrame();
	g_theAudio->EndFrame();
	g_theDevConsole->EndFrame();
	g_theJobSystem->EndFrame();
	g_theInput->EndFrame();
	g_theGame->EndFrame();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
