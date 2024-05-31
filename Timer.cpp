#include "Timer.h"


Timer::Timer(const char* name) : m_name{ name }
{
}
Timer::Timer() : m_start{ std::chrono::high_resolution_clock::now()}, m_name{ "N/A" }
{}
Timer::~Timer()
{
	m_end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(m_end- m_start);
	float duration_ms = duration.count() * 1000.0f;
	std::cout << m_name << ": " << duration.count() << "s" << " (" << duration_ms << "ms)\n";
}
