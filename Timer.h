
#include <chrono>
#include<iostream>

class Timer
{
	std::chrono::steady_clock::time_point m_start ,m_end;
	const char* m_name;
public:
	Timer(const char* name);
	Timer();
	~Timer();
};
