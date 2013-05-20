/*
 * simple_timer.h
 *
 *  Created on: Jul 19, 2011
 *      Author: martin
 */

#ifndef SIMPLE_TIMER_H_
#define SIMPLE_TIMER_H_

#ifdef SGL_TIMING

class SimpleTimer {

private:
	int total;
	int s;

	int times;

	const std::string func;
	const std::string file;
	const int line;

public:

	SimpleTimer(std::string func, std::string file, int line) : total(0), s(0), times(0), func(func), file(file), line(line) {}

	~SimpleTimer() {
		rout << func << " " << static_cast<double>(total) / CLOCKS_PER_SEC <<  " seconds - x"<< times << "." << " (in " << file << " at line "<< line << ") " << std::endl;
	}

	void start() {
		++times;
		s = clock();
	}

	void end() {
		total += clock() - s;
	}
};

class TimerScope {

private:
	SimpleTimer & timer;
public:

	TimerScope(SimpleTimer & timer) : timer(timer) {
		timer.start();
	}

	~TimerScope() {
		timer.end();
	}
};

#define TIMER_START static SimpleTimer timer(__func__, __FILE__, __LINE__); TimerScope timer_scope(timer);
#else
#define TIMER_START
#endif

#endif /* SIMPLE_TIMER_H_ */
