/*
 * interrupt_handling.h
 *
 *  Created on: Jul 21, 2011
 *      Author: martin
 */

#ifndef INTERRUPT_HANDLING_H_
#define INTERRUPT_HANDLING_H_

#include <csignal>
#include <stdexcept>

//Used for signaling
static bool interrupt_process;

void interrupt(int s) {
	interrupt_process = true;
}

void init_interrupt_handler() {

	//Reset signalling
	interrupt_process = false;

	signal(SIGINT, interrupt);
}

#define SGL_INTERRUPT_INIT init_interrupt_handler();
#define SGL_INTERRUPT_CHECK if (interrupt_process) throw std::runtime_error("Interrupted");
#define SGL_INTERRUPT_RESET interrupt_process = false;
#define SGL_INTERRUPT interrupt(0);

#endif /* INTERRUPT_HANDLING_H_ */
