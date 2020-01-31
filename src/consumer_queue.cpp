#include "consumer_queue.h"
#include "sample.h"
#include "send_buffer.h"
#include "common.h"
#include <boost/chrono/duration.hpp>

// === implementation of the consumer_queue class ===

using namespace lsl;

/**
* Create a new queue with a given capacity.
* @param max_capacity The maximum number of samples that can be held by the queue. Beyond that, the oldest samples are dropped.
* @param registry Optionally a pointer to a registration facility, to dispatch samples to all consumers.
*/
consumer_queue::consumer_queue(std::size_t max_capacity, send_buffer_p registry): registry_(registry), buffer_(max_capacity)  {
	if (registry_)
		registry_->register_consumer(this);
}

/**
* Destructor.
* Unregisters from the send buffer, if any.
*/
consumer_queue::~consumer_queue() {
	try {
		if (registry_)
			registry_->unregister_consumer(this);
	} catch(std::exception &e) {
		LOG_F(ERROR,
			"Unexpected error while trying to unregister a consumer queue from its registry: %s",
			e.what());
	}
}

/**
* Push a new sample onto the queue.
*/
void consumer_queue::push_sample(const sample_p &sample) {
	while (!buffer_.push(sample)) {
		sample_p dummy;
		buffer_.pop(dummy);
	}
}

/**
* Pop a sample from the queue.
* Blocks if empty.
* @param timeout Timeout for the blocking, in seconds. If expired, an empty sample is returned.
*/
sample_p consumer_queue::pop_sample(double timeout) {
	sample_p result;
	if (timeout <= 0.0) {
		buffer_.pop(result);
	} else {
		if (!buffer_.pop(result)) {
			// turn timeout into the point in time at which we give up
			timeout += lsl::lsl_clock();
			do {
				if (lsl::lsl_clock() >= timeout)
					break;
				lslboost::this_thread::sleep_for(lslboost::chrono::milliseconds(1));
			} while (!buffer_.pop(result));
		}
	}
	return result;
}

bool consumer_queue::empty() {
	return buffer_.empty();
}
