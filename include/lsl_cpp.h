#ifndef LSL_CPP_H
#define LSL_CPP_H

/**
 * The lab streaming layer provides a set of functions to make instrument data accessible
 * in real time within a lab network. From there, streams can be picked up by recording programs,
 * viewing programs or custom experiment applications that access data streams in real time.
 *
 * The API covers two areas:
 * - The "push API" allows to create stream outlets and to push data (regular or irregular
 * measurement time series, event data, coded audio/video frames, etc.) into them.
 * - The "pull API" allows to create stream inlets and read time-synched experiment data from them
 *   (for recording, viewing or experiment control).
 */

namespace lsl {
/**
 * The stream_info object stores the declaration of a data stream.
 *
 * Represents the following information:
 *  a) stream data format (number of channels, channel format)
 *  b) core information (stream name, content type, sampling rate)
 *  c) optional meta-data about the stream content (channel labels, measurement units, etc.)
 *
 * Whenever a program wants to provide a new stream on the lab network it will typically first
 * create a stream_info to describe its properties and then construct a stream_outlet with it to
 * create the stream on the network. Recipients who discover the outlet can query the stream_info;
 * it is also written to disk when recording the stream (playing a similar role as a file header).
 */
class stream_info {
public:
	/// Default contructor.
	stream_info(): stream_info("untitled", "", 0, 0, cf_undefined, ""){}

private:
	std::shared_ptr<lsl_streaminfo_struct_> obj;
};


// =======================
// ==== Stream Outlet ====
// =======================

class stream_outlet {
public:
	stream_outlet(const stream_info &info, int32_t chunk_size = 0, int32_t max_buffered = 360,
		lsl_transport_options_t flags = transp_default)
		: channel_count(info.channel_count()), sample_rate(info.nominal_srate()),
		  obj(lsl_create_outlet_ex(info.handle().get(), chunk_size, max_buffered, flags),
			  &lsl_destroy_outlet) {}

	// ========================================
	// === Pushing a sample into the outlet ===
	// ========================================

	/** 
	 * @param timestamp Optionally the capture time of the sample, in agreement with
	 * lsl::local_clock(); if omitted, the current time is used.
	 * @param pushthrough Whether to push the sample through to the receivers instead of
	 * buffering it with subsequent samples.
	 * Note that the chunk_size, if specified at outlet construction, takes precedence over the
	 * pushthrough flag.
	 */
	template <class T, int32_t N>
	void push_sample(const T data[N], double timestamp = 0.0, bool pushthrough = true) {
		check_numchan(N);
		push_sample(&data[0], timestamp, pushthrough);
	}

	/** Push a std vector of values as a sample into the outlet.
	 * Each entry in the vector corresponds to one channel. The function handles type checking &
	 * conversion.
	 * @param data A vector of values to push (one for each channel).
	 * @param timestamp Optionally the capture time of the sample, in agreement with local_clock();
	 * if omitted, the current time is used.
	 * @param pushthrough Whether to push the sample through to the receivers instead of buffering
	 * it with subsequent samples. Note that the chunk_size, if specified at outlet construction,
	 * takes precedence over the pushthrough flag.
	 */
	template<typename T>
	void push_sample(
		const std::vector<T> &data, double timestamp = 0.0, bool pushthrough = true) {
		check_numchan(data.size());
		push_sample(data.data(), timestamp, pushthrough);
	}

	/** Push a pointer to some values as a sample into the outlet.
	 * This is a lower-level function for cases where data is available in some buffer.
	 * Handles type checking & conversion.
	 * @param data A pointer to values to push. The number of values pointed to must not be less
	 * than the number of channels in the sample.
	 * @param timestamp Optionally the capture time of the sample, in agreement with local_clock();
	 * if omitted, the current time is used.
	 * @param pushthrough Whether to push the sample through to the receivers instead of buffering
	 * it with subsequent samples. Note that the chunk_size, if specified at outlet construction,
	 * takes precedence over the pushthrough flag.
	 */
	void push_sample(const float *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_ftp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const double *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_dtp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const int64_t *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_ltp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const int32_t *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_itp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const int16_t *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_stp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const char *data, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_ctp(obj.get(), (data), timestamp, pushthrough);
	}
	void push_sample(const std::string *data, double timestamp = 0.0, bool pushthrough = true) {
		std::vector<uint32_t> lengths(channel_count);
		std::vector<const char *> pointers(channel_count);
		for (int32_t k = 0; k < channel_count; k++) {
			pointers[k] = data[k].c_str();
			lengths[k] = (uint32_t)data[k].size();
		}
		lsl_push_sample_buftp(obj.get(), pointers.data(), lengths.data(), timestamp, pushthrough);
	}

	/** Push a packed C struct (of numeric data) as one sample into the outlet (search for
	 * [`#``pragma pack`](https://stackoverflow.com/a/3318475/73299) for information on packing
	 * structs appropriately).<br>
	 * Overall size checking but no type checking or conversion are done.<br>
	 * Can not be used forvariable-size / string-formatted data.
	 * @param sample The sample struct to push.
	 * @param timestamp Optionally the capture time of the sample, in agreement with
	 * local_clock(); if omitted, the current time is used.
	 * @param pushthrough Whether to push the sample through to the receivers instead of
	 * buffering it with subsequent samples. Note that the chunk_size, if specified at outlet
	 * construction, takes precedence over the pushthrough flag.
	 */
	template <class T>
	void push_numeric_struct(const T &sample, double timestamp = 0.0, bool pushthrough = true) {
		if (info().sample_bytes() != sizeof(T))
			throw std::runtime_error(
				"Provided object size does not match the stream's sample size.");
		push_numeric_raw((void *)&sample, timestamp, pushthrough);
	}

	/** Push a pointer to raw numeric data as one sample into the outlet.
	 * This is the lowest-level function; performs no checking whatsoever. Cannot be used for
	 * variable-size / string-formatted channels.
	 * @param sample A pointer to the raw sample data to push.
	 * @param timestamp Optionally the capture time of the sample, in agreement with local_clock();
	 * if omitted, the current time is used.
	 * @param pushthrough Whether to push the sample through to the receivers instead of buffering
	 * it with subsequent samples. Note that the chunk_size, if specified at outlet construction,
	 * takes precedence over the pushthrough flag.
	 */
	void push_numeric_raw(const void *sample, double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_sample_vtp(obj.get(), (sample), timestamp, pushthrough);
	}


	// ===================================================
	// === Pushing an chunk of samples into the outlet ===
	// ===================================================

	/** Push a chunk of samples (batched into an STL vector) into the outlet.
	 * @param samples A vector of samples in some supported format (each sample can be a data
	 * pointer, data array, or std vector of data).
	 * @param timestamp Optionally the capture time of the most recent sample, in agreement with
	 * local_clock(); if omitted, the current time is used. The time stamps of other samples are
	 * automatically derived according to the sampling rate of the stream.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template <class T>
	void push_chunk(
		const std::vector<T> &samples, double timestamp = 0.0, bool pushthrough = true) {
		if (!samples.empty()) {
			if (timestamp == 0.0) timestamp = local_clock();
			if (sample_rate != IRREGULAR_RATE)
				timestamp = timestamp - (samples.size() - 1) / sample_rate;
			push_sample(samples[0], timestamp, pushthrough && samples.size() == 1);
			for (std::size_t k = 1; k < samples.size(); k++)
				push_sample(samples[k], DEDUCED_TIMESTAMP, pushthrough && k == samples.size() - 1);
		}
	}

	/** Push a chunk of samples (batched into an STL vector) into the outlet.
	 * Allows to specify a separate time stamp for each sample (for irregular-rate streams).
	 * @param samples A vector of samples in some supported format (each sample can be a data
	 * pointer, data array, or std vector of data).
	 * @param timestamps A vector of capture times for each sample, in agreement with local_clock().
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template <class T>
	void push_chunk(const std::vector<T> &samples, const std::vector<double> &timestamps,
		bool pushthrough = true) {
		for (unsigned k = 0; k < samples.size() - 1; k++)
			push_sample(samples[k], timestamps[k], false);
		if (!samples.empty()) push_sample(samples.back(), timestamps.back(), pushthrough);
	}

	/** Push a chunk of numeric data as C-style structs (batched into an STL vector) into the
	 * outlet. This performs some size checking but no type checking. Can not be used for
	 * variable-size / string-formatted data.
	 * @param samples A vector of samples, as C structs.
	 * @param timestamp Optionally the capture time of the sample, in agreement with local_clock();
	 * if omitted, the current time is used.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template <class T>
	void push_chunk_numeric_structs(
		const std::vector<T> &samples, double timestamp = 0.0, bool pushthrough = true) {
		if (!samples.empty()) {
			if (timestamp == 0.0) timestamp = local_clock();
			if (sample_rate != IRREGULAR_RATE)
				timestamp = timestamp - (samples.size() - 1) / sample_rate;
			push_numeric_struct(samples[0], timestamp, pushthrough && samples.size() == 1);
			for (std::size_t k = 1; k < samples.size(); k++)
				push_numeric_struct(
					samples[k], DEDUCED_TIMESTAMP, pushthrough && k == samples.size() - 1);
		}
	}

	/** Push a chunk of numeric data from C-style structs (batched into an STL vector), into the
	 * outlet. This performs some size checking but no type checking. Can not be used for
	 * variable-size / string-formatted data.
	 * @param samples A vector of samples, as C structs.
	 * @param timestamps A vector of capture times for each sample, in agreement with local_clock().
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template <class T>
	void push_chunk_numeric_structs(const std::vector<T> &samples,
		const std::vector<double> &timestamps, bool pushthrough = true) {
		for (unsigned k = 0; k < samples.size() - 1; k++)
			push_numeric_struct(samples[k], timestamps[k], false);
		if (!samples.empty()) push_numeric_struct(samples.back(), timestamps.back(), pushthrough);
	}

	/** Push a chunk of multiplexed data into the outlet.
	 * @name Push functions
	 * @param buffer A buffer of channel values holding the data for zero or more successive samples
	 * to send.
	 * @param timestamp Optionally the capture time of the most recent sample, in agreement with
	 * local_clock(); if omitted, the current time is used. The time stamps of other samples are
	 * automatically derived according to the sampling rate of the stream.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template<typename T>
	void push_chunk_multiplexed(
		const std::vector<T> &buffer, double timestamp = 0.0, bool pushthrough = true) {
		if (!buffer.empty())
			push_chunk_multiplexed(
				buffer.data(), static_cast<unsigned long>(buffer.size()), timestamp, pushthrough);
	}

	/** Push a chunk of multiplexed data into the outlet. One timestamp per sample is provided.
	 * Allows to specify a separate time stamp for each sample (for irregular-rate streams).
	 * @param buffer A buffer of channel values holding the data for zero or more successive samples
	 * to send.
	 * @param timestamps A buffer of timestamp values holding time stamps for each sample in the
	 * data buffer.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	template<typename T>
	void push_chunk_multiplexed(const std::vector<T> &buffer,
		const std::vector<double> &timestamps, bool pushthrough = true) {
		if (!buffer.empty() && !timestamps.empty())
			push_chunk_multiplexed(
				buffer.data(), static_cast<unsigned long>(buffer.size()), timestamps.data(), pushthrough);
	}

	/** Push a chunk of multiplexed samples into the outlet. Single timestamp provided.
	 * @warning The provided buffer size is measured in channel values (e.g., floats), not samples.
	 * @param buffer A buffer of channel values holding the data for zero or more successive samples
	 * to send.
	 * @param buffer_elements The number of channel values (of type T) in the buffer. Must be a
	 * multiple of the channel count.
	 * @param timestamp Optionally the capture time of the most recent sample, in agreement with
	 * local_clock(); if omitted, the current time is used. The time stamps of other samples are
	 * automatically derived based on the sampling rate of the stream.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the stream_outlet() constructur parameter @p chunk_size,
	 * if specified at outlet construction, takes precedence over the pushthrough flag.
	 */
	void push_chunk_multiplexed(const float *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_ftp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const double *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_dtp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const int64_t *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_ltp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const int32_t *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_itp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const int16_t *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_stp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const char *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		lsl_push_chunk_ctp(obj.get(), buffer, static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
	}
	void push_chunk_multiplexed(const std::string *buffer, std::size_t buffer_elements,
		double timestamp = 0.0, bool pushthrough = true) {
		if (buffer_elements) {
			std::vector<uint32_t> lengths(buffer_elements);
			std::vector<const char *> pointers(buffer_elements);
			for (std::size_t k = 0; k < buffer_elements; k++) {
				pointers[k] = buffer[k].c_str();
				lengths[k] = (uint32_t)buffer[k].size();
			}
			lsl_push_chunk_buftp(obj.get(), pointers.data(), lengths.data(),
				static_cast<unsigned long>(buffer_elements), timestamp, pushthrough);
		}
	}

	/** Push a chunk of multiplexed samples into the outlet. One timestamp per sample is provided.
	 * @warning Note that the provided buffer size is measured in channel values (e.g., floats)
	 * rather than in samples.
	 * @param data_buffer A buffer of channel values holding the data for zero or more successive
	 * samples to send.
	 * @param timestamp_buffer A buffer of timestamp values holding time stamps for each sample in
	 * the data buffer.
	 * @param data_buffer_elements The number of data values (of type T) in the data buffer. Must be
	 * a multiple of the channel count.
	 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
	 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
	 * precedence over the pushthrough flag.
	 */
	void push_chunk_multiplexed(const float *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_ftnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}
	void push_chunk_multiplexed(const double *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_dtnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}
	void push_chunk_multiplexed(const int64_t *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_ltnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}
	void push_chunk_multiplexed(const int32_t *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_itnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}
	void push_chunk_multiplexed(const int16_t *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_stnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}
	void push_chunk_multiplexed(const char *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		lsl_push_chunk_ctnp(obj.get(), data_buffer, static_cast<unsigned long>(data_buffer_elements),
			(timestamp_buffer), pushthrough);
	}

	void push_chunk_multiplexed(const std::string *data_buffer, const double *timestamp_buffer,
		std::size_t data_buffer_elements, bool pushthrough = true) {
		if (data_buffer_elements) {
			std::vector<uint32_t> lengths(data_buffer_elements);
			std::vector<const char *> pointers(data_buffer_elements);
			for (std::size_t k = 0; k < data_buffer_elements; k++) {
				pointers[k] = data_buffer[k].c_str();
				lengths[k] = (uint32_t)data_buffer[k].size();
			}
			lsl_push_chunk_buftnp(obj.get(), pointers.data(), lengths.data(),
				static_cast<unsigned long>(data_buffer_elements), timestamp_buffer, pushthrough);
		}
	}


	// ===============================
	// === Miscellaneous Functions ===
	// ===============================

	/** Check whether consumers are currently registered.
	 * While it does not hurt, there is technically no reason to push samples if there is no
	 * consumer.
	 */
	bool have_consumers() { return lsl_have_consumers(obj.get()) != 0; }

	/** Wait until some consumer shows up (without wasting resources).
	 * @return True if the wait was successful, false if the timeout expired.
	 */
	bool wait_for_consumers(double timeout) { return lsl_wait_for_consumers(obj.get(), timeout) != 0; }

	/** Retrieve the stream info provided by this outlet.
	 * This is what was used to create the stream (and also has the Additional Network Information
	 * fields assigned).
	 */
	stream_info info() const { return stream_info(lsl_get_info(obj.get())); }

	/// Return a shared pointer to pass to C-API functions that aren't wrapped yet
	///
	/// Example: @code lsl_push_chunk_buft(outlet.handle().get(), data, …); @endcode
	std::shared_ptr<lsl_outlet_struct_> handle() { return obj; }

	/** Destructor.
	 * The stream will no longer be discoverable after destruction and all paired inlets will stop
	 * delivering data.
	 */
	~stream_outlet() = default;

	/// stream_outlet move constructor
	stream_outlet(stream_outlet &&res) noexcept  = default;

	stream_outlet &operator=(stream_outlet &&rhs) noexcept = default;


private:
	// The outlet is a non-copyable object.
	stream_outlet(const stream_outlet &rhs);
	stream_outlet &operator=(const stream_outlet &rhs);

	/// Check whether a given data length matches the number of channels; throw if not
	void check_numchan(std::size_t N) const {
		if (N != static_cast<std::size_t>(channel_count))
			throw std::runtime_error("Provided element count (" + std::to_string(N) +
									 ") does not match the stream's channel count (" +
									 std::to_string(channel_count) + '.');
	}

	int32_t channel_count;
	double sample_rate;
	std::shared_ptr<lsl_outlet_struct_> obj;
};


// ===========================
// ==== Resolve Functions ====
// ===========================

/** Resolve all streams on the network.
 * This function returns all currently available streams from any outlet on the network.
 * The network is usually the subnet specified at the local router, but may also include
 * a multicast group of machines (given that the network supports it), or list of hostnames.
 * These details may optionally be customized by the experimenter in a configuration file
 * (see Network Connectivity in the LSL wiki).
 * This is the default mechanism used by the browsing programs and the recording program.
 * @param wait_time The waiting time for the operation, in seconds, to search for streams.
 * If this is too short (<0.5s) only a subset (or none) of the outlets that are present on the
 * network may be returned.
 * @return A vector of stream info objects (excluding their desc field), any of which can
 *         subsequently be used to open an inlet. The full info can be retrieve from the inlet.
 */
inline std::vector<stream_info> resolve_streams(double wait_time = 1.0) {
	lsl_streaminfo buffer[1024];
	int nres = check_error(lsl_resolve_all(buffer, sizeof(buffer) / sizeof(lsl_streaminfo), wait_time));
	return std::vector<stream_info>(&buffer[0], &buffer[nres]);
}

/** Resolve all streams with a specific value for a given property.
 * If the goal is to resolve a specific stream, this method is preferred over resolving all streams
 * and then selecting the desired one.
 * @param prop The stream_info property that should have a specific value (e.g., "name", "type",
 * "source_id", or "desc/manufaturer").
 * @param value The string value that the property should have (e.g., "EEG" as the type property).
 * @param minimum Return at least this number of streams.
 * @param timeout Optionally a timeout of the operation, in seconds (default: no timeout).
 *                 If the timeout expires, less than the desired number of streams (possibly none)
 * will be returned.
 * @return A vector of matching stream info objects (excluding their meta-data), any of
 *         which can subsequently be used to open an inlet.
 */
inline std::vector<stream_info> resolve_stream(const std::string &prop, const std::string &value,
	int32_t minimum = 1, double timeout = FOREVER) {
	lsl_streaminfo buffer[1024];
	int nres = check_error(
		lsl_resolve_byprop(buffer, sizeof(buffer) / sizeof(lsl_streaminfo), prop.c_str(), value.c_str(), minimum, timeout));
	return std::vector<stream_info>(&buffer[0], &buffer[nres]);
}

/** Resolve all streams that match a given predicate.
 *
 * Advanced query that allows to impose more conditions on the retrieved streams; the given
 * string is an [XPath 1.0](http://en.wikipedia.org/w/index.php?title=XPath_1.0) predicate for
 * the `<info>` node (omitting the surrounding []'s)
 * @param pred The predicate string, e.g. `name='BioSemi'` or
 * `type='EEG' and starts-with(name,'BioSemi') and count(info/desc/channel)=32`
 * @param minimum Return at least this number of streams.
 * @param timeout Optionally a timeout of the operation, in seconds (default: no timeout).
 *                 If the timeout expires, less than the desired number of streams (possibly
 * none) will be returned.
 * @return A vector of matching stream info objects (excluding their meta-data), any of
 *         which can subsequently be used to open an inlet.
 */
inline std::vector<stream_info> resolve_stream(
	const std::string &pred, int32_t minimum = 1, double timeout = FOREVER) {
	lsl_streaminfo buffer[1024];
	int nres =
		check_error(lsl_resolve_bypred(buffer, sizeof(buffer) / sizeof(lsl_streaminfo), pred.c_str(), minimum, timeout));
	return std::vector<stream_info>(&buffer[0], &buffer[nres]);
}


// ======================
// ==== Stream Inlet ====
// ======================

/** A stream inlet.
 * Inlets are used to receive streaming data (and meta-data) from the lab network.
 */
class stream_inlet {
public:
	/**
	 * @param max_buflen Optionally the maximum amount of data to buffer (in seconds if there is a
	 * nominal sampling rate, otherwise x100 in samples). Recording applications want to use a
	 * fairly large buffer size here, while real-time applications would only buffer as much as
	 * they need to perform their next calculation.
	 * @param max_chunklen Optionally the maximum size, in samples, at which chunks are transmitted
	 * (the default corresponds to the chunk sizes used by the sender).
	 * Recording applications can use a generous size here (leaving it to the network how to pack
	 * things), while real-time applications may want a finer (perhaps 1-sample) granularity.
	 * If left unspecified (=0), the sender determines the chunk granularity.
	 */

	double pull_sample(std::string *buffer, int32_t buffer_elements, double timeout = FOREVER) {
		int32_t ec = 0;
		if (buffer_elements) {
			std::vector<char *> result_strings(buffer_elements);
			std::vector<uint32_t> result_lengths(buffer_elements);
			double res = lsl_pull_sample_buf(
				obj.get(), result_strings.data(), result_lengths.data(), buffer_elements, timeout, &ec);
			check_error(ec);
			for (int32_t k = 0; k < buffer_elements; k++) {
				buffer[k].assign(result_strings[k], result_lengths[k]);
				lsl_destroy_string(result_strings[k]);
			}
			return res;
		} else
			throw std::runtime_error(
				"Provided element count does not match the stream's channel count.");
	}

	/**
	 * Pull a sample from the inlet and read it into a custom C-style struct.
	 *
	 * Overall size checking but no type checking or conversion are done.
	 * Do not use for variable-size/string-formatted streams.
	 * @param sample The raw sample object to hold the data (packed C-style struct).
	 * Search for [`#``pragma pack`](https://stackoverflow.com/a/3318475/73299) for information
	 * on how to pack structs correctly.
	 * @param timeout The timeout for this operation, if any. Use 0.0 to make the function
	 * non-blocking.
	 * @return The capture time of the sample on the remote machine, or 0.0 if no new sample was
	 * available. To remap this time stamp to the local clock, add the value returned by
	 * .time_correction() to it.
	 * @throws lost_error (if the stream source has been lost).
	 */
	template <class T> double pull_numeric_struct(T &sample, double timeout = FOREVER) {
		return pull_numeric_raw((void *)&sample, sizeof(T), timeout);
	}

	/**
	 * Pull a sample from the inlet and read it into a pointer to raw data.
	 *
	 * No type checking or conversions are done (not recommended!).<br>
	 * Do not use for variable-size/string-formatted streams.
	 * @param sample A pointer to hold the resulting raw sample data.
	 * @param buffer_bytes The number of bytes allocated in the buffer.<br>
	 * Note: it is the responsibility of the user to allocate enough memory.
	 * @param timeout The timeout for this operation, if any. Use 0.0 to make the function
	 * non-blocking.
	 * @return The capture time of the sample on the remote machine, or 0.0 if no new sample was
	 * available. To remap this time stamp to the local clock, add the value returned by
	 * .time_correction() to it.
	 * @throws lost_error (if the stream source has been lost).
	 */
	double pull_numeric_raw(void *sample, int32_t buffer_bytes, double timeout = FOREVER) {
		int32_t ec = 0;
		double res = lsl_pull_sample_v(obj.get(), sample, buffer_bytes, timeout, &ec);
		check_error(ec);
		return res;
	}


	// =================================================
	// === Pulling a chunk of samples from the inlet ===
	// =================================================

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This is the most complete version, returning both the data and a timestamp for each sample.
	 * @param chunk A vector of vectors to hold the samples.
	 * @param timestamps A vector to hold the time stamps.
	 * @return True if some data was obtained.
	 * @throws lost_error (if the stream source has been lost).
	 */
	template <class T>
	bool pull_chunk(std::vector<std::vector<T>> &chunk, std::vector<double> &timestamps) {
		std::vector<T> sample;
		chunk.clear();
		timestamps.clear();
		while (double ts = pull_sample(sample, 0.0)) {
			chunk.push_back(sample);
			timestamps.push_back(ts);
		}
		return !chunk.empty();
	}

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This version returns only the most recent sample's time stamp.
	 * @param chunk A vector of vectors to hold the samples.
	 * @return The time when the most recent sample was captured
	 *         on the remote machine, or 0.0 if no new sample was available.
	 * @throws lost_error (if the stream source has been lost)
	 */
	template <class T> double pull_chunk(std::vector<std::vector<T>> &chunk) {
		double timestamp = 0.0;
		std::vector<T> sample;
		chunk.clear();
		while (double ts = pull_sample(sample, 0.0)) {
			chunk.push_back(sample);
			timestamp = ts;
		}
		return timestamp;
	}

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This function does not return time stamps for the samples. Invoked as: mychunk =
	 * pull_chunk<float>();
	 * @return A vector of vectors containing the obtained samples; may be empty.
	 * @throws lost_error (if the stream source has been lost)
	 */
	template <class T> std::vector<std::vector<T>> pull_chunk() {
		std::vector<std::vector<T>> result;
		std::vector<T> sample;
		while (pull_sample(sample, 0.0)) result.push_back(sample);
		return result;
	}

	/**
	 * Pull a chunk of data from the inlet into a pre-allocated buffer.
	 *
	 * This is a high-performance function that performs no memory allocations
	 * (useful for very high data rates or on low-powered devices).
	 * @warning The provided buffer size is measured in channel values (e.g., floats), not samples.
	 * @param data_buffer A pointer to a buffer of data values where the results shall be stored.
	 * @param timestamp_buffer A pointer to a buffer of timestamp values where time stamps shall be
	 * stored. If this is NULL, no time stamps will be returned.
	 * @param data_buffer_elements The size of the data buffer, in channel data elements (of type
	 * T). Must be a multiple of the stream's channel count.
	 * @param timestamp_buffer_elements The size of the timestamp buffer. If a timestamp buffer is
	 * provided then this must correspond to the same number of samples as data_buffer_elements.
	 * @param timeout The timeout for this operation, if any. When the timeout expires, the function
	 * may return before the entire buffer is filled. The default value of 0.0 will retrieve only
	 * data available for immediate pickup.
	 * @return data_elements_written Number of channel data elements written to the data buffer.
	 * @throws lost_error (if the stream source has been lost).
	 */
	std::size_t pull_chunk_multiplexed(float *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_f(obj.get(), data_buffer, timestamp_buffer,
			(unsigned long)data_buffer_elements, (unsigned long)timestamp_buffer_elements, timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(double *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_d(obj.get(), data_buffer, timestamp_buffer,
			(unsigned long)data_buffer_elements, (unsigned long)timestamp_buffer_elements, timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(int64_t *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_l(obj.get(), data_buffer, timestamp_buffer,
			(unsigned long)data_buffer_elements, (unsigned long)timestamp_buffer_elements, timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(int32_t *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_i(obj.get(), data_buffer, timestamp_buffer,
			(unsigned long)data_buffer_elements, (unsigned long)timestamp_buffer_elements, timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(int16_t *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_s(obj.get(), data_buffer, timestamp_buffer,
			(unsigned long)data_buffer_elements, (unsigned long)timestamp_buffer_elements, timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(char *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		std::size_t res = lsl_pull_chunk_c(obj.get(), data_buffer, timestamp_buffer,
			static_cast<unsigned long>(data_buffer_elements), static_cast<unsigned long>(timestamp_buffer_elements), timeout,
			&ec);
		check_error(ec);
		return res;
	}
	std::size_t pull_chunk_multiplexed(std::string *data_buffer, double *timestamp_buffer,
		std::size_t data_buffer_elements, std::size_t timestamp_buffer_elements,
		double timeout = 0.0) {
		int32_t ec = 0;
		if (data_buffer_elements) {
			std::vector<char *> result_strings(data_buffer_elements);
			std::vector<uint32_t> result_lengths(data_buffer_elements);
			std::size_t num = lsl_pull_chunk_buf(obj.get(), result_strings.data(), result_lengths.data(),
				timestamp_buffer, static_cast<unsigned long>(data_buffer_elements),
				static_cast<unsigned long>(timestamp_buffer_elements), timeout, &ec);
			check_error(ec);
			for (std::size_t k = 0; k < num; k++) {
				data_buffer[k].assign(result_strings[k], result_lengths[k]);
				lsl_destroy_string(result_strings[k]);
			}
			return num;
		};
		return 0;
	}

	/**
	 * Pull a multiplexed chunk of samples and optionally the sample timestamps from the inlet.
	 *
	 * @param chunk A vector to hold the multiplexed (Sample 1 Channel 1,
	 * S1C2, S2C1, S2C2, S3C1, S3C2, ...) samples
	 * @param timestamps A vector to hold the timestamps or nullptr
	 * @param timeout Time to wait for the first sample. The default value of 0.0 will not wait
	 * for data to arrive, pulling only samples already received.
	 * @param append (True:) Append data or (false:) clear them first
	 * @return True if some data was obtained.
	 * @throws lost_error (if the stream source has been lost).
	 */
	template <typename T>
	bool pull_chunk_multiplexed(std::vector<T> &chunk, std::vector<double> *timestamps = nullptr,
		double timeout = 0.0, bool append = false) {
		if (!append) {
			chunk.clear();
			if (timestamps) timestamps->clear();
		}
		std::vector<T> sample;
		double ts;
		if ((ts = pull_sample(sample, timeout)) == 0.0) return false;
		chunk.insert(chunk.end(), sample.begin(), sample.end());
		if (timestamps) timestamps->push_back(ts);
		const auto target = samples_available();
		chunk.reserve(chunk.size() + target * this->channel_count);
		if (timestamps) timestamps->reserve(timestamps->size() + target);
		while ((ts = pull_sample(sample, 0.0)) != 0.0) {
#if LSL_CPP11
			chunk.insert(chunk.end(), std::make_move_iterator(sample.begin()),
				std::make_move_iterator(sample.end()));
#else
			chunk.insert(chunk.end(), sample.begin(), sample.end());
#endif
			if (timestamps) timestamps->push_back(ts);
		}
		return true;
	}

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This is the most complete version, returning both the data and a timestamp for each sample.
	 * @param chunk A vector of C-style structs to hold the samples.
	 * @param timestamps A vector to hold the time stamps.
	 * @return True if some data was obtained.
	 * @throws lost_error (if the stream source has been lost)
	 */
	template <class T>
	bool pull_chunk_numeric_structs(std::vector<T> &chunk, std::vector<double> &timestamps) {
		T sample;
		chunk.clear();
		timestamps.clear();
		while (double ts = pull_numeric_struct(sample, 0.0)) {
			chunk.push_back(sample);
			timestamps.push_back(ts);
		}
		return !chunk.empty();
	}

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This version returns only the most recent sample's time stamp.
	 * @param chunk A vector of C-style structs to hold the samples.
	 * @return The time when the most recent sample was captured
	 *         on the remote machine, or 0.0 if no new sample was available.
	 * @throws lost_error (if the stream source has been lost)
	 */
	template <class T> double pull_chunk_numeric_structs(std::vector<T> &chunk) {
		double timestamp = 0.0;
		T sample;
		chunk.clear();
		while (double ts = pull_numeric_struct(sample, 0.0)) {
			chunk.push_back(sample);
			timestamp = ts;
		}
		return timestamp;
	}

	/**
	 * Pull a chunk of samples from the inlet.
	 *
	 * This function does not return time stamps. Invoked as: mychunk = pull_chunk<mystruct>();
	 * @return A vector of C-style structs containing the obtained samples; may be empty.
	 * @throws lost_error (if the stream source has been lost)
	 */
	template <class T> std::vector<T> pull_chunk_numeric_structs() {
		std::vector<T> result;
		T sample;
		while (pull_numeric_struct(sample, 0.0)) result.push_back(sample);
		return result;
	}

	std::shared_ptr<lsl_inlet_struct_> obj;
};

inline xml_element stream_info::desc() { return lsl_get_desc(obj.get()); }


// =============================
// ==== Continuous Resolver ====
// =============================

/**
 * A convenience class that resolves streams continuously in the background throughout
 * its lifetime and which can be queried at any time for the set of streams that are currently
 * visible on the network.
 */
class continuous_resolver {
public:
	/**
	 * Construct a new continuous_resolver that resolves all streams on the network.
	 *
	 * This is analogous to the functionality offered by the free function resolve_streams().
	 * @param forget_after When a stream is no longer visible on the network (e.g., because it was
	 * shut down), this is the time in seconds after which it is no longer reported by the resolver.
	 */
	continuous_resolver(double forget_after = 5.0)
		: obj(lsl_create_continuous_resolver(forget_after), &lsl_destroy_continuous_resolver) {}

	/**
	 * Construct a new continuous_resolver that resolves all streams with a specific value for a
	 * given property.
	 *
	 * This is analogous to the functionality provided by the free function resolve_stream(prop,value).
	 * @param prop The stream_info property that should have a specific value (e.g., "name", "type",
	 * "source_id", or "desc/manufaturer").
	 * @param value The string value that the property should have (e.g., "EEG" as the type
	 * property).
	 * @param forget_after When a stream is no longer visible on the network (e.g., because it was
	 * shut down), this is the time in seconds after which it is no longer reported by the resolver.
	 */
	continuous_resolver(
		const std::string &prop, const std::string &value, double forget_after = 5.0)
		: obj(lsl_create_continuous_resolver_byprop((prop.c_str()), (value.c_str()), forget_after),
			  &lsl_destroy_continuous_resolver) {}

	/**
	 * Construct a new continuous_resolver that resolves all streams that match a given XPath 1.0
	 * predicate.
	 *
	 * This is analogous to the functionality provided by the free function resolve_stream(pred).
	 * @param pred The predicate string, e.g.
	 * `name='BioSemi'` or
	 * `type='EEG' and starts-with(name,'BioSemi') and count(info/desc/channel)=32`
	 * @param forget_after When a stream is no longer visible on the network (e.g., because it was
	 * shut down), this is the time in seconds after which it is no longer reported by the resolver.
	 */
	continuous_resolver(const std::string &pred, double forget_after = 5.0)
		: obj(lsl_create_continuous_resolver_bypred((pred.c_str()), forget_after), &lsl_destroy_continuous_resolver) {}

	/**
	 * Obtain the set of currently present streams on the network (i.e. resolve result).
	 * @return A vector of matching stream info objects (excluding their meta-data), any of
	 *         which can subsequently be used to open an inlet.
	 */
	std::vector<stream_info> results() {
		lsl_streaminfo buffer[1024];
		return std::vector<stream_info>(
			buffer, buffer + check_error(lsl_resolver_results(obj.get(), buffer, sizeof(buffer))));
	}

	/// Move constructor for stream_inlet
	continuous_resolver(continuous_resolver &&rhs) noexcept = default;
	continuous_resolver &operator=(continuous_resolver &&rhs) noexcept = default;

private:
	std::unique_ptr<lsl_continuous_resolver_, void(*)(lsl_continuous_resolver_*)> obj;
};

} // namespace lsl

#endif // LSL_CPP_H
