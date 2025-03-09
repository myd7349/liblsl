
/**
 * @note If the timeout expires before a new sample was received the function returns 0.0;
 * ec is *not* set to #lsl_timeout_error (because this case is not considered an error condition).
 * @return The capture time of the sample on the remote machine, or 0.0 if no new sample was
 * available. To remap this time stamp to the local clock, add the value returned by
 * lsl_time_correction() to it.
 * @{
 */

extern LIBLSL_C_API double lsl_pull_sample_str(lsl_inlet in, char **buffer, int32_t buffer_elements, double timeout, int32_t *ec);
///@}

/** @copydoc lsl_pull_sample_f
 * These strings may contains 0's, therefore the lengths are read into the buffer_lengths array.
 * @param buffer_lengths
 * A pointer to an array that holds the resulting lengths for each returned binary string.*/
extern LIBLSL_C_API double lsl_pull_sample_buf(lsl_inlet in, char **buffer, uint32_t *buffer_lengths, int32_t buffer_elements, double timeout, int32_t *ec);

/**
 * Pull a sample from the inlet and read it into a custom struct or buffer.
 *
 * Overall size checking but no type checking or conversion are done.
 * Do not use for variable-size/string-formatted streams.
 * @param in The #lsl_inlet object to act on.
 * @param[out] buffer A pointer to hold the resulting values.
 * @param buffer_bytes Length of the array held by buffer in bytes, not items
 * @param timeout The timeout for this operation, if any.
 * Use #LSL_FOREVER to effectively disable it. It is also permitted to use 0.0 here;
 * in this case a sample is only returned if one is currently buffered.
 * @param[out] ec Error code: can be either no error or #lsl_lost_error
 * (if the stream source has been lost).<br>
 * @note If the timeout expires before a new sample was received the function returns 0.0;
 * ec is *not* set to #lsl_timeout_error (because this case is not considered an error condition).
 * @return The capture time of the sample on the remote machine, or 0.0 if no new sample was
 * available. To remap this time stamp to the local clock, add the value returned by
 * lsl_time_correction() to it.
 */
extern LIBLSL_C_API double lsl_pull_sample_v(lsl_inlet in, void *buffer, int32_t buffer_bytes, double timeout, int32_t *ec);


extern LIBLSL_C_API unsigned long lsl_pull_chunk_str(lsl_inlet in, char **data_buffer, double *timestamp_buffer, unsigned long data_buffer_elements, unsigned long timestamp_buffer_elements, double timeout, int32_t *ec);

///@}

/**
 * Pull a chunk of data from the inlet and read it into an array of binary strings.
 *
 * These strings may contains 0's, therefore the lengths are read into the lengths_buffer array.
 * Handles type checking & conversion.
 * IMPORTANT: Note that the provided data buffer size is measured in channel values (e.g., floats)
 * rather than in samples.
 * @param in The lsl_inlet object to act on.
 * @param[out] data_buffer A pointer to a buffer of data values where the results shall be stored.
 * @param[out] lengths_buffer A pointer to an array that holds the resulting lengths for each
 * returned binary string.
 * @param timestamp_buffer A pointer to a buffer of timestamp values where time stamps shall be
 * stored. If this is NULL, no time stamps will be returned.
 * @param data_buffer_elements The size of the data buffer, in channel data elements (of type T).
 * Must be a multiple of the stream's channel count.
 * @param timestamp_buffer_elements The size of the timestamp buffer. If a timestamp buffer is
 * provided then this must correspond to the same number of samples as data_buffer_elements.
 * @param timeout The timeout for this operation, if any.
 *
 * When the timeout expires, the function may return before the entire buffer is filled.
 *
 * The default value of 0.0 will retrieve only data available for immediate pickup.
 * @param[out] ec Error code: can be either no error or #lsl_lost_error (if the stream source has
 * been lost).
 * @note If the timeout expires before a new sample was received the function returns 0.0; ec is
 * *not* set to #lsl_timeout_error (because this case is not considered an error condition).
 * @return data_elements_written Number of channel data elements written to the data buffer.
 */

extern LIBLSL_C_API unsigned long lsl_pull_chunk_buf(lsl_inlet in, char **data_buffer, uint32_t *lengths_buffer, double *timestamp_buffer, unsigned long data_buffer_elements, unsigned long timestamp_buffer_elements, double timeout, int32_t *ec);


