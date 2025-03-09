#pragma once
#include "./common.h"
#include "types.h"

extern LIBLSL_C_API lsl_outlet lsl_create_outlet(lsl_streaminfo info, int32_t chunk_size, int32_t max_buffered);
/** @copydoc lsl_create_outlet()
 * @param flags An integer that is the result of bitwise OR'ing one or more options from
 * #lsl_transport_options_t together (e.g., #transp_bufsize_samples|#transp_bufsize_thousandths)
 */
extern LIBLSL_C_API lsl_outlet lsl_create_outlet_ex(
	lsl_streaminfo info, int32_t chunk_size, int32_t max_buffered, lsl_transport_options_t flags);

/**
 * @param data A pointer to values to push. The number of values pointed to must be no less than the
 * number of channels in the sample.
 * #lsl_local_clock(); if omitted, the current time is used.
 * @return Error code of the operation or lsl_no_error if successful (usually attributed to the
 * wrong data type).
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_sample_str(lsl_outlet out, const char **data);
extern LIBLSL_C_API int32_t lsl_push_sample_v(lsl_outlet out, const void *data);
/// @}
/** @copydoc lsl_push_sample_f
 * @param timestamp Optionally the capture time of the sample, in agreement with
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_sample_strt(lsl_outlet out, const char **data, double timestamp);
extern LIBLSL_C_API int32_t lsl_push_sample_vt(lsl_outlet out, const void *data, double timestamp);
/// @}
/** @copydoc lsl_push_sample_ft
 * @param pushthrough Whether to push the sample through to the receivers instead of buffering it
 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
 * precedence over the pushthrough flag.
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_sample_strtp(lsl_outlet out, const char **data, double timestamp, int32_t pushthrough);
extern LIBLSL_C_API int32_t lsl_push_sample_vtp(lsl_outlet out, const void *data, double timestamp, int32_t pushthrough);
///@}

/** @copybrief lsl_push_sample_ftp
 * @see lsl_push_sample_ftp
 * @param out The lsl_outlet object through which to push the data.
 * @param data A pointer to values to push. The number of values pointed to must be no less than the
 * number of channels in the sample.
 * @param lengths A pointer the number of elements to push for each channel (string lengths).
 */
extern LIBLSL_C_API int32_t lsl_push_sample_buf(lsl_outlet out, const char **data, const uint32_t *lengths);
/** @copydoc lsl_push_sample_buf
 * @param timestamp @see lsl_push_sample_ftp */
extern LIBLSL_C_API int32_t lsl_push_sample_buft(lsl_outlet out, const char **data, const uint32_t *lengths, double timestamp);
/** @copydoc lsl_push_sample_buft
 * @param pushthrough @see lsl_push_sample_ftp */
extern LIBLSL_C_API int32_t lsl_push_sample_buftp(lsl_outlet out, const char **data, const uint32_t *lengths, double timestamp, int32_t pushthrough);

/** One timestamp per sample is provided.
 *
 * @attention Note that the provided buffer size is measured in channel values (e.g. floats) rather
 * than in samples.
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_chunk_str(lsl_outlet out, const char **data, unsigned long data_elements);
/// @}

/** @copydoc lsl_push_chunk_f
 * @param timestamp Optionally the capture time of the most recent sample, in agreement with
 * lsl_local_clock(); if omitted, the current time is used.
 * The time stamps of other samples are automatically derived based on the sampling rate of the
 * stream.
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_chunk_strt(lsl_outlet out, const char **data, unsigned long data_elements, double timestamp);
/// @}

/** @copydoc lsl_push_chunk_ft
 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
 * precedence over the pushthrough flag.
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_chunk_strtp(lsl_outlet out, const char **data, unsigned long data_elements, double timestamp, int32_t pushthrough);

extern LIBLSL_C_API int32_t lsl_push_chunk_strtn(lsl_outlet out, const char **data, unsigned long data_elements, const double *timestamps);
/// @}

/** @copydoc lsl_push_chunk_ftn
 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
 * precedence over the pushthrough flag.
 * @{
 */
extern LIBLSL_C_API int32_t lsl_push_chunk_strtnp(lsl_outlet out, const char **data, unsigned long data_elements, const double *timestamps, int32_t pushthrough);
///@}

/** @copybrief lsl_push_chunk_ftp
 * @sa lsl_push_chunk_ftp
 * @param out The lsl_outlet object through which to push the data.
 * @param data An array of channel values holding the data to push.
 * @param lengths Pointer the number of elements to push for each value (string lengths) so that
 * `size(data[i])==lengths[i]`.
 * @param data_elements The number of data values in the data buffer.
 * Must be a multiple of the channel count.
 */
extern LIBLSL_C_API int32_t lsl_push_chunk_buf(lsl_outlet out, const char **data, const uint32_t *lengths, unsigned long data_elements);

/** @copydoc lsl_push_chunk_buf @sa lsl_push_chunk_ftp @sa lsl_push_chunk_buf
 * @param timestamp Optionally the capture time of the most recent sample, in agreement with
 * lsl_local_clock(); if omitted, the current time is used.
 * The time stamps of other samples are automatically derived based on the sampling rate of the
 * stream. */
extern LIBLSL_C_API int32_t lsl_push_chunk_buft(lsl_outlet out, const char **data, const uint32_t *lengths, unsigned long data_elements, double timestamp);

/** @copydoc lsl_push_chunk_buft @sa lsl_push_chunk_ftp @sa lsl_push_chunk_buf
 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
 * precedence over the pushthrough flag. */
extern LIBLSL_C_API int32_t lsl_push_chunk_buftp(lsl_outlet out, const char **data, const uint32_t *lengths, unsigned long data_elements, double timestamp, int32_t pushthrough);

/** @copydoc lsl_push_chunk_buf @sa lsl_push_chunk_ftp @sa lsl_push_chunk_buf
 * @param timestamps Buffer holding one time stamp for each sample in the data buffer. */
extern LIBLSL_C_API int32_t lsl_push_chunk_buftn(lsl_outlet out, const char **data, const uint32_t *lengths, unsigned long data_elements, const double *timestamps);

/** @copydoc lsl_push_chunk_buftn @sa lsl_push_chunk_ftp @sa lsl_push_chunk_buf
 * @param pushthrough Whether to push the chunk through to the receivers instead of buffering it
 * with subsequent samples. Note that the chunk_size, if specified at outlet construction, takes
 * precedence over the pushthrough flag. */
extern LIBLSL_C_API int32_t lsl_push_chunk_buftnp(lsl_outlet out, const char **data, const uint32_t *lengths, unsigned long data_elements, const double *timestamps, int32_t pushthrough);

/**
 * Retrieve a handle to the stream info provided by this outlet.
 * This is what was used to create the stream (and also has the Additional Network Information
 * fields assigned).
 * @return A copy of the streaminfo of the outlet or NULL in the event that an error occurred.
 * @note It is the user's responsibility to destroy it when it is no longer needed.
 * @sa lsl_destroy_string()
 */
extern LIBLSL_C_API lsl_streaminfo lsl_get_info(lsl_outlet out);

///@}
