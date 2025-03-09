
/**
 * @attention It is the user's responsibility to either destroy the resulting streaminfo objects or
 * to pass them back to the LSL during during creation of an inlet.
 */
extern LIBLSL_C_API int32_t lsl_resolver_results(lsl_continuous_resolver res, lsl_streaminfo *buffer, uint32_t buffer_elements);


/// @}

/** @defgroup resolve Resolving streams on the network
 * @{*/

/**
 * @param[out] buffer A user-allocated buffer to hold the resolve results.
 * @attention It is the user's responsibility to either destroy the resulting streaminfo objects or
 * to pass them back to the LSL during during creation of an inlet.
 */
extern LIBLSL_C_API int32_t lsl_resolve_all(lsl_streaminfo *buffer, uint32_t buffer_elements, double wait_time);

/**
 * @param[out] buffer A user-allocated buffer to hold the resolve results.
 * @attention It is the user's responsibility to either destroy the resulting streaminfo objects or
 * to pass them back to the LSL during during creation of an inlet.
 */
extern LIBLSL_C_API int32_t lsl_resolve_byprop(lsl_streaminfo *buffer, uint32_t buffer_elements, const char *prop, const char *value, int32_t minimum, double timeout);

/**
 * Advanced query that allows to impose more conditions on the retrieved streams;
 * the given string is an [XPath 1.0 predicate](http://en.wikipedia.org/w/index.php?title=XPath_1.0)
 * for the `<info>` node (omitting the surrounding []'s)
 * @param[out] buffer A user-allocated buffer to hold the resolve results.
 * @attention It is the user's responsibility to either destroy the resulting streaminfo objects or
 * to pass them back to the LSL during during creation of an inlet.
 */
extern LIBLSL_C_API int32_t lsl_resolve_bypred(lsl_streaminfo *buffer, uint32_t buffer_elements, const char *pred, int32_t minimum, double timeout);

/// @}
