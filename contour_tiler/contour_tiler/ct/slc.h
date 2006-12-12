#ifndef SLC_H
#define SLC_H

#define SLC_MAGIC 11111

/*      Data types supported by these algorithms.  C_SCANNED_DATA       */
/*      is used for medical, biological, and scan converted  data,      */
/*      while C_GEOMETRIC_DATA is used for data with geometric          */
/*      descriptions.                                                   */
typedef enum
{
        C_SCALAR_DATA_8BIT,
        C_SCALAR_DATA_16BIT,
        C_SCALAR_DATA_32BIT,
        C_GEOMETRIC_DATA
} C_DataType;

typedef enum
{
        C_BIORAD_CONFOCAL_DATA,
        C_MAGNETIC_RESONANCE_DATA,
        C_COMPUTED_TOMOGRAPHY_DATA,
        C_SIMULATION_DATA,
        C_BINARY_VOXELIZED_DATA,
        C_FUZZY_VOXELIZED_DATA,
        C_FUN_VOXELIZED_DATA,
        C_OTHER_DATA_ORIGIN
} C_DataOrigin;

typedef enum
{
        C_ORIGINAL_DATA,
        C_RESAMPLED_DATA,
        C_RESAMPLED_FILTERED_DATA,
        C_FILTERED_DATA,
        C_OTHER_DATA_MODIFICATION
} C_DataModification;

/*      File types which are supported by the load and save functions   */
typedef enum
{
        C_SLICE_FILE,
        C_TEXTURE_FILE,
        C_IMAGE_FILE,
        C_ANIMATION_FILE,
        C_FUNCTION_FILE,
        C_GEOMETRIC_FILE,
        C_ENVIRONMENT_FILE
} C_FileType;

/*      Unit Types Supported For Data Input.                            */
/*      C_METER_UNIT            - A Unit Refers To A Meter              */
/*      C_MILLIMETER_UNIT       - A Unit Refers To A Millimeter         */
/*      C_MICRON_UNIT           - A Unit Refers To A Micron Distance    */
/*      C_FOOT_UNIT             - A Unit Refers To A Foot Distance      */
/*      C_INCH_UNIT             - A Unit Refers To An Inch Distance     */
typedef enum
{
        C_METER_UNIT,
        C_MILLIMETER_UNIT,
        C_MICRON_UNIT,
        C_FOOT_UNIT,
        C_INCH_UNIT

} C_UnitType;

/*      Interpolation Types Supported During Scanned Data Input.        */
/*      C_NO_INTERP             - Do Not Perform Interpolation          */
/*      C_LINEAR_INTERP         - Linearly Interpolate Input Data       */
typedef enum
{
        C_NO_INTERP,
        C_LINEAR_INTERP
} C_InterpType;

typedef enum
{
        C_NO_COMPRESSION,
        C_RUN_LENGTH_ENCODE
} C_CompressionType;

#endif
