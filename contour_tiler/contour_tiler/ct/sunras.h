/* info on sun rasterfiles
 * taken from rasterfile(5) man page
 */

#define   RAS_MAGIC 0x59a66a95


struct rasterfile
{
     int  ras_magic;
     int  ras_width;
     int  ras_height;
     int  ras_depth;
     int  ras_length;
     int  ras_type;
     int  ras_maptype;
     int  ras_maplength;
};

#define RT_OLD          0       /* Raw pixrect image in 68000 byte order */
#define RT_STANDARD     1       /* Raw pixrect image in 68000 byte order */
#define RT_BYTE_ENCODED 2       /* Run-length compression of bytes */
#define RT_FORMAT_RGB   3       /* XRGB or RGB instead of XBGR or BGR */


#define RMT_RAW		2
#define RMT_NONE	0
#define RMT_EQUAL_RGB	1


#define RAS_RLE 0x80
