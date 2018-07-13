#!/usr/bin/env python
from __future__ import division, print_function
'''
*ReadMarCCDFrame: Read Mar Files*
---------------------------------

'''
""" 
  from /opt/marccd/documentation/header.txt

   MarCCD Header Documentataion 

   from C code in frame.h  and types.h

   Documentation updated by R. Doyle Mon Mar 22 15:04:00 CDT 2010
   Documentation updated by M. Blum Tue Oct 11 11:49:20 CDT 2005

   Description documents marccd v0.20.0


   Summary of file structure:
   |-- 1024 bytes TIFF HEADER -------------|
   |-- 3072 byte frame_header structure ---|
   |-- nfast*nslow*depth byte image -------|

   The full header, as written to the file, is a TIFF header.
   The initial 1024 bytes are a minimal TIFF header with a standard
   TIFF TAG pointing to the image data and a private TIFF TAG
   pointing to this header structure.  As written by mmx/marccd, the
   frame_header structure always begins at byte 1024 and is 3072 bytes long
   making the full header 4096 bytes.

   immediately following the header is the image - it is of arbitrary size defined
   by the header fields nfast, nslow and depth. The total size is
   nfast * nslow * depth bytes.

   The meanings of the data types should be self evident:
   (example:  UINT32 is an unsigned 32 bit integer)
   The exact C language definition is machine dependent but these
   are the most common definitions on a 32bit architecture cpu.
#define UINT16	unsigned short
#define INT16	short
#define UINT32	unsigned int
#define INT32	int

   Currently frames are always written as defined below:
	 origin=UPPER_LEFT
	 orientation=HFAST
	 view_direction=FROM_SOURCE


/* This number is  written into the byte_order fields in the
   native byte order of the machine writing the file */
#define LITTLE_ENDIAN	1234
#define BIG_ENDIAN	4321

/* possible orientations of frame data (stored in orienation field) */
#define HFAST			0	 /* Horizontal axis is fast */
#define VFAST			1	 /* Vertical axis is fast */

/* possible origins of frame data (stored in origin field) */
#define UPPER_LEFT		0
#define LOWER_LEFT		1
#define UPPER_RIGHT		2
#define LOWER_RIGHT		3

/* possible view directions of frame data for
   the given orientation and origin (stored in view_direction field) */
#define FROM_SOURCE		0
#define TOWARD_SOURCE		1

/* possible types of data (in data_type field) */
#define DATA_UNSIGNED_INTEGER   0
#define DATA_SIGNED_INTEGER     1
#define DATA_FLOAT              2

#define MAXIMAGES 9
#define MAXSUBIMAGES 4096
#define MAXFRAMEDIMENSION       8192

typedef struct frame_header_type {
	/* File/header format parameters (256 bytes) */
	UINT32        header_type;	/* flag for header type  (can be used as magic number) */
	char header_name[16];		/* header name (MARCCD) */
	UINT32        header_major_version;	/* header_major_version (n.) */
	UINT32        header_minor_version; 	/* header_minor_version (.n) */
	UINT32        header_byte_order;/* BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel) */
	UINT32        data_byte_order;	/* BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel) */
	UINT32        header_size;	/* in bytes			*/
	UINT32        frame_type;	/* flag for frame type */
	UINT32        magic_number;	/* to be used as a flag - usually to indicate new file */
	UINT32        compression_type;	/* type of image compression    */
	UINT32        compression1;	/* compression parameter 1 */
	UINT32        compression2;	/* compression parameter 2 */
	UINT32        compression3;	/* compression parameter 3 */
	UINT32        compression4;	/* compression parameter 4 */
	UINT32        compression5;	/* compression parameter 4 */
	UINT32        compression6;	/* compression parameter 4 */
	UINT32        nheaders;	        /* total number of headers 	*/
 	UINT32        nfast;		/* number of pixels in one line */
 	UINT32        nslow;		/* number of lines in image     */
 	UINT32        depth;		/* number of bytes per pixel    */
 	UINT32        record_length;	/* number of pixels between succesive rows */
 	UINT32        signif_bits;	/* true depth of data, in bits  */
 	UINT32        data_type;	/* (signed,unsigned,float...) */
 	UINT32        saturated_value;	/* value marks pixel as saturated */
	UINT32        sequence;	        /* TRUE or FALSE */
	UINT32        nimages;	        /* total number of images - size of each is nfast*(nslow/nimages) */
 	UINT32        origin;		/* corner of origin 		*/
 	UINT32        orientation;	/* direction of fast axis 	*/
        UINT32        view_direction;   /* direction to view frame      */
	UINT32        overflow_location;/* FOLLOWING_HEADER, FOLLOWING_DATA */
	UINT32        over_8_bits;	/* # of pixels with counts > 255 */
	UINT32        over_16_bits;	/* # of pixels with count > 65535 */
	UINT32        multiplexed;	/* multiplex flag */
	UINT32        nfastimages;	/* # of images in fast direction */
	UINT32        nslowimages;	/* # of images in slow direction */
	UINT32        darkcurrent_applied; /* flags correction has been applied - hold magic number ? */
	UINT32        bias_applied;	  /* flags correction has been applied - hold magic number ? */
	UINT32        flatfield_applied;  /* flags correction has been applied - hold magic number ? */
	UINT32        distortion_applied; /* flags correction has been applied - hold magic number ? */
	UINT32        original_header_type;	/* Header/frame type from file that frame is read from */
	UINT32        file_saved;         /* Flag that file has been saved, should be zeroed if modified */
	UINT32        n_valid_pixels;     /* Number of pixels holding valid data - first N pixels */
	UINT32        defectmap_applied; /* flags correction has been applied - hold magic number ? */
	UINT32        subimage_nfast; 	    /* when divided into subimages (eg. frameshifted) */
	UINT32        subimage_nslow;       /* when divided into subimages (eg. frameshifted) */
	UINT32        subimage_origin_fast; /* when divided into subimages (eg. frameshifted) */
	UINT32        subimage_origin_slow; /* when divided into subimages (eg. frameshifted) */
	UINT32        readout_pattern;      /* BIT Code - 1 = A, 2 = B, 4 = C, 8 = D */
 	UINT32        saturation_level;	    /* at this value and above, data are not reliable */
 	UINT32        orientation_code;	    /* Describes how this frame needs to be rotated to make it "right" */
 	UINT32        frameshift_multiplexed;  /* frameshift multiplex flag */
	UINT32        prescan_nfast;            /* Number of non-image pixels preceeding imaging pixels - fast direction */
	UINT32        prescan_nslow;            /* Number of non-image pixels preceeding imaging pixels - slow direction */
	UINT32        postscan_nfast;           /* Number of non-image pixels followng imaging pixels - fast direction */
	UINT32        postscan_nslow;           /* Number of non-image pixels followng imaging pixels - slow direction */
	UINT32        prepost_trimmed;          /* trimmed==1 means pre and post scan pixels have been removed */
	char reserve1[(64-55)*sizeof(INT32)-16];

	/* Data statistics (128) */
	UINT32        total_counts[2];	/* 64 bit integer range = 1.85E19*/
	UINT32        special_counts1[2];
	UINT32        special_counts2[2];
	UINT32        min;
	UINT32        max;
	INT32        mean;			/* mean * 1000 */
	UINT32        rms;			/* rms * 1000 */
	UINT32        n_zeros;			/* number of pixels with 0 value  - not included in stats in unsigned data */
	UINT32        n_saturated;		/* number of pixels with saturated value - not included in stats */
	UINT32        stats_uptodate;		/* Flag that stats OK - ie data not changed since last calculation */
        UINT32        pixel_noise[MAXIMAGES];		/* 1000*base noise value (ADUs) */
	char reserve2[(32-13-MAXIMAGES)*sizeof(INT32)];

	/* Sample Changer info */
	char          barcode[16];
	UINT32        barcode_angle;
	UINT32        barcode_status;
	/* Pad to 256 bytes */
	char reserve2a[(64-6)*sizeof(INT32)];

	/* Goniostat parameters (128 bytes) */
        INT32 xtal_to_detector;		/* 1000*distance in millimeters */
        INT32 beam_x;			/* 1000*x beam position (pixels) */
        INT32 beam_y;			/* 1000*y beam position (pixels) */
        INT32 integration_time;		/* integration time in milliseconds */
        INT32 exposure_time;		/* exposure time in milliseconds */
        INT32 readout_time;		/* readout time in milliseconds */
        INT32 nreads;			/* number of readouts to get this image */
        INT32 start_twotheta;		/* 1000*two_theta angle */
        INT32 start_omega;		/* 1000*omega angle */
        INT32 start_chi;			/* 1000*chi angle */
        INT32 start_kappa;		/* 1000*kappa angle */
        INT32 start_phi;			/* 1000*phi angle */
        INT32 start_delta;		/* 1000*delta angle */
        INT32 start_gamma;		/* 1000*gamma angle */
        INT32 start_xtal_to_detector;	/* 1000*distance in mm (dist in um)*/
        INT32 end_twotheta;		/* 1000*two_theta angle */
        INT32 end_omega;			/* 1000*omega angle */
        INT32 end_chi;			/* 1000*chi angle */
        INT32 end_kappa;			/* 1000*kappa angle */
        INT32 end_phi;			/* 1000*phi angle */
        INT32 end_delta;			/* 1000*delta angle */
        INT32 end_gamma;			/* 1000*gamma angle */
        INT32 end_xtal_to_detector;	/* 1000*distance in mm (dist in um)*/
        INT32 rotation_axis;		/* active rotation axis (index into above ie. 0=twotheta,1=omega...) */
        INT32 rotation_range;		/* 1000*rotation angle */
        INT32 detector_rotx;		/* 1000*rotation of detector around X */
        INT32 detector_roty;		/* 1000*rotation of detector around Y */
        INT32 detector_rotz;		/* 1000*rotation of detector around Z */
        INT32 total_dose;		/* Hz-sec (counts) integrated over full exposure */
	char reserve3[(32-29)*sizeof(INT32)]; /* Pad Gonisotat parameters to 128 bytes */

	/* Detector parameters (128 bytes) */
	INT32 detector_type;		/* detector type */
	INT32 pixelsize_x;		/* pixel size (nanometers) */
	INT32 pixelsize_y;		/* pixel size (nanometers) */
        INT32 mean_bias;			/* 1000*mean bias value */
        INT32 photons_per_100adu;	/* photons / 100 ADUs */
        INT32 measured_bias[MAXIMAGES];	/* 1000*mean bias value for each image*/
        INT32 measured_temperature[MAXIMAGES];	/* Temperature of each detector in milliKelvins */
        INT32 measured_pressure[MAXIMAGES];	/* Pressure of each chamber in microTorr */
	/* Retired reserve4 when MAXIMAGES set to 9 from 16 and two fields removed, and temp and pressure added
	char reserve4[(32-(5+3*MAXIMAGES))*sizeof(INT32)];
	*/

	/* X-ray source and optics parameters (128 bytes) */
	/* X-ray source parameters (14*4 bytes) */
        INT32 source_type;		/* (code) - target, synch. etc */
        INT32 source_dx;			/* Optics param. - (size microns) */
        INT32 source_dy;			/* Optics param. - (size microns) */
        INT32 source_wavelength;		/* wavelength (femtoMeters) */
        INT32 source_power;		/* (Watts) */
        INT32 source_voltage;		/* (Volts) */
        INT32 source_current;		/* (microAmps) */
        INT32 source_bias;		/* (Volts) */
        INT32 source_polarization_x;	/* () */
        INT32 source_polarization_y;	/* () */
        INT32 source_intensity_0;	/* (arbitrary units) */
        INT32 source_intensity_1;	/* (arbitrary units) */
	char reserve_source[2*sizeof(INT32)];

	/* X-ray optics_parameters (8*4 bytes) */
        INT32 optics_type;		/* Optics type (code)*/
        INT32 optics_dx;			/* Optics param. - (size microns) */
        INT32 optics_dy;			/* Optics param. - (size microns) */
        INT32 optics_wavelength;		/* Optics param. - (size microns) */
        INT32 optics_dispersion;		/* Optics param. - (*10E6) */
        INT32 optics_crossfire_x;	/* Optics param. - (microRadians) */
        INT32 optics_crossfire_y;	/* Optics param. - (microRadians) */
        INT32 optics_angle;		/* Optics param. - (monoch. 2theta - microradians) */
        INT32 optics_polarization_x;	/* () */
        INT32 optics_polarization_y;	/* () */
	char reserve_optics[4*sizeof(INT32)];

	char reserve5[((32-28)*sizeof(INT32))]; /* Pad X-ray parameters to 128 bytes */

	/* File parameters (1024 bytes) */
	char filetitle[128];		/* Title 				*/
	char filepath[128];		/* path name for data file		*/
	char filename[64];		/* name of data file			*/
        char acquire_timestamp[32];	/* date and time of acquisition		*/
        char header_timestamp[32];	/* date and time of header update	*/
        char save_timestamp[32];	/* date and time file saved 		*/
        char file_comment[512];	/* comments  - can be used as desired 	*/
	char reserve6[1024-(128+128+64+(3*32)+512)]; /* Pad File parameters to 1024 bytes */

	/* Dataset parameters (512 bytes) */
        char dataset_comment[512];	/* comments  - can be used as desired 	*/

	/* Reserved for user definable data - will not be used by Mar! */
	char user_data[512];

	/* char pad[----] USED UP! */     /* pad out to 3072 bytes */

	} frame_header;
"""

import struct as st
import array as ar

MAXIMAGES=9

class marFrame():
    '''A class to extract correct mar header and image info from a MarCCD file

    :param str File: file object [from open()]
    :param byteOrd: '<' (default) or '>'  
    :param dict IFD: ?
    '''
    def __init__(self,File,byteOrd='<',IFD={}):
        # simple TIFF header info
        self.TIFFsizeX = IFD[256][2][0]
        self.TIFFsizeY = IFD[257][2][0]
        self.TIFFbitDepth = IFD[258][2][0]
        self.TIFFcompression = IFD[259][2][0] # 1 = no compression
        self.TIFFphotometricInterpretation = IFD[262][2][0] # 1 = bilevel or grayscale where 0 is imaged as black
        self.TIFFstripOffsets = IFD[273][2][0] # seems to be 4096 for marCCD
        self.TIFForientation = IFD[274][2][0] # 1 = 0th row it top, 0th column is left
        self.TIFFrowsPerStrip = IFD[278][2][0] # varies based on image size
        self.TIFFstripByteCounts = IFD[279][2][0] # number of bytes in a strip also varies based on size
        self.TIFFxResolution = IFD[282][2][0] # pixels per resolutionUnit in X direction (ImageWidth direction)
        self.TIFFyResolution = IFD[283][2][0] # pixels per resolutionUnit in Y direction (ImageLength direction
        self.TIFFresolutionUnit = IFD[296][2][0] # 3 = centimeter
        self.byteDepth = self.TIFFbitDepth//8
        self.arrayTypeCode = ['','B','H','I','I'][self.byteDepth]
        # MarCCD specific header info
        File.seek(IFD[34710][2][0])
        self.headerType = st.unpack(byteOrd+'I',File.read(4))[0] #/* flag for header type  (can be used as magic number) */
        self.headerName = b''.join(st.unpack(byteOrd+16*'s',File.read(16))).replace(b'\x00',b'')
        self.headerMajorVersion = st.unpack(byteOrd+'I',File.read(4))[0] #/* header_major_version (n.) */
        self.headerMinorVersion = st.unpack(byteOrd+'I',File.read(4))[0] #/* header_minor_version (.n) */
        self.headerByteOrder = st.unpack(byteOrd+'I',File.read(4))[0] #/* BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel) */
        self.dataByteOrder = st.unpack(byteOrd+'I',File.read(4))[0] #/* BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel) */
        self.headerSize = st.unpack(byteOrd+'I',File.read(4))[0] #/* in bytes			*/
        self.frameType = st.unpack(byteOrd+'I',File.read(4))[0] #/* flag for frame type */
        self.magicNumber = st.unpack(byteOrd+'I',File.read(4))[0] #/* to be used as a flag - usually to indicate new file */
        self.compressionType = st.unpack(byteOrd+'I',File.read(4))[0] #/* type of image compression    */
        self.compression1 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 1 */
        self.compression2 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 2 */
        self.compression3 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 3 */
        self.compression4 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 4 */
        self.compression5 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 4 */
        self.compression6 = st.unpack(byteOrd+'I',File.read(4))[0] #/* compression parameter 4 */
        self.nheaders = st.unpack(byteOrd+'I',File.read(4))[0] #/* total number of headers 	*/
        self.nfast = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of pixels in one line */
        self.nslow = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of lines in image     */
        self.depth = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of bytes per pixel    */
        self.recordLength = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of pixels between succesive rows */
        self.signifBits = st.unpack(byteOrd+'I',File.read(4))[0] #/* true depth of data, in bits  */
        self.dataType = st.unpack(byteOrd+'I',File.read(4))[0] #/* (signed,unsigned,float...) */
        self.saturatedValue = st.unpack(byteOrd+'I',File.read(4))[0] #/* value marks pixel as saturated */
        self.sequence = st.unpack(byteOrd+'I',File.read(4))[0] #/* TRUE or FALSE */
        self.nimages = st.unpack(byteOrd+'I',File.read(4))[0] #/* total number of images - size of each is nfast*(nslow/nimages) */
        self.origin = st.unpack(byteOrd+'I',File.read(4))[0] #/* corner of origin 		*/
        self.orientation = st.unpack(byteOrd+'I',File.read(4))[0] #/* direction of fast axis 	*/
        self.viewDirection = st.unpack(byteOrd+'I',File.read(4))[0] #/* direction to view frame      */
        self.overflowLocation = st.unpack(byteOrd+'I',File.read(4))[0] #/* FOLLOWING_HEADER, FOLLOWING_DATA */
        self.over8Bits = st.unpack(byteOrd+'I',File.read(4))[0] #/* # of pixels with counts > 255 */
        self.over16Bits = st.unpack(byteOrd+'I',File.read(4))[0] #/* # of pixels with count > 65535 */
        self.multiplexed = st.unpack(byteOrd+'I',File.read(4))[0] #/* multiplex flag */
        self.nfastimages = st.unpack(byteOrd+'I',File.read(4))[0] #/* # of images in fast direction */
        self.nslowimages = st.unpack(byteOrd+'I',File.read(4))[0] #/* # of images in slow direction */
        self.darkcurrentApplied = st.unpack(byteOrd+'I',File.read(4))[0] #/* flags correction has been applied - hold magic number ? */
        self.biasApplied = st.unpack(byteOrd+'I',File.read(4))[0] #/* flags correction has been applied - hold magic number ? */
        self.flatfieldApplied = st.unpack(byteOrd+'I',File.read(4))[0] #/* flags correction has been applied - hold magic number ? */
        self.distortionApplied = st.unpack(byteOrd+'I',File.read(4))[0] #/* flags correction has been applied - hold magic number ? */
        self.originalHeaderType = st.unpack(byteOrd+'I',File.read(4))[0] #/* Header/frame type from file that frame is read from */
        self.fileSaved = st.unpack(byteOrd+'I',File.read(4))[0] #/* Flag that file has been saved, should be zeroed if modified */
        self.nValidPixels = st.unpack(byteOrd+'I',File.read(4))[0] #/* Number of pixels holding valid data - first N pixels */
        self.defectmapApplied = st.unpack(byteOrd+'I',File.read(4))[0] #/* flags correction has been applied - hold magic number ? */
        self.subimageNfast = st.unpack(byteOrd+'I',File.read(4))[0] #/* when divided into subimages (eg. frameshifted) */
        self.subimageNslow = st.unpack(byteOrd+'I',File.read(4))[0] #/* when divided into subimages (eg. frameshifted) */
        self.subimageOriginFast = st.unpack(byteOrd+'I',File.read(4))[0] #/* when divided into subimages (eg. frameshifted) */
        self.subimageOriginSlow = st.unpack(byteOrd+'I',File.read(4))[0] #/* when divided into subimages (eg. frameshifted) */
        self.readoutPattern = st.unpack(byteOrd+'I',File.read(4))[0] #/* BIT Code - 1 = A, 2 = B, 4 = C, 8 = D */
        self.saturationLevel = st.unpack(byteOrd+'I',File.read(4))[0] #/* at this value and above, data are not reliable */
        self.orientationCode = st.unpack(byteOrd+'I',File.read(4))[0] #/* Describes how this frame needs to be rotated to make it "right" */
        self.frameshiftMultiplexed = st.unpack(byteOrd+'I',File.read(4))[0] #/* frameshift multiplex flag */
        self.prescanNfast = st.unpack(byteOrd+'I',File.read(4))[0] #/* Number of non-image pixels preceeding imaging pixels - fast direction */
        self.prescanNslow = st.unpack(byteOrd+'I',File.read(4))[0] #/* Number of non-image pixels preceeding imaging pixels - slow direction */
        self.postscanNfast = st.unpack(byteOrd+'I',File.read(4))[0] #/* Number of non-image pixels followng imaging pixels - fast direction */
        self.postscanNslow = st.unpack(byteOrd+'I',File.read(4))[0] #/* Number of non-image pixels followng imaging pixels - slow direction */
        self.prepostTrimmed = st.unpack(byteOrd+'I',File.read(4))[0] #/* trimmed==1 means pre and post scan pixels have been removed */

        File.seek(IFD[34710][2][0]+256)
        #self.totalCounts = st.unpack(byteOrd+'Q',File.read(8))[0] # /* 64 bit integer range = 1.85E19*/
        #self.specialCounts1 = st.unpack(byteOrd+'Q',File.read(8))[0]
        #self.specialCounts2 = st.unpack(byteOrd+'Q',File.read(8))[0]
        self.totalCounts = st.unpack(byteOrd+'II',File.read(8))
        self.specialCounts1 = st.unpack(byteOrd+'II',File.read(8))
        self.specialCounts2 = st.unpack(byteOrd+'II',File.read(8))
        self.min = st.unpack(byteOrd+'I',File.read(4))[0]
        self.max = st.unpack(byteOrd+'I',File.read(4))[0]
        self.mean = st.unpack(byteOrd+'i',File.read(4))[0] # /* mean * 1000 */
        self.rms = st.unpack(byteOrd+'I',File.read(4))[0] #/* rms * 1000 */
        self.nZeros = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of pixels with 0 value  - not included in stats in unsigned data */
        self.nSaturated = st.unpack(byteOrd+'I',File.read(4))[0] #/* number of pixels with saturated value - not included in stats */
        self.statsUptodate = st.unpack(byteOrd+'I',File.read(4))[0] #/* Flag that stats OK - ie data not changed since last calculation */
        self.pixelNoise = st.unpack(byteOrd+'I'*MAXIMAGES,File.read(4*MAXIMAGES)) # /* 1000*base noise value (ADUs) */

        File.seek(IFD[34710][2][0]+256+128)
        self.barcode = b''.join(st.unpack(byteOrd+16*'s',File.read(16))).replace(b'\x00',b'')
        self.barcodeAngle = st.unpack(byteOrd+'I',File.read(4))[0]
        self.barcodeStatus = st.unpack(byteOrd+'I',File.read(4))[0]

        File.seek(IFD[34710][2][0]+256+128+256)
        self.xtalToDetector = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*distance in millimeters */
        self.beamX = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*x beam position (pixels) */
        self.beamY = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*y beam position (pixels) */
        self.integrationTime = st.unpack(byteOrd+'i',File.read(4))[0] #/* integration time in milliseconds */
        self.exposureTime = st.unpack(byteOrd+'i',File.read(4))[0] #/* exposure time in milliseconds */
        self.readoutTime = st.unpack(byteOrd+'i',File.read(4))[0] #/* readout time in milliseconds */
        self.nreads = st.unpack(byteOrd+'i',File.read(4))[0] #/* number of readouts to get this image */
        self.startTwotheta = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*two_theta angle */
        self.startOmega = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*omega angle */
        self.startChi = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*chi angle */
        self.startKappa = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*kappa angle */
        self.startPhi = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*phi angle */
        self.startDelta = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*delta angle */
        self.startGamma = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*gamma angle */
        self.startXtalToDetector = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*distance in mm (dist in um)*/
        self.endTwotheta = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*two_theta angle */
        self.endOmega = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*omega angle */
        self.endChi = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*chi angle */
        self.endKappa = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*kappa angle */
        self.endPhi = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*phi angle */
        self.endDelta = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*delta angle */
        self.endGamma = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*gamma angle */
        self.endXtalToDetector = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*distance in mm (dist in um)*/
        self.rotationAxis = st.unpack(byteOrd+'i',File.read(4))[0] #/* active rotation axis (index into above ie. 0=twotheta,1=omega...) */
        self.rotationRange = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*rotation angle */
        self.detectorRotx = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*rotation of detector around X */
        self.detectorRoty = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*rotation of detector around Y */
        self.detectorRotz = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*rotation of detector around Z */
        self.totalDose = st.unpack(byteOrd+'i',File.read(4))[0] #/* Hz-sec (counts) integrated over full exposure */

        File.seek(IFD[34710][2][0]+256+128+256+128)
        self.detectorType = st.unpack(byteOrd+'i',File.read(4))[0] #/* detector type */
        self.pixelsizeX = st.unpack(byteOrd+'i',File.read(4))[0] #/* pixel size (nanometers) */
        self.pixelsizeY = st.unpack(byteOrd+'i',File.read(4))[0] #/* pixel size (nanometers) */
        self.meanBias = st.unpack(byteOrd+'i',File.read(4))[0] #/* 1000*mean bias value */
        self.photonsPer100adu = st.unpack(byteOrd+'i',File.read(4))[0] #/* photons / 100 ADUs */
        self.measuredBias = st.unpack(byteOrd+'i'*MAXIMAGES,File.read(4*MAXIMAGES)) # /* 1000*mean bias value for each image*/
        self.measuredTemperature = st.unpack(byteOrd+'i'*MAXIMAGES,File.read(4*MAXIMAGES)) # /* Temperature of each detector in milliKelvins */
        self.measuredPressure = st.unpack(byteOrd+'i'*MAXIMAGES,File.read(4*MAXIMAGES)) # /* Pressure of each chamber in microTorr */

        File.seek(IFD[34710][2][0]+256+128+256+128+128)
        self.sourceType = st.unpack(byteOrd+'i',File.read(4))[0] #/* (code) - target, synch. etc */
        self.sourceDx = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (size microns) */
        self.sourceDy = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (size microns) */
        self.sourceWavelength = st.unpack(byteOrd+'i',File.read(4))[0] #/* wavelength (femtoMeters) */
        self.sourcePower = st.unpack(byteOrd+'i',File.read(4))[0] #/* (Watts) */
        self.sourceVoltage = st.unpack(byteOrd+'i',File.read(4))[0] #/* (Volts) */
        self.sourceCurrent = st.unpack(byteOrd+'i',File.read(4))[0] #/* (microAmps) */
        self.sourceBias = st.unpack(byteOrd+'i',File.read(4))[0] #/* (Volts) */
        self.sourcePolarizationX = st.unpack(byteOrd+'i',File.read(4))[0] #/* () */
        self.sourcePolarizationY = st.unpack(byteOrd+'i',File.read(4))[0] #/* () */
        self.sourceIntensity0 = st.unpack(byteOrd+'i',File.read(4))[0] #/* (arbitrary units) */
        self.sourceIntensity1 = st.unpack(byteOrd+'i',File.read(4))[0] #/* (arbitrary units) */

        File.seek(IFD[34710][2][0]+256+128+256+128+128+14*4)
        self.opticsType = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics type (code)*/
        self.opticsDx = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (size microns) */
        self.opticsDy = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (size microns) */
        self.opticsWavelength = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (size microns) */
        self.opticsDispersion = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (*10E6) */
        self.opticsCrossfireX = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (microRadians) */
        self.opticsCrossfireY = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (microRadians) */
        self.opticsAngle = st.unpack(byteOrd+'i',File.read(4))[0] #/* Optics param. - (monoch. 2theta - microradians) */
        self.opticsPolarizationX = st.unpack(byteOrd+'i',File.read(4))[0] #/* () */
        self.opticsPolarizationY = st.unpack(byteOrd+'i',File.read(4))[0] #/* () */

        File.seek(IFD[34710][2][0]+256+128+256+128+128+128)
        self.filetitle = b''.join(st.unpack(byteOrd+128*'s',File.read(128))).replace(b'\x00',b'')
        self.filepath = b''.join(st.unpack(byteOrd+'s'*128,File.read(128))).replace(b'\x00',b'') #/* path name for data file*/
        self.filename = b''.join(st.unpack(byteOrd+'s'*64,File.read(64))).replace(b'\x00',b'') #/* name of data file*/
        self.acquireTimestamp = b''.join(st.unpack(byteOrd+'s'*32,File.read(32))).replace(b'\x00',b'') #/* date and time of acquisition*/
        self.headerTimestamp = b''.join(st.unpack(byteOrd+'s'*32,File.read(32))).replace(b'\x00',b'') #/* date and time of header update*/
        self.saveTimestamp = b''.join(st.unpack(byteOrd+'s'*32,File.read(32))).replace(b'\x00',b'') #/* date and time file saved */
        self.fileComment = b''.join(st.unpack(byteOrd+'s'*512,File.read(512))).replace(b'\x00',b'') #/* comments  - can be used as desired */
        self.datasetComment = b''.join(st.unpack(byteOrd+'s'*512,File.read(512))).replace(b'\x00',b'') #/* comments  - can be used as desired */

        self.userData = b''.join(st.unpack(byteOrd+'s'*512,File.read(512))).replace(b'\x00',b'')

        File.seek(4096)
        self.image = ar.array(self.arrayTypeCode,File.read(self.byteDepth*self.TIFFsizeX*self.TIFFsizeY))
    # reverse the array so if can have the same view as is read in marccd
    # also switch the view direction 
        self.image.reverse()
        self.viewDirection = abs(self.viewDirection - 1)

    def outputHead(self):
        myHead = []
        for curAttr in dir(self):
            if ( curAttr != '__doc__' and \
                 curAttr != '__init__' and \
                 curAttr != '__module__' and \
                 curAttr != 'outputHead' and \
                 curAttr != 'image' ):
                myHead.append(" %s = %s" % (curAttr,getattr(self,curAttr)))
        return myHead
