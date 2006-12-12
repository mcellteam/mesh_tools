/****************************************************
 *                                                  *
 *    DO NOT COPY OR DISTRIBUTE                     *
 *                                                  *
 *   Copyright (c) 1998 by C. Bajaj, G. Zhuang      *
 *                                                  *
 ****************************************************/



/***********************************************************************************
 * COPYRIGHT STATEMENT                                                             *
 *                                                                                 *
 * NOTE: This copyright statement applies to all individual files or               *
 *       modules in which the developers listed below are mentioned.               *
 *                                                                                 *
 * This software module was originally developed by NuTech Consultants Inc.,       *
 * West Lafayette, Indiana, USA. Contributors to the software include but are      *
 * not limited to: Chandrajit Bajaj, Guozhong Zhuang, Valerio Pascucci.            *
 * This software module is an implementation of a part of one or more              *
 * MPEG-4 Visual 3D Model Coding (ISO/IEC JTC1/SC29) tools as specified            *
 * by MPEG-4 Visual 3D Model Coding (ISO/IEC JTC1/SC29). ISO/IEC gives users       *
 * of the MPEG-4 Visual 3D Model Coding standard free license to this software     *
 * module or modifications thereof for use in hardware or software products        *
 * claiming conformance to the MPEG-4 Visual 3D Model Coding (ISO/IEC JTC1/SC29)   *
 * standard. Those intending to use this software module in hardware or software   *
 * products are advised that its use may infringe existing patents. The original   *
 * developers of this software and his/her company, the subsequent editors and     *
 * their companies, and ISO/IEC have no liability for use of this software module  *
 * or modifications thereof. Copyright is not released for non MPEG-4 Visual 3D    *
 * Model Coding standard conforming products. The original developers retain full  *
 * right to use the code for its own purpose, assign or donate the code for        *
 * non MPEG-4 Visual 3D Model Coding conforming products.                          *
 * This copyright notice must be included in all copies or derivative works.       *
 *                                                                                 *
 * Copyright (c) 1998                                                              *
 ***********************************************************************************/


#ifndef  _WRL2RAW_H_
#define  _WRL2RAW_H_

#include "header.h"

typedef struct _obj_vrml2 {
     int        num_positions;
     int        size_positions;
     float     *positions;
     int        num_triangles;
     int        size_triangles;
     int	   *triangles;
     /*int        num_colors;
     int        size_colors;
     float     *colors;
     int        num_ccolors;
     int        size_ccolors;
     int       *ccolors;
     int        num_normals;
     int        size_normals;
     float     *normals;
     int        num_cnormals;
     int        size_cnormals;
     int       *cnormals; 
     int        num_textures;
     int        size_textures;
     float     *textures;
     int        num_ctextures;
     int        size_ctextures;
     int       *ctextures;*/
} ObjVRML2;

#endif
