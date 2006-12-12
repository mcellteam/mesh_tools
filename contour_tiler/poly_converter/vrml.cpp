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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header.h"
#include "readwrl.h"

#define WRLCHUNK    2000
#define VRML_HEADER "#VRML V2.0 utf8"

extern tempTri* tris;
extern tempVert* verts;
extern int no_of_tris, orig_verts;

int Get_A_Token(FILE* fp, char token[])
{
    char c;
	int ns, len_token;
	
    do 
	{
        len_token = 0;
        token[0] = 0;
        while(1) 
		{
            ns = fscanf(fp,"%c",&c);
            if (ns==0) 
                break;
            if (!(c==' '||c=='\t'||c=='\n'||c==','||c==13)) 
			{
                token[len_token] = c;
                len_token++;
                break;
            } 
        }
		
        while(1) 
		{
            ns = fscanf(fp,"%c",&c);
            if (ns==0) 
                break;
			if (c!=' '&&c!='\t'&&c!='\n'&&c!=','&&c!=13) 
			{
                token[len_token] = c;
                len_token++;
            } 
			else 
				break;
        }
		
        if (len_token>0 && token[0]=='#') 
		{	
			while(c!='\n')
				fscanf(fp,"%c",&c);
        }
        token[len_token] = '\0';
    } while(token[0]=='#');
    return len_token;
}


void ObjVRML2_Init (ObjVRML2 *obj)
{
	obj->num_positions   = 0;
    obj->size_positions  = 0;
	obj->positions = NULL;
    obj->num_triangles   = 0;
    obj->size_triangles  = 0;
	obj->triangles = NULL;
    /*obj->num_colors      = 0;
    obj->size_colors     = 0;
	obj->colors = NULL;
    obj->num_ccolors     = 0;
    obj->size_ccolors    = 0;
	obj->ccolors = NULL;
    obj->num_normals     = 0;
    obj->size_normals    = 0;
	obj->normals = NULL;
    obj->num_cnormals    = 0;
    obj->size_cnormals   = 0;
	obj->cnormals = NULL;
    obj->num_textures    = 0;
    obj->size_textures   = 0;
	obj->textures = NULL;
    obj->num_ctextures   = 0;
    obj->size_ctextures  = 0;
	obj->ctextures = NULL;*/
}

void ErrMessage()
{
	printf("Out of memory\n");
	exit(1);
}

void Wrl_Add_Color(float fval, ObjVRML2 *obj)
{
    /*if (obj->num_colors==0) 
	{
		obj->colors = (float *) malloc(WRLCHUNK*sizeof(float));
		if (obj->colors==NULL)
			ErrMessage();
		else
			obj->size_colors = WRLCHUNK;
	} 
	else 
	{
		if (obj->num_colors==obj->size_colors) 
		{
			obj->colors = (float *) realloc(obj->colors, 
						(obj->size_colors + WRLCHUNK)*sizeof(float));
			if (obj->colors==NULL)
				ErrMessage();
			obj->size_colors += WRLCHUNK;
		}
	}

    obj->colors[obj->num_colors] = fval;
    obj->num_colors += 1;*/
}

void Wrl_Add_Coord(float fval, ObjVRML2 *obj)
{
    if (obj->num_positions==0) 
	{
		obj->positions = (float *) malloc(WRLCHUNK*sizeof(float));
		if (obj->positions==NULL)
			ErrMessage();
		else
			obj->size_positions = WRLCHUNK;
	}
	else 
	{
		if (obj->num_positions==obj->size_positions) 
		{
			obj->positions = (float *) realloc(obj->positions, 
						(obj->size_positions + WRLCHUNK)*sizeof(float));
			if (obj->positions==NULL)
				ErrMessage();
			obj->size_positions += WRLCHUNK;
		}
    }

    obj->positions[obj->num_positions] = fval;
    obj->num_positions += 1;
}               

void Wrl_Add_Normal(float fval, ObjVRML2 *obj)
{
    /*if (obj->num_normals==0) 
	{
		obj->normals = (float *) malloc(WRLCHUNK*sizeof(float));
		if (obj->normals==NULL)
			ErrMessage();
		else
			obj->size_normals = WRLCHUNK;
    } 
	else 
	{
		if (obj->num_normals==obj->size_normals) 
		{
			obj->normals = (float *) realloc(obj->normals, 
						(obj->size_normals + WRLCHUNK)*sizeof(float));
			if (obj->normals==NULL)
				ErrMessage();
			obj->size_normals += WRLCHUNK;
		}
    }

    obj->normals[obj->num_normals] = fval;
    obj->num_normals += 1;*/                     
}


void Wrl_Add_TexCoord(float fval, ObjVRML2 *obj)
{
    /*if (obj->num_textures==0) 
	{
		obj->textures = (float *) malloc (WRLCHUNK*sizeof(float));
		if (obj->textures==NULL)
			ErrMessage();
		else
			obj->size_textures = WRLCHUNK;
    } 
	else 
	{
		if (obj->num_textures==obj->size_textures) 
		{
			obj->textures = (float *) realloc(obj->textures, 
						(obj->size_textures + WRLCHUNK)*sizeof(float));
			if (obj->textures==NULL)
				ErrMessage();
			obj->size_textures += WRLCHUNK;
		}
    }

    obj->textures[obj->num_textures] = fval;
    obj->num_textures += 1;*/
}

void Wrl_Add_ColorIndex(int ival, ObjVRML2 *obj)
{
    /*if (obj->num_ccolors==0) 
	{
		obj->ccolors = (int *) malloc(WRLCHUNK*sizeof(int));
		if (obj->ccolors==NULL)
			ErrMessage();
		else
			obj->size_ccolors = WRLCHUNK;
    } 
	else 
	{
		if (obj->num_ccolors==obj->size_ccolors) 
		{
			obj->ccolors = (int *) realloc(obj->ccolors,
						(obj->size_ccolors + WRLCHUNK)*sizeof(int));
			if (obj->ccolors==NULL)
				ErrMessage();
			obj->size_ccolors += WRLCHUNK;
		}
    }

    obj->ccolors[obj->num_ccolors] = ival;
    obj->num_ccolors += 1;*/                                      
}

void Wrl_Add_CoordIndex(int ival, ObjVRML2 *obj)
{
    /*if (obj->num_triangles==0) 
	{
		obj->triangles = (int *) malloc(WRLCHUNK*sizeof(int));
		if (obj->triangles==NULL)
			ErrMessage();
		else
			obj->size_triangles = WRLCHUNK;
    }
	else 
	{
		if (obj->num_triangles==obj->size_triangles) 
		{
			obj->triangles = (int *) realloc(obj->triangles,
						(obj->size_triangles + WRLCHUNK)*sizeof(int));
			if (obj->triangles==NULL)
				ErrMessage();
			obj->size_triangles += WRLCHUNK;
		}
    }

    obj->triangles[obj->num_triangles] = ival;
    obj->num_triangles += 1;*/
}                              


void Wrl_Add_NormalIndex(int ival, ObjVRML2 *obj)
{
    /*if (obj->num_cnormals==0) 
	{
		obj->cnormals = (int *) malloc(WRLCHUNK*sizeof(int));
		if (obj->cnormals==NULL)
			ErrMessage();
		else
			obj->size_cnormals = WRLCHUNK;
    } 
	else 
	{
		if (obj->num_cnormals==obj->size_cnormals) 
		{
			obj->cnormals = (int *) realloc(obj->cnormals,
						(obj->size_cnormals + WRLCHUNK)*sizeof(int));
			if (obj->cnormals==NULL)
				ErrMessage();
			obj->size_cnormals += WRLCHUNK;
		}
    }

    obj->cnormals[obj->num_cnormals] = ival;
    obj->num_cnormals += 1;*/
}                                     

void Wrl_Add_TexCoordIndex(int ival, ObjVRML2 *obj)
{
    /*if (obj->num_ctextures==0) 
	{
		obj->ctextures = (int *) malloc(WRLCHUNK*sizeof(int));
		if (obj->ctextures==NULL)
			ErrMessage();
		else
			obj->size_ctextures = WRLCHUNK;
    }
	else 
	{
		if (obj->num_ctextures==obj->size_ctextures) 
		{
			obj->ctextures = (int *) realloc(obj->ctextures,
						(obj->size_ctextures + WRLCHUNK)*sizeof(int));
			if (obj->ctextures==NULL)
				ErrMessage();
			obj->size_ctextures += WRLCHUNK;
		}
    }

    obj->ctextures[obj->num_ctextures] = ival;
    obj->num_ctextures += 1;*/
}



//Parse some kinds of VRML 2.0 files
int WRL_Parser(FILE *fp, ObjVRML2 *obj)
{
	ObjVRML2_Init(obj);

    char header[16];
    memset(header,'\0',16);
    fscanf(fp,"%15c",header);
    //if (strcmp(header,VRML_HEADER)!=0)
     //   return 0;

    int  len_token, ns;
    char token[128];
      
    while (1) 
	{
        len_token = Get_A_Token(fp,token);
        if (len_token<=0) return 0;
		if (strcmp(token, "IndexedFaceSet")!=0)
            continue;

        len_token = Get_A_Token(fp,token);
        if (len_token<=0 || strcmp(token,"{"))
            return 0;

        while ((len_token=Get_A_Token(fp,token))>0) 
		{
            if (!strcmp(token,"creaseAngle") || strcmp(token, "solid")==0) 
                len_token = Get_A_Token(fp,token);
			else if (!strcmp(token,"color") || !strcmp(token,"coord")
                     || !strcmp(token,"normal") || !strcmp(token,"texCoord")) 
			{
				int type = 0;
                if (!strcmp(token,"color"))     type = 1;
                if (!strcmp(token,"coord"))     type = 2;
                if (!strcmp(token,"normal"))    type = 3;
                if (!strcmp(token,"texCoord"))  type = 4;

                len_token = Get_A_Token(fp,token);
                if (len_token<=0) return 0;
				if(!strcmp(token, "DEF"))
				{
					len_token = Get_A_Token(fp,token);
					len_token = Get_A_Token(fp,token);
				}
                if (type==1 && strcmp(token,"Color")) return 0;
                if (type==2 && strcmp(token,"Coordinate")) return 0;
                if (type==3 && strcmp(token,"Normal")) return 0;
                if (type==4 && strcmp(token,"TextureCoordinate")) return 0;

                len_token = Get_A_Token(fp,token);
                if (len_token<=0 || strcmp(token,"{")) return 0;

                len_token = Get_A_Token(fp,token);
                if (len_token <= 0 ) return 0; 
                if (type==1 && strcmp(token,"color")) return 0;
                if (type==2 && strcmp(token,"point")) return 0;
                if (type==3 && strcmp(token,"vector")) return 0; 
                if (type==4 && strcmp(token,"point")) return 0;

                len_token = Get_A_Token(fp,token);
                if (len_token<=0 || strcmp(token,"["))
					break;

                float fvalue;
                while(1) 
				{
                    len_token = Get_A_Token(fp,token); 
                    if (len_token==0) return 0;
                    ns = sscanf(token,"%f",&fvalue);
                    if (ns<=0)  break;
                    switch(type) 
					{
                        case 1: 
                            Wrl_Add_Color(fvalue, obj);
                            break;
                        case 2: 
                            Wrl_Add_Coord(fvalue, obj);
                            break;
                        case 3: 
                            Wrl_Add_Normal(fvalue, obj);
                            break;
                        case 4: 
                            Wrl_Add_TexCoord(fvalue, obj);
                            break;
					}
				}

                if (len_token<=0 || strcmp(token,"]")) 
				{
                    return 0;
                }

				len_token = Get_A_Token(fp,token);
				if (len_token<=0 || strcmp(token,"}"))
                    return 0;
			}  else if (!strcmp(token,"colorIndex") ||
                        !strcmp(token,"coordIndex") ||
                        !strcmp(token,"normalIndex") ||
                        !strcmp(token,"texCoordIndex")) 
			{
				int type = 0;
				if (!strcmp(token,"colorIndex"))     type = 1;
				if (!strcmp(token,"coordIndex"))     type = 2;
				if (!strcmp(token,"normalIndex"))    type = 3;
				if (!strcmp(token,"texCoordIndex"))  type = 4;
					len_token = Get_A_Token(fp,token);
				if (len_token<=0 || strcmp(token,"["))
					return 0;
				
				int ivalue;
				while (1) 
				{
					len_token = Get_A_Token(fp,token);
					if (len_token<=0)  break;
					ns = sscanf(token,"%d",&ivalue);
					if (ns<=0)  break;
					switch(type) {
						case 1:
							Wrl_Add_ColorIndex(ivalue, obj);
							break;
						case 2:
							Wrl_Add_CoordIndex(ivalue, obj);
							break;
						case 3:
							Wrl_Add_NormalIndex(ivalue, obj);
							break;
						case 4:
							Wrl_Add_TexCoordIndex(ivalue, obj);
							break;
					}
				}
				if ( len_token<=0 || strcmp(token,"]")) 
					return 0;
			} 
			else if(strcmp(token, "ccw") == 0 ||
					strcmp(token, "solid") == 0 ||
					strcmp(token, "convex") == 0 ||
					strcmp(token, "colorPerVertex") == 0 ||
					strcmp(token, "normalPerVertex") == 0) {
				len_token = Get_A_Token(fp,token);
			}
			else {
				break;
			}
		}
		break;
	}
	return 1;
}


void Wrl_Add_A_Polygon(int v0, int v1, int v2) 
{
	if ((no_of_tris%WRLCHUNK)==0) 
	{
        tris= (tempTri*) realloc(tris, (no_of_tris+WRLCHUNK)*sizeof(tempTri));
        if (tris==NULL) 
		{
            fprintf(stderr, "Wrl_Add_A_Polygon(): out of space\n");
            exit(1);
        }
    }

    tris[no_of_tris].v1 = v0;
    tris[no_of_tris].v2 = v1;
    tris[no_of_tris].v3 = v2;
    no_of_tris++;
}


void Wrl_Add_Polygons(int pidx, int polygons[]) 
{
    if (pidx<3)  return;

    if (pidx==3) 
	{
        Wrl_Add_A_Polygon(polygons[0], polygons[1], polygons[2]);
        return;
    }

    for (int i=1; i<pidx-1; i++) 
        Wrl_Add_A_Polygon(polygons[0], polygons[i], polygons[i+1]);
}



void ReadMeshFromWrl(FILE *fp)
{
    ObjVRML2  obj;
    
    int  flag, pidx, i;
    int  polygons[500];

    flag = WRL_Parser(fp, &obj);
    if (flag==0) 
	{
        fprintf(stderr, "ReadMeshFromWrl(): parse error of input data \n");
        fprintf(stderr, "Make sure it is either RAW or WRL format \n");
        exit(1);
    }

    no_of_tris=0;
	orig_verts=0;
    tris = (tempTri*) malloc(obj.num_triangles*sizeof(tempTri));
    if (tris==NULL) 
	{
		fprintf(stderr, "Out of space\n");
		exit(1);
    }

    i=0;
    while (i<obj.num_triangles) 
	{
        int  index;
        pidx = 0; 
        while (1) 
		{
            index = obj.triangles[i];
            i++;
            if (index==-1) 
			{
				if (pidx>2) 
					Wrl_Add_Polygons(pidx, polygons);
				break;
            } 
			else 
			{
				polygons[pidx] = index;
				pidx++;
            }
        }
    }

    free (obj.triangles);

    orig_verts = obj.num_positions/3;
    verts = (tempVert*) malloc(orig_verts*sizeof(tempVert));
    if (verts==NULL) 
	{
        fprintf(stderr, "ReadMeshFromWrl(): out of space\n");
        exit(1);
    }

    for (i=0; i<obj.num_positions; i=i+3) 
	{
        pidx = i/3;
        verts[pidx].x = obj.positions[i];
        verts[pidx].y = obj.positions[i+1];
        verts[pidx].z = obj.positions[i+2];
    }
    free (obj.positions);
}
